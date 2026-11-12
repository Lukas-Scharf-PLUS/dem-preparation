"""
Command-line entry point: build the Vienna DEM for ISCE2 and COMPASS.

Steps (run one, or ``all`` in this order):

    check     find the bursts that touch the AOI and check that the coarse DEM covers them
    download  download the BEV ALS DSM for the AOI (and Copernicus GLO-30 if
              coarse_dem.source is 'download')
    build     merge both DEMs on the output grid and write the outputs

Examples
--------
    python src/prepare_dem.py --config configs/vienna.yaml all
    python src/prepare_dem.py --config configs/vienna.yaml check
    python src/prepare_dem.py --config configs/vienna.yaml --resolution_m 5 build
"""

import argparse
import logging
import sys
from datetime import datetime

from burst_coverage import (CoverageError, check_dem_covers, find_bursts, format_bounds,
                            log_burst_summary, required_dem_bounds)
from bev_download import MOSAIC_FILE_NAME, BevDownloadError, aoi_window_3035, download_bev_dsm, window_lonlat_bounds
from config import ConfigError, load_settings
from copdem_download import COPDEM_VERTICAL_DATUM, CopDemDownloadError, copdem_output_path, download_copdem
from dem_merge import (build_mixed_resolution_vrt, check_difference, count_nodata, prepare_coarse_part,
                       prepare_fine_part, write_final_geotiff, write_report)
from isce_export import export_isce_dem
from target_grid import output_grid
from vertical_datum import VerticalDatumError, bev_height_reference, enable_proj_network, etrs89_to_itrf_shift

logger = logging.getLogger("prepare_dem")
STEPS = ("check", "download", "build", "all")
# Errors with a clear message for the user; anything else is a bug and shows a traceback.
EXPECTED_ERRORS = (CoverageError, BevDownloadError, CopDemDownloadError, VerticalDatumError)


def parse_arguments(argv=None):
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=STEPS, help="what to run")
    parser.add_argument("--config", required=True, help="YAML config file, e.g. configs/vienna.yaml")
    parser.add_argument("--resolution_m", type=float, help="override merge.resolution_m")
    parser.add_argument("--edition", help="override fine_dem.edition (BEV reference date, e.g. 20250915)")
    parser.add_argument("--coarse_source", choices=("file", "download"), help="override coarse_dem.source")
    return parser.parse_args(argv)


def setup_logging(log_file):
    """Log to the console and to ``log_file``."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, encoding="utf-8")],
    )


def aoi_center(bbox):
    """Return (lon, lat) of the AOI center."""
    west, south, east, north = bbox
    return (west + east) / 2, (south + north) / 2


def run_check(settings):
    """Find the bursts and check the coarse DEM. Returns (bursts_per_orbit, required_bounds).

    Raises
    ------
    CoverageError
        If coarse_dem.source is 'file' and that DEM does not cover the bursts.
    """
    bursts_per_orbit = find_bursts(settings.bursts.burst_db, settings.bursts.relative_orbits,
                                   settings.aoi.bbox)
    required_bounds = required_dem_bounds(bursts_per_orbit, settings.bursts.margin_deg)
    log_burst_summary(bursts_per_orbit, required_bounds)

    if settings.coarse_dem.source == "download":
        logger.info("coarse_dem.source is 'download': Copernicus GLO-30 is downloaded for this area.")
        return bursts_per_orbit, required_bounds

    covers, dem_bounds = check_dem_covers(settings.coarse_dem.file, required_bounds)
    logger.info("Coarse DEM %s covers %s", settings.coarse_dem.file, format_bounds(dem_bounds))
    if not covers:
        raise CoverageError(
            f"The coarse DEM does NOT cover all bursts.\n  DEM:    {format_bounds(dem_bounds)}\n"
            f"  needed: {format_bounds(required_bounds)}\n"
            "Set 'coarse_dem.source: download' (or --coarse_source download) to download "
            "Copernicus GLO-30 for the needed area from AWS (no account needed)."
        )
    logger.info("OK: the coarse DEM covers all bursts of orbits %s.", settings.bursts.relative_orbits)
    return bursts_per_orbit, required_bounds


def coarse_dem_source(settings, required_bounds, download_if_missing):
    """Return (coarse DEM file, its height reference) according to coarse_dem.source."""
    if settings.coarse_dem.source == "file":
        return settings.coarse_dem.file, settings.coarse_dem.vertical_datum
    path = copdem_output_path(required_bounds, settings.work_dir)
    if not path.is_file():
        if not download_if_missing:
            raise CoverageError(f"Copernicus DEM not downloaded yet ({path}); run the 'download' step.")
        path = download_copdem(required_bounds, settings.work_dir)
    return path, COPDEM_VERTICAL_DATUM


def run_download(settings, required_bounds):
    """Download the BEV DSM window (and GLO-30 if configured). Returns the BEV mosaic path."""
    coarse_dem_source(settings, required_bounds, download_if_missing=True)
    return download_bev_dsm(settings.aoi.bbox, settings.fine_dem.buffer_m, settings.fine_dem.edition,
                            settings.fine_dem.overview_level, settings.work_dir)


def run_build(settings, bursts_per_orbit, required_bounds):
    """Merge both DEMs and write all configured outputs."""
    coarse_file, coarse_datum = coarse_dem_source(settings, required_bounds, download_if_missing=False)
    bev_mosaic = settings.work_dir / "bev" / MOSAIC_FILE_NAME
    if not bev_mosaic.is_file():
        raise CoverageError(f"BEV DSM not downloaded yet ({bev_mosaic}); run the 'download' step.")

    center_lon, center_lat = aoi_center(settings.aoi.bbox)
    grid = output_grid(required_bounds, settings.merge.resolution_m, center_lat)
    bev_window = aoi_window_3035(settings.aoi.bbox, settings.fine_dem.buffer_m)
    col_offset, row_offset, fine_grid = grid.subgrid(window_lonlat_bounds(bev_window))
    logger.info("Output grid: %d x %d pixels, %.3e x %.3e deg (~%g m), %s", grid.width, grid.height,
                grid.delta_lon, grid.delta_lat, settings.merge.resolution_m, format_bounds(grid.bounds))

    base = settings.output_base_name
    parts_dir = settings.output_dir / f"{base}_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    coarse_part = parts_dir / "coarse_ellipsoid.tif"
    fine_part = parts_dir / f"fine_{settings.merge.resolution_m:g}m_ellipsoid.tif"

    prepare_coarse_part(coarse_file, coarse_datum, grid, coarse_part)

    fine_reference = bev_height_reference(settings.fine_dem.edition, settings.fine_dem.vertical_datum)
    frame_shift = etrs89_to_itrf_shift(settings.merge.horizontal_frame, settings.merge.frame_epoch,
                                       center_lon, center_lat)
    logger.info("Horizontal frame ETRS89 -> %s @ %.1f: shift %.2e deg lon, %.2e deg lat, %.3f m height",
                settings.merge.horizontal_frame, settings.merge.frame_epoch, frame_shift.delta_lon_deg,
                frame_shift.delta_lat_deg, frame_shift.delta_height_m)
    statistics = prepare_fine_part(bev_mosaic, fine_reference, frame_shift, fine_grid, coarse_part,
                                   settings.aoi.bbox, settings.merge, fine_part)
    check_difference(statistics)

    mixed_vrt = settings.output_dir / f"{base}_mixed_resolution.vrt"
    build_mixed_resolution_vrt(coarse_part, fine_part, grid, mixed_vrt)
    metadata = {
        "HEIGHT_REFERENCE": "WGS84 ellipsoid (EPSG:4979 heights)",
        "COARSE_SOURCE": f"{coarse_file} ({coarse_datum})",
        "FINE_SOURCE": f"BEV ALS DSM edition {settings.fine_dem.edition} ({fine_reference})",
        "HORIZONTAL_FRAME": f"{settings.merge.horizontal_frame} @ {settings.merge.frame_epoch}",
        "CREATED_BY": "dem_preparation/src/prepare_dem.py",
    }

    final_dem = mixed_vrt
    if settings.output.write_geotiff:
        final_dem = settings.output_dir / f"{base}.tif"
        write_final_geotiff(mixed_vrt, final_dem, metadata)
        logger.info("Final DEM: %s", final_dem)
    nodata_pixels = count_nodata(final_dem)
    if nodata_pixels:
        logger.warning("The final DEM has %d no-data pixels.", nodata_pixels)

    isce_dem = None
    if settings.output.write_isce:
        isce_dem = export_isce_dem(final_dem, settings.output_dir / base)
    if not settings.output.write_mixed_resolution_vrt and settings.output.write_geotiff:
        mixed_vrt.unlink()
        logger.info("Mixed-resolution VRT removed (output.write_mixed_resolution_vrt is false); "
                    "the parts folder is kept as intermediate data.")

    write_report(settings.output_dir / f"{base}_report.json", {
        "created": datetime.now().isoformat(timespec="seconds"),
        "config_file": settings.config_file,
        "aoi_bbox": settings.aoi.bbox,
        "bursts": {orbit: [burst.burst_id for burst in bursts] for orbit, bursts in bursts_per_orbit.items()},
        "required_bounds": required_bounds,
        "grid": grid,
        "fine_grid": {"grid": fine_grid, "col_offset": col_offset, "row_offset": row_offset},
        "coarse_dem": {"file": coarse_file, "height_reference": coarse_datum},
        "fine_dem": {"edition": settings.fine_dem.edition, "height_reference": fine_reference,
                     "overview_level": settings.fine_dem.overview_level},
        "frame_shift": frame_shift,
        "merge": settings.merge,
        "bev_minus_coarse_in_aoi_m": statistics,
        "nodata_pixels": nodata_pixels,
        "outputs": {"geotiff": final_dem if settings.output.write_geotiff else None,
                    "mixed_resolution_vrt": mixed_vrt if settings.output.write_mixed_resolution_vrt else None,
                    "isce": isce_dem},
    })


def main(argv=None):
    """Run the requested step(s). Returns the process exit code."""
    arguments = parse_arguments(argv)
    try:
        settings = load_settings(arguments.config, {
            "merge.resolution_m": arguments.resolution_m,
            "fine_dem.edition": arguments.edition,
            "coarse_dem.source": arguments.coarse_source,
        })
    except ConfigError as error:
        print(f"Config error: {error}", file=sys.stderr)
        return 2

    setup_logging(settings.output_dir / f"{settings.output_base_name}.log")
    logger.info("Step '%s' with config %s", arguments.step, settings.config_file)
    enable_proj_network()
    try:
        bursts_per_orbit, required_bounds = run_check(settings)
        if arguments.step in ("download", "all"):
            run_download(settings, required_bounds)
        if arguments.step in ("build", "all"):
            run_build(settings, bursts_per_orbit, required_bounds)
    except EXPECTED_ERRORS as error:
        logger.error("%s", error)
        return 1
    logger.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
