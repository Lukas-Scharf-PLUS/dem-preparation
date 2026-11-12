"""
Merge the BEV ALS DSM (fine, inside the AOI) and the Copernicus DEM (coarse,
everywhere else) into one DEM with WGS84 ellipsoid heights.

Steps
-----
1. Coarse part: cut the coarse DEM to the output area and convert it to
   ellipsoid heights, at its own 30 m resolution.
2. Fine part (only the BEV area = AOI + buffer, on the output grid):
   a. warp the 1 m BEV data (EPSG:3035) to the lat/lon grid by *averaging*
      all 1 m pixels in each output pixel, shifted from ETRS89 to ITRF,
   b. add the height offset of the BEV height reference (and the small
      ETRS89 -> ITRF height change),
   c. interpolate small no-data holes from the surrounding BEV heights,
   d. blend into the coarse heights over ``blend_width_m`` at the edge of the
      BEV area, so there is no height step between the two DEMs.
3. Mixed-resolution VRT: coarse part below, fine part on top, presented on
   the output grid (GDAL resamples the coarse part bilinearly on read).
4. Final GeoTIFF: the VRT written out as one compressed file.

Why one grid: ISCE2 and COMPASS both read a single regular raster. A raster
cannot have two pixel sizes, so outside the AOI the 30 m heights are just
interpolated to the fine grid. The VRT stores them at 30 m instead, which
saves disk space; GDAL-based readers (COMPASS, isce2-dolphin-mintpy) see the
same grid either way.
"""

import json
import logging
from dataclasses import asdict, dataclass

import numpy as np
from osgeo import gdal
from scipy import ndimage

from target_grid import GeoGrid, meters_per_degree
from vertical_datum import ELLIPSOID, height_offset_function, offset_grid

logger = logging.getLogger(__name__)
gdal.UseExceptions()

OUTPUT_NODATA = -32768.0
BEV_NODATA = -9999.0
GEOTIFF_OPTIONS = ["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER",
                   "NUM_THREADS=ALL_CPUS"]
WARP_MEMORY_MB = 1024
ROWS_PER_BLOCK = 1024
# Spacing of the grid on which height offsets are evaluated before interpolation.
# The BEV height grids have structure at the km scale, EGM2008 is smoother.
FINE_OFFSET_STEP_DEG = 0.005
COARSE_OFFSET_STEP_DEG = 0.02
# Coarse pixels added around the output area so bilinear resampling has neighbors.
COARSE_MARGIN_PIXELS = 2
# If BEV and Copernicus differ more than this (median, meters) inside the AOI,
# the height reference of one of them is probably wrong.
MAX_PLAUSIBLE_MEDIAN_DIFFERENCE_M = 3.0
WGS84_SRS = "EPSG:4326"
ETRS89_SRS = "EPSG:4258"
BEV_SRS = "EPSG:3035"


@dataclass
class DifferenceStatistics:
    """BEV minus coarse heights (meters) inside the AOI, as a plausibility check."""

    median: float
    percentile_5: float
    percentile_95: float
    pixel_count: int


def write_array(path, array, grid, nodata=OUTPUT_NODATA, metadata=None):
    """Write a float32 array on ``grid`` as a compressed GeoTIFF (EPSG:4326)."""
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(str(path), grid.width, grid.height, 1, gdal.GDT_Float32,
                            options=GEOTIFF_OPTIONS)
    dataset.SetGeoTransform(grid.geotransform)
    dataset.SetProjection(WGS84_SRS)
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(np.where(np.isfinite(array), array, nodata).astype(np.float32))
    if metadata:
        dataset.SetMetadata(metadata)
    dataset.FlushCache()
    dataset = None


def read_valid(dataset, first_row, rows, nodata):
    """Read rows of band 1 as float32 with no-data pixels set to NaN."""
    block = dataset.GetRasterBand(1).ReadAsArray(0, first_row, dataset.RasterXSize, rows).astype(np.float32)
    if nodata is not None:
        block[block == nodata] = np.nan
    return block


# --------------------------------------------------------------------------- coarse part

def prepare_coarse_part(coarse_file, vertical_datum, grid, out_path):
    """Cut the coarse DEM to the output area and convert it to ellipsoid heights.

    The result keeps the coarse resolution; no-data is written as OUTPUT_NODATA.

    Parameters
    ----------
    coarse_file : Path
        Lat/lon DEM covering ``grid``.
    vertical_datum : str
        'ellipsoid', 'egm2008' or 'egm96'.
    grid : GeoGrid
        Output grid.
    out_path : Path
    """
    source = gdal.Open(str(coarse_file))
    source_nodata = source.GetRasterBand(1).GetNoDataValue()
    _, pixel_lon, _, _, _, pixel_lat = source.GetGeoTransform()
    margin_lon = COARSE_MARGIN_PIXELS * pixel_lon
    margin_lat = COARSE_MARGIN_PIXELS * abs(pixel_lat)
    cropped = gdal.Translate(
        "/vsimem/coarse_crop.tif", source,
        projWin=[grid.west - margin_lon, grid.north + margin_lat,
                 grid.east + margin_lon, grid.south - margin_lat],
    )
    width, height = cropped.RasterXSize, cropped.RasterYSize
    geotransform = cropped.GetGeoTransform()
    offset_function = None if vertical_datum == ELLIPSOID else height_offset_function(vertical_datum)

    output = gdal.GetDriverByName("GTiff").Create(str(out_path), width, height, 1, gdal.GDT_Float32,
                                                  options=GEOTIFF_OPTIONS)
    output.SetGeoTransform(geotransform)
    output.SetProjection(WGS84_SRS)
    out_band = output.GetRasterBand(1)
    out_band.SetNoDataValue(OUTPUT_NODATA)

    longitudes = geotransform[0] + (np.arange(width) + 0.5) * geotransform[1]
    for first_row in range(0, height, ROWS_PER_BLOCK):
        rows = min(ROWS_PER_BLOCK, height - first_row)
        heights = read_valid(cropped, first_row, rows, source_nodata)
        if offset_function is not None:
            latitudes = geotransform[3] + (first_row + np.arange(rows) + 0.5) * geotransform[5]
            heights += offset_grid(latitudes, longitudes, offset_function, COARSE_OFFSET_STEP_DEG)
        out_band.WriteArray(np.where(np.isfinite(heights), heights, OUTPUT_NODATA), 0, first_row)

    output.SetMetadata({"HEIGHT_REFERENCE": "WGS84 ellipsoid",
                        "SOURCE": str(coarse_file), "SOURCE_HEIGHT_REFERENCE": vertical_datum})
    output.FlushCache()
    output = None
    gdal.Unlink("/vsimem/coarse_crop.tif")
    logger.info("Coarse part (%s -> ellipsoid, %d x %d pixels): %s", vertical_datum, width, height, out_path)


# --------------------------------------------------------------------------- fine part

def warp_to_grid(source, grid, source_srs, resample, source_nodata, shift_lon=0.0, shift_lat=0.0,
                 target_srs=WGS84_SRS):
    """Warp a raster onto ``grid`` and return a float32 array (no-data = NaN).

    ``shift_lon``/``shift_lat`` move the *source* positions: a source point at
    (lon, lat) ends up at (lon + shift_lon, lat + shift_lat) on the grid.
    """
    west, south, east, north = grid.bounds
    warped = gdal.Warp(
        "", str(source), format="MEM",
        srcSRS=source_srs, dstSRS=target_srs,
        outputBounds=(west - shift_lon, south - shift_lat, east - shift_lon, north - shift_lat),
        width=grid.width, height=grid.height,
        resampleAlg=resample, srcNodata=source_nodata, dstNodata=np.nan,
        outputType=gdal.GDT_Float32, warpMemoryLimit=WARP_MEMORY_MB, multithread=True,
    )
    return warped.GetRasterBand(1).ReadAsArray().astype(np.float32)


def fill_small_holes(heights, max_distance_px):
    """Interpolate no-data holes inside the data that are at most ``max_distance_px`` wide.

    Only holes surrounded by data are filled; no-data touching the array edge
    (outside the BEV area) is left alone, so the data is not extrapolated.

    Returns
    -------
    np.ndarray
        Copy of ``heights`` with small holes filled (NaN elsewhere unchanged).
    """
    invalid = ~np.isfinite(heights)
    if not invalid.any() or max_distance_px <= 0:
        return heights.copy()
    labels, _ = ndimage.label(invalid)
    edge_labels = np.unique(np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]]))
    interior_holes = invalid & ~np.isin(labels, edge_labels)
    if not interior_holes.any():
        return heights.copy()

    dataset = gdal.GetDriverByName("MEM").Create("", heights.shape[1], heights.shape[0], 1, gdal.GDT_Float32)
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(OUTPUT_NODATA)
    band.WriteArray(np.where(invalid, OUTPUT_NODATA, heights))
    gdal.FillNodata(band, None, maxSearchDist=float(max_distance_px), smoothingIterations=0)
    filled = band.ReadAsArray().astype(np.float32)
    filled[filled == OUTPUT_NODATA] = np.nan

    result = heights.copy()
    result[interior_holes] = filled[interior_holes]
    logger.info("Filled %d of %d BEV no-data pixels inside the data.",
                int(np.isfinite(result[interior_holes]).sum()), int(interior_holes.sum()))
    return result


def blend_weights(valid, blend_width_m, pixel_height_m, pixel_width_m):
    """Return BEV weights 0..1: rising linearly from the edge of the valid area over ``blend_width_m``.

    The array border counts as edge, so the weight also drops towards the
    border of the BEV area.
    """
    padded = np.pad(valid, 1, constant_values=False)
    distance = ndimage.distance_transform_edt(padded, sampling=(pixel_height_m, pixel_width_m))[1:-1, 1:-1]
    return np.clip(distance / blend_width_m, 0.0, 1.0).astype(np.float32)


def blend(fine, coarse, weights):
    """Return ``weights * fine + (1 - weights) * coarse``, using coarse where fine is NaN."""
    merged = np.where(np.isfinite(fine), weights * fine + (1.0 - weights) * coarse, coarse)
    return merged.astype(np.float32)


def aoi_mask(grid, aoi_bbox):
    """Return a boolean array marking the pixels of ``grid`` whose center is inside the AOI box."""
    latitudes, longitudes = (np.asarray(values) for values in grid.pixel_center_coordinates())
    west, south, east, north = aoi_bbox
    rows = (latitudes >= south) & (latitudes <= north)
    cols = (longitudes >= west) & (longitudes <= east)
    return rows[:, None] & cols[None, :]


def difference_statistics(fine, coarse, mask):
    """Return DifferenceStatistics of ``fine - coarse`` where both are valid and ``mask`` is True."""
    difference = (fine - coarse)[mask & np.isfinite(fine) & np.isfinite(coarse)]
    if difference.size == 0:
        return DifferenceStatistics(float("nan"), float("nan"), float("nan"), 0)
    p5, median, p95 = np.percentile(difference, [5, 50, 95])
    return DifferenceStatistics(float(median), float(p5), float(p95), int(difference.size))


def prepare_fine_part(bev_mosaic, bev_height_reference, frame_shift, fine_grid, coarse_part,
                      aoi_bbox, merge_settings, out_path):
    """Create the fine part: BEV heights on ``fine_grid``, ellipsoidal, blended into the coarse DEM.

    Parameters
    ----------
    bev_mosaic : Path
        VRT of the downloaded BEV tiles (EPSG:3035).
    bev_height_reference : str
        'evrf2000_austria' or 'gha'.
    frame_shift : FrameShift
        ETRS89 -> output frame shift.
    fine_grid : GeoGrid
        Sub-grid of the output grid covering the BEV area.
    coarse_part : Path
        Output of :func:`prepare_coarse_part`.
    aoi_bbox : tuple
    merge_settings : MergeSettings
    out_path : Path

    Returns
    -------
    DifferenceStatistics
        BEV minus coarse inside the AOI (before blending).
    """
    logger.info("Warping BEV DSM to the %d x %d output grid (averaging all source pixels per output pixel) ...",
                fine_grid.width, fine_grid.height)
    bev = warp_to_grid(bev_mosaic, fine_grid, BEV_SRS, "average", BEV_NODATA,
                       shift_lon=frame_shift.delta_lon_deg, shift_lat=frame_shift.delta_lat_deg,
                       target_srs=ETRS89_SRS)
    latitudes, longitudes = fine_grid.pixel_center_coordinates()
    offset = offset_grid(latitudes, longitudes, height_offset_function(bev_height_reference),
                         FINE_OFFSET_STEP_DEG)
    logger.info("BEV height offset (%s -> ellipsoid): %.3f .. %.3f m, frame shift dh %.3f m",
                bev_height_reference, float(offset.min()), float(offset.max()), frame_shift.delta_height_m)
    bev = bev + offset + frame_shift.delta_height_m

    coarse = warp_to_grid(coarse_part, fine_grid, WGS84_SRS, "bilinear", OUTPUT_NODATA)
    statistics = difference_statistics(bev, coarse, aoi_mask(fine_grid, aoi_bbox))

    meters_per_deg_lat, meters_per_deg_lon = meters_per_degree(fine_grid.north - fine_grid.height
                                                               * fine_grid.delta_lat / 2)
    pixel_height_m = fine_grid.delta_lat * meters_per_deg_lat
    pixel_width_m = fine_grid.delta_lon * meters_per_deg_lon
    hole_fill_px = merge_settings.hole_fill_max_m / min(pixel_height_m, pixel_width_m)
    bev = fill_small_holes(bev, hole_fill_px)

    weights = blend_weights(np.isfinite(bev), merge_settings.blend_width_m, pixel_height_m, pixel_width_m)
    merged = blend(bev, coarse, weights)
    write_array(out_path, merged, fine_grid, metadata={
        "HEIGHT_REFERENCE": "WGS84 ellipsoid",
        "SOURCE": f"BEV ALS DSM ({bev_mosaic}), blended into the coarse DEM",
        "SOURCE_HEIGHT_REFERENCE": bev_height_reference,
        "HORIZONTAL_FRAME": merge_settings.horizontal_frame,
        "FRAME_EPOCH": str(merge_settings.frame_epoch),
    })
    logger.info("Fine part (%d x %d pixels): %s", fine_grid.width, fine_grid.height, out_path)
    return statistics


# --------------------------------------------------------------------------- final outputs

def build_mixed_resolution_vrt(coarse_part, fine_part, grid, out_path):
    """Write a VRT on ``grid``: coarse part below (bilinear), fine part on top.

    Source paths are stored relative to the VRT, so the output folder can be moved as a whole.
    """
    gdal.BuildVRT(
        str(out_path), [str(coarse_part), str(fine_part)],  # later files are drawn on top
        outputBounds=grid.bounds, xRes=grid.delta_lon, yRes=grid.delta_lat,
        resampleAlg="bilinear", srcNodata=OUTPUT_NODATA, VRTNodata=OUTPUT_NODATA,
    )
    logger.info("Mixed-resolution VRT: %s", out_path)


def write_final_geotiff(vrt_path, out_path, metadata):
    """Write the merged VRT as one compressed GeoTIFF on the output grid."""
    logger.info("Writing the final GeoTIFF %s ...", out_path)
    part_path = out_path.with_name(out_path.name + ".part")
    gdal.Translate(str(part_path), str(vrt_path), format="GTiff", creationOptions=GEOTIFF_OPTIONS,
                   callback=gdal.TermProgress_nocb)
    dataset = gdal.Open(str(part_path), gdal.GA_Update)
    dataset.SetMetadata(metadata)
    dataset = None
    part_path.rename(out_path)


def count_nodata(raster_path):
    """Return the number of no-data pixels in a raster (read block by block)."""
    dataset = gdal.Open(str(raster_path))
    band = dataset.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    count = 0
    for first_row in range(0, dataset.RasterYSize, ROWS_PER_BLOCK):
        rows = min(ROWS_PER_BLOCK, dataset.RasterYSize - first_row)
        block = band.ReadAsArray(0, first_row, dataset.RasterXSize, rows)
        count += int((~np.isfinite(block) | (block == nodata)).sum())
    return count


def write_report(path, report):
    """Write the run report (dict) as pretty JSON."""
    def default(value):
        if isinstance(value, GeoGrid) or hasattr(value, "__dataclass_fields__"):
            return asdict(value)
        return str(value)

    path.write_text(json.dumps(report, indent=2, default=default), encoding="utf-8")
    logger.info("Report: %s", path)


def check_difference(statistics):
    """Log the BEV - coarse difference and warn if it suggests a wrong height reference."""
    logger.info("BEV - coarse DEM inside the AOI: median %.2f m, 5%% %.2f m, 95%% %.2f m (%d pixels)",
                statistics.median, statistics.percentile_5, statistics.percentile_95,
                statistics.pixel_count)
    if not abs(statistics.median) <= MAX_PLAUSIBLE_MEDIAN_DIFFERENCE_M:
        logger.warning(
            "The median difference is larger than %.0f m. Check coarse_dem.vertical_datum and "
            "fine_dem.vertical_datum: one of the height references is probably wrong.",
            MAX_PLAUSIBLE_MEDIAN_DIFFERENCE_M)
