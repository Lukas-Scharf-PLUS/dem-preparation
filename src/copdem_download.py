"""
Download the Copernicus DEM GLO-30 for an area, as one GeoTIFF.

Only used when ``coarse_dem.source: download``, i.e. when the local
``cop_dem.tif`` does not cover all bursts. The tiles come from the public AWS
Open Data bucket (no account or credentials needed):

    https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N48_00_E016_00_DEM/...

Each tile covers 1 x 1 degree and is a Cloud-Optimized GeoTIFF; GDAL reads
only the requested area. Heights are above the EGM2008 geoid; the merge step
converts them to ellipsoid heights.
"""

import logging
import math

from osgeo import gdal

logger = logging.getLogger(__name__)
gdal.UseExceptions()

COPDEM_BASE_URL = "https://copernicus-dem-30m.s3.amazonaws.com"
COPDEM_VERTICAL_DATUM = "egm2008"
REMOTE_READ_OPTIONS = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "3",
}
LOCAL_GEOTIFF_OPTIONS = ["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER"]


class CopDemDownloadError(RuntimeError):
    """Raised when no Copernicus tile could be found for the area."""


def tile_name(latitude, longitude):
    """Return the Copernicus tile name for the 1x1 degree tile with lower-left corner (lat, lon)."""
    lat_tag = f"{'N' if latitude >= 0 else 'S'}{abs(latitude):02d}_00"
    lon_tag = f"{'E' if longitude >= 0 else 'W'}{abs(longitude):03d}_00"
    return f"Copernicus_DSM_COG_10_{lat_tag}_{lon_tag}_DEM"


def tiles_for_bounds(bounds):
    """Return (lat, lon) lower-left corners of all 1x1 degree tiles touching ``bounds`` (W, S, E, N)."""
    west, south, east, north = bounds
    return [(lat, lon)
            for lat in range(math.floor(south), math.ceil(north))
            for lon in range(math.floor(west), math.ceil(east))]


def tile_url(latitude, longitude):
    """Return the /vsicurl/ path of one Copernicus tile."""
    name = tile_name(latitude, longitude)
    return f"/vsicurl/{COPDEM_BASE_URL}/{name}/{name}.tif"


def copdem_output_path(bounds, work_dir):
    """Return where the GLO-30 download for ``bounds`` is stored (one file per area)."""
    west, south, east, north = bounds
    return work_dir / "copdem" / f"copdem_glo30_egm2008_W{west:.2f}_S{south:.2f}_E{east:.2f}_N{north:.2f}.tif"


def download_copdem(bounds, work_dir):
    """Download GLO-30 for ``bounds`` into one GeoTIFF and return its path.

    Tiles that do not exist (open sea) are skipped. An existing file is reused.

    Raises
    ------
    CopDemDownloadError
        If no tile exists for the area.
    """
    out_path = copdem_output_path(bounds, work_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    west, south, east, north = bounds
    if out_path.is_file():
        logger.info("Copernicus DEM already downloaded: %s", out_path)
        return out_path

    with gdal.config_options(REMOTE_READ_OPTIONS):
        available = []
        for latitude, longitude in tiles_for_bounds(bounds):
            url = tile_url(latitude, longitude)
            if gdal.VSIStatL(url) is None:
                logger.warning("No Copernicus tile %s (open sea?), skipped.", tile_name(latitude, longitude))
                continue
            available.append(url)
        if not available:
            raise CopDemDownloadError(f"No Copernicus GLO-30 tile found for {bounds}.")

        logger.info("Downloading Copernicus GLO-30 (%d tiles) for %s", len(available), bounds)
        # resolution="highest": tiles above 50 deg N have a coarser longitude spacing.
        mosaic = gdal.BuildVRT("/vsimem/copdem_mosaic.vrt", available, resolution="highest")
        part_path = out_path.with_name(out_path.name + ".part")
        gdal.Translate(str(part_path), mosaic, format="GTiff", projWin=[west, north, east, south],
                       creationOptions=LOCAL_GEOTIFF_OPTIONS, callback=gdal.TermProgress_nocb)
        mosaic = None
        gdal.Unlink("/vsimem/copdem_mosaic.vrt")
    part_path.rename(out_path)
    logger.info("Copernicus DEM written: %s", out_path)
    return out_path
