"""
Download the BEV ALS DSM (Digitales Oberflächenmodell, 1 m) for the AOI.

The BEV publishes the DSM as free open data in 50 km x 50 km tiles in
EPSG:3035 (ETRS89 / LAEA Europe), one folder per edition (reference date):

    https://data.bev.gv.at/download/ALS/DSM/<edition>/ALS_DSM_CRS3035RES50000mN<n>E<e>.tif

``N<n>E<e>`` is the lower-left corner of the tile in meters. Each tile is a
6-12 GB Cloud-Optimized GeoTIFF (COG). Instead of downloading whole tiles,
GDAL reads only the part that covers the AOI (plus a buffer) through HTTP range
requests (``/vsicurl/``), so ~2-3 GB are transferred for Vienna instead of
~35 GB. No account or credentials are needed.

The result is one local GeoTIFF per tile (only the AOI part) and a VRT that
mosaics them. Finished files are kept, so a second run downloads nothing.
"""

import logging
import math

import numpy as np
from osgeo import gdal
from pyproj import Transformer

logger = logging.getLogger(__name__)
gdal.UseExceptions()

BEV_DSM_BASE_URL = "https://data.bev.gv.at/download/ALS/DSM"
# Newer editions use the "ALS_DSM_" prefix, the 2019 edition does not.
TILE_FILE_PATTERNS = (
    "ALS_DSM_CRS3035RES50000mN{north}E{east}.tif",
    "CRS3035RES50000mN{north}E{east}.tif",
)
TILE_SIZE_M = 50_000
BEV_EPSG = 3035
WGS84_EPSG = 4326
MOSAIC_FILE_NAME = "bev_dsm_mosaic.vrt"

# GDAL settings for reading remote COGs robustly and efficiently.
REMOTE_READ_OPTIONS = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",  # do not list the remote folder
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "GDAL_HTTP_MULTIRANGE": "YES",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "3",
    "GDAL_HTTP_TIMEOUT": "120",
}
LOCAL_GEOTIFF_OPTIONS = ["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER",
                         "NUM_THREADS=ALL_CPUS"]
# Points per edge when reprojecting the AOI box to EPSG:3035 (edges become curves).
EDGE_SAMPLES = 50


class BevDownloadError(RuntimeError):
    """Raised when a BEV tile cannot be found or read."""


def aoi_window_3035(aoi_bbox, buffer_m):
    """Return the AOI box in EPSG:3035 meters, grown by ``buffer_m``, on whole meters.

    Returns
    -------
    tuple
        (xmin, ymin, xmax, ymax) in meters.
    """
    west, south, east, north = aoi_bbox
    transformer = Transformer.from_crs(WGS84_EPSG, BEV_EPSG, always_xy=True)
    edge = np.linspace(0.0, 1.0, EDGE_SAMPLES)
    longitudes = np.concatenate([west + edge * (east - west), np.full(EDGE_SAMPLES, east),
                                 east - edge * (east - west), np.full(EDGE_SAMPLES, west)])
    latitudes = np.concatenate([np.full(EDGE_SAMPLES, south), south + edge * (north - south),
                                np.full(EDGE_SAMPLES, north), north - edge * (north - south)])
    x, y = transformer.transform(longitudes, latitudes)
    return (
        math.floor(x.min() - buffer_m),
        math.floor(y.min() - buffer_m),
        math.ceil(x.max() + buffer_m),
        math.ceil(y.max() + buffer_m),
    )


def window_lonlat_bounds(window):
    """Return the (west, south, east, north) box enclosing an EPSG:3035 window."""
    xmin, ymin, xmax, ymax = window
    transformer = Transformer.from_crs(BEV_EPSG, WGS84_EPSG, always_xy=True)
    edge = np.linspace(0.0, 1.0, EDGE_SAMPLES)
    x = np.concatenate([xmin + edge * (xmax - xmin), np.full(EDGE_SAMPLES, xmax),
                        xmax - edge * (xmax - xmin), np.full(EDGE_SAMPLES, xmin)])
    y = np.concatenate([np.full(EDGE_SAMPLES, ymin), ymin + edge * (ymax - ymin),
                        np.full(EDGE_SAMPLES, ymax), ymax - edge * (ymax - ymin)])
    longitudes, latitudes = transformer.transform(x, y)
    return (longitudes.min(), latitudes.min(), longitudes.max(), latitudes.max())


def tiles_for_window(window):
    """Return the lower-left corners (north, east) of all 50 km tiles touching the window."""
    xmin, ymin, xmax, ymax = window
    easts = range(math.floor(xmin / TILE_SIZE_M) * TILE_SIZE_M, xmax, TILE_SIZE_M)
    norths = range(math.floor(ymin / TILE_SIZE_M) * TILE_SIZE_M, ymax, TILE_SIZE_M)
    return [(north, east) for north in norths for east in easts]


def clip_window_to_tile(window, tile):
    """Return the part of ``window`` inside tile (north, east) as (xmin, ymin, xmax, ymax)."""
    north, east = tile
    xmin, ymin, xmax, ymax = window
    return (max(xmin, east), max(ymin, north),
            min(xmax, east + TILE_SIZE_M), min(ymax, north + TILE_SIZE_M))


def tile_url_candidates(edition, tile):
    """Return the possible URLs of one tile (file naming changed between editions)."""
    north, east = tile
    return [f"{BEV_DSM_BASE_URL}/{edition}/{pattern.format(north=north, east=east)}"
            for pattern in TILE_FILE_PATTERNS]


def find_tile_url(edition, tile):
    """Return the /vsicurl/ path of the tile that exists on the server.

    Raises
    ------
    BevDownloadError
        If none of the candidate URLs exists.
    """
    candidates = tile_url_candidates(edition, tile)
    with gdal.config_options(REMOTE_READ_OPTIONS):
        for url in candidates:
            if gdal.VSIStatL(f"/vsicurl/{url}") is not None:
                return f"/vsicurl/{url}"
    raise BevDownloadError(
        f"BEV DSM tile N{tile[0]}E{tile[1]} of edition {edition} not found. Tried:\n  "
        + "\n  ".join(candidates)
        + "\nCheck fine_dem.edition (e.g. 20250915) and whether the BEV server is reachable."
    )


def local_tile_path(directory, edition, tile, overview_level):
    """Return the file name of the downloaded AOI part of one tile."""
    resolution = "1m" if overview_level is None else f"ovr{overview_level}"
    return directory / f"ALS_DSM_{edition}_N{tile[0]}E{tile[1]}_{resolution}_aoi.tif"


def download_tile_window(remote_path, window, out_path, overview_level):
    """Copy one window of a remote BEV tile to a local compressed GeoTIFF.

    Writes ``<out_path>.part`` first and renames it at the end, so an
    interrupted download never leaves a file that looks complete.
    """
    xmin, ymin, xmax, ymax = window
    open_options = [] if overview_level is None else [f"OVERVIEW_LEVEL={overview_level}"]
    part_path = out_path.with_name(out_path.name + ".part")
    logger.info("Downloading %s window x %d..%d, y %d..%d -> %s",
                remote_path.rsplit("/", 1)[-1], xmin, xmax, ymin, ymax, out_path.name)
    with gdal.config_options(REMOTE_READ_OPTIONS):
        source = gdal.OpenEx(remote_path, gdal.OF_RASTER, open_options=open_options)
        gdal.Translate(
            str(part_path), source,
            format="GTiff",
            projWin=[xmin, ymax, xmax, ymin],
            creationOptions=LOCAL_GEOTIFF_OPTIONS,
            callback=gdal.TermProgress_nocb,
        )
    part_path.rename(out_path)


def download_bev_dsm(aoi_bbox, buffer_m, edition, overview_level, work_dir):
    """Download the BEV DSM covering the AOI plus buffer and return the mosaic VRT path.

    Parameters
    ----------
    aoi_bbox : tuple
        (west, south, east, north) in degrees.
    buffer_m : float
        Buffer around the AOI in meters.
    edition : str
        BEV edition folder, e.g. "20250915".
    overview_level : int or None
        None = full 1 m data; 0 = 2 m overview, 1 = 4 m, ...
    work_dir : Path
        Files are written to ``work_dir / "bev"``.

    Returns
    -------
    Path
        VRT mosaicking all downloaded tile windows (EPSG:3035).
    """
    directory = work_dir / "bev"
    directory.mkdir(parents=True, exist_ok=True)
    window = aoi_window_3035(aoi_bbox, buffer_m)
    tiles = tiles_for_window(window)
    logger.info("BEV DSM edition %s: AOI + %g m buffer needs %d tiles: %s", edition, buffer_m,
                len(tiles), ", ".join(f"N{n}E{e}" for n, e in tiles))

    local_files = []
    for tile in tiles:
        out_path = local_tile_path(directory, edition, tile, overview_level)
        if out_path.is_file():
            logger.info("Already downloaded: %s", out_path.name)
        else:
            download_tile_window(find_tile_url(edition, tile), clip_window_to_tile(window, tile),
                                 out_path, overview_level)
        local_files.append(out_path)

    mosaic_path = directory / MOSAIC_FILE_NAME
    gdal.BuildVRT(str(mosaic_path), [str(path) for path in local_files])
    logger.info("BEV DSM mosaic: %s", mosaic_path)
    return mosaic_path
