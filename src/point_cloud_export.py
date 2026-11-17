"""
Export part of the downloaded 1 m BEV ALS DSM as a LAZ point cloud, e.g. for CloudCompare.

Every valid 1 m DSM pixel inside a lon/lat box becomes one point at its pixel
center. The points use the same height reference and frame as the DEM built
for ISCE2, so they can be compared directly with InSAR points geocoded by it:

* heights are converted to WGS84 *ellipsoid* heights (like the DEM build),
* positions are shifted from ETRS89 to the configured ITRF frame and epoch
  (``merge.horizontal_frame`` / ``merge.frame_epoch``),
* x/y are written in the WGS84 UTM zone of the box (meters), so all three
  axes have the same unit in a 3D viewer. The CRS is stored in the file.

The DSM itself is not resampled: one output point per 1 m input pixel. A
1 km x 1 km box gives about one million points; the export refuses boxes
above ``MAX_POINTS`` because CloudCompare would need too much memory.

Why LAZ: it is the compressed lidar standard (about 5-10x smaller than an
uncompressed format) and opens in CloudCompare, QGIS and PDAL.
"""

import logging
import math

import numpy as np
from osgeo import gdal
from pyproj import CRS, Transformer

from bev_download import BEV_EPSG, aoi_window_3035
from vertical_datum import height_offset_function

logger = logging.getLogger(__name__)
gdal.UseExceptions()

ETRS89_EPSG = 4258
WGS84_EPSG = 4326
# EPSG codes of the WGS84 UTM zones: 326<zone> north, 327<zone> south.
UTM_NORTH_EPSG_BASE = 32600
UTM_SOUTH_EPSG_BASE = 32700
UTM_ZONE_WIDTH_DEG = 6
# LAS stores coordinates as integers times a scale; 1 mm is far below the DSM accuracy.
COORDINATE_SCALE_M = 0.001
# LAS 1.4 point format 6 is the first one that stores the CRS as WKT.
LAS_POINT_FORMAT = 6
LAS_VERSION = "1.4"
# Upper limit of points in one file (~1.2 GB as LAZ, several GB in CloudCompare's memory).
MAX_POINTS = 100_000_000
# DSM rows converted and written at a time, to keep memory use small.
ROWS_PER_BLOCK = 500
# A downloaded overview (overview_level >= 0) has pixels of 2 m or more.
FULL_RESOLUTION_M = 1.0


class PointCloudError(RuntimeError):
    """Raised when the requested box cannot be exported."""


def utm_epsg(longitude, latitude):
    """Return the EPSG code of the WGS84 UTM zone containing (longitude, latitude)."""
    zone = int(math.floor((longitude + 180.0) / UTM_ZONE_WIDTH_DEG)) % 60 + 1
    base = UTM_NORTH_EPSG_BASE if latitude >= 0 else UTM_SOUTH_EPSG_BASE
    return base + zone


def default_output_name(aoi_name, edition, bbox):
    """Return a file name like ``vienna_bev_dsm_20250915_16.3700_48.2000_16.3900_48.2150.laz``."""
    corners = "_".join(f"{value:.4f}" for value in bbox)
    return f"{aoi_name}_bev_dsm_{edition}_{corners}.laz"


def read_dsm_window(bev_mosaic, bbox):
    """Read the DSM pixels covering a lon/lat box.

    Parameters
    ----------
    bev_mosaic : Path
        VRT of the downloaded BEV tiles (EPSG:3035).
    bbox : tuple
        (west, south, east, north) in degrees.

    Returns
    -------
    heights : np.ndarray
        2D float32 heights in the BEV height reference, NaN where no data.
    x_centers, y_centers : np.ndarray
        1D EPSG:3035 coordinates of the column and row pixel centers.

    Raises
    ------
    PointCloudError
        If the box is not completely inside the downloaded DSM, or would give
        more than ``MAX_POINTS`` points.
    """
    xmin, ymin, xmax, ymax = aoi_window_3035(bbox, buffer_m=0)
    source = gdal.Open(str(bev_mosaic))
    origin_x, pixel_width, _, origin_y, _, pixel_height = source.GetGeoTransform()
    pixel_height = -pixel_height  # north-up raster: the geotransform stores it negative
    if pixel_width > FULL_RESOLUTION_M:
        logger.warning("The downloaded DSM has %g m pixels, not 1 m (fine_dem.overview_level is set).",
                       pixel_width)

    first_col = math.floor((xmin - origin_x) / pixel_width)
    last_col = math.ceil((xmax - origin_x) / pixel_width)
    first_row = math.floor((origin_y - ymax) / pixel_height)
    last_row = math.ceil((origin_y - ymin) / pixel_height)
    inside = (first_col >= 0 and first_row >= 0
              and last_col <= source.RasterXSize and last_row <= source.RasterYSize)
    if not inside:
        raise PointCloudError(
            f"The box {bbox} is not completely inside the downloaded BEV DSM ({bev_mosaic}). "
            "Choose a box inside aoi.bbox (plus fine_dem.buffer_m), or download a larger area."
        )

    width, height = last_col - first_col, last_row - first_row
    if width * height > MAX_POINTS:
        raise PointCloudError(
            f"The box has {width} x {height} = {width * height / 1e6:.0f} million DSM pixels, more than "
            f"the limit of {MAX_POINTS / 1e6:.0f} million. Choose a smaller box "
            "(1 km x 1 km is about one million points)."
        )

    band = source.GetRasterBand(1)
    heights = band.ReadAsArray(first_col, first_row, width, height).astype(np.float32)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        heights[heights == nodata] = np.nan
    x_centers = origin_x + (first_col + np.arange(width) + 0.5) * pixel_width
    y_centers = origin_y - (first_row + np.arange(height) + 0.5) * pixel_height
    return heights, x_centers, y_centers


def to_utm_points(x_3035, y_3035, heights, height_offset, frame_shift, utm_crs_epsg):
    """Convert DSM pixel centers to UTM x/y and WGS84 ellipsoid heights.

    Parameters
    ----------
    x_3035, y_3035, heights : np.ndarray
        1D arrays of the same length: EPSG:3035 position and BEV height of each point.
    height_offset : callable
        offset(lon, lat) -> ellipsoid height minus BEV height, in meters
        (from :func:`vertical_datum.height_offset_function`).
    frame_shift : FrameShift
        ETRS89 -> output frame shift (from :func:`vertical_datum.etrs89_to_itrf_shift`).
    utm_crs_epsg : int
        EPSG code of the output UTM zone.

    Returns
    -------
    tuple of np.ndarray
        (x, y, z): UTM easting and northing in meters, ellipsoid height in meters.
    """
    to_lonlat = Transformer.from_crs(BEV_EPSG, ETRS89_EPSG, always_xy=True)
    longitudes, latitudes = to_lonlat.transform(x_3035, y_3035)
    ellipsoid_heights = (heights.astype(np.float64) + height_offset(longitudes, latitudes)
                         + frame_shift.delta_height_m)

    # The shifted lon/lat are in ITRF, which WGS84 (and its UTM zones) follows within centimeters.
    to_utm = Transformer.from_crs(WGS84_EPSG, utm_crs_epsg, always_xy=True)
    x, y = to_utm.transform(longitudes + frame_shift.delta_lon_deg, latitudes + frame_shift.delta_lat_deg)
    return np.asarray(x), np.asarray(y), ellipsoid_heights


def import_laspy():
    """Return the laspy module, or raise PointCloudError explaining how to install it.

    laspy is imported only for this export, so the DEM steps work without it.
    """
    try:
        import laspy
    except ImportError as error:
        raise PointCloudError(
            "The pointcloud step needs laspy with LAZ support. Install it with:\n"
            "  conda install -n dem-preparation -c conda-forge laspy lazrs-python"
        ) from error
    return laspy


def las_header(utm_crs_epsg, x_offset, y_offset):
    """Return a LAS 1.4 header with millimeter scales and the UTM CRS stored as WKT."""
    laspy = import_laspy()
    header = laspy.LasHeader(point_format=LAS_POINT_FORMAT, version=LAS_VERSION)
    header.scales = np.full(3, COORDINATE_SCALE_M)
    # Offsets near the data keep the stored integers small.
    header.offsets = np.array([math.floor(x_offset), math.floor(y_offset), 0.0])
    header.add_crs(CRS.from_epsg(utm_crs_epsg))
    return header


def export_point_cloud(bev_mosaic, bbox, bev_height_reference, frame_shift, out_path):
    """Write the BEV DSM pixels inside a lon/lat box as a LAZ point cloud.

    Parameters
    ----------
    bev_mosaic : Path
        VRT of the downloaded BEV tiles (EPSG:3035).
    bbox : tuple
        (west, south, east, north) in degrees.
    bev_height_reference : str
        'evrf2000_austria' or 'gha'.
    frame_shift : FrameShift
        ETRS89 -> output frame shift.
    out_path : Path
        LAZ file to write. Written as ``<out_path>.part`` first and renamed at
        the end, so an interrupted export never leaves a file that looks complete.

    Returns
    -------
    int
        Number of points written.

    Raises
    ------
    PointCloudError
        If the box is outside the downloaded DSM, too large, or has no valid pixels.
    """
    laspy = import_laspy()
    heights, x_centers, y_centers = read_dsm_window(bev_mosaic, bbox)
    valid_pixels = int(np.isfinite(heights).sum())
    if valid_pixels == 0:
        raise PointCloudError(f"The BEV DSM has no valid pixels inside the box {bbox}.")

    west, south, east, north = bbox
    utm_crs_epsg = utm_epsg((west + east) / 2, (south + north) / 2)
    height_offset = height_offset_function(bev_height_reference)
    corner_x, corner_y = Transformer.from_crs(WGS84_EPSG, utm_crs_epsg, always_xy=True).transform(west, south)
    header = las_header(utm_crs_epsg, corner_x, corner_y)
    logger.info("Exporting %d DSM points (%d x %d pixels) to EPSG:%d, ellipsoid heights: %s",
                valid_pixels, len(x_centers), len(y_centers), utm_crs_epsg, out_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = out_path.with_name(out_path.name + ".part")
    written = 0
    with laspy.open(str(part_path), mode="w", header=header, do_compress=True) as writer:
        for first_row in range(0, len(y_centers), ROWS_PER_BLOCK):
            block = heights[first_row:first_row + ROWS_PER_BLOCK]
            rows, cols = np.nonzero(np.isfinite(block))
            if rows.size == 0:
                continue
            x, y, z = to_utm_points(x_centers[cols], y_centers[first_row + rows], block[rows, cols],
                                    height_offset, frame_shift, utm_crs_epsg)
            points = laspy.ScaleAwarePointRecord.zeros(rows.size, header=header)
            points.x, points.y, points.z = x, y, z
            writer.write_points(points)
            written += rows.size
    part_path.rename(out_path)
    logger.info("Point cloud written: %d points, %.1f MB: %s", written,
                out_path.stat().st_size / 1e6, out_path)
    return written
