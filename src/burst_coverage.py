"""
Find the area the DEM must cover: all Sentinel-1 bursts of the configured
relative orbits that touch the AOI.

Burst footprints come from the OPERA burst database
(``opera-burst-bbox-only.sqlite3``, the file COMPASS uses). It stores, for
every burst, a bounding box in the burst's UTM zone that already includes a
5 km margin. COMPASS geocodes each burst onto exactly this box, so a DEM that
covers it is enough for COMPASS.

ISCE2 (stackSentinel.py) cuts its DEM differently: it uses lat/lon rectangles
computed for heights between -500 m and 9000 m and adds 0.1 deg. These reach
~0.3 deg beyond the real footprint in range direction, so a configurable
margin (``bursts.margin_deg``) is added to the union of the boxes.
"""

import logging
import sqlite3
from dataclasses import dataclass

import numpy as np
import shapely
from osgeo import gdal
from pyproj import Transformer
from shapely.geometry import box

logger = logging.getLogger(__name__)
gdal.UseExceptions()

WGS84_EPSG = 4326
# Densify box edges every 1 km before reprojecting, so curved edges stay accurate.
DENSIFY_STEP_M = 1000.0


class CoverageError(RuntimeError):
    """Raised when the burst database is unusable or a DEM does not cover the bursts."""


@dataclass
class Burst:
    """One Sentinel-1 burst: its JPL burst ID and lat/lon bounding box (W, S, E, N)."""

    burst_id: str
    bounds: tuple


def utm_box_to_lonlat_polygon(epsg, xmin, ymin, xmax, ymax):
    """Return the UTM box (densified) as a lon/lat shapely polygon."""
    transformer = Transformer.from_crs(epsg, WGS84_EPSG, always_xy=True)
    polygon = box(xmin, ymin, xmax, ymax).segmentize(DENSIFY_STEP_M)
    return shapely.transform(polygon, lambda xy: np.column_stack(transformer.transform(xy[:, 0], xy[:, 1])))


def read_bursts_of_orbit(connection, relative_orbit):
    """Return all rows (burst_id, epsg, xmin, ymin, xmax, ymax) of one relative orbit."""
    # JPL burst IDs look like "t073_154933_iw2"; the track number is the relative orbit.
    pattern = f"t{relative_orbit:03d}\\_%"
    return connection.execute(
        "SELECT burst_id_jpl, epsg, xmin, ymin, xmax, ymax FROM burst_id_map "
        "WHERE burst_id_jpl LIKE ? ESCAPE '\\'",
        (pattern,),
    ).fetchall()


def find_bursts(burst_db, relative_orbits, aoi_bbox):
    """Return the bursts of the given relative orbits whose box intersects the AOI.

    Parameters
    ----------
    burst_db : Path
        OPERA ``opera-burst-bbox-only.sqlite3`` file.
    relative_orbits : list[int]
    aoi_bbox : tuple
        (west, south, east, north) in degrees.

    Returns
    -------
    dict[int, list[Burst]]
        Bursts per relative orbit.

    Raises
    ------
    CoverageError
        If the database is missing or an orbit has no burst touching the AOI.
    """
    if not burst_db.is_file():
        raise CoverageError(
            f"Burst database not found: {burst_db}\n"
            "Download 'opera-burst-bbox-only.sqlite3' from "
            "https://github.com/opera-adt/burst_db/releases and set bursts.burst_db."
        )
    aoi = box(*aoi_bbox)
    bursts_per_orbit = {}
    with sqlite3.connect(f"file:{burst_db}?mode=ro", uri=True) as connection:
        for orbit in relative_orbits:
            selected = []
            for burst_id, epsg, xmin, ymin, xmax, ymax in read_bursts_of_orbit(connection, orbit):
                footprint = utm_box_to_lonlat_polygon(epsg, xmin, ymin, xmax, ymax)
                if footprint.intersects(aoi):
                    selected.append(Burst(burst_id, footprint.bounds))
            if not selected:
                raise CoverageError(f"No burst of relative orbit {orbit} touches the AOI {aoi_bbox}.")
            bursts_per_orbit[orbit] = sorted(selected, key=lambda burst: burst.burst_id)
    return bursts_per_orbit


def union_bounds(bounds_list):
    """Return the (west, south, east, north) box enclosing all given boxes."""
    return (
        min(b[0] for b in bounds_list),
        min(b[1] for b in bounds_list),
        max(b[2] for b in bounds_list),
        max(b[3] for b in bounds_list),
    )


def expand_bounds(bounds, margin_deg):
    """Return ``bounds`` grown by ``margin_deg`` on every side."""
    west, south, east, north = bounds
    return (west - margin_deg, south - margin_deg, east + margin_deg, north + margin_deg)


def required_dem_bounds(bursts_per_orbit, margin_deg):
    """Return the area the DEM must cover: union of all burst boxes plus the margin."""
    all_bounds = [burst.bounds for bursts in bursts_per_orbit.values() for burst in bursts]
    return expand_bounds(union_bounds(all_bounds), margin_deg)


def raster_bounds(raster_file):
    """Return (west, south, east, north) of a north-up lat/lon raster.

    Raises
    ------
    CoverageError
        If the raster is not in geographic coordinates or is rotated.
    """
    dataset = gdal.Open(str(raster_file))
    spatial_ref = dataset.GetSpatialRef()
    if spatial_ref is None or not spatial_ref.IsGeographic():
        raise CoverageError(f"DEM must be in geographic lat/lon coordinates: {raster_file}")
    west, delta_x, rotation_x, north, rotation_y, delta_y = dataset.GetGeoTransform()
    if rotation_x or rotation_y or delta_y >= 0:
        raise CoverageError(f"DEM must be a north-up raster without rotation: {raster_file}")
    east = west + dataset.RasterXSize * delta_x
    south = north + dataset.RasterYSize * delta_y
    return (west, south, east, north)


def bounds_contain(outer, inner):
    """Return True if box ``outer`` fully contains box ``inner`` (both W, S, E, N)."""
    return outer[0] <= inner[0] and outer[1] <= inner[1] and outer[2] >= inner[2] and outer[3] >= inner[3]


def format_bounds(bounds):
    """Return a box as a readable 'W .. E, S .. N' string."""
    west, south, east, north = bounds
    return f"lon {west:.3f} .. {east:.3f}, lat {south:.3f} .. {north:.3f}"


def check_dem_covers(raster_file, required_bounds):
    """Return (covers, dem_bounds) for a lat/lon DEM file and the required box."""
    dem_bounds = raster_bounds(raster_file)
    return bounds_contain(dem_bounds, required_bounds), dem_bounds


def log_burst_summary(bursts_per_orbit, required_bounds):
    """Log the selected bursts per orbit and the resulting DEM area."""
    for orbit, bursts in bursts_per_orbit.items():
        orbit_bounds = union_bounds([burst.bounds for burst in bursts])
        logger.info("Relative orbit %3d: %2d bursts (%s), %s", orbit, len(bursts),
                    ", ".join(burst.burst_id for burst in bursts), format_bounds(orbit_bounds))
    logger.info("DEM must cover (bursts + margin): %s", format_bounds(required_bounds))
