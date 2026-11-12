"""
Convert heights to WGS84 ellipsoid heights, and positions from ETRS89 to ITRF.

Why this is needed
------------------
ISCE2 and COMPASS treat DEM values as heights above the WGS84 *ellipsoid*.
The input DEMs store other heights:

* Copernicus DEM: heights above the EGM2008 geoid.
* BEV ALS DSM: Austrian heights, EVRF2000 Austria (EPSG:9274); the edition
  20210915 uses Gebrauchshöhen Adria (GHA, EPSG:5778).

ellipsoid height = height + offset(lat, lon). In Vienna the offset is about
44.6 m (EVRF2000 Austria), 44.2 m (GHA) and 44.7 m (EGM2008). Using the wrong
reference shifts the radar geolocation by roughly offset / tan(incidence).

The offsets come from the official PROJ grids (BEV: ``at_bev_*.tif``, NGA:
``us_nga_egm08_25.tif``). PROJ downloads them once from cdn.proj.org. PROJ
silently uses a zero shift when a grid is missing, so transformations are
created with ``only_best=True``, which raises an error instead. Offsets are
smooth, so they are computed on a coarse grid and interpolated bilinearly.

Horizontal frame: BEV coordinates are ETRS89, which is fixed to the Eurasian
plate and moves ~2.5 cm/year relative to ITRF. Sentinel-1 orbits (and WGS84
today) are realized in ITRF. In 2024 the difference in Vienna is about
0.7 m east and 0.6 m north. It is applied as one constant shift, because it
varies by less than a millimeter across the area.
"""

import logging
from dataclasses import dataclass

import numpy as np
from pyproj import Transformer, network
from pyproj.exceptions import ProjError

logger = logging.getLogger(__name__)

# Source CRS (2D position + vertical height) and target 3D CRS (ellipsoid heights)
# for every height reference this tool converts.
HEIGHT_REFERENCES = {
    "egm2008": ("EPSG:4326+3855", "EPSG:4979"),
    "egm96": ("EPSG:4326+5773", "EPSG:4979"),
    "evrf2000_austria": ("EPSG:4258+9274", "EPSG:4937"),
    "gha": ("EPSG:4258+5778", "EPSG:4937"),
}
ELLIPSOID = "ellipsoid"

# BEV: "Grundsätzlich: EPSG:9274. Ausnahme: alle Kacheln der Serie mit Stichtag
# 15.09.2021: EPSG:5778" (BEV ALS-DOM product page).
BEV_GHA_EDITIONS = ("20210915",)

# ETRF2000 (the ETRS89 realization used in Austria) and the ITRF 3D CRS.
ETRF2000_3D_EPSG = "EPSG:7931"
ITRF_3D_EPSG = {"itrf2014": "EPSG:7912", "itrf2020": "EPSG:9989"}
ETRS89 = "etrs89"

# Height offsets on Earth lie between about -107 m and +86 m.
OFFSET_RANGE_M = (-110.0, 90.0)


class VerticalDatumError(RuntimeError):
    """Raised when PROJ cannot convert heights or positions exactly."""


@dataclass(frozen=True)
class FrameShift:
    """Constant shift from ETRS89 to the target frame, at the AOI center."""

    delta_lon_deg: float
    delta_lat_deg: float
    delta_height_m: float


def enable_proj_network():
    """Let PROJ download missing transformation grids from cdn.proj.org (cached locally)."""
    network.set_network_enabled(True)


def bev_height_reference(edition, configured):
    """Return the height reference of a BEV edition ('evrf2000_austria' or 'gha').

    ``configured`` is the config value; 'auto' picks the reference documented by the BEV.
    """
    if configured != "auto":
        return configured
    return "gha" if edition in BEV_GHA_EDITIONS else "evrf2000_austria"


def height_offset_function(height_reference):
    """Return a function offset(lon, lat) -> ellipsoid height minus height, in meters.

    Raises
    ------
    VerticalDatumError
        If PROJ cannot use the exact grid (not available and not downloadable).
    """
    source, target = HEIGHT_REFERENCES[height_reference]
    enable_proj_network()
    try:
        transformer = Transformer.from_crs(source, target, always_xy=True, only_best=True)
    except ProjError as error:
        raise VerticalDatumError(
            f"PROJ cannot convert {height_reference} heights ({source} -> {target}) exactly: "
            f"{error}\nIs cdn.proj.org reachable? Grids can also be fetched once with 'projsync'."
        ) from error

    def offset(longitudes, latitudes):
        longitudes = np.asarray(longitudes, dtype=np.float64)
        _, _, ellipsoid_height = transformer.transform(
            longitudes, np.asarray(latitudes, dtype=np.float64), np.zeros_like(longitudes)
        )
        return np.asarray(ellipsoid_height)

    return offset


def offset_grid(latitudes, longitudes, offset_function, step_deg):
    """Return the height offset for every pixel of a lat/lon grid (rows x columns).

    The offset is evaluated on a coarse grid with ``step_deg`` spacing and
    interpolated bilinearly to the pixel centers.

    Parameters
    ----------
    latitudes, longitudes : 1D arrays
        Pixel-center coordinates of the rows and columns.
    offset_function : callable
        offset(lon, lat) -> meters, e.g. from :func:`height_offset_function`.
    step_deg : float
        Spacing of the coarse evaluation grid.

    Raises
    ------
    VerticalDatumError
        If an offset is not finite or outside the physical range.
    """
    from scipy.interpolate import RegularGridInterpolator

    latitudes = np.asarray(latitudes, dtype=np.float64)
    longitudes = np.asarray(longitudes, dtype=np.float64)
    coarse_lat = np.arange(latitudes.min() - step_deg, latitudes.max() + 2 * step_deg, step_deg)
    coarse_lon = np.arange(longitudes.min() - step_deg, longitudes.max() + 2 * step_deg, step_deg)
    lon_mesh, lat_mesh = np.meshgrid(coarse_lon, coarse_lat)
    coarse_offset = offset_function(lon_mesh.ravel(), lat_mesh.ravel()).reshape(lat_mesh.shape)

    low, high = OFFSET_RANGE_M
    if not np.all(np.isfinite(coarse_offset)) or coarse_offset.min() < low or coarse_offset.max() > high:
        raise VerticalDatumError(
            f"Height offset is not finite or outside {low} .. {high} m "
            f"(min {np.nanmin(coarse_offset)}, max {np.nanmax(coarse_offset)}). "
            "Is the area outside the grid of this height reference?"
        )

    interpolator = RegularGridInterpolator((coarse_lat, coarse_lon), coarse_offset)
    lat_pixels, lon_pixels = np.meshgrid(latitudes, longitudes, indexing="ij")
    return interpolator((lat_pixels, lon_pixels)).astype(np.float32)


def etrs89_to_itrf_shift(frame, epoch, longitude, latitude):
    """Return the FrameShift from ETRS89 (ETRF2000) to ``frame`` at ``epoch`` and one point.

    Parameters
    ----------
    frame : str
        'etrs89' (no shift), 'itrf2014' or 'itrf2020'.
    epoch : float
        Decimal year, e.g. 2024.0.
    longitude, latitude : float
        Where the shift is evaluated (AOI center).

    Raises
    ------
    VerticalDatumError
        If PROJ has no exact transformation.
    """
    if frame == ETRS89:
        return FrameShift(0.0, 0.0, 0.0)
    try:
        transformer = Transformer.from_crs(ETRF2000_3D_EPSG, ITRF_3D_EPSG[frame],
                                           always_xy=True, only_best=True)
    except ProjError as error:
        raise VerticalDatumError(f"PROJ cannot transform ETRF2000 -> {frame}: {error}") from error
    reference_height_m = 200.0  # any height; the shift does not depend on it noticeably
    lon, lat, height, _ = transformer.transform(longitude, latitude, reference_height_m, epoch)
    return FrameShift(lon - longitude, lat - latitude, height - reference_height_m)
