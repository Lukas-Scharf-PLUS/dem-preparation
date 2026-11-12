"""
The lat/lon grids the DEM is written on.

ISCE2 and COMPASS both expect one regular grid in geographic WGS84
coordinates. A pixel size given in meters is converted to degrees at the
latitude of the AOI, separately for latitude and longitude, so pixels are
about square on the ground (at 48.2 deg N, 10 m = 0.32 arcsec in latitude
but 0.48 arcsec in longitude).

All grids are anchored at integer multiples of their pixel size (like
``gdalwarp -tap``). A sub-grid with the same pixel size therefore lines up
exactly with the main grid, and its pixels can be copied into it by an
integer row/column offset.
"""

import math
from dataclasses import dataclass

# WGS84 ellipsoid
WGS84_SEMI_MAJOR_AXIS_M = 6378137.0
WGS84_ECCENTRICITY_SQUARED = 0.00669437999014


def meters_per_degree(latitude_deg):
    """Return (meters per degree latitude, meters per degree longitude) on the WGS84 ellipsoid.

    Uses the meridian and prime-vertical radii of curvature at ``latitude_deg``.
    """
    phi = math.radians(latitude_deg)
    sin_squared = math.sin(phi) ** 2
    denominator = 1.0 - WGS84_ECCENTRICITY_SQUARED * sin_squared
    meridian_radius = WGS84_SEMI_MAJOR_AXIS_M * (1.0 - WGS84_ECCENTRICITY_SQUARED) / denominator ** 1.5
    prime_vertical_radius = WGS84_SEMI_MAJOR_AXIS_M / math.sqrt(denominator)
    to_radians = math.pi / 180.0
    return meridian_radius * to_radians, prime_vertical_radius * math.cos(phi) * to_radians


def pixel_size_degrees(resolution_m, latitude_deg):
    """Return (delta_lon, delta_lat) in degrees for a pixel of ``resolution_m`` meters."""
    meters_per_deg_lat, meters_per_deg_lon = meters_per_degree(latitude_deg)
    return resolution_m / meters_per_deg_lon, resolution_m / meters_per_deg_lat


@dataclass(frozen=True)
class GeoGrid:
    """A north-up lat/lon grid.

    Attributes
    ----------
    west, north : float
        Outer corner of the upper-left pixel (degrees).
    delta_lon, delta_lat : float
        Pixel size in degrees, both positive.
    width, height : int
        Number of columns and rows.
    """

    west: float
    north: float
    delta_lon: float
    delta_lat: float
    width: int
    height: int

    @property
    def east(self):
        """Outer east edge of the grid (degrees)."""
        return self.west + self.width * self.delta_lon

    @property
    def south(self):
        """Outer south edge of the grid (degrees)."""
        return self.north - self.height * self.delta_lat

    @property
    def bounds(self):
        """(west, south, east, north) of the outer pixel edges."""
        return (self.west, self.south, self.east, self.north)

    @property
    def geotransform(self):
        """GDAL geotransform (corner convention)."""
        return (self.west, self.delta_lon, 0.0, self.north, 0.0, -self.delta_lat)

    def pixel_center_coordinates(self):
        """Return (latitudes of the rows, longitudes of the columns) at pixel centers, as lists."""
        latitudes = [self.north - (row + 0.5) * self.delta_lat for row in range(self.height)]
        longitudes = [self.west + (col + 0.5) * self.delta_lon for col in range(self.width)]
        return latitudes, longitudes

    def subgrid(self, bounds):
        """Return the part of this grid that covers ``bounds``, and its offset.

        Parameters
        ----------
        bounds : tuple
            (west, south, east, north) in degrees; clipped to this grid.

        Returns
        -------
        (col_offset, row_offset, GeoGrid)
            The sub-grid uses the same pixels as this grid; its upper-left pixel
            is pixel (row_offset, col_offset) of this grid.

        Raises
        ------
        ValueError
            If ``bounds`` does not overlap this grid.
        """
        west, south, east, north = bounds
        first_col = max(0, math.floor((west - self.west) / self.delta_lon))
        last_col = min(self.width, math.ceil((east - self.west) / self.delta_lon))
        first_row = max(0, math.floor((self.north - north) / self.delta_lat))
        last_row = min(self.height, math.ceil((self.north - south) / self.delta_lat))
        if first_col >= last_col or first_row >= last_row:
            raise ValueError(f"Bounds {bounds} do not overlap the grid {self.bounds}.")
        sub = GeoGrid(
            west=self.west + first_col * self.delta_lon,
            north=self.north - first_row * self.delta_lat,
            delta_lon=self.delta_lon,
            delta_lat=self.delta_lat,
            width=last_col - first_col,
            height=last_row - first_row,
        )
        return first_col, first_row, sub


def aligned_grid(bounds, delta_lon, delta_lat):
    """Return the smallest grid anchored at multiples of the pixel size that covers ``bounds``.

    Parameters
    ----------
    bounds : tuple
        (west, south, east, north) in degrees.
    delta_lon, delta_lat : float
        Pixel size in degrees (positive).
    """
    west, south, east, north = bounds
    first_col = math.floor(west / delta_lon)
    last_col = math.ceil(east / delta_lon)
    first_row = math.floor(-north / delta_lat)  # rows count southwards from the equator * -1
    last_row = math.ceil(-south / delta_lat)
    return GeoGrid(
        west=first_col * delta_lon,
        north=-first_row * delta_lat,
        delta_lon=delta_lon,
        delta_lat=delta_lat,
        width=last_col - first_col,
        height=last_row - first_row,
    )


def output_grid(bounds, resolution_m, reference_latitude_deg):
    """Return the final DEM grid covering ``bounds`` with ``resolution_m`` pixels.

    The pixel size in degrees is computed at ``reference_latitude_deg``
    (the AOI center), so pixels are ~square meters inside the AOI.
    """
    delta_lon, delta_lat = pixel_size_degrees(resolution_m, reference_latitude_deg)
    return aligned_grid(bounds, delta_lon, delta_lat)
