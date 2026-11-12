"""Tests for burst selection (with a small fake burst database) and the coverage check."""

import sqlite3

import numpy as np
import pytest
from osgeo import gdal

from burst_coverage import (CoverageError, bounds_contain, check_dem_covers, find_bursts,
                            required_dem_bounds)

AOI = (16.18, 48.12, 16.59, 48.32)
UTM_33N = 32633


def make_burst_db(path, rows):
    """Write a minimal burst_id_map table like the OPERA bbox-only database."""
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE burst_id_map (burst_id_jpl TEXT, epsg INTEGER, "
                           "xmin INTEGER, ymin INTEGER, xmax INTEGER, ymax INTEGER)")
        connection.executemany("INSERT INTO burst_id_map VALUES (?, ?, ?, ?, ?, ?)", rows)
    return path


@pytest.fixture
def burst_db(tmp_path):
    # Vienna is at roughly x 600 km, y 5340 km in UTM 33N.
    return make_burst_db(tmp_path / "bursts.sqlite3", [
        ("t073_000001_iw1", UTM_33N, 560_000, 5_320_000, 650_000, 5_345_000),  # touches the AOI
        ("t073_000002_iw1", UTM_33N, 560_000, 5_400_000, 650_000, 5_425_000),  # far north
        ("t022_000001_iw2", UTM_33N, 590_000, 5_330_000, 680_000, 5_355_000),  # touches the AOI
        ("t173_000001_iw1", UTM_33N, 560_000, 5_320_000, 650_000, 5_345_000),  # other orbit (173), must not match 73
    ])


def test_only_bursts_touching_the_aoi_are_selected(burst_db):
    bursts = find_bursts(burst_db, [73, 22], AOI)
    assert [burst.burst_id for burst in bursts[73]] == ["t073_000001_iw1"]
    assert [burst.burst_id for burst in bursts[22]] == ["t022_000001_iw2"]


def test_orbit_without_bursts_raises(burst_db):
    with pytest.raises(CoverageError, match="146"):
        find_bursts(burst_db, [146], AOI)


def test_missing_database_raises(tmp_path):
    with pytest.raises(CoverageError, match="not found"):
        find_bursts(tmp_path / "missing.sqlite3", [73], AOI)


def test_required_bounds_is_union_plus_margin(burst_db):
    bursts = find_bursts(burst_db, [73, 22], AOI)
    west, south, east, north = required_dem_bounds(bursts, 0.3)
    union_west = min(b.bounds[0] for orbit in bursts.values() for b in orbit)
    union_north = max(b.bounds[3] for orbit in bursts.values() for b in orbit)
    assert west == pytest.approx(union_west - 0.3)
    assert north == pytest.approx(union_north + 0.3)
    assert bounds_contain((west, south, east, north), AOI)


def write_latlon_raster(path, west, north, width, height, pixel=0.01):
    """Write a small zero-valued EPSG:4326 GeoTIFF."""
    dataset = gdal.GetDriverByName("GTiff").Create(str(path), width, height, 1, gdal.GDT_Float32)
    dataset.SetGeoTransform((west, pixel, 0, north, 0, -pixel))
    dataset.SetProjection("EPSG:4326")
    dataset.GetRasterBand(1).WriteArray(np.zeros((height, width), np.float32))
    dataset = None
    return path


def test_dem_coverage_check(tmp_path):
    dem = write_latlon_raster(tmp_path / "dem.tif", 14.0, 50.0, 500, 300)  # 14-19 E, 47-50 N
    covers, dem_bounds = check_dem_covers(dem, (15.1, 47.6, 18.0, 48.8))
    assert covers
    assert dem_bounds == pytest.approx((14.0, 47.0, 19.0, 50.0))
    covers, _ = check_dem_covers(dem, (15.1, 46.9, 18.0, 48.8))
    assert not covers
