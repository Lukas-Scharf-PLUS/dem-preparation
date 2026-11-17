"""Tests for the LAZ point cloud export on a small synthetic 1 m DSM near Vienna."""

import numpy as np
import pytest
from osgeo import gdal, osr
from pyproj import Transformer

import point_cloud_export
from point_cloud_export import PointCloudError, export_point_cloud, utm_epsg
from prepare_dem import parse_arguments
from vertical_datum import FrameShift

laspy = pytest.importorskip("laspy")

# Lower-left corner of the synthetic DSM in EPSG:3035 (central Vienna), 1 m pixels.
ORIGIN_X, ORIGIN_Y = 4_796_000.0, 2_808_000.0
WIDTH, HEIGHT = 40, 30
BEV_NODATA = -9999.0
CONSTANT_OFFSET_M = 44.0
FRAME_SHIFT = FrameShift(delta_lon_deg=0.0, delta_lat_deg=0.0, delta_height_m=0.5)


@pytest.fixture
def synthetic_dsm(tmp_path):
    """Write a 40 x 30 m DSM (EPSG:3035) whose height = 200 + column, with one no-data pixel."""
    heights = 200.0 + np.tile(np.arange(WIDTH, dtype=np.float32), (HEIGHT, 1))
    heights[15, 20] = BEV_NODATA
    path = tmp_path / "dsm.tif"
    dataset = gdal.GetDriverByName("GTiff").Create(str(path), WIDTH, HEIGHT, 1, gdal.GDT_Float32)
    dataset.SetGeoTransform((ORIGIN_X, 1.0, 0.0, ORIGIN_Y + HEIGHT, 0.0, -1.0))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3035)
    dataset.SetProjection(srs.ExportToWkt())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(BEV_NODATA)
    band.WriteArray(heights)
    dataset = None
    return path


@pytest.fixture(autouse=True)
def constant_height_offset(monkeypatch):
    """Replace the PROJ geoid grid (network download) with a constant offset."""
    monkeypatch.setattr(point_cloud_export, "height_offset_function",
                        lambda reference: lambda lon, lat: np.full(np.shape(lon), CONSTANT_OFFSET_M))


def lonlat_box(x_min, y_min, x_max, y_max):
    """Return the (west, south, east, north) box of an EPSG:3035 rectangle's corners."""
    to_lonlat = Transformer.from_crs(3035, 4258, always_xy=True)
    lon, lat = to_lonlat.transform([x_min, x_max, x_min, x_max], [y_min, y_min, y_max, y_max])
    return (min(lon), min(lat), max(lon), max(lat))


def test_vienna_is_in_utm_zone_33_north():
    assert utm_epsg(16.37, 48.21) == 32633


def test_points_have_ellipsoid_heights_and_utm_positions(synthetic_dsm, tmp_path):
    bbox = lonlat_box(ORIGIN_X + 10, ORIGIN_Y + 8, ORIGIN_X + 30, ORIGIN_Y + 22)
    out_path = tmp_path / "cloud.laz"
    written = export_point_cloud(synthetic_dsm, bbox, "evrf2000_austria", FRAME_SHIFT, out_path)

    cloud = laspy.read(str(out_path))
    assert len(cloud.points) == written > 0
    assert cloud.header.parse_crs().to_epsg() == 32633
    assert not out_path.with_name(out_path.name + ".part").exists()

    # Back to EPSG:3035 to find each point's DSM pixel; height = 200 + column + offsets.
    to_3035 = Transformer.from_crs(32633, 3035, always_xy=True)
    x_3035, y_3035 = to_3035.transform(np.asarray(cloud.x), np.asarray(cloud.y))
    columns = np.floor(x_3035 - ORIGIN_X)
    np.testing.assert_allclose(x_3035 - ORIGIN_X - columns, 0.5, atol=0.01)  # pixel centers
    expected_z = 200.0 + columns + CONSTANT_OFFSET_M + FRAME_SHIFT.delta_height_m
    np.testing.assert_allclose(np.asarray(cloud.z), expected_z, atol=0.002)


def test_no_data_pixels_are_left_out(synthetic_dsm, tmp_path):
    full_box = lonlat_box(ORIGIN_X + 5, ORIGIN_Y + 5, ORIGIN_X + 35, ORIGIN_Y + 25)
    cloud_path = tmp_path / "cloud.laz"
    export_point_cloud(synthetic_dsm, full_box, "evrf2000_austria", FRAME_SHIFT, cloud_path)
    z = np.asarray(laspy.read(str(cloud_path)).z)
    assert np.all(z > 200.0)  # the -9999 no-data pixel is not written


def test_box_outside_the_download_is_refused(synthetic_dsm, tmp_path):
    bbox = lonlat_box(ORIGIN_X - 100, ORIGIN_Y, ORIGIN_X + 10, ORIGIN_Y + 10)
    with pytest.raises(PointCloudError, match="not completely inside"):
        export_point_cloud(synthetic_dsm, bbox, "evrf2000_austria", FRAME_SHIFT, tmp_path / "x.laz")


def test_too_large_box_is_refused(synthetic_dsm, tmp_path, monkeypatch):
    monkeypatch.setattr(point_cloud_export, "MAX_POINTS", 10)
    bbox = lonlat_box(ORIGIN_X + 10, ORIGIN_Y + 8, ORIGIN_X + 30, ORIGIN_Y + 22)
    with pytest.raises(PointCloudError, match="smaller box"):
        export_point_cloud(synthetic_dsm, bbox, "evrf2000_austria", FRAME_SHIFT, tmp_path / "x.laz")


def test_pointcloud_step_requires_a_valid_bbox():
    with pytest.raises(SystemExit):
        parse_arguments(["pointcloud", "--config", "c.yaml"])
    with pytest.raises(SystemExit):
        parse_arguments(["pointcloud", "--config", "c.yaml", "--bbox", "16.4", "48.2", "16.3", "48.3"])
    arguments = parse_arguments(["pointcloud", "--config", "c.yaml", "--bbox", "16.3", "48.2", "16.4", "48.3"])
    assert arguments.bbox == [16.3, 48.2, 16.4, 48.3]
