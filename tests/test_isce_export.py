"""Tests for the ISCE2 DEM export (raw file, XML header, VRT)."""

import numpy as np
import pytest
from osgeo import gdal

from dem_merge import OUTPUT_NODATA, write_array
from isce_export import export_isce_dem, isce_header_values
from target_grid import aligned_grid


def test_header_uses_pixel_centers():
    values = isce_header_values((16.0, 0.1, 0.0, 49.0, 0.0, -0.1), width=4, length=3)
    assert values["first_lon"] == pytest.approx(16.05)
    assert values["first_lat"] == pytest.approx(48.95)
    assert values["end_lon"] == pytest.approx(16.45)


def test_export_writes_raw_xml_and_vrt(tmp_path):
    grid = aligned_grid((16.0, 48.0, 16.04, 48.03), 0.01, 0.01)
    heights = np.arange(grid.width * grid.height, dtype=np.float32).reshape(grid.height, grid.width) + 200
    heights[0, 0] = np.nan  # becomes OUTPUT_NODATA in the GeoTIFF, 0 in the ISCE DEM
    source = tmp_path / "dem.tif"
    write_array(source, heights, grid)

    raw = export_isce_dem(source, tmp_path / "out")
    assert raw.name == "out.dem.wgs84"
    data = np.fromfile(raw, dtype="<f4").reshape(grid.height, grid.width)
    assert data[0, 0] == 0.0
    np.testing.assert_allclose(data[1:], heights[1:])

    # The VRT must give the same pixels and the correct (corner) geotransform.
    vrt = gdal.Open(str(raw) + ".vrt")
    np.testing.assert_allclose(vrt.GetGeoTransform(), grid.geotransform)
    np.testing.assert_allclose(vrt.ReadAsArray(), data)
    assert OUTPUT_NODATA not in data


def test_isce2_reads_the_header(tmp_path):
    pytest.importorskip("isce")  # must be imported before isceobj
    isceobj = pytest.importorskip("isceobj")
    grid = aligned_grid((16.0, 48.0, 16.04, 48.03), 0.01, 0.01)
    source = tmp_path / "dem.tif"
    write_array(source, np.full((grid.height, grid.width), 250.0, np.float32), grid)
    raw = export_isce_dem(source, tmp_path / "out")

    image = isceobj.createDemImage()
    image.load(str(raw) + ".xml")
    assert image.getWidth() == grid.width
    assert image.getLength() == grid.height
    assert image.getFirstLongitude() == pytest.approx(16.005)
    assert image.getFirstLatitude() == pytest.approx(48.025)
    assert image.getDeltaLatitude() == pytest.approx(-0.01)
