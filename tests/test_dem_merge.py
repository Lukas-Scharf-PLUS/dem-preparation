"""Tests for the merge helpers and the coarse/fine/VRT chain on small synthetic rasters."""

import numpy as np
import pytest
from osgeo import gdal

from dem_merge import (OUTPUT_NODATA, aoi_mask, blend, blend_weights, build_mixed_resolution_vrt,
                       difference_statistics, fill_small_holes, prepare_coarse_part, write_array,
                       write_final_geotiff)
from target_grid import aligned_grid


def test_blend_weights_rise_from_the_edge():
    valid = np.ones((41, 41), bool)
    weights = blend_weights(valid, blend_width_m=100.0, pixel_height_m=10.0, pixel_width_m=10.0)
    # distance to the array border: 1 pixel (10 m) at the edge, 21 pixels (210 m) in the middle
    assert weights[20, 0] == pytest.approx(0.1)
    assert weights[20, 20] == pytest.approx(1.0)
    assert np.all(np.diff(weights[20, :10]) > 0)


def test_blend_uses_coarse_where_fine_is_missing():
    fine = np.array([[10.0, np.nan]], np.float32)
    coarse = np.array([[0.0, 5.0]], np.float32)
    merged = blend(fine, coarse, np.array([[0.5, 1.0]], np.float32))
    np.testing.assert_allclose(merged, [[5.0, 5.0]])


def test_only_interior_holes_are_filled():
    heights = np.full((20, 20), 100.0, np.float32)
    heights[9:11, 9:11] = np.nan   # interior hole
    heights[:, :3] = np.nan        # no-data at the edge (outside the BEV area)
    filled = fill_small_holes(heights, max_distance_px=5)
    assert np.all(np.isfinite(filled[9:11, 9:11]))
    np.testing.assert_allclose(filled[9:11, 9:11], 100.0)
    assert np.all(np.isnan(filled[:, :3]))


def test_difference_statistics_and_aoi_mask():
    grid = aligned_grid((16.0, 48.0, 16.1, 48.1), 0.01, 0.01)
    mask = aoi_mask(grid, (16.02, 48.02, 16.08, 48.08))
    assert mask.sum() == 36
    fine = np.full((grid.height, grid.width), 12.0, np.float32)
    coarse = np.full_like(fine, 10.0)
    statistics = difference_statistics(fine, coarse, mask)
    assert statistics.median == pytest.approx(2.0)
    assert statistics.pixel_count == 36


def test_coarse_part_conversion_and_mixed_vrt(tmp_path):
    # coarse "DEM": 0.01 deg pixels, constant 100 m, with no-data 0 in one pixel
    coarse_grid = aligned_grid((16.0, 48.0, 16.5, 48.5), 0.01, 0.01)
    heights = np.full((coarse_grid.height, coarse_grid.width), 100.0, np.float32)
    heights[0, 0] = 0.0
    coarse_file = tmp_path / "coarse.tif"
    write_array(coarse_file, heights, coarse_grid, nodata=0.0)

    output = aligned_grid((16.1, 48.1, 16.4, 48.4), 0.0025, 0.0025)
    coarse_part = tmp_path / "coarse_part.tif"
    prepare_coarse_part(coarse_file, "ellipsoid", output, coarse_part)
    part = gdal.Open(str(coarse_part))
    assert part.GetRasterBand(1).GetNoDataValue() == OUTPUT_NODATA
    assert part.RasterXSize == 34  # 30 pixels + 2 margin pixels per side

    # fine part: a sub-grid of the output grid with 120 m heights
    _, _, fine_grid = output.subgrid((16.2, 48.2, 16.3, 48.3))
    fine_part = tmp_path / "fine_part.tif"
    write_array(fine_part, np.full((fine_grid.height, fine_grid.width), 120.0, np.float32), fine_grid)

    vrt = tmp_path / "mixed.vrt"
    build_mixed_resolution_vrt(coarse_part, fine_part, output, vrt)
    final = tmp_path / "final.tif"
    write_final_geotiff(vrt, final, {"HEIGHT_REFERENCE": "WGS84 ellipsoid"})

    dataset = gdal.Open(str(final))
    assert (dataset.RasterXSize, dataset.RasterYSize) == (output.width, output.height)
    np.testing.assert_allclose(dataset.GetGeoTransform(), output.geotransform)
    values = dataset.ReadAsArray()
    row = int((output.north - 48.25) / output.delta_lat)
    assert values[row, int((16.25 - output.west) / output.delta_lon)] == pytest.approx(120.0)
    assert values[row, int((16.12 - output.west) / output.delta_lon)] == pytest.approx(100.0)
    assert dataset.GetMetadataItem("HEIGHT_REFERENCE") == "WGS84 ellipsoid"
