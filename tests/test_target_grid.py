"""Tests for the output grid: pixel size in degrees and aligned sub-grids."""

import pytest

from target_grid import aligned_grid, meters_per_degree, output_grid, pixel_size_degrees


def test_meters_per_degree_in_vienna():
    per_lat, per_lon = meters_per_degree(48.2)
    assert per_lat == pytest.approx(111_200, abs=100)
    assert per_lon == pytest.approx(74_260, abs=100)


def test_pixel_size_is_about_ten_meters_on_the_ground():
    delta_lon, delta_lat = pixel_size_degrees(10.0, 48.2)
    per_lat, per_lon = meters_per_degree(48.2)
    assert delta_lon * per_lon == pytest.approx(10.0)
    assert delta_lat * per_lat == pytest.approx(10.0)


def test_grid_covers_the_requested_bounds():
    bounds = (15.1, 47.6, 18.0, 48.8)
    grid = output_grid(bounds, 10.0, 48.2)
    assert grid.west <= bounds[0] and grid.south <= bounds[1]
    assert grid.east >= bounds[2] and grid.north >= bounds[3]
    # less than one pixel too large on each side
    assert bounds[0] - grid.west < grid.delta_lon
    assert grid.north - bounds[3] < grid.delta_lat


def test_subgrid_lines_up_with_the_main_grid():
    grid = aligned_grid((16.0, 48.0, 17.0, 49.0), 0.01, 0.01)
    col, row, sub = grid.subgrid((16.123, 48.2, 16.4, 48.555))
    assert sub.west == pytest.approx(grid.west + col * grid.delta_lon)
    assert sub.north == pytest.approx(grid.north - row * grid.delta_lat)
    assert sub.west <= 16.123 and sub.east >= 16.4
    assert sub.south <= 48.2 and sub.north >= 48.555


def test_subgrid_is_clipped_to_the_main_grid():
    grid = aligned_grid((16.0, 48.0, 17.0, 49.0), 0.01, 0.01)
    col, row, sub = grid.subgrid((15.5, 48.5, 16.5, 49.5))
    assert (col, row) == (0, 0)
    assert sub.west == pytest.approx(16.0) and sub.north == pytest.approx(49.0)


def test_subgrid_outside_raises():
    grid = aligned_grid((16.0, 48.0, 17.0, 49.0), 0.01, 0.01)
    with pytest.raises(ValueError):
        grid.subgrid((18.0, 48.0, 19.0, 49.0))
