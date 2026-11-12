"""Tests for height offsets and the ETRS89 -> ITRF shift.

Tests marked ``network`` use real PROJ grids, which PROJ may download from
cdn.proj.org on first use; they are skipped when that fails.
"""

import numpy as np
import pytest

from vertical_datum import (FrameShift, VerticalDatumError, bev_height_reference, etrs89_to_itrf_shift,
                            height_offset_function, offset_grid)

STEPHANSPLATZ = (16.3725, 48.2085)


def test_bev_height_reference_auto():
    assert bev_height_reference("20250915", "auto") == "evrf2000_austria"
    assert bev_height_reference("20210915", "auto") == "gha"
    assert bev_height_reference("20250915", "gha") == "gha"


def test_offset_grid_interpolates_a_smooth_function():
    def plane(longitudes, latitudes):
        return 40.0 + 2.0 * (np.asarray(longitudes) - 16.0) + 3.0 * (np.asarray(latitudes) - 48.0)

    latitudes = np.linspace(48.3, 48.1, 21)
    longitudes = np.linspace(16.1, 16.6, 51)
    grid = offset_grid(latitudes, longitudes, plane, step_deg=0.05)
    expected = plane(longitudes[None, :], latitudes[:, None])
    assert grid.shape == (21, 51)
    np.testing.assert_allclose(grid, expected, atol=1e-4)


def test_offset_grid_rejects_unphysical_values():
    with pytest.raises(VerticalDatumError):
        offset_grid(np.array([48.2]), np.array([16.3]), lambda lon, lat: np.full(np.shape(lon), 500.0), 0.05)


def test_no_frame_shift_for_etrs89():
    assert etrs89_to_itrf_shift("etrs89", 2024.0, *STEPHANSPLATZ) == FrameShift(0.0, 0.0, 0.0)


def proj_or_skip(function, *args):
    """Call ``function``; skip the test if PROJ cannot get its grid (no network)."""
    try:
        return function(*args)
    except VerticalDatumError as error:
        pytest.skip(f"PROJ grid not available: {error}")


@pytest.mark.parametrize("reference, expected", [
    ("evrf2000_austria", 44.64), ("gha", 44.17), ("egm2008", 44.74)])
def test_real_offsets_in_vienna(reference, expected):
    offset = proj_or_skip(height_offset_function, reference)
    value = float(offset(np.array([STEPHANSPLATZ[0]]), np.array([STEPHANSPLATZ[1]]))[0])
    if value == 0.0:
        pytest.skip("PROJ returned no offset (grid not downloadable)")
    assert value == pytest.approx(expected, abs=0.05)


def test_itrf2020_shift_in_vienna():
    shift = proj_or_skip(etrs89_to_itrf_shift, "itrf2020", 2024.0, *STEPHANSPLATZ)
    east_m = shift.delta_lon_deg * 74_260
    north_m = shift.delta_lat_deg * 111_200
    # ETRF2000 vs ITRF2020 in 2024: ~0.7 m east, ~0.6 m north
    assert east_m == pytest.approx(0.71, abs=0.05)
    assert north_m == pytest.approx(0.60, abs=0.05)
    assert abs(shift.delta_height_m) < 0.05
