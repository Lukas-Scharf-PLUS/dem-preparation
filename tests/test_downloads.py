"""Tests for the tile logic of the BEV and Copernicus downloaders (no network access)."""

from bev_download import (aoi_window_3035, clip_window_to_tile, local_tile_path, tile_url_candidates,
                          tiles_for_window, window_lonlat_bounds)
from copdem_download import tile_name, tile_url, tiles_for_bounds

VIENNA = (16.18099875016676, 48.11744311620298, 16.587587299348154, 48.324307364563005)


def test_vienna_needs_four_bev_tiles():
    window = aoi_window_3035(VIENNA, 1000)
    tiles = tiles_for_window(window)
    assert sorted(tiles) == [(2750000, 4750000), (2750000, 4800000),
                             (2800000, 4750000), (2800000, 4800000)]


def test_window_contains_the_aoi_plus_buffer():
    window = aoi_window_3035(VIENNA, 1000)
    west, south, east, north = window_lonlat_bounds(window)
    assert west < VIENNA[0] and south < VIENNA[1] and east > VIENNA[2] and north > VIENNA[3]


def test_window_is_clipped_to_the_tile():
    window = (4_778_000, 2_796_000, 4_812_000, 2_824_000)
    assert clip_window_to_tile(window, (2_800_000, 4_750_000)) == (4_778_000, 2_800_000, 4_800_000, 2_824_000)


def test_tile_urls_cover_both_naming_schemes():
    urls = tile_url_candidates("20250915", (2800000, 4800000))
    assert urls[0] == ("https://data.bev.gv.at/download/ALS/DSM/20250915/"
                       "ALS_DSM_CRS3035RES50000mN2800000E4800000.tif")
    assert urls[1].endswith("/20250915/CRS3035RES50000mN2800000E4800000.tif")


def test_local_tile_name_contains_edition_and_resolution(tmp_path):
    assert local_tile_path(tmp_path, "20250915", (2800000, 4800000), None).name == \
        "ALS_DSM_20250915_N2800000E4800000_1m_aoi.tif"
    assert local_tile_path(tmp_path, "20250915", (2800000, 4800000), 0).name.endswith("_ovr0_aoi.tif")


def test_copernicus_tile_names():
    assert tile_name(48, 16) == "Copernicus_DSM_COG_10_N48_00_E016_00_DEM"
    assert tile_name(-1, -2) == "Copernicus_DSM_COG_10_S01_00_W002_00_DEM"
    assert tile_url(48, 16).startswith("/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/")


def test_copernicus_tiles_for_bounds():
    assert tiles_for_bounds((15.1, 47.6, 17.9, 48.8)) == [
        (47, 15), (47, 16), (47, 17), (48, 15), (48, 16), (48, 17)]
