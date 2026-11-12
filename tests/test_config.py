"""Tests for loading and validating the YAML config."""

from pathlib import Path

import pytest
import yaml

from config import ConfigError, build_settings, load_settings

REPO_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "vienna.yaml"


def valid_values():
    """Return the repository's Vienna config as a dict (a known-good starting point)."""
    return yaml.safe_load(REPO_CONFIG.read_text(encoding="utf-8"))


def test_repository_config_is_valid():
    settings = build_settings(valid_values())
    assert settings.aoi.bbox[0] < settings.aoi.bbox[2]
    assert settings.bursts.relative_orbits == [22, 73, 124, 146]
    assert settings.output_base_name == "vienna_bevdsm_cop30_ellipsoid_10m"


def test_command_line_override_wins(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text(yaml.safe_dump(valid_values()), encoding="utf-8")
    settings = load_settings(config_file, {"merge.resolution_m": 5, "fine_dem.edition": None})
    assert settings.merge.resolution_m == 5
    assert settings.fine_dem.edition == "20250915"
    assert settings.output_base_name.endswith("_5m")


@pytest.mark.parametrize("bbox", [[16.5, 48.1, 16.2, 48.3], [16.2, 48.3, 16.5, 48.1], [16.2, 48.1, 16.5]])
def test_invalid_bbox_is_rejected(bbox):
    values = valid_values()
    values["aoi"]["bbox"] = bbox
    with pytest.raises(ConfigError):
        build_settings(values)


def test_blend_zone_must_fit_into_buffer():
    values = valid_values()
    values["merge"]["blend_width_m"] = 2000
    with pytest.raises(ConfigError, match="buffer_m"):
        build_settings(values)


def test_unknown_vertical_datum_is_rejected():
    values = valid_values()
    values["coarse_dem"]["vertical_datum"] = "nn"
    with pytest.raises(ConfigError, match="vertical_datum"):
        build_settings(values)


def test_missing_section_is_reported():
    values = valid_values()
    del values["merge"]
    with pytest.raises(ConfigError, match="merge"):
        build_settings(values)
