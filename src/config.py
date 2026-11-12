"""
Configuration handling for the DEM preparation.

All settings come from one YAML file (see ``configs/vienna.yaml``). A few
values can be overridden on the command line (see ``prepare_dem.py``). After
loading, every value is validated so that mistakes are reported before
anything is downloaded or computed.
"""

from dataclasses import dataclass, field
from pathlib import Path

import yaml

# Height references the tool can convert to WGS84 ellipsoid heights.
COARSE_VERTICAL_DATUMS = ("ellipsoid", "egm2008", "egm96")
FINE_VERTICAL_DATUMS = ("auto", "evrf2000_austria", "gha")
COARSE_SOURCES = ("file", "download")
HORIZONTAL_FRAMES = ("etrs89", "itrf2014", "itrf2020")

MAX_RELATIVE_ORBIT = 175  # Sentinel-1 has 175 relative orbits per 12-day cycle
# BEV editions are folders named by their reference date, e.g. 20250915.
EDITION_LENGTH = 8


class ConfigError(ValueError):
    """Raised when the configuration is missing values or is invalid."""


@dataclass
class AoiSettings:
    """Area of interest: where the fine (BEV) DEM is used."""

    name: str
    bbox: tuple  # (west, south, east, north) in degrees


@dataclass
class BurstSettings:
    """Which Sentinel-1 bursts the DEM must cover."""

    burst_db: Path
    relative_orbits: list
    margin_deg: float


@dataclass
class CoarseDemSettings:
    """The 30 m DEM used outside the AOI."""

    source: str  # "file" or "download"
    file: Path | None
    vertical_datum: str


@dataclass
class FineDemSettings:
    """The BEV ALS DSM used inside the AOI."""

    edition: str
    vertical_datum: str
    buffer_m: float
    overview_level: int | None


@dataclass
class MergeSettings:
    """How both DEMs are combined on the final grid."""

    resolution_m: float
    blend_width_m: float
    hole_fill_max_m: float
    horizontal_frame: str
    frame_epoch: float


@dataclass
class OutputSettings:
    """Which output files are written."""

    name: str
    write_geotiff: bool = True
    write_mixed_resolution_vrt: bool = True
    write_isce: bool = True


@dataclass
class Settings:
    """All settings of one DEM preparation run."""

    aoi: AoiSettings
    bursts: BurstSettings
    coarse_dem: CoarseDemSettings
    fine_dem: FineDemSettings
    merge: MergeSettings
    output: OutputSettings
    work_dir: Path
    output_dir: Path
    config_file: Path | None = field(default=None)

    @property
    def output_base_name(self):
        """Base name of the output files, e.g. ``vienna_bevdsm_cop30_ellipsoid_10m``."""
        resolution = f"{self.merge.resolution_m:g}".replace(".", "p")
        return f"{self.output.name}_{resolution}m"


def load_yaml_config(config_file):
    """Read a YAML config file and return its content as a dict.

    Raises
    ------
    ConfigError
        If the file does not exist or does not contain a YAML mapping.
    """
    config_file = Path(config_file)
    if not config_file.is_file():
        raise ConfigError(f"Config file not found: {config_file}")
    with open(config_file, encoding="utf-8") as stream:
        content = yaml.safe_load(stream)
    if not isinstance(content, dict):
        raise ConfigError(f"Config file {config_file} must contain a YAML mapping.")
    return content


def required_section(values, name):
    """Return the sub-dict ``values[name]`` or raise a ConfigError naming the missing section."""
    section = values.get(name)
    if not isinstance(section, dict):
        raise ConfigError(f"Config section '{name}' is missing.")
    return section


def required_value(section, section_name, key):
    """Return ``section[key]`` or raise a ConfigError naming the missing key."""
    if section.get(key) is None:
        raise ConfigError(f"Config value '{section_name}.{key}' is missing.")
    return section[key]


def parse_bbox(value):
    """Validate a [west, south, east, north] box in degrees and return it as a tuple.

    Raises
    ------
    ConfigError
        If the box does not have four numbers or west/east, south/north are swapped.
    """
    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise ConfigError(f"aoi.bbox must be [west, south, east, north], got {value!r}")
    west, south, east, north = (float(v) for v in value)
    if not (-180 <= west < east <= 180) or not (-90 <= south < north <= 90):
        raise ConfigError(
            f"aoi.bbox must be [west, south, east, north] with west < east and "
            f"south < north, got {value!r}"
        )
    return (west, south, east, north)


def parse_relative_orbits(value):
    """Validate the list of relative orbits (1..175) and return it as a list of ints."""
    if not isinstance(value, (list, tuple)) or not value:
        raise ConfigError("bursts.relative_orbits must be a non-empty list, e.g. [22, 73].")
    orbits = [int(orbit) for orbit in value]
    invalid = [orbit for orbit in orbits if not 1 <= orbit <= MAX_RELATIVE_ORBIT]
    if invalid:
        raise ConfigError(f"Relative orbits must be between 1 and {MAX_RELATIVE_ORBIT}: {invalid}")
    return orbits


def parse_choice(value, allowed, name):
    """Return ``value`` (lower case) if it is one of ``allowed``, else raise a ConfigError."""
    value = str(value).lower()
    if value not in allowed:
        raise ConfigError(f"{name} must be one of {', '.join(allowed)}, got {value!r}")
    return value


def parse_positive(value, name):
    """Return ``value`` as a float if it is > 0, else raise a ConfigError."""
    number = float(value)
    if number <= 0:
        raise ConfigError(f"{name} must be > 0, got {value!r}")
    return number


def build_settings(values, config_file=None):
    """Turn the raw config dict into validated Settings.

    Parameters
    ----------
    values : dict
        Content of the YAML config (after command-line overrides).
    config_file : Path, optional
        Where the values came from; only stored for logging.

    Raises
    ------
    ConfigError
        If a value is missing or invalid.
    """
    aoi = required_section(values, "aoi")
    bursts = required_section(values, "bursts")
    paths = required_section(values, "paths")
    coarse = required_section(values, "coarse_dem")
    fine = required_section(values, "fine_dem")
    merge = required_section(values, "merge")
    output = required_section(values, "output")

    edition = str(required_value(fine, "fine_dem", "edition"))
    if len(edition) != EDITION_LENGTH or not edition.isdigit():
        raise ConfigError(f"fine_dem.edition must look like 20250915, got {edition!r}")

    overview_level = fine.get("overview_level")
    if overview_level is not None and int(overview_level) < 0:
        raise ConfigError(f"fine_dem.overview_level must be null or >= 0, got {overview_level!r}")

    settings = Settings(
        aoi=AoiSettings(
            name=str(required_value(aoi, "aoi", "name")),
            bbox=parse_bbox(required_value(aoi, "aoi", "bbox")),
        ),
        bursts=BurstSettings(
            burst_db=Path(required_value(bursts, "bursts", "burst_db")).expanduser(),
            relative_orbits=parse_relative_orbits(required_value(bursts, "bursts", "relative_orbits")),
            margin_deg=float(bursts.get("margin_deg", 0.3)),
        ),
        coarse_dem=CoarseDemSettings(
            source=parse_choice(coarse.get("source", "file"), COARSE_SOURCES, "coarse_dem.source"),
            file=Path(coarse["file"]).expanduser() if coarse.get("file") else None,
            vertical_datum=parse_choice(required_value(coarse, "coarse_dem", "vertical_datum"),
                                        COARSE_VERTICAL_DATUMS, "coarse_dem.vertical_datum"),
        ),
        fine_dem=FineDemSettings(
            edition=edition,
            vertical_datum=parse_choice(fine.get("vertical_datum", "auto"),
                                        FINE_VERTICAL_DATUMS, "fine_dem.vertical_datum"),
            buffer_m=parse_positive(fine.get("buffer_m", 1000), "fine_dem.buffer_m"),
            overview_level=None if overview_level is None else int(overview_level),
        ),
        merge=MergeSettings(
            resolution_m=parse_positive(required_value(merge, "merge", "resolution_m"),
                                        "merge.resolution_m"),
            blend_width_m=parse_positive(merge.get("blend_width_m", 300), "merge.blend_width_m"),
            hole_fill_max_m=float(merge.get("hole_fill_max_m", 100)),
            horizontal_frame=parse_choice(merge.get("horizontal_frame", "itrf2020"),
                                          HORIZONTAL_FRAMES, "merge.horizontal_frame"),
            frame_epoch=float(merge.get("frame_epoch", 2024.0)),
        ),
        output=OutputSettings(
            name=str(required_value(output, "output", "name")),
            write_geotiff=bool(output.get("write_geotiff", True)),
            write_mixed_resolution_vrt=bool(output.get("write_mixed_resolution_vrt", True)),
            write_isce=bool(output.get("write_isce", True)),
        ),
        work_dir=Path(required_value(paths, "paths", "work_dir")).expanduser(),
        output_dir=Path(required_value(paths, "paths", "output_dir")).expanduser(),
        config_file=config_file,
    )
    validate_settings(settings)
    return settings


def validate_settings(settings):
    """Check rules that involve more than one value.

    Raises
    ------
    ConfigError
        If the settings contradict each other.
    """
    if settings.fine_dem.buffer_m <= settings.merge.blend_width_m:
        # Otherwise the blend zone reaches into the AOI and parts of Vienna
        # would get a mix of BEV and Copernicus heights.
        raise ConfigError(
            f"fine_dem.buffer_m ({settings.fine_dem.buffer_m:g}) must be larger than "
            f"merge.blend_width_m ({settings.merge.blend_width_m:g})."
        )
    if settings.coarse_dem.source == "file" and settings.coarse_dem.file is None:
        raise ConfigError("coarse_dem.file is required when coarse_dem.source is 'file'.")
    if settings.merge.resolution_m > 30:
        raise ConfigError(
            f"merge.resolution_m ({settings.merge.resolution_m:g}) is coarser than the "
            "30 m Copernicus DEM; that would throw away detail everywhere."
        )


def load_settings(config_file, overrides=None):
    """Load, override and validate the config.

    Parameters
    ----------
    config_file : str or Path
    overrides : dict, optional
        ``{"section.key": value}`` pairs from the command line; ``None`` values are ignored.

    Returns
    -------
    Settings
    """
    values = load_yaml_config(config_file)
    for dotted_key, value in (overrides or {}).items():
        if value is None:
            continue
        section, key = dotted_key.split(".", 1)
        values.setdefault(section, {})[key] = value
    return build_settings(values, Path(config_file))
