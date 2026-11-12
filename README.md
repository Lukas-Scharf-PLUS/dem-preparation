# dem_preparation

Builds **one DEM for the Sentinel-1 InSAR pipelines of Vienna** (ISCE2 and COMPASS):

- **inside Vienna:** the BEV airborne-laser-scanning **surface model (ALS DSM, 1 m)**,
  averaged to the output resolution (default **10 m**),
- **around Vienna:** the **Copernicus DEM GLO-30** (your `cop_dem.tif`), so the DEM covers
  every burst of the relative orbits 22, 73, 124 and 146 that touches Vienna,
- **heights above the WGS84 ellipsoid**, which ISCE2 and COMPASS expect,
- written as a **GeoTIFF**, a space-saving **mixed-resolution VRT** and an **ISCE2 DEM**
  (`.dem.wgs84` + `.xml` + `.vrt`).

The goal is better geolocation of scatterers: a 10 m ALS surface model in the city is much
closer to the real roof and street heights than the 30 m Copernicus DEM.

> Work in progress. Check the outputs (see [Checking the result](#checking-the-result))
> before you use them in production runs.

---

## Credentials

**None are needed.** All data is public and is downloaded over HTTPS:

| Host | What | When |
| --- | --- | --- |
| `data.bev.gv.at` | BEV ALS DSM tiles (open data, free of charge) | `download` step |
| `cdn.proj.org` | official height/geoid grids for PROJ (downloaded once, cached) | `build` step |
| `copernicus-dem-30m.s3.amazonaws.com` | Copernicus GLO-30 (AWS Open Data) | only with `coarse_dem.source: download` |

---

## Setup

```bash
conda env create -f environment.yml                                  # first time
conda env update -n dem-preparation -f environment.yml --prune      # after changes
conda activate dem-preparation
```

The `isce2-dolphin-mintpy` environment also contains every package needed, so it works too.

You also need the **OPERA burst database** `opera-burst-bbox-only.sqlite3` (the file COMPASS
uses; already at `data/input_data/s1-burst-db/`). If it is missing, download it from the
[opera-adt/burst_db releases](https://github.com/opera-adt/burst_db/releases) and set
`bursts.burst_db`.

---

## Usage

```bash
# everything: coverage check -> download -> build
python src/prepare_dem.py --config configs/vienna.yaml all

# or step by step
python src/prepare_dem.py --config configs/vienna.yaml check      # bursts + coverage only, no download
python src/prepare_dem.py --config configs/vienna.yaml download   # BEV (+ GLO-30 if configured)
python src/prepare_dem.py --config configs/vienna.yaml build      # merge and write the outputs
```

Command-line options override the config:

| Option | Overrides | Example |
| --- | --- | --- |
| `--resolution_m` | `merge.resolution_m` | `--resolution_m 5` |
| `--edition` | `fine_dem.edition` | `--edition 20240915` |
| `--coarse_source` | `coarse_dem.source` | `--coarse_source download` |

Every step runs the coverage check first. `download` and `build` reuse files that already
exist, so you can repeat them, or change the resolution and run `build` again without
downloading anything.

**Time and disk space (Vienna, measured):** with the 2 m overview (`fine_dem.overview_level: 0`)
the download is 0.54 GB and the whole run takes about 1 minute. The full 1 m data is about
4× more download (~2–2.5 GB, a few minutes). Outputs at 10 m: 1.07 GB GeoTIFF, 1.9 GB
ISCE2 raw file, 0.21 GB for the mixed-resolution parts.

### Using the result

All outputs are in `paths.output_dir`, named `<output.name>_<resolution>m*`:

| File | Use it for |
| --- | --- |
| `…_10m.tif` | **COMPASS** (`dem_file`), **isce2-dolphin-mintpy** (`paths.dem_file` with `paths.dem_vertical_datum: ellipsoid`), the Argo pipelines (upload instead of `cop_dem.tif`) |
| `…_10m_mixed_resolution.vrt` + `…_10m_parts/` | same as the GeoTIFF, but 30 m outside the BEV area on disk (see below). Keep the folder next to the VRT. |
| `…_10m.dem.wgs84` (+ `.xml`, `.vrt`) | **ISCE2 directly**: `stackSentinel.py -d …_10m.dem.wgs84` |
| `…_10m_report.json` | what was done: bursts, grids, sources, frame shift, BEV − Copernicus statistics |
| `…_10m.log` | the log of all runs |

All files hold **WGS84 ellipsoid heights** on an EPSG:4326 grid. Do not convert them again:
in pipelines that ask for the height reference, set it to `ellipsoid`.

---

## Your question: does the resolution have to be the same everywhere?

**For the files that ISCE2 and COMPASS read: yes.** Both read one regular raster. A raster
has one pixel size, so the final DEM has the fine resolution everywhere. Outside Vienna,
the 30 m Copernicus heights are only interpolated to the fine grid. This adds no
information and no harm, it only costs disk space.

**Disk space can still be saved**, in two ways this repository uses:

1. **The DEM only covers what the bursts need.** It does not cover all of `cop_dem.tif`
   (14–19° E, 47–50° N), only the bursts + margin (≈ 14.8–18.3° E, 47.4–49.1° N).
2. **Mixed-resolution VRT.** `…_mixed_resolution.vrt` stores the area outside the BEV
   data at 30 m (`_parts/coarse_ellipsoid.tif`) and only the city at 10 m
   (`_parts/fine_10m_ellipsoid.tif`). GDAL builds the fine grid on the fly when the file
   is read. COMPASS and isce2-dolphin-mintpy read DEMs through GDAL, so they can use the
   VRT directly. That needs 0.21 GB instead of 1.07 GB. The ISCE2 raw format cannot do
   this; `stackSentinel.py` needs the full raw file.

Does a finer resolution outside the city matter? Only where there is real fine-scale
information. Outside Vienna the source is the 30 m Copernicus DEM. Its effect on ISCE2 is
also small: ISCE2 interpolates the DEM onto the radar grid (~2.3 m × 14 m), so the DEM
pixel size changes only how detailed the heights are, not the size of the products.

The resolution inside the city is `merge.resolution_m` (or `--resolution_m`). 10 m is a
good match for Sentinel-1 IW (range ~2.3 m × azimuth ~14 m). Values below ~5 m need a lot
of memory in the `build` step (the whole fine part is processed in memory), and for 2–4 m
set `fine_dem.overview_level: null` (full 1 m data).

---

## How it works

```
check     OPERA burst DB ──► bursts of orbits 22/73/124/146 touching the AOI
                         ──► union of their boxes + margin = area the DEM must cover
                         ──► does coarse_dem.file cover it?  (no → source: download)

download  BEV ALS DSM COG tiles (EPSG:3035) ──► only the AOI + buffer window, via HTTP ranges
          [GLO-30 from AWS, only if coarse_dem.source: download]

build     coarse DEM ──► cut to the area, convert to ellipsoid heights (30 m)  = coarse part
          BEV window ──► average to the 10 m lat/lon grid, ETRS89 → ITRF shift,
                         + height offset (EVRF2000 Austria → ellipsoid), fill small holes,
                         blend into the coarse heights at the edge               = fine part
          coarse part + fine part ──► mixed-resolution VRT ──► GeoTIFF ──► ISCE2 DEM
```

### 1. Which bursts, which area (`src/burst_coverage.py`)

The OPERA burst database stores a box for every Sentinel-1 burst (in UTM, already with a
5 km margin; COMPASS geocodes each burst onto exactly this box). All bursts of the
configured orbits whose box touches the AOI are selected. The DEM must cover the union of
their boxes plus `bursts.margin_deg` (0.3°). The margin is needed for ISCE2.
`stackSentinel.py` cuts the DEM with lat/lon rectangles computed for heights up to 9000 m
and adds 0.1°, and these rectangles reach ~0.3° beyond the real burst.

For Vienna (bbox from `configs/vienna.yaml`):

| Orbit | Bursts touching Vienna | Burst area |
| --- | --- | --- |
| 22 | 4 (IW1) | 15.56–17.05° E, 47.76–48.70° N |
| 73 | 8 (IW1, IW2) | 15.14–17.86° E, 47.66–48.75° N |
| 124 | 7 (IW2, IW3) | 15.39–17.96° E, 47.79–48.71° N |
| 146 | 4 (IW3) | 15.46–16.84° E, 47.75–48.67° N |
| **DEM must cover** (+0.3°) | | **14.84–18.26° E, 47.36–49.05° N** |

`cop_dem.tif` covers 14–19° E, 47–50° N, so **your existing Copernicus DEM covers all
bursts** and nothing has to be downloaded from AWS. If you add orbits or a larger AOI and
the check fails, set `coarse_dem.source: download` (or `--coarse_source download`). The
tool then downloads GLO-30 for exactly the needed area.

### 2. BEV download (`src/bev_download.py`)

The BEV publishes the ALS DSM as 50 km × 50 km tiles in EPSG:3035, one folder per edition
(reference date, *Stichtag*), e.g.
`https://data.bev.gv.at/download/ALS/DSM/20250915/ALS_DSM_CRS3035RES50000mN2800000E4750000.tif`
(`N…E…` = lower-left corner). Vienna touches **4 tiles** (N2750000/N2800000 ×
E4750000/E4800000). Each tile is a 6–12 GB Cloud-Optimized GeoTIFF. GDAL reads only the
part that covers Vienna + `fine_dem.buffer_m` through HTTP range requests, so about 2–3 GB
are transferred, not 35 GB. Downloads are written as `.part` files and renamed when
complete, so an interrupted download is redone.

Editions on the server: 20190915, 20220915 – 20250915 (each is the newest data of all of
Austria at that date). `fine_dem.edition` chooses one.

### 3. Heights and positions (`src/vertical_datum.py`)

ISCE2 and COMPASS expect **heights above the WGS84 ellipsoid**. The inputs use other heights:

| Source | Height reference | Offset to ellipsoid in Vienna |
| --- | --- | --- |
| BEV ALS DSM (all editions except 2021) | EVRF2000 Austria, EPSG:9274 | +44.6 m |
| BEV ALS DSM edition 20210915 | Gebrauchshöhen Adria (GHA), EPSG:5778 | +44.2 m |
| Copernicus GLO-30 (original tiles) | EGM2008 geoid | +44.7 m |
| `cop_dem.tif` | already WGS84 ellipsoid | 0 |

The offsets come from the official PROJ grids: the BEV's own Austrian grid
(`at_bev_GEOID_GRS80_Oesterreich.tif`, 5 cm accuracy) and NGA's EGM2008 grid. PROJ fetches
them once from `cdn.proj.org`. If a grid cannot be loaded, the tool **stops with an error**
instead of silently using no offset.

Why this matters: a height error *Δh* moves a point in the radar geometry by about
*Δh / tan(incidence)*. A wrong reference (~45 m) gives errors of 50–80 m. Mixing up GHA and
EVRF2000 (0.5 m) gives ~0.6 m.

**Horizontal frame.** BEV positions are in ETRS89, which moves with the Eurasian plate.
The Sentinel-1 orbits, and WGS84 today, are in ITRF. In 2024 the two differ by about
**0.7 m east and 0.6 m north** in Vienna. With `merge.horizontal_frame: itrf2020`
(default), the BEV data is shifted by this amount at `merge.frame_epoch` (set it to the
middle of your time series). The shift is the same everywhere in the area to well below
1 mm, so one constant shift is used. The Copernicus DEM is not shifted, because its
absolute accuracy is several meters anyway.

### 4. Merging (`src/dem_merge.py`)

- **Fine part:** the 1 m BEV pixels inside each output pixel are **averaged**
  (no nearest-neighbour picking of single roof or street pixels), converted as described
  above, and holes up to `merge.hole_fill_max_m` are interpolated from the BEV data around
  them.
- **Blending:** BEV and Copernicus differ by a few meters (different dates, sensors and
  resolution). Over `merge.blend_width_m` (300 m) at the edge of the BEV area, the weight
  of the BEV heights rises linearly from 0 to 1, so there is no height step. The edge lies
  `fine_dem.buffer_m` (1000 m) outside the AOI, so **all of Vienna has pure BEV heights**.
- **Plausibility check:** the median of *BEV − Copernicus* inside the AOI is logged and
  stored in the report. It should be within ±1–2 m. A value of tens of meters means a
  height reference is set wrong, and the tool then prints a warning.

### 5. ISCE2 export (`src/isce_export.py`)

Writes the raw float32 file and the same XML header ISCE2's `renderHdr()` writes, without
needing ISCE2 installed. ISCE2 counts the first coordinate from the *center* of the first
pixel, GDAL from its *corner*, so the half-pixel difference is applied. (ISCE2's own VRT
writer and `gdal2isce_xml.py` skip this; the `.vrt` written here is correct.) The test
suite checks that ISCE2 reads the header correctly when ISCE2 is installed.

---

## Checking the result

1. Read `…_report.json` → `bev_minus_coarse_in_aoi_m.median` should be within ±1–2 m
   (test run: +0.24 m; 5 % / 95 %: −6.5 / +7.5 m, which is normal for a 30 m vs. ALS DSM).
2. Compare with a terrain model on open ground, where DSM = DTM. In the test run, the most
   common value of *DSM − older BEV ALS DTM* (`data/dem/vienna_bev_als_dtm_…dem.wgs84`) was
   0.00–0.25 m. Single points such as squares are less reliable: a 10 m DSM pixel averages
   everything in it (statues, trees, building edges).
3. Look at a hillshade in QGIS: no visible step at the edge of the BEV area, no holes.
   `…_report.json` → `nodata_pixels` should be 0.

---

## Configuration

See [`configs/vienna.yaml`](configs/vienna.yaml); every value is commented there.

| Section | Key | Meaning |
| --- | --- | --- |
| `aoi` | `bbox` | `[west, south, east, north]` in degrees, area with BEV heights |
| `bursts` | `burst_db`, `relative_orbits`, `margin_deg` | which bursts the DEM must cover |
| `paths` | `work_dir`, `output_dir` | downloads/intermediate files, final files |
| `coarse_dem` | `source` (`file`/`download`), `file`, `vertical_datum` (`ellipsoid`/`egm2008`/`egm96`) | the 30 m DEM |
| `fine_dem` | `edition`, `vertical_datum` (`auto`/`evrf2000_austria`/`gha`), `buffer_m`, `overview_level` | the BEV DSM; `overview_level: 0` downloads the 2 m overview (4× less data) |
| `merge` | `resolution_m`, `blend_width_m`, `hole_fill_max_m`, `horizontal_frame`, `frame_epoch` | how the DEMs are combined |
| `output` | `name`, `write_geotiff`, `write_mixed_resolution_vrt`, `write_isce` | which files are written |

---

## Project structure

```
dem_preparation/
├── CLAUDE.md              instructions for Claude (same as the other ADUCAT repositories)
├── README.md              this file
├── environment.yml        conda environment "dem-preparation"
├── configs/
│   └── vienna.yaml        config for Vienna, orbits 22/73/124/146
├── src/
│   ├── prepare_dem.py     command line: check / download / build / all
│   ├── config.py          loading and validating the YAML config
│   ├── burst_coverage.py  bursts touching the AOI, area the DEM must cover, coverage check
│   ├── bev_download.py    BEV ALS DSM: tile selection, windowed COG download
│   ├── copdem_download.py Copernicus GLO-30 download from AWS (optional)
│   ├── vertical_datum.py  height offsets to the ellipsoid, ETRS89 → ITRF shift (PROJ)
│   ├── target_grid.py     output lat/lon grid from a resolution in meters
│   ├── dem_merge.py       coarse part, fine part, blending, VRT, GeoTIFF, report
│   └── isce_export.py     ISCE2 DEM (.dem.wgs84 + .xml + .vrt)
└── tests/                 pytest tests (synthetic data; a few use PROJ grids if reachable)
```

---

## Tests

```bash
conda activate dem-preparation
python -m pytest -q
```

The tests use small synthetic rasters and a fake burst database and do not download DEM
data. The tests for the real Vienna height offsets and the ITRF shift need PROJ grids from
`cdn.proj.org` and are skipped when they cannot be loaded. The ISCE2 header test runs only
where ISCE2 is installed (e.g. in `isce2-dolphin-mintpy`).

---

## Notes and limitations

- **DSM, not DTM.** The DSM contains buildings and trees. That is what the radar sees and
  what persistent scatterers on roofs sit on. Remember that heights estimated by the
  InSAR time series are then relative to the surface, not the ground.
- **Data dates differ.** Copernicus is from 2011–2015, the BEV edition is newer. Buildings
  built or torn down in between create local differences. These only affect the blend zone
  outside Vienna.
- **BEV terms of use.** The product is free of charge. The BEV product page links its
  terms of use. Cite the BEV as data source in publications.
- **The coverage check is an estimate.** It uses the OPERA burst boxes plus a margin.
  `isce2-dolphin-mintpy` checks again with the real burst geometry and reports the area it
  needs if the DEM is too small.
