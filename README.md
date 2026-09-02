# CONUS-30m-LAI: Monthly 30-m Leaf Area Index for the Contiguous United States (2000–2022)

This repository contains the production, gap-filling, and validation code for the **CONUS-30m-LAI**
dataset — a monthly, 30-meter Leaf Area Index (LAI) product covering the contiguous United States
(CONUS) from 2000 to 2022, together with a **gap-filled** version, a per-pixel **QA** layer, and a
per-year count of valid months (**nvalid**).

## Dataset Overview

| Attribute | Value |
|---|---|
| Variable | Leaf Area Index (LAI) |
| Spatial coverage | Contiguous United States (CONUS) |
| Spatial resolution | 30 m |
| Temporal coverage | 2000–2022 |
| Temporal resolution | Monthly (276 months) |
| Projection | EPSG:5070 (NAD83 / Conus Albers) |
| Storage / scale | **int16**, **scale factor 0.001** (real LAI = stored value × 0.001; valid LAI > 0) |
| No-data | retrieved product: gaps **and** non-vegetation are stored as **0** (valid LAI ≥ 0.1, so 0 is unambiguous). Gap-filled product: pixels that remain unfilled (non-vegetation, permanent water, no climatological basis) are **−32768** |
| Format | tiled GeoTIFF, DEFLATE compression |

Two LAI products are provided:
- **Retrieved LAI** — the observed Random-Forest retrieval; cloud / SLC-off / no-observation pixels, and non-vegetated surfaces, are stored as 0.
- **Gap-filled LAI** — the retrieved values preserved exactly, with gaps filled by a climatology-anchored
  temporal method (`gapfill/`), every filled pixel tagged in the QA layer.

### QA codes (gap-filled product)
| Code | Meaning |
|---|---|
| 0 | original observation (unchanged) |
| 1 | within-year linear interpolation |
| 2 | climatology fill (≥ 3 observed years) |
| 3 | climatology fill (1–2 observed years) |
| 4 | spatial-fallback climatology |
| 255 | not filled — no climatological basis (e.g. permanent water/ocean) |

**nvalid** gives, per pixel per year, the number of months with a valid *observed* LAI (0–12).

## Repository Structure

```
CONUS-30m-LAI-mapping/
├── gee/                                  # Google Earth Engine (retrieval, re-export, validation refs)
│   ├── lai_state_monthly_export.js       # main LAI retrieval, exported per state × year
│   ├── lai_state_JUNE2022_export.js      # re-export of June 2022 (fixes one corrupt month)
│   ├── obs_count_export.js               # per-pixel observation-count export
│   ├── obs_count_analysis.js             # observation-count analysis
│   ├── export_MODIS_VIIRS_median.js      # monthly-median MODIS/VIIRS LAI at 500 m (validation reference)
│   └── export_PFT_yearly.js              # year-matched NLCD biome maps at 500 m (validation strata)
├── postprocessing/                       # mosaic state tiles onto the CONUS grid
│   ├── merge_lai_to_ref_grid.py          # warp + merge state tiles → CONUS monthly mosaics
│   └── build_june2022_conus.py           # reconcile + mosaic the re-exported June 2022
├── gapfill/                              # gap-filling pipeline
│   ├── GAPFILL_PLAN.md                   # method description
│   └── gapfill_lai.py                    # climatology-anchored temporal gap-fill + QA + nvalid
├── validation/                          # independent cross-comparison with MODIS
│   ├── validate_local.py                 # PFT × QA-stratified sampling (retrieved / filled / overall)
│   ├── agg_validate.py                   # scale-matched 30 m → 500 m aggregation vs MODIS
│   ├── period_breakdown.py               # accuracy by period (2000–05 / 06–18 / 19–22)
│   └── period_pft_breakdown.py           # accuracy by period × biome
├── download/
│   └── download_conus_lai.py             # download the published GEE asset to local GeoTIFFs
├── LICENSE
└── README.md
```

## Methods Summary

### 1 — LAI retrieval (Google Earth Engine) — `gee/lai_state_monthly_export.js`
Monthly LAI is retrieved state-by-state with the **biome-stratified Random Forest** algorithm of
[Kang et al. (2021)](https://doi.org/10.1016/j.rse.2021.112383):
1. Landsat C02 L2 surface reflectance (Landsat 5/7/8/9) is loaded, scaled, and cloud/shadow/water masked via `QA_PIXEL`.
2. Spectral indices (NDVI, NDWI) and solar geometry are computed.
3. Each pixel is assigned to one of 9 biome types using **year-matched NLCD** (2001/2004/2006/2008/2011/2013/2016/2019/2021 epochs).
4. A Random Forest (100 trees, trained on 2006–2018 MODIS LAI) is applied per sensor × biome.
5. Valid pixels are composited by **monthly median** and exported as **int16 (LAI × 1000)**.

Each run processes one state × one year → 12 monthly GeoTIFFs.

### 2 — CONUS mosaicking (Python) — `postprocessing/merge_lai_to_ref_grid.py`
State tiles are warped (nearest-neighbor) onto the common EPSG:5070 CONUS reference grid and merged
(valid-wins) into 12 monthly mosaics per year. June 2022 was re-exported and reconciled separately
via `build_june2022_conus.py`.

### 3 — Gap-filling (Python) — `gapfill/gapfill_lai.py`
Gaps are filled with a **climatology-anchored temporal** method (Kandasamy et al., 2013): a robust
per-pixel, per-month climatology is built across years; within-year interior gaps are linearly
interpolated; remaining gaps are filled from the amplitude-scaled climatology, with a coarse spatial
fallback where a pixel's own climatology is undefined. **Original observations are preserved exactly**
(`valid = LAI > 0`), and every filled pixel is flagged in the QA layer. Pixels with no climatological
basis (permanent water/ocean) are left unfilled (QA 255). See `gapfill/GAPFILL_PLAN.md`.

### 4 — Validation (Python) — `validation/`
Retrieved and gap-filled LAI are cross-compared against **monthly-median MODIS (MOD15A2H)** at 500 m,
stratified by biome and period. The scale-matched comparison (`agg_validate.py`; 30 m aggregated to the
500 m MODIS grid) gives, over all biomes and years:

| | R² | RMSE (LAI) |
|---|---|---|
| Retrieved vs MODIS | 0.77 | 0.83 |
| Gap-filled vs MODIS | 0.76 | 0.83 |

Agreement is **stable across 2000–2005, 2006–2018, and 2019–2022** (R² ≈ 0.76–0.77) — no degradation
outside the RF training window — and gap-filled pixels agree with MODIS as well as the observations.

## Requirements
- **GEE scripts:** a Google Earth Engine account.
- **Python:** `pip install rasterio numpy pandas` (gap-fill / validation); `geedim` for `download/`.
- Input/output paths are defined at the top of each Python script — adjust to your environment before running.

## Products
- **GEE asset (ImageCollection):** `projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m` — 276 monthly gap-filled LAI images.
  - https://code.earthengine.google.com/?asset=projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m
- Companion **QA** and **nvalid** collections: `…/CONUS_Monthly_LAI_30m_QA`, `…/CONUS_Monthly_LAI_30m_nvalid`.
- Archived dataset (DOI): *to be added on deposit.*

## Citation
If you use this dataset or code, please cite:

> You, H., et al. (2026). A 30-m monthly leaf area index dataset for the contiguous United States from 2000 to 2022. *Scientific Data*. [DOI TBD]

Underlying LAI algorithm:

> Kang, Y., Ozdogan, M., Gao, F., Anderson, M. C., White, W. A., Yang, Y., Yang, Y., & Erickson, T. A. (2021).
> A data-driven approach to estimate leaf area index for Landsat images over the contiguous US.
> *Remote Sensing of Environment*, 258, 112383. https://doi.org/10.1016/j.rse.2021.112383

Gap-filling method:

> Kandasamy, S., Baret, F., Verger, A., Neveux, P., & Weiss, M. (2013). A comparison of methods for smoothing
> and gap filling time series of remote sensing observations. *Biogeosciences*, 10, 4055–4071.

## License
Code: [MIT License](LICENSE) · Dataset: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
