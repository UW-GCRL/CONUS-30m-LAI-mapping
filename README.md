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
| Projection / grid | Albers Equal Area Conic on the WGS 84 datum, with the EPSG:5070 (NAD83 / Conus Albers) projection parameters (standard parallels 29.5° N and 45.5° N, latitude of origin 23° N, central meridian 96° W); 30 m pixels, grid origin (−2417805, 3321225), 160,169 × 105,432 pixels |
| Storage / scale | **int16**, **scale factor 0.001** (real LAI = stored value × 0.001; valid LAI > 0) |
| No-data | **−32768** in both LAI products: pixels without a value (unfilled pixel-years, which include non-vegetated pixels and water, and the area outside CONUS) are masked in Earth Engine and are −32768 (declared NoData) in GeoTIFF exports made as shown under "Exporting to GeoTIFF". QA: 255 (masked). nvalid: 0 (masked) |
| Format | Earth Engine ImageCollections (export to GeoTIFF from Earth Engine) |

Two LAI products:
- **Retrieved LAI**: the Random Forest retrieval, monthly median of valid retrievals. Pixel-months without
  a valid retrieval (no clear observation, or retrievals rejected by the quality checks) and
  non-vegetated pixels have no value. It equals
  the gap-filled LAI where QA = 0, so it is obtained from the published collections as
  `LAI.updateMask(QA.eq(0))` (−32768 elsewhere when exported as shown below).
- **Gap-filled LAI**: the retrieved values preserved exactly, with gaps filled by a climatology-anchored
  temporal method (`gapfill/`); every value is labelled in the QA layer.

### QA codes (gap-filled product)
| Code | Meaning |
|---|---|
| 0 | retrieved value (not filled) |
| 1 | linear interpolation between two valid months of the same year |
| 2 | climatology fill; the pixel-year has at least three valid months |
| 3 | climatology fill; the pixel-year has one or two valid months |
| 4 | climatology fill using the spatial-fallback (block) climatology, because fewer than three years were valid for that calendar month |
| 255 | not filled: no valid month in the pixel-year (includes non-vegetated pixels and water) |

**nvalid** gives, per pixel per year, the number of months with a valid retrieval (1–12; pixel-years
with no valid month are masked, stored 0). In the released files, LAI = −32768, QA = 255, and nvalid = 0
mark exactly the same pixels.

### Exporting to GeoTIFF
Earth Engine GeoTIFF downloads write masked pixels with their stored values, so a mask added with
`updateMask` is lost unless the image is unmasked, and 8-bit downloads declare NoData = 0, which would
hide QA code 0. Unmask to the NoData value, declare it, and export on the native grid (no resampling):

```javascript
var LAI = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m')
  .filter(ee.Filter.eq('year', 2010)).filter(ee.Filter.eq('month', 7)).first();
var QA = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m_QA')
  .filter(ee.Filter.eq('year', 2010)).filter(ee.Filter.eq('month', 7)).first();
var grid = LAI.projection().getInfo();   // native grid: crs = grid.wkt, crsTransform = grid.transform

// gap-filled LAI (int16, LAI = value * 0.001; -32768 = no value)
Export.image.toDrive({image: LAI.unmask(-32768), description: 'LAI_2010_Jul', region: geometry,
  crs: grid.wkt, crsTransform: grid.transform, maxPixels: 1e13,
  formatOptions: {cloudOptimized: true, noData: -32768}});
// retrieved LAI only (-32768 where QA != 0)
Export.image.toDrive({image: LAI.updateMask(QA.eq(0)).unmask(-32768), description: 'LAIretrieved_2010_Jul',
  region: geometry, crs: grid.wkt, crsTransform: grid.transform, maxPixels: 1e13,
  formatOptions: {cloudOptimized: true, noData: -32768}});
// QA (uint8; 255 = not filled)
Export.image.toDrive({image: QA.unmask(255), description: 'QA_2010_Jul', region: geometry,
  crs: grid.wkt, crsTransform: grid.transform, maxPixels: 1e13,
  formatOptions: {cloudOptimized: true, noData: 255}});
```

## Repository Structure

```
CONUS-30m-LAI-mapping/
├── gee/                                  # Google Earth Engine (retrieval, re-export, validation refs)
│   ├── lai_state_monthly_export.js       # main LAI retrieval, exported per state × year
│   ├── lai_state_JUNE2022_export.js      # re-export of June 2022 (missing from the first export)
│   ├── obs_count_export.js               # observation-density statistics and 12 multi-year mean maps (500 m)
│   ├── obs_count_analysis.js             # observation-density analysis
│   ├── export_MODIS_VIIRS_median.js      # monthly-median MODIS/VIIRS LAI at 500 m (validation reference)
│   ├── export_PFT_yearly.js              # year-matched NLCD biome maps at 500 m (validation strata)
│   ├── conus_lai_viewer.js               # source of the public Earth Engine viewer app
│   └── set_asset_metadata.py             # descriptions / QA codes / DOI on the GEE collections
├── postprocessing/                       # mosaic state tiles onto the CONUS grid
│   ├── merge_lai_to_ref_grid.py          # warp + merge state tiles → CONUS monthly mosaics
│   └── build_june2022_conus.py           # reconcile + mosaic the re-exported June 2022
├── gapfill/                              # gap-filling pipeline
│   ├── GAPFILL_PLAN.md                   # method description
│   └── gapfill_lai.py                    # climatology-anchored temporal gap-fill + QA + nvalid
├── validation/                           # field, 20 m map, withholding, and MODIS checks
│   ├── gbov_download.py                  # download GBOV reference data (NEON sites)
│   ├── gbov_validate.py                  # product vs GBOV RM-7 field PAI (plot-level)
│   ├── gbov_lp3_validate.py              # product vs GBOV LP-3 Sentinel-2 20 m maps
│   ├── gapfill_holdout.py                # withholding test of the gap filling, by QA code
│   ├── qa_frequency.py                   # frequency of QA codes in the released product
│   ├── agg_validate.py                   # scale-matched 30 m → 500 m comparison with MODIS (R1)
│   ├── agg_validate_r2.py                # same, plus common cell set and valid-area sensitivity
│   ├── figures_validation.py             # validation figures
│   ├── validate_local.py                 # PFT × QA-stratified sampling (retrieved / filled / overall)
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
5. Valid pixels are composited by **monthly median** and exported per state as **uint16 (LAI × 100)** with a companion obs_count band. The CONUS mosaics are re-encoded to **int16 (LAI × 1000)** during merging (`postprocessing/`).

Each run processes one state × one year → 12 monthly GeoTIFFs.

### 2 — CONUS mosaicking (Python) — `postprocessing/merge_lai_to_ref_grid.py`
State tiles (EPSG:5070, 30 m) are warped (nearest-neighbor) onto the CONUS reference grid (see
Projection / grid above) and merged (valid-wins) into 12 monthly mosaics per year; pixels without a
valid retrieval are −32768. June 2022 was re-exported and reconciled separately
via `build_june2022_conus.py`.

### 3 — Gap-filling (Python) — `gapfill/gapfill_lai.py`
Gaps are filled with a **climatology-anchored temporal** method (Kandasamy et al., 2013), with
values > 0 treated as valid:
- **Climatology:** for each pixel and calendar month, the mean of the valid values across the 23
  years. Where fewer than three years are valid, the pixel climatology is replaced by a coarse
  spatial climatology: the mean climatology of the 160 × 160-pixel block (~4.8 km) containing the
  pixel, with empty blocks filled by interpolation from surrounding blocks. No land-cover
  constraint is applied.
- **Filling, within each pixel-year with at least one valid month:** gaps between two valid months
  are linearly interpolated (QA 1); the remaining gaps (before the first or after the last valid
  month) are filled with the climatology times an amplitude factor, the ratio of the sum of valid
  values to the sum of the corresponding climatology values in that year, clamped to 0.3–3.0
  (QA 2, 3, or 4; see table). If no climatology were available after the fallback, the nearest valid
  month of the same year would be copied; this did not occur in the released product.
- **Retrieved values are preserved exactly**, and every value is labelled in the QA layer.
  Pixel-years without any valid month (including non-vegetated pixels and water) are not filled
  (QA 255).

See `gapfill/GAPFILL_PLAN.md`.

### 4 — Validation (Python) — `validation/`
R² is the squared Pearson correlation; RMSE and bias are in m² m⁻² (bias = product minus reference).

| Check | Reference | N | R² | RMSE | Bias |
|---|---|---|---|---|---|
| Field validation (`gbov_validate.py`) | GBOV RM-7 field PAI, 27 NEON sites, 2013–2022 | 6,195 plot-dates | 0.78 | 1.06 | +0.29 |
| Spatial detail (`gbov_lp3_validate.py`) | GBOV LP-3 Sentinel-2 20 m maps | 888,940 pixels (1,781 maps) | 0.79 | 1.10 | −0.06 |
| Consistency with MODIS (`agg_validate_r2.py`), retrieved | MOD15A2H monthly median, 500 m | 3,482,166 cell-months | 0.77 | 0.83 | +0.38 |
| Same cells, gap-filled | MOD15A2H monthly median, 500 m | 3,482,166 cell-months | 0.77 | 0.83 | +0.38 |

- MODIS also supplied the training targets (MCD15A3H), so the MODIS comparison measures consistency,
  not independent accuracy. Agreement is similar for 2000–2005, 2006–2018, and 2019–2022
  (retrieved R² 0.78, 0.77, 0.76).
- Gap filling (`gapfill_holdout.py`, `qa_frequency.py`): withheld retrievals reconstructed by the
  production routine give an approximate RMSE of 0.69 m² m⁻² across filled values (QA codes weighted
  by their frequency, short-gap result for QA 1), with negligible bias.
- Known limitations: retrieved LAI levels off near 5–6 m² m⁻² in dense canopies and cannot exceed
  6.9 m² m⁻², the largest MODIS LAI in the training samples; field errors are largest in evergreen
  needleleaf forests (RMSE 1.48 m² m⁻²); and climatology-based filling damps interannual anomalies
  (use QA 0 values for trend, disturbance, and extreme-event analyses).

## Requirements
- **GEE scripts:** a Google Earth Engine account.
- **Python:** `pip install rasterio numpy pandas pyproj requests matplotlib tqdm` (mosaicking, gap-fill,
  validation); `earthengine-api` for `gee/set_asset_metadata.py`; `earthengine-api` and `geedim` for `download/`.
- Input/output paths are defined at the top of each Python script — adjust to your environment before running.

## Products
- **GEE asset (ImageCollection):** `projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m` — 276 monthly gap-filled LAI images.
  - https://code.earthengine.google.com/?asset=projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m
- Companion **QA** and **nvalid** collections: `…/CONUS_Monthly_LAI_30m_QA`, `…/CONUS_Monthly_LAI_30m_nvalid`.
- **Retrieved LAI:** the gap-filled LAI where QA = 0 (−32768 elsewhere when exported); see "Two LAI
  products" and "Exporting to GeoTIFF" above.
- **Observation count:** `…/CONUS_Monthly_LAI_30m_obscount`, 276 monthly images of the number of
  cloud-free Landsat observations, on a 500 m grid (float, EPSG:5070).
- Dataset record (DOI): https://doi.org/10.5281/zenodo.22292385 (Zenodo; metadata and access information).

## Citation
If you use this dataset or code, please cite the data record:

> You, H., Kang, Y., & Chen, M. (2026). Monthly 30-m Leaf Area Index Dataset for the Contiguous United
> States from 2000 to 2022. Zenodo. https://doi.org/10.5281/zenodo.22292385

Data descriptor (in review): You, H., Kang, Y., & Chen, M. Monthly 30-m Leaf Area Index Dataset for the
Contiguous United States from 2000 to 2022. *Scientific Data*.

Underlying LAI algorithm:

> Kang, Y., Ozdogan, M., Gao, F., Anderson, M. C., White, W. A., Yang, Y., Yang, Y., & Erickson, T. A. (2021).
> A data-driven approach to estimate leaf area index for Landsat images over the contiguous US.
> *Remote Sensing of Environment*, 258, 112383. https://doi.org/10.1016/j.rse.2021.112383

Gap-filling method:

> Kandasamy, S., Baret, F., Verger, A., Neveux, P., & Weiss, M. (2013). A comparison of methods for smoothing
> and gap filling time series of remote sensing observations – application to MODIS LAI products.
> *Biogeosciences*, 10, 4055–4071.

## License
Code: [MIT License](LICENSE) · Dataset: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
