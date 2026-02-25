# CONUS-30m-LAI: Monthly 30-m Leaf Area Index for the Contiguous United States (2000–2022)

This repository contains the production code for generating the **CONUS-30m-LAI** dataset — a monthly, 30-meter resolution Leaf Area Index (LAI) product covering the contiguous United States (CONUS) from 2000 to 2022.

## Dataset Overview

| Attribute | Value |
|---|---|
| Variable | Leaf Area Index (LAI) |
| Spatial coverage | Contiguous United States (CONUS) |
| Spatial resolution | 30 m |
| Temporal coverage | 2000–2022 |
| Temporal resolution | Monthly |
| Projection | EPSG:5070 (Albers Equal Area Conic) |
| Data type | float32 (real LAI units, 0–8) |
| Format | Cloud-Optimized GeoTIFF |

## Repository Structure

```
CONUS-30m-LAI-mapping/
├── gee/
│   └── lai_state_monthly_export.js     # GEE script: exports monthly LAI per state
├── postprocessing/
│   └── merge_lai_to_ref_grid.py        # Python: merges state tiles into CONUS mosaics
├── LICENSE
└── README.md
```

## Methods Summary

### Step 1 — LAI retrieval (Google Earth Engine)
Monthly LAI composites are generated state-by-state using the **biome-stratified Random Forest** algorithm of [Kang et al. (2021)](https://doi.org/10.1016/j.rse.2021.112383):

1. Landsat C02 L2 surface reflectance (Landsat 5/7/8/9) is loaded, band-renamed, scaled, and cloud/shadow/water masked via `QA_PIXEL`.
2. Spectral indices (NDVI, NDWI) and solar geometry are computed.
3. Each pixel is assigned to one of 9 biome types using the year-matched [NLCD](https://www.usgs.gov/centers/eros/science/national-land-cover-database).
4. A separate Random Forest model (100 trees, trained against MODIS LAI) is applied per sensor × biome combination.
5. A QA band flags pixels with out-of-range inputs (bit 0), out-of-range LAI (bit 1), or non-vegetation biome (bit 2).
6. Valid pixels are composited (median) across the month and exported as uint16 (scale factor = 0.01).

Each GEE run processes one state × one year → 12 monthly GeoTIFFs.

### Step 2 — CONUS mosaicking (Python)
State-level tiles are warped to a common CONUS reference grid and merged into 12 monthly mosaics per year using `merge_lai_to_ref_grid.py`:

- Warping: nearest-neighbor resampling onto EPSG:5070 reference grid
- Merge policy: last-one-wins (valid pixels overwrite)
- Scale factor 0.01 applied to convert uint16 → real LAI (float32)
- Output: LZW-compressed, tiled GeoTIFF with overviews

## Requirements

### GEE script
- Google Earth Engine account
- Access to `projects/ee-yanghuikang/assets/LAI/` (public training assets)

### Python postprocessing
```bash
pip install rasterio numpy tqdm
```

## Citation

If you use this dataset or code, please cite:

> You, H., et al. (2025). A 30-m monthly leaf area index dataset for the contiguous United States from 2000 to 2022. *Scientific Data*. [DOI TBD]

The underlying LAI algorithm should also be cited:

> Kang, Y., Ozdogan, M., Gao, F., Anderson, M. C., White, W. A., Yang, Y., Yang, Y., & Erickson, T. A. (2021). A data-driven approach to estimate leaf area index for Landsat images over the contiguous US. *Remote Sensing of Environment*, 258, 112383. https://doi.org/10.1016/j.rse.2021.112383

## License

Code: [MIT License](LICENSE)
Dataset: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
