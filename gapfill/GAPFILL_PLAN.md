# CONUS Monthly LAI — Gap-Filling Plan (robust / everywhere-applicable)

## Goal
Per pixel, within each year: if the pixel has **≥1 valid** monthly LAI that year, fill the
missing months of that year. Preserve original valid values exactly. Emit a **confidence-tiered
QA layer** so heavy-gap fills are visible, not hidden. Method must degrade gracefully from
data-rich to chronically-gappy regions.

## Inputs
- 276 monthly GeoTIFFs `YYYY_Mon_LAI.tif` (2000–2022 × 12), int16, nodata −32768,
  160169 × 105432, EPSG:5070, in `K:\Hangkai\CONUS_LAI`.
- `2022_Jun` currently empty at source → being regenerated (see JUNE2022 script). Gap-fill
  treats it as nodata until the corrected file is dropped in.

## Method — robust climatology-anchored temporal infill (preserves observations)
Smoothers (Savitzky–Golay / Whittaker) rejected: they overwrite real values. Interpolation/
infill touches only gaps. Climatology is the literature-recommended backbone for long/dense
gaps (Kandasamy et al. 2013).

### Stage 1 — robust per-pixel climatology (per calendar month)
1. `clim_m`  = mean of month *m* across all years, ignoring nodata; `cnt_m` = # valid years.
2. `coarse_m` = block-mean of `clim_m` to a coarse grid (~5 km); coarse holes filled by
   nearest — a dense, low-res seasonal field that exists everywhere.
3. **Robust climatology** `rclim_m` = `clim_m` where `cnt_m ≥ MIN_YEARS` (default 3),
   else the upsampled `coarse_m` (spatial fallback). `csrc_m` records native(0)/fallback(1).

### Stage 2 — per pixel × year fill (12 monthly values)
- `n_valid` = # valid months that pixel-year.
- **n_valid = 0** → leave all nodata, QA = `255` (not filled).
- **n_valid ≥ 1** → keep valid months (QA `0`); for each missing month:
  - **Interior gap** (valid months both sides within the year): **linear interpolation**
    between them. QA = `1` (high confidence).
  - **Edge / exterior gap** (or n_valid=1): **climatology `rclim_m` × robust scale**, where
    `scale = median(obs/rclim over valid months)`, clamped to [0.3, 3] (one outlier can't blow
    up the fill). If `rclim_m` unavailable, hold nearest valid month.
    QA = `2` if n_valid ≥ 3; `3` if n_valid ∈ {1,2}; **`4`** if the climatology used came from
    the spatial fallback (`csrc_m = 1`) — lowest confidence.

## Outputs (`K:\Hangkai\CONUS_LAI_gapfilled\`)
- `YYYY_Mon_LAI_gf.tif` — gap-filled LAI, int16, same grid/scale/nodata as inputs.
- `YYYY_Mon_LAI_qa.tif` — uint8 status/confidence:
  `0` original · `1` linear-interp · `2` climatology (year ≥3 valid) · `3` climatology (year 1–2 valid)
  · `4` spatial-fallback climatology · `255` not filled (no valid that pixel-year).
  → filled = {1,2,3,4}; original = 0; empty = 255; and 2→3→4 = decreasing confidence.
- `YYYY_nvalid.tif` — uint8, # valid months that pixel-year (0–12): direct confidence measure.
- `_climatology\` — intermediate `rclim_Mon.tif`, `cnt_Mon.tif`, `csrc_Mon.tif`.

## Robustness guarantees
- Every pixel-year with ≥1 observation is filled; the fill always has a defined seasonal
  basis (own climatology, else spatial fallback), so no region is left ragged.
- Original observations are never altered.
- QA + n_valid expose exactly where/how each value was produced → you can trust, down-weight,
  or mask heavy-gap fills.

## Tunables
`MIN_YEARS=3` (climatology trust threshold), `COARSE≈160` px (~4.8 km fallback grid),
`SCALE_CLAMP=(0.3,3)`, `TILE=2048`.

## Scale / cost
~2.7 TB read (climatology) + ~2.7 TB read + ~2.7 TB write (fill) + small QA/nvalid.
~8–9 TB I/O, ~2.7 TB new disk (K: has ~19 TB free). Many hours to ~2 days, tiled/parallelizable,
resumable. **Prototype with `--tile-test` before the full run.**
