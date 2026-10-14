# Gap filling of the CONUS monthly LAI (method as implemented in `gapfill_lai.py`)

## Input
The 276 retrieved monthly mosaics (2000–2022, int16, LAI × 1000) on the CONUS grid. Values > 0 are
valid retrievals; any value ≤ 0 (0 or −32768) is treated as missing.

## Stage 1: climatology (`--climatology`)
For each calendar month:
1. Per-pixel climatology = mean of the valid values across the 23 years, with the number of valid
   years (`cnt`).
2. Coarse climatology = mean of the per-pixel climatology over 160 × 160-pixel blocks (~4.8 km);
   empty blocks are filled by interpolation from surrounding blocks (`rasterio.fill.fillnodata`).
3. Climatology used for filling = the per-pixel climatology where `cnt` ≥ 3 (`MIN_YEARS`), otherwise
   the coarse climatology of the block containing the pixel (spatial fallback).

## Stage 2: filling (`--fill`), per pixel and year
- No valid month in the year: not filled (QA 255, LAI −32768).
- Otherwise valid months are kept unchanged (QA 0), and each missing month is filled:
  - between two valid months of the same year: linear interpolation (QA 1);
  - before the first or after the last valid month: climatology × amplitude factor, where the factor
    is the sum of the valid values divided by the sum of the corresponding climatology values in
    that year, clamped to 0.3–3.0. QA 2 if the year has ≥ 3 valid months, QA 3 if it has 1–2, and
    QA 4 if the climatology came from the spatial fallback. If no climatology were available, the
    nearest valid month of the year would be copied; this did not occur in the released product.

## Outputs
Per month: gap-filled LAI (int16, LAI × 1000, −32768 = not filled) and QA (uint8: 0, 1, 2, 3, 4,
255). Per year: nvalid (uint8, number of valid months, 0–12).

## Parameters
`MIN_YEARS = 3`, `COARSE = 160` pixels, `SCALE_CLAMP = (0.3, 3.0)`, `TILE = 1024` (processing tile).
