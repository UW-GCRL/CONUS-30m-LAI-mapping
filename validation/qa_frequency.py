"""
Frequency of gap-fill QA codes in the released gap-filled product (Reviewer 2, Major 2).

Sample: N_WIN random 512 x 512 windows across the CONUS grid (same windows for every file).
Population: "land pixels with retrieval history" = pixels whose climatology count (cnt_<Mon>)
is > 0 in at least one calendar month, i.e. pixels with at least one valid retrieval in
2000-2022. Water and never-vegetated pixels are excluded.

Reported per calendar month and overall:
  * share of pixel-months by QA code (0, 1, 2, 3, 4, 255);
  * QA 1 split by the pixel-year's nvalid (<= 2: long interpolation spans; >= 3: other);
  * share of pixel-months whose climatology is undefined after the spatial fallback (rclim NaN),
    where the production code copies the nearest valid month instead.
Run: python qa_frequency.py   (writes qa_frequency.csv)
"""
import os, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window

QA_DIR = r"K:\Hangkai\CONUS_LAI_gapfilled\_up_QA"
NV_DIR = r"K:\Hangkai\CONUS_LAI_gapfilled\_up_nvalid"
CLIM = r"K:\Hangkai\CONUS_LAI_gapfilled\_climatology"
OUT = r"K:\Hangkai\CONUS_LAI_validation\qa_frequency.csv"
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
YEARS = range(2000, 2023); N_WIN = 120; WS = 512
rng = np.random.default_rng(7)


def main():
    with rasterio.open(os.path.join(CLIM, "cnt_Jul.tif")) as ds:
        W, H = ds.width, ds.height
    wins = [Window(int(rng.integers(0, W - WS)), int(rng.integers(0, H - WS)), WS, WS) for _ in range(N_WIN)]
    # population mask + climatology-undefined mask per month
    hist = [np.zeros((WS, WS), bool) for _ in wins]
    nan_clim = {m: [] for m in range(12)}
    for m, mon in enumerate(MON):
        with rasterio.open(os.path.join(CLIM, f"cnt_{mon}.tif")) as c, \
             rasterio.open(os.path.join(CLIM, f"rclim_{mon}.tif")) as r:
            for i, w in enumerate(wins):
                hist[i] |= c.read(1, window=w) > 0
                nan_clim[m].append(~np.isfinite(r.read(1, window=w)))
        print(f"climatology {mon} read", flush=True)
    npop = sum(h.sum() for h in hist)
    print(f"population pixels in sample: {npop:,}", flush=True)
    rows = []
    for y in YEARS:
        with rasterio.open(os.path.join(NV_DIR, f"{y}_nvalid.tif")) as d:
            nv = [d.read(1, window=w) for w in wins]
        for m, mon in enumerate(MON):
            p = os.path.join(QA_DIR, f"{y}_{mon}_QA.tif")
            if not os.path.exists(p):
                continue
            cnt = {k: 0 for k in ['qa0', 'qa1', 'qa1_nv_le2', 'qa2', 'qa3', 'qa4', 'qa255', 'other', 'clim_nan_fill']}
            with rasterio.open(p) as d:
                for i, w in enumerate(wins):
                    q = d.read(1, window=w)[hist[i]]
                    v = nv[i][hist[i]]
                    for k in [0, 1, 2, 3, 4, 255]:
                        cnt[f'qa{k}'] += int((q == k).sum())
                    cnt['other'] += int((~np.isin(q, [0, 1, 2, 3, 4, 255])).sum())
                    cnt['qa1_nv_le2'] += int(((q == 1) & (v <= 2)).sum())
                    cnt['clim_nan_fill'] += int((np.isin(q, [2, 3, 4]) & nan_clim[m][i][hist[i]]).sum())
            rows.append(dict(year=y, month=m + 1, pop=int(npop), **cnt))
        print(f"{y} done", flush=True)
    df = pd.DataFrame(rows); df.to_csv(OUT, index=False)
    tot = df.drop(columns=['year', 'month']).sum()
    print("\n== Share of land pixel-months (%), 2000-2022")
    for k in ['qa0', 'qa1', 'qa1_nv_le2', 'qa2', 'qa3', 'qa4', 'qa255', 'other', 'clim_nan_fill']:
        print(f"  {k:<14}{100 * tot[k] / tot['pop']:8.3f}")
    filled = tot[['qa1', 'qa2', 'qa3', 'qa4']].sum()
    print(f"\n== Share of FILLED pixel-months (%): " + ", ".join(
        f"{k} {100 * tot[k] / filled:.1f}" for k in ['qa1', 'qa1_nv_le2', 'qa2', 'qa3', 'qa4', 'clim_nan_fill']))
    mm = df.groupby('month').sum()
    print("\n== By calendar month (% of land pixel-months)")
    print((100 * mm[['qa0', 'qa1', 'qa2', 'qa3', 'qa4', 'qa255']].div(mm['pop'], axis=0)).round(2).to_string())


if __name__ == "__main__":
    main()
