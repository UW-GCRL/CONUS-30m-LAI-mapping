"""
Hold-out evaluation of the gap-filling procedure, stratified by the QA code each fill receives.

Field data rarely fall in gap months (field visits are in the growing season), so the
climatology-based fills (QA 2-4) cannot be checked against field LAI. Here real retrieved
values are withheld, refilled with the production routine (gapfill.gapfill_lai.fill_block), and
compared with the withheld values.

Leakage control: for each target year the climatology is rebuilt from the other 22 years only,
reproducing build_climatology(): pixel mean where >= MIN_YEARS valid years, otherwise the mean
of the 160 x 160-pixel block (computed within the tile, aligned to the production block grid).

Scenarios applied to each pixel-year (retrieved values, LAI > 0):
  random25  withhold a random 25% of valid months (pixels with >= 4 valid months)
  keep2     keep 2 random valid months, withhold the rest (pixels with >= 6 valid months)
  keep1     keep 1 random valid month, withhold the rest (pixels with >= 6 valid months)

  python gapfill_holdout.py --run      # tiles x years -> holdout_samples.csv
  python gapfill_holdout.py --stats
"""
import os, sys, argparse
import numpy as np, pandas as pd, rasterio
from rasterio.windows import Window

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gapfill'))
from gapfill_lai import fill_block, MIN_YEARS, COARSE, NODATA          # production routine

RET_DIR = r"K:\Hangkai\CONUS_LAI"
OUT = r"K:\Hangkai\CONUS_LAI_validation\gapfill_holdout_samples.csv"
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
YEARS = list(range(2000, 2023))
TILE = 4 * COARSE            # 640 px, aligned to the 160-px fallback grid
N_TILES, YEARS_PER_TILE, SAMPLES_PER_CASE = 30, 4, 3000
rng = np.random.default_rng(42)


def read_stack(win):
    """(23, 12, h, w) int16 retrieved LAI x1000; gaps and -32768 -> 0."""
    st = np.zeros((len(YEARS), 12, win.height, win.width), np.int16)
    for yi, y in enumerate(YEARS):
        for m in range(12):
            p = os.path.join(RET_DIR, f"{y}_{MON[m]}_LAI.tif")
            if os.path.exists(p):
                with rasterio.open(p) as ds:
                    a = ds.read(1, window=win)
                st[yi, m] = np.where(a > 0, a, 0)
    return st


def climatology_excluding(st, yi):
    """rclim (12,h,w) float32 and csrc (12,h,w) uint8 built without year index yi."""
    v = st > 0
    S = np.where(v, st, 0).astype(np.float64).sum(0) - np.where(v[yi], st[yi], 0)
    C = v.sum(0) - v[yi]
    clim = np.where(C > 0, S / np.maximum(C, 1), np.nan)
    h, w = clim.shape[1:]
    coarse = np.full_like(clim, np.nan)
    for r0 in range(0, h, COARSE):
        for c0 in range(0, w, COARSE):
            blk = clim[:, r0:r0 + COARSE, c0:c0 + COARSE]
            with np.errstate(invalid='ignore'):
                mu = np.nanmean(blk.reshape(12, -1), axis=1)
            coarse[:, r0:r0 + COARSE, c0:c0 + COARSE] = mu[:, None, None]
    native = C >= MIN_YEARS
    rclim = np.where(native, clim, coarse).astype(np.float32)
    return rclim, np.where(native, 0, 1).astype(np.uint8)


def withhold(valid, scenario):
    """boolean (12,h,w) mask of valid months to withhold."""
    nv = valid.sum(0)
    u = rng.random(valid.shape)
    u[~valid] = np.inf
    rank = u.argsort(0).argsort(0)                 # random rank among valid months
    if scenario == 'random25':
        k = np.floor(0.25 * nv).astype(int)
        return valid & (rank < k[None]) & (nv[None] >= 4)
    keep = 2 if scenario == 'keep2' else 1
    return valid & (rank >= keep) & (nv[None] >= 6)


def run():
    with rasterio.open(os.path.join(RET_DIR, '2015_Jul_LAI.tif')) as ds:
        W, H = ds.width, ds.height
    rows, done = [], 0
    while done < N_TILES:
        c0 = int(rng.integers(0, (W - TILE) // COARSE)) * COARSE
        r0 = int(rng.integers(0, (H - TILE) // COARSE)) * COARSE
        win = Window(c0, r0, TILE, TILE)
        with rasterio.open(os.path.join(RET_DIR, '2015_Jul_LAI.tif')) as ds:
            if (ds.read(1, window=win) > 0).mean() < 0.3:
                continue                                      # mostly water / non-vegetated
        st = read_stack(win)
        done += 1
        for yi in rng.choice(len(YEARS), YEARS_PER_TILE, replace=False):
            rclim, csrc = climatology_excluding(st, yi)
            truth = st[yi]
            valid = truth > 0
            for sc in ['random25', 'keep2', 'keep1']:
                held = withhold(valid, sc)
                if not held.any():
                    continue
                stack = np.where(held, 0, truth).astype(np.int16)
                gf, qa, _ = fill_block(stack, rclim, csrc)
                idx = np.flatnonzero(held)
                idx = rng.choice(idx, min(SAMPLES_PER_CASE, idx.size), replace=False)
                m, rr, cc = np.unravel_index(idx, held.shape)
                g = gf.reshape(-1)[idx].astype(float)
                rows.append(pd.DataFrame(dict(
                    tile=done, year=YEARS[yi], scenario=sc, month=m + 1,
                    truth=truth.reshape(-1)[idx] / 1000.0,
                    filled=np.where(g == NODATA, np.nan, g / 1000.0),
                    qa=qa.reshape(-1)[idx], n_kept=(stack > 0).sum(0)[rr, cc])))
        print(f"tile {done}/{N_TILES} ({c0},{r0}) done", flush=True)
    pd.concat(rows).to_csv(OUT, index=False)
    print(f"wrote {OUT}")


def metrics(t, f):
    m = np.isfinite(t) & np.isfinite(f); t, f = t[m], f[m]
    return dict(N=len(t), R2_pearson=np.corrcoef(t, f)[0, 1] ** 2 if len(t) > 2 else np.nan,
                R2_1to1=1 - np.sum((f - t) ** 2) / np.sum((t - t.mean()) ** 2),
                RMSE=np.sqrt(np.mean((f - t) ** 2)), MAE=np.mean(np.abs(f - t)),
                bias=np.mean(f - t))


def stats():
    d = pd.read_csv(OUT)
    pd.set_option('display.width', 200)
    print("== By QA code assigned to the refilled value")
    print(pd.DataFrame([dict(qa=int(q), **metrics(g['truth'].values, g['filled'].values))
                        for q, g in d.groupby('qa')]).round(3).to_string())
    print("\n== By scenario and QA code")
    print(pd.DataFrame([dict(scenario=s, qa=int(q), **metrics(g['truth'].values, g['filled'].values))
                        for (s, q), g in d.groupby(['scenario', 'qa'])]).round(3).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--stats', action='store_true')
    a = ap.parse_args()
    if a.run:
        run()
    if a.stats:
        stats()
