"""
Local PFT x QA-stratified validation of retrieved / gap-filled / overall LAI against
monthly-MEDIAN MODIS (MOD15A2H). MODIS only (VIIRS skipped).

Inputs already on disk:
  - retrieved LAI : SRC_DIR\{YYYY}_{Mon}_LAI.tif        (gaps = 0, valid LAI > 0)
  - gap-filled    : GF_DIR\{YYYY}_{Mon}_LAI_gf.tif      (valid/filled > 0, unfilled = -32768)
  - QA            : GF_DIR\{YYYY}_{Mon}_LAI_qa.tif       (0 orig, 1-4 filled, 255 unfilled)
  - MODIS median  : VAL_DIR\MODISmed_{YYYY}.tif          (500 m, 12 monthly bands, x1000, -32768 = gap)
  - year PFT      : VAL_DIR\PFT_{YYYY}.tif               (500 m, 1-8 biomes, 0 = other; year-matched)

Method: random 30 m windows/month; within each, classify pixels by QA and pair each with its
co-located 500 m MODIS-median cell + PFT. Report accuracy per PFT for:
   RETRIEVED (QA=0, value=retrieved) | GAP-FILLED (QA 1-4, value=gf) | OVERALL (QA 0-4, value=gf)
Run: python validate_local.py
"""
import os, csv
import numpy as np
import rasterio
from rasterio.windows import Window

SRC_DIR = r"K:\Hangkai\CONUS_LAI"
GF_DIR  = r"K:\Hangkai\CONUS_LAI_gapfilled"
VAL_DIR = r"K:\Hangkai\CONUS_LAI_gapfilled\MODIS"     # MODISmed_YYYY.tif + PFT_YYYY.tif
OUT_CSV = r"K:\Hangkai\CONUS_LAI\lai_val_samples.csv"
MON     = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
YEARS   = range(2000, 2023)
NODATA  = -32768
SCALE   = 0.001                       # int16 x1000 -> LAI (product & MODIS both)
WIN     = 256                         # 30 m sampling window edge
WINS_PER_MONTH = 60
RET_PER_WIN = 8                       # retrieved samples kept per window (abundant -> cap)
PFT_NAMES = {1:'decid',2:'evergreen',3:'mixed',4:'shrub',5:'grass',6:'crop',7:'woodywet',8:'herbwet'}

rng = np.random.default_rng(0)
def period(y): return '2000-2005' if y <= 2005 else ('2019-2022' if y >= 2019 else '2006-2018')
def read_win(ds, win, band=1): return ds.read(band, window=win, boundless=True, fill_value=NODATA)

def coarse_at(ds, xs, ys, band=1):
    """value from a 500 m raster at arrays of EPSG:5070 x,y (nearest cell), given band."""
    inv = ~ds.transform
    cc = (inv * (xs, ys))
    col = np.floor(cc[0]).astype(int); row = np.floor(cc[1]).astype(int)
    ok = (row >= 0) & (row < ds.height) & (col >= 0) & (col < ds.width)
    out = np.full(xs.shape, NODATA, np.int16)
    if ok.any():
        r0, r1 = row[ok].min(), row[ok].max()+1; c0, c1 = col[ok].min(), col[ok].max()+1
        block = ds.read(band, window=Window(c0, r0, c1-c0, r1-r0))
        out[ok] = block[row[ok]-r0, col[ok]-c0]
    return out

def main():
    rows = []
    with rasterio.open(os.path.join(SRC_DIR, "2015_Jul_LAI.tif")) as ref:
        W, H, tr = ref.width, ref.height, ref.transform
    for y in YEARS:
        mpath = os.path.join(VAL_DIR, f"MODISmed_{y}.tif")
        ppath = os.path.join(VAL_DIR, f"PFT_{y}.tif")
        if not (os.path.exists(mpath) and os.path.exists(ppath)):
            print(f"{y}: MODIS/PFT missing, skip"); continue
        mds = rasterio.open(mpath); pds = rasterio.open(ppath)
        for mi, mon in enumerate(MON):
            rp = os.path.join(SRC_DIR, f"{y}_{mon}_LAI.tif")
            gp = os.path.join(GF_DIR,  f"{y}_{mon}_LAI_gf.tif")
            qp = os.path.join(GF_DIR,  f"{y}_{mon}_LAI_qa.tif")
            if not (os.path.exists(rp) and os.path.exists(gp) and os.path.exists(qp)): continue
            rds, gds, qds = rasterio.open(rp), rasterio.open(gp), rasterio.open(qp)
            for _ in range(WINS_PER_MONTH):
                c0 = int(rng.integers(0, W - WIN)); r0 = int(rng.integers(0, H - WIN))
                win = Window(c0, r0, WIN, WIN)
                ret = read_win(rds, win); gf = read_win(gds, win); qa = read_win(qds, win)
                jj, ii = np.meshgrid(np.arange(WIN), np.arange(WIN))
                xs = (tr.c + (c0 + jj + 0.5) * tr.a).ravel()
                ys = (tr.f + (r0 + ii + 0.5) * tr.e).ravel()
                mod = coarse_at(mds, xs, ys, band=mi+1).reshape(WIN, WIN)   # MODIS month band
                pf  = coarse_at(pds, xs, ys, band=1).reshape(WIN, WIN)
                valid_mod = mod != NODATA
                fm = (qa >= 1) & (qa <= 4) & (gf > 0)  & valid_mod          # filled (rare) -> keep all
                rm = (qa == 0)            & (ret > 0) & valid_mod           # retrieved -> subsample
                idx_r = np.argwhere(rm)
                if len(idx_r) > RET_PER_WIN:
                    idx_r = idx_r[rng.choice(len(idx_r), RET_PER_WIN, replace=False)]
                idx_f = np.argwhere(fm)
                for (a, b) in list(idx_r) + list(idx_f):
                    rows.append((y, mi+1, period(y), int(pf[a,b]), int(qa[a,b]),
                                 int(ret[a,b]), int(gf[a,b]), int(mod[a,b])))
            rds.close(); gds.close(); qds.close()
        mds.close(); pds.close()
        print(f"{y}: {len(rows)} samples so far", flush=True)

    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['year','month','period','pft','qa','retrieved','overall','modis'])
        w.writerows(rows)
    print("wrote", OUT_CSV, "n=", len(rows))
    report(np.array(rows, dtype=object))

def metrics(a, b):
    a = a.astype(float); b = b.astype(float)
    if len(a) < 5: return (len(a), np.nan, np.nan, np.nan, np.nan)
    rmse = np.sqrt(np.mean((a-b)**2)); mae = np.mean(np.abs(a-b)); bias = np.mean(a-b)
    r2 = np.corrcoef(a, b)[0,1]**2 if a.std()>0 and b.std()>0 else np.nan
    return (len(a), r2, rmse, mae, bias)

def report(R):
    if len(R) == 0: print("no samples"); return
    y,mo,per,pft,qa,ret,ov,mod = [R[:,k] for k in range(8)]
    ret=ret.astype(float)*SCALE; ov=ov.astype(float)*SCALE; mod=mod.astype(float)*SCALE
    qa=qa.astype(int); pft=pft.astype(int)
    def show(name, mask, ours):
        print(f"\n=== {name} vs MODIS-median (LAI units) ===")
        print(f"{'PFT':<10}{'N':>8}{'R2':>7}{'RMSE':>7}{'MAE':>7}{'bias':>7}")
        for p in range(1,9):
            m = mask & (pft==p) & (mod>=0)
            n,r2,rm,ma,bi = metrics(ours[m], mod[m])
            print(f"{PFT_NAMES[p]:<10}{n:>8}{r2:>7.2f}{rm:>7.2f}{ma:>7.2f}{bi:>7.2f}")
        n,r2,rm,ma,bi = metrics(ours[mask & (mod>=0)], mod[mask & (mod>=0)])
        print(f"{'ALL':<10}{n:>8}{r2:>7.2f}{rm:>7.2f}{ma:>7.2f}{bi:>7.2f}")
    show("RETRIEVED (QA=0)",         qa==0,                 ret)
    show("GAP-FILLED (QA 1-4)",      (qa>=1)&(qa<=4),       ov)
    show("OVERALL (retrieved+fill)", (qa>=0)&(qa<=4),       ov)

if __name__ == "__main__":
    main()
