"""
Scale-matched MODIS consistency check, revised for Scientific Data round 2 (Reviewer 2, Minor 4).
Differences from agg_validate.py (R1):
  * minimum valid area is a FRACTION of the 30 m pixels falling in each 500 m cell, evaluated at
    0.5 (main) and 50 px ~ 0.18 (R1 setting, sensitivity);
  * retrieved and gap-filled products are also compared on a COMMON cell set (cells passing the
    threshold for both);
  * per-PFT sample sizes are reported for every case;
  * file locations updated to the released gap-filled and QA files.
Aggregation (unchanged): mean of 30 m pixels whose centres fall in the 500 m cell (exact
coordinate binning); only LAI > 0 contributes (non-vegetated and unobserved pixels excluded);
retrieved = QA 0 pixels; gap-filled = all filled or retrieved pixels. R2 = squared Pearson r.
Run: python agg_validate_r2.py   (writes agg_validate_r2_results.csv)
"""
import os, numpy as np, pandas as pd, rasterio, collections
from rasterio.windows import Window

SRC = r"K:\Hangkai\CONUS_LAI"
GF_DIRS = [r"K:\Hangkai\CONUS_LAI_gapfilled\_up_LAI", r"K:\Hangkai\CONUS_LAI_gapfilled\_up_LAI_done"]
QA_DIR = r"K:\Hangkai\CONUS_LAI_gapfilled\_up_QA"
VAL = r"K:\Hangkai\CONUS_LAI_gapfilled\MODIS"
OUT = r"K:\Hangkai\CONUS_LAI_validation\agg_validate_r2_results.csv"
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
YEARS = range(2000, 2023); SCALE = 0.001
WIN5 = 48; WINS_PER_MONTH = 12
THRESH = {'frac50': ('frac', 0.5), 'px50_R1': ('px', 50)}
PFT = {1:'decid',2:'evergreen',3:'mixed',4:'shrub',5:'grass',6:'crop',7:'woodywet',8:'herbwet'}
rng = np.random.default_rng(0)
def period(y): return '2000-2005' if y <= 2005 else ('2019-2022' if y >= 2019 else '2006-2018')
acc = collections.defaultdict(lambda: np.zeros(7))
def add(key, a, b):
    if a.size == 0: return
    v = acc[key]; v[0]+=a.size; v[1]+=a.sum(); v[2]+=b.sum(); v[3]+=(a*a).sum(); v[4]+=(b*b).sum(); v[5]+=(a*b).sum(); v[6]+=np.abs(a-b).sum()

with rasterio.open(os.path.join(SRC, "2015_Jul_LAI.tif")) as r: tr30 = r.transform; W30 = r.width; H30 = r.height

def find(dirs, name):
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.exists(p): return p

def cellstats(gd, qd, sd, mds, c5, r5, w5, h5):
    """per 500 m cell: sums and counts for gap-filled and retrieved, plus total 30 m pixels."""
    t5 = mds.transform
    xs = [t5.c + c5*t5.a, t5.c + (c5+w5)*t5.a]; ys = [t5.f + r5*t5.e, t5.f + (r5+h5)*t5.e]
    inv30 = ~tr30
    cols = [(inv30*(x, y))[0] for x in xs for y in ys]; rows = [(inv30*(x, y))[1] for x in xs for y in ys]
    cc0 = max(0, int(np.floor(min(cols)))); cc1 = min(W30, int(np.ceil(max(cols))))
    rr0 = max(0, int(np.floor(min(rows)))); rr1 = min(H30, int(np.ceil(max(rows))))
    if cc1 <= cc0 or rr1 <= rr0: return None
    w30 = Window(cc0, rr0, cc1-cc0, rr1-rr0)
    gf = gd.read(1, window=w30); qa = qd.read(1, window=w30); ret = sd.read(1, window=w30)
    jj, ii = np.meshgrid(np.arange(cc1-cc0), np.arange(rr1-rr0))
    x = tr30.c + (cc0+jj+0.5)*tr30.a; y = tr30.f + (rr0+ii+0.5)*tr30.e
    inv5 = ~t5
    c5f = np.floor(inv5.a*x + inv5.b*y + inv5.c).astype(int) - c5
    r5f = np.floor(inv5.d*x + inv5.e*y + inv5.f).astype(int) - r5
    inside = (c5f >= 0) & (c5f < w5) & (r5f >= 0) & (r5f < h5)
    lin = (r5f*w5 + c5f); N = w5*h5
    def sc(mask, val):
        m = inside & mask
        return (np.bincount(lin[m], weights=val[m].astype(float), minlength=N),
                np.bincount(lin[m], minlength=N))
    s_gf, c_gf = sc(gf > 0, gf)
    s_rt, c_rt = sc((qa == 0) & (ret > 0), ret)
    tot = np.bincount(lin[inside], minlength=N)
    sh = (h5, w5)
    return s_gf.reshape(sh), c_gf.reshape(sh), s_rt.reshape(sh), c_rt.reshape(sh), tot.reshape(sh)

def main():
    for y in YEARS:
        mp = os.path.join(VAL, f"MODISmed_{y}.tif"); pp = os.path.join(VAL, f"PFT_{y}.tif")
        if not (os.path.exists(mp) and os.path.exists(pp)): print(f"{y}: skip"); continue
        mds = rasterio.open(mp); pds = rasterio.open(pp); W5 = mds.width; H5 = mds.height
        for mi, mon in enumerate(MON):
            rp = os.path.join(SRC, f"{y}_{mon}_LAI.tif"); gp = find(GF_DIRS, f"{y}_{mon}_LAI.tif")
            qp = os.path.join(QA_DIR, f"{y}_{mon}_QA.tif")
            if not (os.path.exists(rp) and gp and os.path.exists(qp)): print(f"  {y} {mon}: missing"); continue
            sd = rasterio.open(rp); gd = rasterio.open(gp); qd = rasterio.open(qp)
            for _ in range(WINS_PER_MONTH):
                c5 = int(rng.integers(0, W5-WIN5)); r5 = int(rng.integers(0, H5-WIN5))
                res = cellstats(gd, qd, sd, mds, c5, r5, WIN5, WIN5)
                if res is None: continue
                s_gf, c_gf, s_rt, c_rt, tot = res
                mod = mds.read(mi+1, window=Window(c5, r5, WIN5, WIN5)).astype(float)*SCALE
                pf = pds.read(1, window=Window(c5, r5, WIN5, WIN5))
                mgf = s_gf/np.maximum(c_gf, 1)*SCALE; mrt = s_rt/np.maximum(c_rt, 1)*SCALE
                for tn, (kind, val) in THRESH.items():
                    need = val*tot if kind == 'frac' else np.full(tot.shape, val)
                    okg = (c_gf >= need) & (c_gf > 0); okr = (c_rt >= need) & (c_rt > 0)
                    for k in PFT:
                        kk = (pf == k) & (mod >= 0)
                        for cls, m, est in [('F', kk & okg, mgf), ('R', kk & okr, mrt),
                                            ('Fc', kk & okg & okr, mgf), ('Rc', kk & okg & okr, mrt)]:
                            add((tn, period(y), k, cls), est[m], mod[m])
            sd.close(); gd.close(); qd.close()
        mds.close(); pds.close()
        print(f"{y} done", flush=True)
    report()

def stats(v):
    n, Sa, Sb, Saa, Sbb, Sab, Sd = v
    if n < 20: return dict(N=int(n))
    num = (n*Sab - Sa*Sb)**2; den = (n*Saa - Sa*Sa)*(n*Sbb - Sb*Sb)
    return dict(N=int(n), R2=num/den if den > 0 else np.nan,
                RMSE=np.sqrt(max(0, (Saa - 2*Sab + Sbb)/n)), MAE=Sd/n, bias=(Sa - Sb)/n)

def report():
    rows = []
    for t in THRESH:
        for c in ['R', 'F', 'Rc', 'Fc']:
            for k in list(PFT) + ['ALL']:
                for p in ['ALL', '2000-2005', '2006-2018', '2019-2022']:
                    tot = np.zeros(7)
                    for (tt, pp, kk, cc), v in acc.items():
                        if tt == t and cc == c and (k == 'ALL' or kk == k) and (p == 'ALL' or pp == p):
                            tot += v
                    rows.append(dict(threshold=t, product=c, pft=PFT.get(k, k), period=p, **stats(tot)))
    df = pd.DataFrame(rows); df.to_csv(OUT, index=False)
    pd.set_option('display.width', 200)
    print(df[(df.period == 'ALL')].round(3).to_string())
    print(df[(df.pft == 'ALL') & (df['product'].isin(['Rc', 'Fc']))].round(3).to_string())

if __name__ == "__main__": main()
