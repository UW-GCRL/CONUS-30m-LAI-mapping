"""
Pixel-level comparison of CONUS-30m-LAI with GBOV LP-3 high-resolution reference LAI maps.

Reference: current GBOV LP-3 versions (V3.3-V3.5.0, 20 m), all derived from Sentinel-2
(S2A/S2B/S2G), produced by calibrating Sentinel-2 SL2P LAI against the NEON DHP plot
measurements of each site. Earlier versions (V2.1-V3.2), largely Landsat-8-based and so sharing
our Landsat input, are not used (their Sentinel-2 maps are excluded too). Pixels flagged
input/output out of range are excluded. Maps cover about 5 km x 5 km.

For each site and calendar month, one LP-3 map is used: the date closest to mid-month whose
valid fraction is >= MIN_VALID. The 20 m map is averaged onto the 30 m EPSG:5070 product grid
(rasterio average resampling), then compared with the retrieved product (LAI > 0) of the same
month. Pooled metrics use up to PIX_PER_MAP random pixels per map; per-map metrics describe the
within-site spatial agreement (Pearson r and the ratio of spatial standard deviations).

  python gbov_lp3_validate.py --extract
  python gbov_lp3_validate.py --stats
"""
import os, re, glob, argparse
import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds, Window

LP3_GLOB = r"K:\Hangkai\CONUS_LAI_gapfilled\GBOV\COPERNICUS_GBOV_LP3_RM7_*\LP3\*.TIF"
RET_DIR = r"K:\Hangkai\CONUS_LAI"
OUT_DIR = r"K:\Hangkai\CONUS_LAI_validation"
PIX_OUT = os.path.join(OUT_DIR, "gbov_lp3_pixels.csv")
MAP_OUT = os.path.join(OUT_DIR, "gbov_lp3_maps.csv")
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
MIN_VALID, PIX_PER_MAP, SCALE = 0.5, 500, 0.001
rng = np.random.default_rng(0)


def inventory():
    rows = []
    for f in glob.glob(LP3_GLOB):
        b = os.path.basename(f).upper()
        m = re.match(r'GBOV_LP03_(S2[ABG])_(\w{4})_(\d{8})_\d{8}_\w+?_UOS_(V[\d.]+?)_20M_W[SG]{2}84\.TIF', b)
        if m:
            rows.append(dict(f=f, sensor=m.group(1), site=m.group(2),
                             date=pd.Timestamp(m.group(3)), ver=m.group(4)))
    d = pd.DataFrame(rows)
    d['vkey'] = d['ver'].str[1:].apply(lambda v: tuple(int(x) for x in v.split('.')))
    d = d.sort_values('vkey').groupby(['site', 'date']).tail(1)        # newest version per site-date
    d = d[d['date'].dt.year.between(2000, 2022)]
    d['year'], d['month'] = d['date'].dt.year, d['date'].dt.month
    d['dmid'] = (d['date'].dt.day - 15).abs()
    return d.sort_values(['site', 'year', 'month', 'dmid'])


def read_lp3(f):
    with rasterio.open(f) as s:
        a = s.read().astype('float64')
        lai = a[0]
        ok = np.isfinite(lai) & (lai >= 0) & (lai <= 10) & (a[1] == 0) & (a[2] == 0)
        crs = s.crs
        if crs is None and -180 <= s.bounds.left <= 180 and -90 <= s.bounds.bottom <= 90:
            crs = rasterio.crs.CRS.from_epsg(4326)                   # *_WGS84 files missing the tag
        return np.where(ok, lai, np.nan), (ok.mean() if crs else 0.0), s.transform, crs, s.bounds


def extract():
    inv = inventory()
    print(f"Sentinel-2 LP-3 maps (newest version per site-date): {len(inv)}")
    with rasterio.open(os.path.join(RET_DIR, '2015_Jul_LAI.tif')) as ds:
        ptr, pcrs = ds.transform, ds.crs
    pix, maps = [], []
    for (site, yr, mo), g in inv.groupby(['site', 'year', 'month']):
        for _, row in g.iterrows():                                     # closest to mid-month first
            try:
                lai, vf, tr, crs, bnd = read_lp3(row['f'])
            except Exception as e:
                print(f"  skip unreadable {os.path.basename(row['f'])}: {e}", flush=True)
                continue
            if vf >= MIN_VALID:
                break
        else:
            continue
        rp = os.path.join(RET_DIR, f"{yr}_{MON[mo - 1]}_LAI.tif")
        if not os.path.exists(rp):
            continue
        b5070 = transform_bounds(crs, pcrs, *bnd)
        with rasterio.open(rp) as ds:
            win = from_bounds(*b5070, transform=ptr).round_offsets().round_lengths()
            prod = ds.read(1, window=win).astype('float64')
            wtr = ds.window_transform(win)
        ref = np.full(prod.shape, np.nan)
        reproject(lai, ref, src_transform=tr, src_crs=crs, src_nodata=np.nan,
                  dst_transform=wtr, dst_crs=pcrs, dst_nodata=np.nan,
                  resampling=Resampling.average)
        ok = np.isfinite(ref) & (prod > 0)
        if ok.sum() < 100:
            continue
        p, r = prod[ok] * SCALE, ref[ok]
        maps.append(dict(site=site, year=yr, month=mo, date=row['date'].date(), sensor=row['sensor'],
                         ver=row['ver'], n=int(ok.sum()), ref_mean=r.mean(), prod_mean=p.mean(),
                         r=np.corrcoef(r, p)[0, 1], sd_ref=r.std(), sd_prod=p.std(),
                         rmse=np.sqrt(np.mean((p - r) ** 2)), bias=np.mean(p - r)))
        k = rng.choice(ok.sum(), min(PIX_PER_MAP, ok.sum()), replace=False)
        pix.append(pd.DataFrame(dict(site=site, year=yr, month=mo, ref=r[k], prod=p[k])))
        if len(maps) % 50 == 0:
            print(f"  {len(maps)} maps ({site} {yr}-{mo:02d})", flush=True)
    pd.DataFrame(maps).to_csv(MAP_OUT, index=False)
    pd.concat(pix).to_csv(PIX_OUT, index=False)
    print(f"wrote {MAP_OUT} ({len(maps)} maps) and {PIX_OUT}")


def metrics(r, p):
    return dict(N=len(r), R2_pearson=np.corrcoef(r, p)[0, 1] ** 2,
                R2_1to1=1 - np.sum((p - r) ** 2) / np.sum((r - r.mean()) ** 2),
                RMSE=np.sqrt(np.mean((p - r) ** 2)), bias=np.mean(p - r))


def stats():
    px = pd.read_csv(PIX_OUT); mp = pd.read_csv(MAP_OUT)
    pd.set_option('display.width', 200)
    print("== Pooled pixels (30 m)"); print(pd.Series(metrics(px['ref'].values, px['prod'].values)).round(3))
    print("\n== Per site, pooled pixels")
    print(pd.DataFrame([dict(site=s, **metrics(g['ref'].values, g['prod'].values))
                        for s, g in px.groupby('site')]).round(3).to_string())
    print("\n== Within-map spatial agreement (median over maps)")
    mp['sd_ratio'] = mp['sd_prod'] / mp['sd_ref']
    print(mp[['r', 'sd_ratio', 'rmse', 'bias']].median().round(3))
    print(mp.groupby('site')[['n', 'r', 'sd_ratio', 'rmse', 'bias']].median().round(2).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--extract', action='store_true')
    ap.add_argument('--stats', action='store_true')
    a = ap.parse_args()
    if a.extract:
        extract()
    if a.stats:
        stats()
