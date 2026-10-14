"""
Field validation of the RELEASED CONUS-30m-LAI product against GBOV RM-7 (NEON DHP plot LAI).

Reference: GBOV RM-7 v2.0, Warren-Wilson clumping-corrected LAI (PAI), upward + downward
photos combined, as in Kang et al. (2021). Forest plots require a valid upward (overstory)
photo; short-vegetation plots use the downward photo (plus upward if present).

Plot locations: GBOV RM-7 files give a single site coordinate for all plots, so plot centres
are taken from the NEON locations API (<plotID>.basePlot.dhp; cached to neon_dhp_plot_coords.csv).

Product: for each plot-date, the 30 m pixel containing the plot centre (and the 3x3 mean as a
sensitivity case) in the product month containing the field date:
  retrieved   K:\\Hangkai\\CONUS_LAI\\YYYY_Mon_LAI.tif           int16 x1000, 0 = no retrieval
  gap-filled  K:\\Hangkai\\CONUS_LAI_gapfilled\\_up_LAI[_done]\\  int16 x1000, -32768 = unfilled
  QA          K:\\Hangkai\\CONUS_LAI_gapfilled\\_up_QA\\YYYY_Mon_QA.tif  uint8 (0 = observed)

Subsets: 'kang' = Kang et al. (2021) sites and years (19 sites, 2013-2018); 'new' = all other
records (2019-2022 at those sites, or any year at the 8 sites Kang did not use).

  python gbov_validate.py --extract   # build matchups CSV (reads NAS rasters)
  python gbov_validate.py --stats     # metrics from the matchups CSV
"""
import os, glob, io, zipfile, argparse, time
import numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from pyproj import Transformer

RM7_ZIP = r"K:\Hangkai\CONUS_LAI_gapfilled\GBOV\COPERNICUS_GBOV_RM7_2026913233546.zip"
RET_DIR = r"K:\Hangkai\CONUS_LAI"
GF_DIRS = [r"K:\Hangkai\CONUS_LAI_gapfilled\_up_LAI", r"K:\Hangkai\CONUS_LAI_gapfilled\_up_LAI_done"]
QA_DIR  = r"K:\Hangkai\CONUS_LAI_gapfilled\_up_QA"
OUT_DIR = r"K:\Hangkai\CONUS_LAI_validation"
MATCH   = os.path.join(OUT_DIR, "gbov_rm7_matchups.csv")
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
SCALE = 0.001
KANG_SITES = {'BART','BLAN','CPER','DSNY','HARV','JERC','JORN','MOAB','NIWO','ONAQ',
              'ORNL','OSBS','SCBI','SERC','STEI','STER','TALL','UNDE','WOOD'}
FOREST = {'Mixed Forest','Evergreen Needleleaf Forest','Deciduous Broadleaf Forest'}


def load_rm7():
    z = zipfile.ZipFile(RM7_ZIP)
    dfs = []
    for n in z.namelist():
        if n.endswith('.csv'):
            d = pd.read_csv(io.BytesIO(z.read(n)), sep=';')
            d['plot'] = '_'.join(os.path.basename(n).split('_')[3:5])
            dfs.append(d)
    d = pd.concat(dfs, ignore_index=True).drop_duplicates(['plot', 'TIME_IS'])
    d['site'] = d['plot'].str[:4]
    d['date'] = pd.to_datetime(d['TIME_IS'], format='%Y%m%dT%H%M%SZ')
    d['year'], d['month'] = d['date'].dt.year, d['date'].dt.month
    for c in ['LAI_Warren_up', 'LAI_Warren_down', 'LAIe_Warren_up', 'LAIe_Warren_down']:
        d[c] = pd.to_numeric(d[c], errors='coerce').where(lambda s: s != -999)
    up_ok = (d['up_flag'] == 0) & d['LAI_Warren_up'].notna()
    dn_ok = (d['down_flag'] == 0) & d['LAI_Warren_down'].notna()
    forest = d['IGBP_class'].isin(FOREST)
    d['usable'] = dn_ok & (up_ok | (~forest & d['LAI_Warren_up'].isna()))
    d['pai']  = d['LAI_Warren_down'] + d['LAI_Warren_up'].where(up_ok, 0)
    d['paie'] = d['LAIe_Warren_down'] + d['LAIe_Warren_up'].where(up_ok, 0)
    d['subset'] = np.where(d['site'].isin(KANG_SITES) & (d['year'] <= 2018), 'kang', 'new')
    return d


def plot_coords(plots):
    """NEON plot-centre coordinates. GBOV RM-7 files carry one site-level coordinate for every
       plot, so plot locations come from the NEON locations API (<plot>.basePlot.dhp), cached."""
    import requests
    p_csv = os.path.join(OUT_DIR, 'neon_dhp_plot_coords.csv')
    have = pd.read_csv(p_csv) if os.path.exists(p_csv) else pd.DataFrame(columns=['plot', 'lat'])
    have = have[have['lat'].notna()]
    todo = [p for p in plots if p not in set(have['plot'])]
    rows = []
    for p in todo:
        rec = dict(plot=p)
        for suf in ['basePlot.dhp', 'towerPlot.dhp', 'basePlot.all', 'towerPlot.all']:
            r = None
            for wait in [0, 2, 10, 30]:
                time.sleep(wait)
                try:
                    r = requests.get(f"https://data.neonscience.org/api/v0/locations/{p}.{suf}", timeout=60)
                except requests.RequestException:
                    continue
                if r.status_code == 429 or r.status_code >= 500:
                    continue
                break
            dd = r.json().get('data') if (r is not None and r.ok) else None
            if dd:
                props = {q['locationPropertyName']: q['locationPropertyValue']
                         for q in dd.get('locationProperties', [])}
                rec.update(loc=dd['locationName'], lat=dd['locationDecimalLatitude'],
                           lon=dd['locationDecimalLongitude'],
                           plotSize=props.get('Value for Plot size'),
                           nlcd=props.get('Value for National Land Cover Database (2001)'),
                           coordUnc=props.get('Value for Coordinate uncertainty'))
                break
        rows.append(rec)
    out = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(p_csv, index=False)
    return out


def find(dirs, name):
    for dd in dirs:
        p = os.path.join(dd, name)
        if os.path.exists(p):
            return p
    return None


def sample(path, rows, cols, nodata_vals):
    """centre pixel + 3x3 mean (valid only) for each (row, col)."""
    c1 = np.full(len(rows), np.nan); c3 = np.full(len(rows), np.nan)
    if path is None:
        return c1, c3
    with rasterio.open(path) as ds:
        for k, (r, c) in enumerate(zip(rows, cols)):
            a = ds.read(1, window=Window(c - 1, r - 1, 3, 3)).astype(float)
            bad = np.isin(a, nodata_vals)
            if not bad[1, 1]:
                c1[k] = a[1, 1]
            if (~bad).any():
                c3[k] = a[~bad].mean()
    return c1, c3


def extract():
    os.makedirs(OUT_DIR, exist_ok=True)
    d = load_rm7()
    d = d[d['usable'] & d['year'].between(2000, 2022)].copy()
    pc = plot_coords(sorted(d['plot'].unique()))
    d = d.merge(pc[['plot', 'lat', 'lon', 'plotSize', 'nlcd']], on='plot', how='left')
    print(f"plots with NEON coordinates: {d.drop_duplicates('plot')['lat'].notna().sum()}"
          f" / {d['plot'].nunique()}")
    d = d[d['lat'].notna()].copy()
    with rasterio.open(os.path.join(RET_DIR, '2015_Jul_LAI.tif')) as ds:
        inv = ~ds.transform
    x, y = Transformer.from_crs('EPSG:4326', 'EPSG:5070', always_xy=True).transform(
        d['lon'].values, d['lat'].values)
    cr = [inv * (xi, yi) for xi, yi in zip(x, y)]
    d['col'] = [int(np.floor(c)) for c, _ in cr]
    d['row'] = [int(np.floor(r)) for _, r in cr]
    cols = {k: np.full(len(d), np.nan) for k in ['ret', 'ret3', 'gf', 'gf3', 'qa']}
    groups = d.groupby(['year', 'month']).indices
    for i, ((yr, mo), idx) in enumerate(sorted(groups.items())):
        name = f"{yr}_{MON[mo - 1]}"
        rr, cc = d['row'].values[idx], d['col'].values[idx]
        a, a3 = sample(find([RET_DIR], f"{name}_LAI.tif"), rr, cc, [0, -32768])
        cols['ret'][idx], cols['ret3'][idx] = a * SCALE, a3 * SCALE
        g, g3 = sample(find(GF_DIRS, f"{name}_LAI.tif"), rr, cc, [-32768])
        cols['gf'][idx], cols['gf3'][idx] = g * SCALE, g3 * SCALE
        q, _ = sample(find([QA_DIR], f"{name}_QA.tif"), rr, cc, [])
        cols['qa'][idx] = q
        print(f"  {name}: {len(idx)} records ({i + 1}/{len(groups)})", flush=True)
    for k, v in cols.items():
        d[k] = v
    keep = ['site', 'plot', 'date', 'year', 'month', 'IGBP_class', 'subset', 'lat', 'lon',
            'plotSize', 'nlcd',
            'pai', 'paie', 'row', 'col', 'ret', 'ret3', 'gf', 'gf3', 'qa']
    d[keep].to_csv(MATCH, index=False)
    print(f"wrote {MATCH}: {len(d)} plot-dates")


def metrics(ref, est):
    m = np.isfinite(ref) & np.isfinite(est)
    r, e = ref[m], est[m]
    n = len(r)
    if n < 10:
        return dict(N=n)
    rmse = np.sqrt(np.mean((e - r) ** 2))
    return dict(N=n,
                R2_pearson=np.corrcoef(r, e)[0, 1] ** 2,
                R2_1to1=1 - np.sum((e - r) ** 2) / np.sum((r - r.mean()) ** 2),
                RMSE=rmse, MAE=np.mean(np.abs(e - r)), bias=np.mean(e - r),
                within=100 * np.mean(np.abs(e - r) <= np.maximum(1.0, 0.2 * r)))


def table(d, by, ref='pai', est='ret'):
    rows = [dict(group='ALL', **metrics(d[ref].values, d[est].values))]
    for g, s in d.groupby(by):
        rows.append(dict(group=g, **metrics(s[ref].values, s[est].values)))
    return pd.DataFrame(rows).round(3)


def stats():
    d = pd.read_csv(MATCH, parse_dates=['date'])
    pd.set_option('display.width', 200)
    print(f"plot-dates: {len(d)}; with retrieval: {d['ret'].notna().sum()}; "
          f"with gap-filled value: {d['gf'].notna().sum()}")
    print("\n== Retrieved product (QA 0), centre pixel, by biome"); print(table(d, 'IGBP_class'))
    print("\n== Retrieved, by subset (Kang records vs new records)"); print(table(d, 'subset'))
    print("\n== Retrieved, 3x3 mean, by biome"); print(table(d, 'IGBP_class', est='ret3'))
    print("\n== Retrieved vs effective PAI"); print(table(d, 'subset', ref='paie'))
    com = d[d['ret'].notna() & d['gf'].notna()]
    print(f"\n== Common set (both retrieved and gap-filled): {len(com)}")
    print(pd.concat([table(com, 'subset', est='ret').assign(product='retrieved'),
                     table(com, 'subset', est='gf').assign(product='gap-filled')]))
    fill = d[d['ret'].isna() & d['gf'].notna()]
    print(f"\n== Gap-filled only (no retrieval that month), by QA code: {len(fill)}")
    print(table(fill, 'qa', est='gf'))
    print("\n== Retrieved, leaf-on (May-Sep) vs other months")
    print(table(d.assign(season=np.where(d['month'].between(5, 9), 'May-Sep', 'other')), 'season'))
    print("\n== Retrieved, by site"); print(table(d, 'site').to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--extract', action='store_true')
    ap.add_argument('--stats', action='store_true')
    a = ap.parse_args()
    if a.extract:
        extract()
    if a.stats:
        stats()
