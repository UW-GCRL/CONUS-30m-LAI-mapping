"""
Robust gap-fill for CONUS monthly LAI (see GAPFILL_PLAN.md).

Per pixel, per year: keep valid months; fill gaps within a year that has >=1 valid month.
Interior gaps -> linear interpolation. Edge/sparse gaps -> per-pixel climatology * robust
scale. Climatology uses a min-valid-years threshold with a coarse spatial fallback so it is
defined everywhere. Confidence-tiered QA. Zero-valid years are left empty (flagged 255).

Stages:
  python gapfill_lai.py --climatology   # build robust climatology (rclim/cnt/csrc per month)
  python gapfill_lai.py --fill          # fill all years -> *_gf.tif, *_qa.tif, YYYY_nvalid.tif
  python gapfill_lai.py --fill --year 2005
  python gapfill_lai.py --fill --tile-test COL ROW     # one 2048 tile, all years (prototype)

Resumable (skips existing outputs). Run --climatology before --fill.
"""
import os, glob, re, argparse, time
import numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.fill import fillnodata
from rasterio.enums import Resampling

def open_retry(path, tries=30):
    """rasterio.open (read) resilient to stale NAS 'file used by other process' locks.
       These clear only after being left UNTOUCHED for ~15-30 min, so use exponential
       backoff with long quiet windows rather than hammering (which keeps the lock alive)."""
    waits = [20, 40, 90, 180, 300, 600, 600, 900, 900] + [900] * 30
    for i in range(tries):
        try:
            return rasterio.open(path)
        except rasterio.errors.RasterioIOError as e:
            if 'used by another process' in str(e) or 'used by other process' in str(e):
                if i == 0:
                    print(f"  [lock] {os.path.basename(path)} locked; backing off (clears when left alone)...", flush=True)
                time.sleep(waits[min(i, len(waits) - 1)]); continue
            raise
    return rasterio.open(path)

SRC_DIR  = r"K:\Hangkai\CONUS_LAI"
OUT_DIR  = r"K:\Hangkai\CONUS_LAI_gapfilled"
CLIM_DIR = os.path.join(OUT_DIR, "_climatology")
NODATA   = -32768
MONTHS   = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
TILE     = 1024              # 1024 (not 2048): 16x smaller per-tile arrays -> avoids long-run memory fragmentation OOM
CKPT_EVERY = 200             # tiles between flush+checkpoint (bounds crash/reboot loss)
MIN_YEARS   = 3               # min valid years to trust a pixel-month climatology
COARSE      = 160             # spatial-fallback block (~4.8 km at 30 m)
SCALE_CLAMP = (0.3, 3.0)
# QA codes
QA_ORIG, QA_LIN, QA_CLIM_HI, QA_CLIM_LO, QA_CLIM_SP, QA_NOFILL = 0, 1, 2, 3, 4, 255

def index_files():
    idx = {}
    for p in glob.glob(os.path.join(SRC_DIR, "*_LAI.tif")):
        m = re.match(r"(\d{4})_([A-Za-z]{3})_LAI\.tif$", os.path.basename(p))
        if m:
            idx.setdefault(int(m.group(1)), {})[MONTHS.index(m.group(2))] = p
    return idx, sorted(idx)

def prof_of(path, dtype, nodata, count=1):
    with rasterio.open(path) as ds:
        pr = ds.profile
    pr.update(count=count, dtype=dtype, nodata=nodata, compress='deflate', predictor=2,
              tiled=True, blockxsize=512, blockysize=512, BIGTIFF='YES')
    return pr

def wins(W, H, tile=TILE):
    for r in range(0, H, tile):
        for c in range(0, W, tile):
            yield Window(c, r, min(tile, W-c), min(tile, H-r))

def strips(W, H, sh):
    """full-width horizontal strips (sequential top-to-bottom reads)."""
    for r in range(0, H, sh):
        yield Window(0, r, W, min(sh, H-r))

# ------------------------------------------------------------ Stage 1
def build_climatology():
    """Sequential full-width-strip climatology: each yearly file is read ONCE, top-to-bottom,
       instead of 23-way random tile access (turns seek-thrashing into streaming reads)."""
    os.makedirs(CLIM_DIR, exist_ok=True)
    idx, years = index_files()
    ref = idx[years[0]][0]
    with rasterio.open(ref) as ds: W, H = ds.width, ds.height
    fprof = prof_of(ref, 'float32', np.nan)
    cprof = prof_of(ref, 'int16', 0)
    bprof = prof_of(ref, 'uint8', 0)
    SH = 2048   # strip height (multiple of the 512 output block); ~2.6 GB working set

    for mi, mon in enumerate(MONTHS):
        rclim_p = os.path.join(CLIM_DIR, f"rclim_{mon}.tif")
        if os.path.exists(rclim_p):
            print(f"clim {mon}: exists, skip"); continue
        srcs = [idx[y][mi] for y in years if mi in idx[y]]
        clim_p = os.path.join(CLIM_DIR, f"clim_{mon}.tif")
        cnt_p  = os.path.join(CLIM_DIR, f"cnt_{mon}.tif")
        print(f"clim {mon}: mean+count over {len(srcs)} years (sequential strips)", flush=True)

        # 1a) per-pixel mean + valid-year count -- full-width strips, each file read in order
        with rasterio.open(clim_p, 'w', **fprof) as cd, rasterio.open(cnt_p, 'w', **cprof) as nd:
            hs = [rasterio.open(s) for s in srcs]
            try:
                for w in strips(W, H, SH):
                    acc = np.zeros((w.height, w.width), np.float32)
                    cnt = np.zeros((w.height, w.width), np.int16)
                    for h in hs:
                        a = h.read(1, window=w); v = a > 0   # gaps are 0 (NOT -32768); valid LAI floors at 100
                        acc[v] += a[v]; cnt[v] += 1
                    mean = np.where(cnt > 0, acc/np.maximum(cnt, 1), np.nan).astype(np.float32)
                    cd.write(mean, 1, window=w); nd.write(cnt, 1, window=w)
            finally:
                for h in hs: h.close()
        print(f"clim {mon}: 1a done, building coarse+robust", flush=True)

        # 1b) coarse seasonal field = decimated average read of the climatology (one fast GDAL pass)
        cH, cW = (H + COARSE - 1)//COARSE, (W + COARSE - 1)//COARSE
        with rasterio.open(clim_p) as cd:
            coarse = cd.read(1, out_shape=(cH, cW), resampling=Resampling.average).astype(np.float32)
        cmask = np.isfinite(coarse).astype(np.uint8)
        coarse = fillnodata(np.where(cmask == 1, coarse, 0).astype(np.float32),
                            mask=cmask, max_search_distance=100)

        # 1c) robust climatology = native where cnt>=MIN_YEARS else coarse-upsampled; + source flag
        with rasterio.open(clim_p) as cd, rasterio.open(cnt_p) as nd, \
             rasterio.open(rclim_p, 'w', **fprof) as rd, \
             rasterio.open(os.path.join(CLIM_DIR, f"csrc_{mon}.tif"), 'w', **bprof) as sd:
            for w in strips(W, H, SH):
                cl = cd.read(1, window=w); cn = nd.read(1, window=w)
                rr = (int(w.row_off) + np.arange(w.height))//COARSE
                cc = (np.arange(w.width))//COARSE
                up = coarse[np.ix_(rr, cc)]
                native = cn >= MIN_YEARS
                out = np.where(native, cl, up).astype(np.float32)
                src = np.where(native, 0, 1).astype(np.uint8)
                rd.write(out, 1, window=w); sd.write(src, 1, window=w)
        print(f"clim {mon}: robust climatology done", flush=True)

# ------------------------------------------------------------ fill core
def fill_block(stack, rclim, csrc):
    """stack:(12,h,w) int16/NODATA; rclim:(12,h,w) float32/nan; csrc:(12,h,w) uint8.
       -> gf int16, qa uint8, nvalid uint8"""
    T = stack.shape[0]; shp = stack.shape[1:]
    valid = stack > 0                 # gaps encoded as 0 (and -32768 in the June-2022 mosaic); valid LAI > 0
    nvalid = valid.sum(0).astype(np.uint8)
    fillable = nvalid >= 1

    # FAST PATH: tile has no fillable gaps (all-gap ocean/water tiles, or fully-complete land).
    # Then output == input (valid kept, gaps -> clean nodata), QA = original/not-filled.
    if not ((~valid) & fillable[None]).any():
        qa = np.where(valid, QA_ORIG, QA_NOFILL).astype(np.uint8)
        return np.where(valid, stack, NODATA).astype(np.int16), qa, nvalid

    f = stack.astype(np.float32); f[~valid] = np.nan

    # prev/next valid value & time index (12-step fwd/bwd fill)
    prev_v = np.full(stack.shape, np.nan, np.float32); prev_i = np.full(stack.shape, -1, np.int16)
    lv = np.full(shp, np.nan, np.float32); li = np.full(shp, -1, np.int16)
    for t in range(T):
        vt = valid[t]; lv = np.where(vt, f[t], lv); li = np.where(vt, t, li)
        prev_v[t] = lv; prev_i[t] = li
    next_v = np.full(stack.shape, np.nan, np.float32); next_i = np.full(stack.shape, -1, np.int16)
    nv = np.full(shp, np.nan, np.float32); ni = np.full(shp, -1, np.int16)
    for t in range(T-1, -1, -1):
        vt = valid[t]; nv = np.where(vt, f[t], nv); ni = np.where(vt, t, ni)
        next_v[t] = nv; next_i[t] = ni

    # amplitude scale: sum(obs)/sum(rclim) over valid months, clamped (fast, robust via clamp)
    with np.errstate(invalid='ignore', divide='ignore'):
        obs_sum  = np.nansum(np.where(valid, f, np.float32(0)), axis=0)      # float32 zero: no float64 upcast
        clim_sum = np.nansum(np.where(valid, rclim, np.float32(0)), axis=0)
        scale = np.where(clim_sum > 0, obs_sum / np.maximum(clim_sum, 1e-6), 1.0)
    scale = np.clip(np.where(np.isfinite(scale), scale, 1.0), *SCALE_CLAMP)

    tt = np.arange(T, dtype=np.float32).reshape(T, 1, 1)   # float32: (tt - prev_i) stays float32, no int64/float64 blowup
    interior = (~valid) & (prev_i >= 0) & (next_i >= 0) & fillable[None]
    edge     = (~valid) & (~interior) & fillable[None]

    gf = f.copy()
    qa = np.where(valid, QA_ORIG, QA_NOFILL).astype(np.uint8)

    # interior: linear interpolation between bracketing valid months
    denom = (next_i - prev_i).astype(np.float32); denom[denom == 0] = 1
    gf = np.where(interior, prev_v + (next_v - prev_v)*(tt - prev_i)/denom, gf)
    qa = np.where(interior, QA_LIN, qa)

    # edge/sparse: climatology * scale ; fallback to nearest valid if rclim missing
    clim_fill = rclim * scale[None]
    nearest = np.where(np.isnan(prev_v), next_v, prev_v)
    edge_val = np.where(np.isnan(clim_fill), nearest, clim_fill)
    gf = np.where(edge, edge_val, gf)
    # QA tier for climatology fills
    base = np.where(nvalid[None] >= 3, QA_CLIM_HI, QA_CLIM_LO)
    tier = np.where(csrc == 1, QA_CLIM_SP, base).astype(np.uint8)
    qa = np.where(edge, tier, qa)

    gf_out = np.where(np.isnan(gf), NODATA, np.rint(gf)).astype(np.int16)
    return gf_out, qa, nvalid

# ------------------------------------------------------------ Stage 2
def fill_year(year, idx, W, H, gfp, qap, nvp, only=None):
    mons = MONTHS
    gf_paths = [os.path.join(OUT_DIR, f"{year}_{m}_LAI_gf.tif") for m in mons]
    qa_paths = [os.path.join(OUT_DIR, f"{year}_{m}_LAI_qa.tif") for m in mons]
    nv_path  = os.path.join(OUT_DIR, f"{year}_nvalid.tif")
    if only is None and all(os.path.exists(p) for p in gf_paths):
        print(f"{year}: exists, skip"); return
    sh = [rasterio.open(idx[year][mi]) if mi in idx[year] else None for mi in range(12)]
    rc = [rasterio.open(os.path.join(CLIM_DIR, f"rclim_{m}.tif")) for m in mons]
    sc = [rasterio.open(os.path.join(CLIM_DIR, f"csrc_{m}.tif"))  for m in mons]
    gd = [rasterio.open(p, 'w', **gfp) for p in gf_paths]
    qd = [rasterio.open(p, 'w', **qap) for p in qa_paths]
    nd = rasterio.open(nv_path, 'w', **nvp)
    try:
        ws = [only] if only else list(wins(W, H))
        for k, w in enumerate(ws):
            stack = np.full((12, w.height, w.width), NODATA, np.int16)
            for mi in range(12):
                if sh[mi] is not None: stack[mi] = sh[mi].read(1, window=w)
            rclim = np.stack([rc[mi].read(1, window=w) for mi in range(12)])
            csrc  = np.stack([sc[mi].read(1, window=w) for mi in range(12)])
            gf, qa, nvalid = fill_block(stack, rclim, csrc)
            for mi in range(12):
                gd[mi].write(gf[mi], 1, window=w); qd[mi].write(qa[mi], 1, window=w)
            nd.write(nvalid, 1, window=w)
            if k % 50 == 0: print(f"  {year}: tile {k+1}/{len(ws)}", flush=True)
    finally:
        for h in sh+rc+sc+gd+qd+[nd]:
            if h: h.close()
    print(f"{year}: done")

def _year_done(y):
    return os.path.exists(os.path.join(OUT_DIR, f"{y}.done"))

def _open_outs(todo, mons, gfp, qap, nvp, resume):
    """open per-year gf/qa/nvalid datasets: 'w' fresh, or 'r+' to resume (preserves written tiles)."""
    def o(path, prof): return rasterio.open(path, 'r+') if resume else rasterio.open(path, 'w', **prof)
    gd = {y: [o(os.path.join(OUT_DIR, f"{y}_{m}_LAI_gf.tif"), gfp) for m in mons] for y in todo}
    qd = {y: [o(os.path.join(OUT_DIR, f"{y}_{m}_LAI_qa.tif"), qap) for m in mons] for y in todo}
    nd = {y:  o(os.path.join(OUT_DIR, f"{y}_nvalid.tif"), nvp) for y in todo}
    return gd, qd, nd

def _close_outs(todo, gd, qd, nd):
    for y in todo:
        for h in gd[y] + qd[y] + [nd[y]]:
            if h: h.close()

def fill_batched(idx, years, W, H, gfp, qap, nvp, batch=6, max_tiles=None, ckpt_tag='single'):
    """Tile-outer / year-inner in batches: read each climatology tile ONCE per batch and reuse
       across the batch's years. Resumable at TILE level: every CKPT_EVERY tiles the outputs are
       flushed (close+reopen 'r+') and a checkpoint tile index is recorded, so a crash/reboot
       loses <= CKPT_EVERY tiles, not the whole batch. Per-year .done marks batch completion."""
    mons = MONTHS
    rc = [open_retry(os.path.join(CLIM_DIR, f"rclim_{m}.tif")) for m in mons]
    sc = [open_retry(os.path.join(CLIM_DIR, f"csrc_{m}.tif"))  for m in mons]
    all_wins = list(wins(W, H))
    if max_tiles: all_wins = all_wins[:max_tiles]
    try:
        for bi in range(0, len(years), batch):
            todo = [y for y in years[bi:bi+batch] if not _year_done(y)]
            if not todo:
                print(f"batch {years[bi:bi+batch]}: all done, skip", flush=True); continue
            ckpt = os.path.join(OUT_DIR, f".ckpt_{ckpt_tag}_b{bi}")
            start = 0
            outs_exist = all(os.path.exists(os.path.join(OUT_DIR, f"{y}_{m}_LAI_gf.tif")) for y in todo for m in mons)
            if os.path.exists(ckpt) and outs_exist:
                try: start = max(0, int(open(ckpt).read().strip()))
                except Exception: start = 0
            resume = start > 0
            print(f"batch years {todo} ({len(all_wins)} tiles)" + (f" RESUME@tile {start}" if resume else ""), flush=True)
            sh = {y: [open_retry(idx[y][mi]) if mi in idx[y] else None for mi in range(12)] for y in todo}
            gd, qd, nd = _open_outs(todo, mons, gfp, qap, nvp, resume)
            try:
                for k in range(start, len(all_wins)):
                    w = all_wins[k]
                    rclim = np.stack([rc[mi].read(1, window=w) for mi in range(12)])   # ONCE per tile
                    csrc  = np.stack([sc[mi].read(1, window=w) for mi in range(12)])
                    for y in todo:
                        stack = np.full((12, w.height, w.width), NODATA, np.int16)
                        for mi in range(12):
                            if sh[y][mi] is not None: stack[mi] = sh[y][mi].read(1, window=w)
                        gf, qa, nvalid = fill_block(stack, rclim, csrc)
                        for mi in range(12):
                            gd[y][mi].write(gf[mi], 1, window=w); qd[y][mi].write(qa[mi], 1, window=w)
                        nd[y].write(nvalid, 1, window=w)
                    if k % 25 == 0: print(f"  tiles {k+1}/{len(all_wins)}", flush=True)
                    if not max_tiles and (k + 1) % CKPT_EVERY == 0:      # flush to disk + record checkpoint
                        _close_outs(todo, gd, qd, nd)
                        with open(ckpt, 'w') as cf: cf.write(str(k + 1))
                        gd, qd, nd = _open_outs(todo, mons, gfp, qap, nvp, resume=True)
            finally:
                for y in todo:
                    for h in sh[y]:
                        if h: h.close()
                _close_outs(todo, gd, qd, nd)
            if not max_tiles:
                for y in todo:
                    open(os.path.join(OUT_DIR, f"{y}.done"), 'w').close()
                if os.path.exists(ckpt):
                    try: os.remove(ckpt)
                    except OSError: pass
            print(f"batch {todo} done", flush=True)
    finally:
        for h in rc + sc: h.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--climatology', action='store_true')
    ap.add_argument('--fill', action='store_true')
    ap.add_argument('--year', type=int)
    ap.add_argument('--tile-test', nargs=2, type=int, metavar=('COL', 'ROW'))
    ap.add_argument('--smoke', type=int, help='fill only first N tiles (test)')
    ap.add_argument('--group', nargs=2, type=int, metavar=('G', 'N'),
                    help='process year-group G of N (for parallel workers)')
    a = ap.parse_args()
    if a.climatology:
        build_climatology(); return
    if a.fill:
        os.makedirs(OUT_DIR, exist_ok=True)
        idx, years = index_files()
        ref = idx[years[0]][0]
        with rasterio.open(ref) as ds: W, H = ds.width, ds.height
        gfp = prof_of(ref, 'int16', NODATA); qap = prof_of(ref, 'uint8', QA_NOFILL)
        nvp = prof_of(ref, 'uint8', 0)
        if a.tile_test or a.year:
            only = Window(a.tile_test[0]*TILE, a.tile_test[1]*TILE, TILE, TILE) if a.tile_test else None
            for y in ([a.year] if a.year else years):
                fill_year(y, idx, W, H, gfp, qap, nvp, only=only)
        else:
            tag = 'single'
            if a.group:
                G, N = a.group
                years = [y for i, y in enumerate(years) if i % N == G]   # interleave years across workers
                tag = f"g{G}"
                print(f"worker group {G}/{N}: years {years}", flush=True)
            fill_batched(idx, years, W, H, gfp, qap, nvp, max_tiles=a.smoke, ckpt_tag=tag)
        return
    ap.print_help()

if __name__ == "__main__":
    main()
