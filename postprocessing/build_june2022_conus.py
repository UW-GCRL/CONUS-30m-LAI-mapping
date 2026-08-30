"""
build_june2022_conus.py

Rebuild the CONUS 2022-June monthly LAI mosaic from regenerated per-state GEE tiles and
make it a drop-in replacement matching the existing CONUS_Monthly_LAI_30m files.

The per-state export tiles are:  uint16, 2 bands (LAI[x100], obs_count), EPSG:5070.
The existing CONUS files are:     int16,  1 band  (LAI[x1000]),  nodata -32768, custom
                                  Albers grid 160169 x 105432.

Reconciliation per pixel:
    LAI_out = round(LAI_x100 * 10)      where obs_count > 0        (x100 -> x1000)
            = -32768 (nodata)           where obs_count == 0 (masked / no observation)

Two stages (memory-safe, tiled; run in order, resumable):
  1) reconcile + reproject each state tile onto the exact reference grid  -> ALIGN_DIR/*.tif
  2) composite the aligned tiles (valid-wins) into OUT_FILE on the reference grid

Run:  python build_june2022_conus.py            # both stages
      python build_june2022_conus.py --align    # stage 1 only
      python build_june2022_conus.py --mosaic   # stage 2 only
"""
import os, glob, argparse
import numpy as np
import rasterio
from rasterio import Affine
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import Window

TILE_DIR  = r"G:\Hangkai\Download\june2022_tiles"
ALIGN_DIR = r"G:\Hangkai\Download\june2022_aligned"
REF_FILE  = r"K:\Hangkai\CONUS_LAI\2021_Jun_LAI.tif"    # any existing file = reference grid
OUT_FILE  = r"K:\Hangkai\CONUS_LAI\2022_Jun_LAI.tif"
SCALE_MULT = 10          # export LAI x100  ->  existing LAI x1000
NODATA     = -32768
BLOCK      = 2048

def ref_grid():
    with rasterio.open(REF_FILE) as ds:
        return ds.crs, ds.transform, ds.width, ds.height, ds.profile.copy()

# ---------------------------------------------------------------- stage 1
def reconcile_and_align(tp, ref_crs, ref_tr, ref_w, ref_h):
    out = os.path.join(ALIGN_DIR, os.path.splitext(os.path.basename(tp))[0] + ".tif")
    if os.path.exists(out):
        return out
    with rasterio.open(tp) as src:
        # tile footprint -> reference pixel bbox
        l, b, r, t = transform_bounds(src.crs, ref_crs, *src.bounds, densify_pts=21)
        inv = ~ref_tr
        cols = [ (inv*(x, y))[0] for x, y in [(l, t), (r, t), (l, b), (r, b)] ]
        rows = [ (inv*(x, y))[1] for x, y in [(l, t), (r, t), (l, b), (r, b)] ]
        c0 = max(0, int(np.floor(min(cols)))); c1 = min(ref_w, int(np.ceil(max(cols))))
        r0 = max(0, int(np.floor(min(rows)))); r1 = min(ref_h, int(np.ceil(max(rows))))
        if c1 <= c0 or r1 <= r0:
            print(f"  {os.path.basename(tp)}: no overlap, skip"); return None
        w, h = c1 - c0, r1 - r0
        dst_tr = ref_tr * Affine.translation(c0, r0)

        lai = src.read(1).astype(np.float32)
        obs = src.read(2)
        val = np.where(obs > 0, np.rint(lai * SCALE_MULT), NODATA).astype(np.int16)

        dst = np.full((h, w), NODATA, np.int16)
        reproject(source=val, destination=dst,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=dst_tr, dst_crs=ref_crs,
                  src_nodata=NODATA, dst_nodata=NODATA,
                  resampling=Resampling.nearest)

        prof = dict(driver='GTiff', height=h, width=w, count=1, dtype='int16',
                    crs=ref_crs, transform=dst_tr, nodata=NODATA,
                    compress='deflate', predictor=2, tiled=True,
                    blockxsize=512, blockysize=512, BIGTIFF='YES')
        with rasterio.open(out, 'w', **prof) as d:
            d.write(dst, 1)
    print(f"  aligned {os.path.basename(out)}  ({w}x{h} @ col{c0},row{r0})")
    return out

def stage_align():
    os.makedirs(ALIGN_DIR, exist_ok=True)
    crs, tr, W, H, _ = ref_grid()
    tiles = sorted(glob.glob(os.path.join(TILE_DIR, "LAI_*_2022_06*.tif")))
    print(f"stage 1: reconcile+align {len(tiles)} tiles")
    for tp in tiles:
        reconcile_and_align(tp, crs, tr, W, H)

# ---------------------------------------------------------------- stage 2
def stage_mosaic():
    crs, tr, W, H, prof = ref_grid()
    prof.update(count=1, dtype='int16', nodata=NODATA, compress='deflate',
                predictor=2, tiled=True, blockxsize=512, blockysize=512, BIGTIFF='YES')
    aligned = sorted(glob.glob(os.path.join(ALIGN_DIR, "LAI_*_2022_06*.tif")))
    print(f"stage 2: composite {len(aligned)} aligned tiles -> {OUT_FILE}")

    # index each aligned tile by its pixel offset/extent on the reference grid
    meta = []
    inv = ~tr
    for p in aligned:
        with rasterio.open(p) as ds:
            c = int(round((inv * (ds.transform.c, ds.transform.f))[0]))
            r = int(round((inv * (ds.transform.c, ds.transform.f))[1]))
            meta.append((p, c, r, ds.width, ds.height))

    with rasterio.open(OUT_FILE, 'w', **prof) as dst:
        n = 0
        for row in range(0, H, BLOCK):
            for col in range(0, W, BLOCK):
                wh, ww = min(BLOCK, H - row), min(BLOCK, W - col)
                acc = np.full((wh, ww), NODATA, np.int16)
                for p, tc, tr_, tw, th in meta:
                    ic0, ic1 = max(col, tc), min(col + ww, tc + tw)
                    ir0, ir1 = max(row, tr_), min(row + wh, tr_ + th)
                    if ic1 <= ic0 or ir1 <= ir0:
                        continue
                    with rasterio.open(p) as ds:
                        sub = ds.read(1, window=Window(ic0 - tc, ir0 - tr_, ic1 - ic0, ir1 - ir0))
                    a0, b0 = ir0 - row, ic0 - col
                    tgt = acc[a0:a0 + sub.shape[0], b0:b0 + sub.shape[1]]
                    m = sub != NODATA
                    tgt[m] = sub[m]
                dst.write(acc, 1, window=Window(col, row, ww, wh))
            n += 1
            if n % 10 == 0:
                print(f"  row {row+BLOCK}/{H}", flush=True)
    print("mosaic done:", OUT_FILE)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--align', action='store_true')
    ap.add_argument('--mosaic', action='store_true')
    a = ap.parse_args()
    if a.align:  stage_align();  return
    if a.mosaic: stage_mosaic(); return
    stage_align(); stage_mosaic()

if __name__ == "__main__":
    main()
