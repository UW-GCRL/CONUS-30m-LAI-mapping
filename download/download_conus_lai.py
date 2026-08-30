"""
Batch-download the CONUS_Monthly_LAI_30m ImageCollection from Earth Engine
to local GeoTIFFs, one file per monthly image -- PARALLEL version.

Speed strategy: the bottleneck is per-tile round-trip latency to Earth Engine,
not local CPU or network bandwidth. So we (a) run several images concurrently
and (b) use many tile-threads per image, keeping lots of requests in flight.

Tunables (env vars, with defaults):
    LAI_IMG_WORKERS   number of images downloaded concurrently   (default 4)
    LAI_TILE_THREADS  tile threads per image                     (default 8)

Run:
    python download_conus_lai.py

Resumable: finished .tif files are skipped; interrupted images leave a
.tif.incomplete that is discarded and retried on the next run.
"""

import os
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import ee
from geedim.download import BaseImage

ASSET   = "projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m"
OUTDIR  = r"K:\Hangkai\CONUS_LAI"
PROJECT = "ee-hyou34"

IMG_WORKERS  = int(os.environ.get("LAI_IMG_WORKERS", "4"))
TILE_THREADS = int(os.environ.get("LAI_TILE_THREADS", "8"))

_print_lock = threading.Lock()
def say(msg):
    with _print_lock:
        print(msg, flush=True)

def download_one(i, n, idx):
    out = os.path.join(OUTDIR, f"{idx}.tif")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        say(f"[{i}/{n}] skip  {idx}  (already downloaded)")
        return ("skip", idx, None)

    img_id = f"{ASSET}/{idx}"
    tmp = out + ".incomplete"
    t0 = time.time()
    say(f"[{i}/{n}] fetch {idx} ...")
    try:
        img = BaseImage.from_id(img_id)
        # full footprint at native 30 m projection; many threads per image
        img.download(tmp, overwrite=True, num_threads=TILE_THREADS)
        os.replace(tmp, out)
        dt = (time.time() - t0) / 60
        mb = os.path.getsize(out) / 1e6
        say(f"[{i}/{n}] DONE  {idx}  {mb:,.0f} MB in {dt:,.1f} min ({mb/max(dt,1e-9):,.0f} MB/min)")
        return ("ok", idx, None)
    except Exception as e:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        say(f"[{i}/{n}] FAIL  {idx}: {e}")
        return ("fail", idx, str(e))

def main():
    ee.Initialize(project=PROJECT)
    os.makedirs(OUTDIR, exist_ok=True)

    coll = ee.ImageCollection(ASSET)
    ids = sorted(coll.aggregate_array("system:index").getInfo())
    n = len(ids)
    say(f"Collection has {n} images. Output dir: {OUTDIR}")
    say(f"Parallelism: {IMG_WORKERS} images x {TILE_THREADS} tile-threads = "
        f"{IMG_WORKERS*TILE_THREADS} concurrent tile requests\n")

    ok = skipped = 0
    failed = []
    with ThreadPoolExecutor(max_workers=IMG_WORKERS) as ex:
        futs = [ex.submit(download_one, i, n, idx) for i, idx in enumerate(ids, 1)]
        for fut in as_completed(futs):
            status, idx, err = fut.result()
            if status == "ok":
                ok += 1
            elif status == "skip":
                skipped += 1
            else:
                failed.append((idx, err))

    say(f"\nDone. downloaded={ok}  skipped={skipped}  failed={len(failed)}")
    if failed:
        say("Failed (re-run to retry): " + ", ".join(idx for idx, _ in failed))
        sys.exit(1)

if __name__ == "__main__":
    main()
