"""
merge_lai_to_ref_grid.py
========================
Merge the per-state monthly LAI GeoTIFFs (exported from GEE) into a single
CONUS-wide monthly mosaic aligned to a reference grid, in the delivered encoding.

Encoding
--------
Per-state export tiles (input):  band 1 = LAI as uint16, scale factor 0.01
    (LAI x100); masked (cloud/shadow/water) and non-vegetated pixels are 0.
    A second uint16 obs_count band may be present but is not required here.
Delivered CONUS mosaic (output): int16, single band, scale factor 0.001
    (LAI x1000). Unobserved and non-vegetated pixels are 0; valid LAI >= 0.1.
The x100 -> x1000 re-encoding is a fixed integer multiply (SCALE_MULT = 10).

These retrieved mosaics are the input to the gap-filling step
(``gapfill/gapfill_lai.py``), which produces the gap-filled product (also int16,
LAI x1000) with unfilled pixels set to -32768.

Usage
-----
    python merge_lai_to_ref_grid.py

Configure the paths in the ``__main__`` block at the bottom.

Algorithm
---------
1. Read the reference raster (CRS, transform, width, height).
2. Pre-fill the int16 output raster with 0.
3. For each per-state tile:
   a. Warp band 1 (LAI x100) to the reference grid via WarpedVRT (nearest).
   b. Process block-by-block (default 512x512 pixels).
   c. Re-encode: value = round(LAI_x100 * 10) as int16 where valid (LAI_x100 > 0),
      else 0.
   d. Write valid pixels to the output (merge policy: last-one-wins).
4. Build overview levels (2, 4, 8, 16, 32) for fast display.

Dependencies
------------
    pip install rasterio numpy tqdm

Reference
---------
Kang, Y., et al. (2021). A data-driven approach to estimate leaf area
index for Landsat images over the contiguous US. Remote Sensing of
Environment, 258, 112383. https://doi.org/10.1016/j.rse.2021.112383
"""

from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window
from rasterio.vrt import WarpedVRT
from tqdm import tqdm

NODATA = -32768          # declared no-data of the delivered files (retrieved gaps are stored as 0)
SCALE_MULT = 10          # per-state LAI x100  ->  delivered LAI x1000


def merge_lai_to_ref_grid(
    in_dir,
    ref_tif,
    out_tif,
    recursive=True,
    blocksize=512,
    prefer="src",
    resampling=Resampling.nearest,
    compress="LZW",
    bigtiff="IF_SAFER",
    scale_mult=SCALE_MULT,
):
    """
    Merge per-state LAI tiles onto the exact grid of ref_tif and write the delivered
    int16 (LAI x1000) CONUS mosaic; unobserved and non-vegetated pixels are 0.

    Parameters
    ----------
    in_dir : str or Path
        Directory of per-state GeoTIFF tiles (band 1 = LAI x100, uint16).
    ref_tif : str or Path
        Reference raster defining the output grid (CRS, transform, size).
    out_tif : str or Path
        Path for the output CONUS mosaic GeoTIFF.
    recursive : bool
        Whether to search subdirectories for .tif files.
    blocksize : int
        Tile size (pixels) for block-by-block processing.
    prefer : str
        Merge policy: 'src' = last-one-wins (overwrites); 'dest' = first-one-wins
        (fills only where the output is still 0).
    resampling : rasterio.enums.Resampling
        Resampling for warping (default: nearest, to preserve 30 m values).
    compress : str
        GDAL compression codec (default: 'LZW').
    bigtiff : str
        BigTIFF mode (default: 'IF_SAFER').
    scale_mult : int
        Integer multiplier from the per-state LAI x100 band to the delivered
        LAI x1000 encoding (default 10).
    """
    in_dir = Path(in_dir)
    tifs = sorted(in_dir.rglob("*.tif") if recursive else in_dir.glob("*.tif"))
    tifs += sorted(in_dir.rglob("*.tiff") if recursive else in_dir.glob("*.tiff"))
    if not tifs:
        raise SystemExit(f"No .tif/.tiff files found in: {in_dir}")

    with rasterio.open(ref_tif) as ref:
        dst_crs       = ref.crs
        dst_transform = ref.transform
        width, height = ref.width, ref.height

    profile = {
        "driver":    "GTiff",
        "width":     width,
        "height":    height,
        "count":     1,
        "dtype":     "int16",
        "crs":       dst_crs,
        "transform": dst_transform,
        "nodata":    NODATA,
        "tiled":     True,
        "blockxsize": blocksize,
        "blockysize": blocksize,
        "compress":  compress,
        "BIGTIFF":   bigtiff,
        "predictor": 2,
    }

    out_tif = Path(out_tif)
    out_tif.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: create output raster pre-filled with 0 (unobserved / non-vegetation)
    with rasterio.open(out_tif, "w", **profile) as dst_w:
        blank = np.zeros((1, blocksize, blocksize), dtype="int16")
        for row in range(0, height, blocksize):
            h = min(blocksize, height - row)
            for col in range(0, width, blocksize):
                w = min(blocksize, width - col)
                dst_w.write(blank[:, :h, :w], window=Window(col, row, w, h))

    # Step 2: warp and merge each state tile
    with rasterio.open(out_tif, "r+") as dst:
        for tif in tqdm(tifs, desc="Merging LAI tiles"):
            with rasterio.open(str(tif)) as src:
                # warp band 1 (LAI x100) to float32 so the x10 scaling is exact
                with WarpedVRT(
                    src,
                    crs=dst_crs,
                    transform=dst_transform,
                    width=width,
                    height=height,
                    resampling=resampling,
                    src_nodata=src.nodata,
                    nodata=np.nan,
                    dtype="float32",
                ) as vrt:
                    for row in range(0, height, blocksize):
                        h = min(blocksize, height - row)
                        for col in range(0, width, blocksize):
                            w = min(blocksize, width - col)
                            win = Window(col, row, w, h)

                            lai = vrt.read(1, window=win)          # LAI x100 (float32)
                            valid = np.isfinite(lai) & (lai > 0)
                            src_i16 = np.where(
                                valid, np.rint(lai * scale_mult), 0
                            ).astype(np.int16)                     # gaps -> 0

                            dest_block = dst.read(1, window=win)
                            out_block = dest_block.copy()
                            if prefer == "src":
                                out_block[valid] = src_i16[valid]
                            elif prefer == "dest":
                                fill = (dest_block == 0) & valid
                                out_block[fill] = src_i16[fill]
                            else:
                                raise ValueError("prefer must be 'src' or 'dest'")

                            dst.write(out_block[np.newaxis, :, :], window=win)

        try:
            dst.build_overviews([2, 4, 8, 16, 32], Resampling.average)
            dst.update_tags(ns="rio_overview", resampling="average")
        except Exception as e:
            print(f"Overview build skipped: {e}")

    print(f"Done. Wrote: {out_tif}")


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # Configure paths and run for all 12 months of a given year           #
    # ------------------------------------------------------------------ #
    root_in  = "/mnt/cephfs-mount/hangkai/CONUS_LAI/2021/"
    ref_tif  = "/mnt/cephfs-mount/hangkai/CONUS_LAI/Sample/CONUS_Disturbance_1985.tif"
    out_root = "/mnt/cephfs-mount/hangkai/CONUS_LAI/backup_all/"

    month_labels = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]

    for m in month_labels:
        in_dir  = str(Path(root_in) / m)
        out_tif = str(Path(out_root) / f"2021_{m}_LAI.tif")

        if not Path(in_dir).exists():
            print(f"[Skip] {in_dir} not found.")
            continue

        print(f"==> Processing {m} ...")
        merge_lai_to_ref_grid(
            in_dir=in_dir,
            ref_tif=ref_tif,
            out_tif=out_tif,
            recursive=True,
            blocksize=512,
            prefer="src",
            resampling=Resampling.nearest,
        )
