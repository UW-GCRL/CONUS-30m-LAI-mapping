"""
merge_lai_to_ref_grid.py
========================
Merge per-state monthly LAI GeoTIFFs (exported from GEE) into a single
CONUS-wide monthly mosaic aligned to a reference grid.

Usage
-----
    python merge_lai_to_ref_grid.py

Configure the paths in the ``__main__`` block at the bottom.

Algorithm
---------
1. Read the reference raster (CRS, transform, width, height).
2. Pre-fill a float32 output raster with nodata (-9999 by default).
3. For each state tile:
   a. Warp to the reference grid via WarpedVRT (nearest resampling).
   b. Process block-by-block (default 512x512 pixels).
   c. Apply scale factor 0.01 (converts uint16 x100 to real LAI units).
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


def merge_continuous_to_ref_grid(
    in_dir,
    ref_tif,
    out_tif,
    recursive=True,
    blocksize=512,
    prefer="src",
    dst_nodata=None,
    src_valid_fn=None,
    resampling=Resampling.nearest,
    compress="LZW",
    bigtiff="IF_SAFER",
    scale_factor=0.01,
):
    """
    Merge continuous rasters (e.g., LAI) onto the exact grid of ref_tif.

    Parameters
    ----------
    in_dir : str or Path
        Directory containing state-level GeoTIFF files.
    ref_tif : str or Path
        Reference raster that defines the output grid
        (CRS, transform, width, height).
    out_tif : str or Path
        Path for the output CONUS mosaic GeoTIFF.
    recursive : bool
        Whether to search subdirectories for .tif files.
    blocksize : int
        Tile size (pixels) for block-by-block processing.
    prefer : str
        Merge policy: ``'src'`` = last-one-wins (overwrites);
        ``'dest'`` = first-one-wins (fills gaps only).
    dst_nodata : float or None
        Output nodata value. Defaults to ref raster nodata or -9999.
    src_valid_fn : callable or None
        Function(array, nodata) -> bool mask defining valid pixels.
        Defaults to: finite and != nodata.
    resampling : rasterio.enums.Resampling
        Resampling method for warping (default: nearest).
    compress : str
        GDAL compression codec (default: 'LZW').
    bigtiff : str
        BigTIFF mode (default: 'IF_SAFER').
    scale_factor : float
        Multiply valid source pixels by this value before writing.
        Use 0.01 to convert LAI stored as uint16 x100 to real units.
    """
    in_dir = Path(in_dir)
    tifs = sorted(in_dir.rglob("*.tif") if recursive else in_dir.glob("*.tif"))
    tifs += sorted(in_dir.rglob("*.tiff") if recursive else in_dir.glob("*.tiff"))
    if not tifs:
        raise SystemExit(f"No .tif/.tiff files found in: {in_dir}")

    # Read reference grid metadata
    with rasterio.open(ref_tif) as ref:
        dst_crs       = ref.crs
        dst_transform = ref.transform
        width, height = ref.width, ref.height
        ref_nodata    = ref.nodata
        out_dtype     = "float32"

    if dst_nodata is None:
        dst_nodata = ref_nodata if ref_nodata is not None else -9999.0

    profile = {
        "driver":    "GTiff",
        "width":     width,
        "height":    height,
        "count":     1,
        "dtype":     out_dtype,
        "crs":       dst_crs,
        "transform": dst_transform,
        "nodata":    dst_nodata,
        "tiled":     True,
        "blockxsize": blocksize,
        "blockysize": blocksize,
        "compress":  compress,
        "BIGTIFF":   bigtiff,
        "predictor": 2,
    }

    if src_valid_fn is None:
        def src_valid_fn(arr, nodata_val):
            if nodata_val is None:
                return np.isfinite(arr)
            return np.isfinite(arr) & (arr != nodata_val)

    out_tif = Path(out_tif)
    out_tif.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Create output raster pre-filled with nodata
    with rasterio.open(out_tif, "w", **profile) as dst_w:
        blank = np.full((1, blocksize, blocksize), dst_nodata, dtype=out_dtype)
        for row in range(0, height, blocksize):
            h = min(blocksize, height - row)
            for col in range(0, width, blocksize):
                w = min(blocksize, width - col)
                dst_w.write(blank[:, :h, :w], window=Window(col, row, w, h))

    # Step 2: Warp and merge each state tile
    with rasterio.open(out_tif, "r+") as dst:
        for tif in tqdm(tifs, desc="Merging LAI tiles"):
            with rasterio.open(str(tif)) as src:
                src_nodata = src.nodata
                with WarpedVRT(
                    src,
                    crs=dst_crs,
                    transform=dst_transform,
                    width=width,
                    height=height,
                    resampling=resampling,
                    src_nodata=src_nodata,
                    nodata=dst_nodata,
                    dtype=out_dtype,
                ) as vrt:
                    for row in range(0, height, blocksize):
                        h = min(blocksize, height - row)
                        for col in range(0, width, blocksize):
                            w = min(blocksize, width - col)
                            win = Window(col, row, w, h)

                            dest_block = dst.read(1, window=win)
                            src_block  = vrt.read(1, window=win)

                            dest_valid = src_valid_fn(dest_block, dst_nodata)
                            src_valid  = src_valid_fn(src_block,  dst_nodata)

                            # Scale valid pixels from uint16 x100 to real LAI
                            if scale_factor is not None and scale_factor != 1.0:
                                src_block_scaled = src_block.astype("float32", copy=False)
                                src_block_scaled[src_valid] = (
                                    src_block_scaled[src_valid] * scale_factor
                                )
                            else:
                                src_block_scaled = src_block

                            out_block = dest_block.copy()
                            if prefer == "src":
                                out_block[src_valid] = src_block_scaled[src_valid]
                            elif prefer == "dest":
                                fill = (~dest_valid) & src_valid
                                out_block[fill] = src_block_scaled[fill]
                            else:
                                raise ValueError("prefer must be 'src' or 'dest'")

                            dst.write(out_block[np.newaxis, :, :], window=win)

        # Build overviews for fast display
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
        merge_continuous_to_ref_grid(
            in_dir=in_dir,
            ref_tif=ref_tif,
            out_tif=out_tif,
            recursive=True,
            blocksize=512,
            prefer="src",
            dst_nodata=None,
            src_valid_fn=None,
            resampling=Resampling.nearest,
            compress="LZW",
            bigtiff="IF_SAFER",
            scale_factor=0.01,
        )
