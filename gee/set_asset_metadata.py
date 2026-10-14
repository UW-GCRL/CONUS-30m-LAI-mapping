"""
Set (and print back) the descriptive metadata of the three CONUS-30m-LAI Earth Engine collections,
so that the asset descriptions match the README and the gap-filling code (gapfill/gapfill_lai.py).

  python set_asset_metadata.py            # show current collection + first-image properties only
  python set_asset_metadata.py --apply    # write the properties below, then show them again

Requires an authenticated Earth Engine session (`earthengine authenticate`).
"""
import sys
import ee

ee.Initialize(project="ee-hyou34")
BASE = "projects/ee-hyou34/assets"
GH = "https://github.com/UW-GCRL/CONUS-30m-LAI-mapping"
DOI = "10.5281/zenodo.22292385"

COMMON = {
    "github_repository": GH,
    "contact": "hyou34@wisc.edu; min.chen@wisc.edu",
    "doi": DOI,
    "doi_url": "https://doi.org/" + DOI,
    "citation": ("You, H., Kang, Y. & Chen, M. (2026). CONUS-30m-LAI: A monthly 30-m leaf area index "
                 "dataset for the contiguous United States, 2000-2022 [Data set]. Zenodo. "
                 "https://doi.org/" + DOI),
}
QA_CODES = ("0 retrieved value (not filled); 1 linear interpolation between two valid months of the "
            "same year; 2 climatology fill, pixel-year with >=3 valid months; 3 climatology fill, "
            "pixel-year with 1-2 valid months; 4 climatology fill using the spatial-fallback (block) "
            "climatology, because <3 years were valid for that calendar month; 255 not filled (no "
            "valid month in the pixel-year, including non-vegetated pixels and water)")
PROPS = {
    "CONUS_Monthly_LAI_30m": dict(COMMON, nodata_value=-32768, projection="EPSG:5070 (NAD83 / Conus Albers)",
        description=(
        "Monthly gap-filled 30-m Leaf Area Index (LAI) over the contiguous United States, 2000-2022. "
        "int16, scale factor 0.001 (LAI = value*0.001); unfilled pixel-years = -32768. Retrieved "
        "(not gap-filled) LAI = this collection where the companion QA = 0. Companion QA, nvalid, "
        "and obscount collections. Code, methods, and documentation: " + GH)),
    "CONUS_Monthly_LAI_30m_QA": dict(COMMON, qa_codes=QA_CODES, description=(
        "Gap-fill quality-assurance layer (uint8) for CONUS_Monthly_LAI_30m. Codes: " + QA_CODES +
        ". Docs: " + GH)),
    "CONUS_Monthly_LAI_30m_nvalid": dict(COMMON, description=(
        "Per-year count of months with a valid retrieval (0-12) per pixel (uint8) for "
        "CONUS_Monthly_LAI_30m. Docs: " + GH)),
    # existing description and title are kept; only the shared fields are added
    "CONUS_Monthly_LAI_30m_obscount": dict(COMMON),
}


def show():
    for name in PROPS:
        aid = f"{BASE}/{name}"
        props = ee.data.getAsset(aid).get("properties", {})
        print(f"\n== {name} (collection properties)")
        for k in sorted(props):
            print(f"  {k}: {props[k]}")
        first = ee.ImageCollection(aid).first()
        info = first.getInfo()
        print(f"  -- first image: {info['id']}")
        print(f"     bands: {[(b['id'], b['data_type'].get('precision'), b['data_type'].get('min'), b['data_type'].get('max')) for b in info['bands']]}")
        print(f"     properties: {info.get('properties', {})}")


if "--apply" in sys.argv:
    for name, pr in PROPS.items():
        ee.data.setAssetProperties(f"{BASE}/{name}", pr)
        print(f"set: {name}")
    # per-image QA property (all 276 images carried the old years-based wording)
    ids = ee.ImageCollection(f"{BASE}/CONUS_Monthly_LAI_30m_QA").aggregate_array("system:index").getInfo()
    for k, idx in enumerate(ids):
        ee.data.setAssetProperties(f"{BASE}/CONUS_Monthly_LAI_30m_QA/{idx}", {"qa_codes": QA_CODES})
        if (k + 1) % 50 == 0:
            print(f"  qa_codes updated on {k + 1}/{len(ids)} QA images", flush=True)
    print(f"  qa_codes updated on {len(ids)} QA images")
show()
