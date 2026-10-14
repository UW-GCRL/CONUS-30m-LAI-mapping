"""
Set (and print back) the descriptive metadata of the three CONUS-30m-LAI Earth Engine collections,
so that the asset descriptions match the README and the gap-filling code (gapfill/gapfill_lai.py).

  python set_asset_metadata.py            # show current collection + first-image properties only
  python set_asset_metadata.py --apply    # write the properties below, then show them again
  python set_asset_metadata.py --apply --collections-only   # skip the 276 per-image QA updates

The descriptions were checked against the data: LAI = -32768, QA = 255, and nvalid = 0 are masked and
mark the same pixels; Earth Engine GeoTIFF downloads write masked pixels with their stored values and
declare NoData = 0 for 8-bit images, hence the unmask + noData export settings.

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
    # must match the Zenodo record exactly (title, creators; Zenodo record type is "Data paper")
    "citation": ("You, H., Kang, Y. & Chen, M. (2026). Monthly 30-m Leaf Area Index Dataset for the "
                 "Contiguous United States from 2000 to 2022. Zenodo. https://doi.org/" + DOI),
}
PROJ = ("Albers Equal Area Conic on the WGS 84 datum, with the EPSG:5070 projection parameters "
        "(standard parallels 29.5N and 45.5N, latitude of origin 23N, central meridian 96W); "
        "30 m grid, origin (-2417805, 3321225)")
QA_CODES = ("0 retrieved value (not filled); 1 linear interpolation between two valid months of the "
            "same year; 2 climatology fill, pixel-year with >=3 valid months; 3 climatology fill, "
            "pixel-year with 1-2 valid months; 4 climatology fill using the spatial-fallback (block) "
            "climatology, because <3 years were valid for that calendar month; 255 not filled (no "
            "valid month in the pixel-year, including non-vegetated pixels and water)")
PROPS = {
    "CONUS_Monthly_LAI_30m": dict(COMMON, nodata_value=-32768, projection=PROJ, description=(
        "Monthly gap-filled 30-m Leaf Area Index (LAI) for the contiguous United States, January 2000 "
        "to December 2022 (276 images; band b1). LAI was retrieved from Landsat 5, 7, 8, and 9 with the "
        "Random Forest algorithm of Kang et al. (2021); gaps were filled with a climatology-anchored "
        "temporal method, and retrieved values are unchanged.\n\n"
        "LAI = value × 0.001 (int16). No-data: -32768 (masked). Retrieved (not gap-filled) LAI: pixels "
        "where CONUS_Monthly_LAI_30m_QA equals 0.\n\n"
        "Documentation and export settings: " + GH + "\n"
        "Citation: " + COMMON["citation"])),
    "CONUS_Monthly_LAI_30m_QA": dict(COMMON, qa_codes=QA_CODES, projection=PROJ, description=(
        "Gap-fill quality-assurance layer (uint8) for CONUS_Monthly_LAI_30m. Codes: " + QA_CODES +
        ". Code 255 is masked. GeoTIFF: 8-bit downloads declare NoData = 0, which would hide code 0, "
        "so export QA.unmask(255) with formatOptions {noData: 255}. Docs: " + GH)),
    "CONUS_Monthly_LAI_30m_nvalid": dict(COMMON, projection=PROJ, description=(
        "Per-year count of months with a valid retrieval (1-12) per pixel (uint8) for "
        "CONUS_Monthly_LAI_30m; pixel-years with no valid month are masked (stored 0). Docs: " + GH)),
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
if "--apply" in sys.argv and "--collections-only" not in sys.argv:
    # per-image QA property (all 276 images carried the old years-based wording)
    ids = ee.ImageCollection(f"{BASE}/CONUS_Monthly_LAI_30m_QA").aggregate_array("system:index").getInfo()
    for k, idx in enumerate(ids):
        ee.data.setAssetProperties(f"{BASE}/CONUS_Monthly_LAI_30m_QA/{idx}", {"qa_codes": QA_CODES})
        if (k + 1) % 50 == 0:
            print(f"  qa_codes updated on {k + 1}/{len(ids)} QA images", flush=True)
    print(f"  qa_codes updated on {len(ids)} QA images")
show()
