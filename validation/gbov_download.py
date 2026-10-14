"""
Download GBOV LAI reference data (RM-7 plot LAI/PAI from NEON DHP, and LP-3 20 m/300 m
upscaled LAI) for the CONUS NEON sites, via the GBOV public order API.

Credentials (GBOV portal account, https://gbov.land.copernicus.eu) are read from the
GBOV_USER / GBOV_PASS environment variables, or from ~/.gbov_credentials (line 1 = user,
line 2 = password). They are never written by this script.

  python gbov_download.py                 # order + poll + download all sites
  python gbov_download.py --sites UKFS LENO
  python gbov_download.py --status        # list existing orders only

Orders one site at a time (small, restartable). Skips sites whose zip already exists.
The full 2013-2022 record is downloaded; the independence split (dates/sites outside the
Kang et al. 2021 validation set) is applied in the analysis step, not here.
"""
import os, sys, json, time, argparse, zipfile
import requests, urllib3

urllib3.disable_warnings()                      # GBOV manual uses curl -k
API   = "https://gbov.land.copernicus.eu/api/v1"
OUT   = r"K:\Hangkai\CONUS_LAI_validation\GBOV"
START, STOP = "2013-01-01T00:00:00", "2023-01-01T00:00:00"
CONUS_EXCLUDE = {"AK", "HI", "PR"}
DONE_STATES = {"done", "completed", "complete", "finished", "ready", "success", "succeeded"}
FAIL_STATES = {"failed", "error", "cancelled", "canceled"}


def creds():
    u, p = os.environ.get("GBOV_USER"), os.environ.get("GBOV_PASS")
    if u and p:
        return u, p
    f = os.path.join(os.path.expanduser("~"), ".gbov_credentials")
    if os.path.exists(f):
        lines = [l.strip() for l in open(f, encoding="utf-8") if l.strip()]
        if len(lines) >= 2:
            return lines[0], lines[1]
    sys.exit("No GBOV credentials: set GBOV_USER/GBOV_PASS or create ~/.gbov_credentials")


def conus_neon_sites():
    v = requests.get(f"{API}/sites", verify=False, timeout=60).json()["value"]
    out = []
    for s in v:
        if s["network"] != "NEON" or s["country"] != "United States":
            continue
        if s["lon"] < -125 or s["lon"] > -66 or s["lat"] < 24 or s["lat"] > 50:   # drop HI/AK/PR
            continue
        out.append(s)
    return sorted(out, key=lambda s: s["codename"])


def state_of(o):
    for k in ("status", "state", "order_status"):
        if k in o and o[k] is not None:
            return str(o[k]).lower()
    return ""


def order_site(sess, site):
    body = {"sites": [site["gbov_site_id"]], "rm_types": ["RM7"], "lp_types": ["LP3"],
            "include_rm": True, "start": START, "stop": STOP}
    r = sess.post(f"{API}/orders", json=body, verify=False, timeout=120)
    r.raise_for_status()
    o = r.json()
    o = o.get("value", o) if isinstance(o, dict) else o
    return o["id"] if isinstance(o, dict) else o


def wait_and_download(sess, oid, dest, poll=60, max_wait=6 * 3600):
    t0 = time.time()
    while time.time() - t0 < max_wait:
        o = sess.get(f"{API}/orders/{oid}", verify=False, timeout=60).json()
        o = o.get("value", o)
        st = state_of(o)
        if st in FAIL_STATES:
            raise RuntimeError(f"order {oid} failed: {json.dumps(o)[:500]}")
        if st in DONE_STATES or not st:
            r = sess.get(f"{API}/orders/{oid}/download", verify=False, timeout=600, stream=True)
            ctype = r.headers.get("Content-Type", "")
            if r.ok and ("zip" in ctype or "octet" in ctype):
                tmp = dest + ".part"
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                zipfile.ZipFile(tmp).testzip()
                os.replace(tmp, dest)
                return
        print(f"    order {oid}: {st or 'pending'} ({int(time.time() - t0)} s)", flush=True)
        time.sleep(poll)
    raise TimeoutError(f"order {oid} not ready after {max_wait} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", nargs="*", help="NEON codes (default: all CONUS GBOV NEON sites)")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    sess = requests.Session()
    sess.auth = creds()
    sess.headers["Content-type"] = "application/json"

    if a.status:
        print(json.dumps(sess.get(f"{API}/orders", verify=False, timeout=60).json(), indent=1)[:5000])
        return

    sites = conus_neon_sites()
    json.dump(sites, open(os.path.join(OUT, "gbov_conus_neon_sites.json"), "w"), indent=1)
    if a.sites:
        sites = [s for s in sites if s["codename"] in set(a.sites)]
    print(f"{len(sites)} CONUS NEON sites: {' '.join(s['codename'] for s in sites)}")

    log_p = os.path.join(OUT, "orders.json")
    log = json.load(open(log_p)) if os.path.exists(log_p) else {}
    for s in sites:
        code = s["codename"]
        dest = os.path.join(OUT, f"GBOV_{code}_RM7_LP3_2013-2022.zip")
        if os.path.exists(dest):
            print(f"{code}: have zip, skip"); continue
        try:
            oid = log.get(code) or order_site(sess, s)
            log[code] = oid; json.dump(log, open(log_p, "w"), indent=1)
            print(f"{code}: order {oid}", flush=True)
            wait_and_download(sess, oid, dest)
            with zipfile.ZipFile(dest) as z:
                z.extractall(os.path.join(OUT, code))
            print(f"{code}: downloaded + extracted", flush=True)
        except Exception as e:
            print(f"{code}: ERROR {e}", flush=True)


if __name__ == "__main__":
    main()
