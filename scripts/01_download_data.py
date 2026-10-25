"""Step 1: download all public, keyless inputs into RAW_DIR.

Sources: CPRA CIMS (CRMS vegetation and hourly hydrology, basin polygons), NOAA NHC HURDAT2,
NOAA NCEI GHCN-Daily and nClimDiv, and the FGDC CSDGM XML schema. The current file names for
HURDAT2 and nClimDiv carry a date stamp, so they are discovered from the directory listings.
MODIS NDVI is read later (step 6) straight from Microsoft Planetary Computer.
Run: python -I scripts/01_download_data.py
"""
import datetime as dt
import json
import re
import sys
import zipfile
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (RAW_DIR, SOURCES, CLIMDIV, GHCN, GHCN_STATIONS, HURDAT_DIR, FGDC_XSD)  # noqa: E402

UA = {"User-Agent": "Mozilla/5.0 (research download script)"}


def fetch(url, dest, force=False):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"  exists: {dest.name}")
        return
    print(f"  downloading {dest.name}")
    with requests.get(url, stream=True, timeout=900, headers=UA) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        tmp.rename(dest)


def latest(listing_url, pattern):
    html = requests.get(listing_url, timeout=60, headers=UA).text
    names = sorted(set(re.findall(pattern, html)))
    if not names:
        raise RuntimeError(f"no file matching {pattern} at {listing_url}")
    return names[-1]


def main():
    today = dt.date.today().isoformat()
    jobs = dict(SOURCES)
    # HURDAT2 Atlantic: pick the newest 1851-YYYY file (sorted by season then revision date)
    def stamp(t):                     # revision stamps are MMDDYY or MMDDYYYY
        m, d, y = t[1][:2], t[1][2:4], t[1][4:]
        return (t[0], (2000 + int(y)) if len(y) == 2 else int(y), m, d)
    hd = sorted(set(re.findall(r"hurdat2-1851-(\d{4})-(\d+)\.txt",
                               requests.get(HURDAT_DIR, timeout=60, headers=UA).text)), key=stamp)
    hname = f"hurdat2-1851-{hd[-1][0]}-{hd[-1][1]}.txt"
    jobs["hurdat2"] = (HURDAT_DIR + hname, f"hurdat/{hname}")
    for var in ("pcpndv", "pdsidv", "sp06dv"):
        n = latest(CLIMDIV, rf"climdiv-{var}-v1\.0\.0-\d{{8}}")
        jobs[f"climdiv_{var}"] = (CLIMDIV + n, f"climdiv/{n}")
    for sid in GHCN_STATIONS:
        jobs[f"ghcn_{sid}"] = (f"{GHCN}/{sid}.csv.gz", f"ghcn/{sid}.csv.gz")
    jobs["ghcn_stations"] = (GHCN.replace("by_station", "ghcnd-stations.txt"), "ghcn/ghcnd-stations.txt")
    for k in ["", "-sect01", "-sect02", "-sect03", "-sect04", "-sect05", "-sect06", "-sect07",
              "-sect08", "-sect09", "-sect10"]:
        n = f"fgdc-std-001-1998{k}.xsd"
        jobs[f"xsd{k}"] = (FGDC_XSD + n, f"schemas/{n}")

    log = {}
    for key, (url, rel) in jobs.items():
        dest = RAW_DIR / rel
        fetch(url, dest)
        log[key] = {"url": url, "file": rel, "bytes": dest.stat().st_size, "accessed": today}
    for z in ["crms/CRMS_Marsh_Vegetation.zip", "gis/Basins.zip", "climdiv/CONUS_CLIMATE_DIVISIONS.shp.zip"]:
        out = RAW_DIR / z.replace(".shp.zip", "").replace(".zip", "")
        if not out.exists():
            with zipfile.ZipFile(RAW_DIR / z) as zf:
                zf.extractall(out)
    old = {}
    lp = RAW_DIR / "download_log.json"
    if lp.exists():
        old = json.loads(lp.read_text())
    for k, v in log.items():          # keep the first access date for files already on disk
        if k in old and old[k]["file"] == v["file"]:
            v["accessed"] = old[k]["accessed"]
    lp.write_text(json.dumps(log, indent=2))
    for k, v in log.items():
        print(f"{k:22s} {v['bytes'] / 1e6:10.2f} MB  {v['accessed']}")


if __name__ == "__main__":
    main()
