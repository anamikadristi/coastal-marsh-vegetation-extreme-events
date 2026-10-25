"""Step 6: MODIS Terra MOD13Q1 (Collection 6.1, 250 m, 16-day) NDVI at CRMS site buffers.

Read anonymously from Microsoft Planetary Computer (STAC collection modis-13Q1-061; the
signing call only adds a public read token, no account or key). For every composite from
2006 to 2025 and each of the four MODIS tiles over coastal Louisiana (h09v05, h09v06,
h10v05, h10v06), only the small window that holds the site buffer pixels is read.

Buffer: MODIS pixels whose centers are within BUFFER_M (1 km) of the site location.
Quality: pixel_reliability 0 (good) or 1 (marginal) kept; others set to missing.
Output: WORK_DIR/modis/ndvi_pixels.npz (pixel time series) and WORK_DIR/modis/pixels.csv
Run: python -I scripts/06_modis_ndvi.py   (about 30-60 minutes; resumable by year)
"""
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import WORK_DIR, BUFFER_M, VEG_YEARS  # noqa: E402

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
SIN = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
X0, Y0, TILE, NPIX = -20015109.354, 10007554.677, 1111950.5196666666, 4800
CELL = TILE / NPIX
OUT = WORK_DIR / "modis"
OUT.mkdir(exist_ok=True)
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_TIMEOUT="60",
                GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="5")
ASSETS = {"ndvi": "250m_16_days_NDVI", "rel": "250m_16_days_pixel_reliability"}


def buffer_pixels(sites):
    tr = Transformer.from_crs("EPSG:4326", SIN, always_xy=True)
    x, y = tr.transform(sites.lon.values, sites.lat.values)
    rows = []
    k = int(np.ceil(BUFFER_M / CELL)) + 1
    for s, sx, sy in zip(sites.site, x, y):
        gc0, gr0 = int((sx - X0) // CELL), int((Y0 - sy) // CELL)
        for dr in range(-k, k + 1):
            for dc in range(-k, k + 1):
                gr, gc = gr0 + dr, gc0 + dc
                cx, cy = X0 + (gc + 0.5) * CELL, Y0 - (gr + 0.5) * CELL
                d = np.hypot(cx - sx, cy - sy)
                if d <= BUFFER_M:
                    rows.append(dict(site=s, grow=gr, gcol=gc, h=gc // NPIX, v=gr // NPIX,
                                     r=gr % NPIX, c=gc % NPIX, dist_m=round(d, 1)))
    px = pd.DataFrame(rows)
    # a pixel can fall in two overlapping buffers; keep one row per site and pixel
    return px.drop_duplicates(["site", "grow", "gcol"]).reset_index(drop=True)


def read_window(href, win):
    for attempt in range(5):
        try:
            with rasterio.Env(**GDAL_ENV):
                with rasterio.open(pc.sign(href)) as r:
                    return r.read(1, window=win)
        except Exception as e:      # retry transient HTTP errors
            err = e
            time.sleep(5 * (attempt + 1))
    print("  FAILED", href.split("/")[-1], repr(err)[:100], flush=True)
    return None


def main():
    sy = pd.read_csv(WORK_DIR / "veg_site_year.csv")
    sites = sy.groupby("site", as_index=False)[["lat", "lon"]].mean()
    px = buffer_pixels(sites)
    px.to_csv(OUT / "pixels.csv", index=False)
    upx = px.drop_duplicates(["grow", "gcol"])[["grow", "gcol", "h", "v", "r", "c"]].reset_index(drop=True)
    tiles = {}
    for (h, v), g in upx.groupby(["h", "v"]):
        r0, r1, c0, c1 = g.r.min(), g.r.max(), g.c.min(), g.c.max()
        tiles[f"h{h:02d}v{v:02d}"] = dict(win=Window(c0, r0, c1 - c0 + 1, r1 - r0 + 1),
                                          idx=g.index.values, rr=g.r.values - r0, cc=g.c.values - c0)
    print(f"sites {len(sites)}, buffer pixels {len(px)}, unique pixels {len(upx)}, tiles "
          + ", ".join(f"{k} {t['win'].height}x{t['win'].width}" for k, t in tiles.items()), flush=True)
    cat = pystac_client.Client.open(STAC)
    ydir = OUT / "years"
    ydir.mkdir(exist_ok=True)
    for year in range(VEG_YEARS[0], VEG_YEARS[1] + 1):
        yf = ydir / f"{year}.npz"
        if yf.exists():
            continue
        t0 = time.time()
        items = [i for i in cat.search(collections=["modis-13Q1-061"], bbox=[-94.0, 28.8, -88.7, 30.7],
                                       datetime=f"{year}-01-01/{year}-12-31").items()
                 if i.id.startswith("MOD13Q1") and i.id.split(".")[2] in tiles
                 and i.id.split(".")[1][1:5] == str(year)]
        dates = sorted({i.id.split(".")[1][1:] for i in items})
        didx = {d: k for k, d in enumerate(dates)}
        ndvi = np.full((len(dates), len(upx)), -3000, dtype=np.int16)
        rel = np.full((len(dates), len(upx)), 255, dtype=np.int16)
        jobs = [(it, a) for it in items for a in ASSETS]

        def job(ia):
            it, a = ia
            t = tiles[it.id.split(".")[2]]
            return it.id, a, read_window(it.assets[ASSETS[a]].href, t["win"])
        failed = []
        with ThreadPoolExecutor(4) as ex:
            for iid, a, arr in ex.map(job, jobs):
                if arr is None:
                    failed.append(f"{iid}:{a}")
                    continue
                t = tiles[iid.split(".")[2]]
                (ndvi if a == "ndvi" else rel)[didx[iid.split(".")[1][1:]], t["idx"]] = arr[t["rr"], t["cc"]]
        np.savez_compressed(yf, dates=np.array(dates), ndvi=ndvi, rel=rel)
        json.dump(dict(year=year, items=len(items), dates=len(dates), reads=len(jobs), failed=failed),
                  open(ydir / f"{year}.json", "w"))
        print(f"{year}: {len(items)} items, {len(dates)} composites, {len(failed)} failed reads, "
              f"{time.time() - t0:.0f} s", flush=True)
    parts = [np.load(ydir / f"{y}.npz") for y in range(VEG_YEARS[0], VEG_YEARS[1] + 1)]
    np.savez_compressed(OUT / "ndvi_pixels.npz", dates=np.concatenate([p["dates"] for p in parts]),
                        ndvi=np.concatenate([p["ndvi"] for p in parts]),
                        rel=np.concatenate([p["rel"] for p in parts]),
                        grow=upx.grow.values, gcol=upx.gcol.values)
    logs = [json.load(open(ydir / f"{y}.json")) for y in range(VEG_YEARS[0], VEG_YEARS[1] + 1)]
    json.dump(dict(items=sum(l["items"] for l in logs), composites=sum(l["dates"] for l in logs),
                   reads=sum(l["reads"] for l in logs), failed=sum((l["failed"] for l in logs), [])),
              open(OUT / "modis_log.json", "w"), indent=1)
    print("wrote ndvi_pixels.npz", flush=True)


if __name__ == "__main__":
    main()
