"""Step 5: extreme-event exposure for each CRMS vegetation site.

Hurricanes (NOAA HURDAT2): 6-hourly best-track points are interpolated to 30-minute steps.
For each site and storm: closest-approach distance (km), time of closest approach,
maximum sustained wind (kt) at that time, and the highest wind-radius class (34, 50, 64 kt)
whose quadrant radius reached the site. A simple wind proxy is also given:
  wind_proxy_kt = Vmax * min(1, (40 km / distance) ** 0.5)
(a modified Rankine decay with a fixed 40 km radius of maximum wind). It is a rough
exposure index only, not modeled wind.

Freeze (NOAA GHCN-Daily): minimum TMIN and the number of days with TMIN <= 0 C and
<= -4 C from 14 to 20 February 2021, assigned from the nearest of five stations.

Drought (NOAA nClimDiv): monthly PDSI, SPI-6 and precipitation for Louisiana climate
divisions 6 (East Central), 7 (Southwest), 8 (South Central) and 9 (Southeast); sites are joined to divisions.

Outputs: WORK_DIR/storm_exposure.csv, WORK_DIR/freeze_exposure.csv,
WORK_DIR/climdiv_monthly.csv, WORK_DIR/site_divisions.csv, outputs/tables/storm_list.csv,
outputs/tables/freeze_stations_feb2021.csv
Run: python -I scripts/05_event_exposure.py
"""
import gzip
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (RAW_DIR, WORK_DIR, TAB_DIR, STORMS, STORM_RADIUS_KM, GHCN_STATIONS,  # noqa: E402
                    FREEZE_START, FREEZE_END, CLIMDIV_STATE)

R_EARTH = 6371.0


def hav(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(a))


def bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    x = np.sin(lon2 - lon1) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def read_hurdat(path):
    rows, cur = [], None
    for line in open(path):
        p = [x.strip() for x in line.split(",")]
        if p[0].startswith("AL") and len(p) >= 3 and p[2].isdigit():
            cur = (p[0], p[1])
            continue
        lat = float(p[4][:-1]) * (1 if p[4][-1] == "N" else -1)
        lon = float(p[5][:-1]) * (-1 if p[5][-1] == "W" else 1)
        rad = [float(v) for v in p[8:20]]   # 34 kt NE SE SW NW, 50 kt ..., 64 kt ...
        rows.append([cur[0], cur[1], pd.Timestamp(p[0] + p[1].zfill(4)), p[2], p[3], lat, lon,
                     float(p[6])] + rad)
    cols = ["atcf", "name", "time", "record", "status", "lat", "lon", "vmax"] + \
        [f"r{k}_{q}" for k in (34, 50, 64) for q in ("ne", "se", "sw", "nw")]
    return pd.DataFrame(rows, columns=cols)


def storm_exposure(sites, hd):
    out, slist = [], []
    for name, yr in STORMS:
        t = hd[(hd.name == name) & (hd.time.dt.year == yr)].sort_values("time")
        if t.empty:
            print("missing storm", name, yr)
            continue
        lf = t[(t.record == "L") & t.lon.between(-94.1, -88.8) & t.lat.between(28.8, 31)]
        slist.append(dict(storm=name, year=yr, atcf=t.atcf.iat[0],
                          la_landfall_time=lf.time.iat[0] if len(lf) else pd.NaT,
                          la_landfall_lon=lf.lon.iat[0] if len(lf) else np.nan,
                          la_landfall_vmax_kt=lf.vmax.iat[0] if len(lf) else np.nan,
                          peak_vmax_kt=t.vmax.max()))
        t = t.drop_duplicates("time").set_index("time")
        num = t[["lat", "lon", "vmax"] + [c for c in t.columns if c[:3] in ("r34", "r50", "r64")]].replace(-999, np.nan)
        fine = num.resample("30min").interpolate("time")
        for s in sites.itertuples():
            dist = hav(s.lat, s.lon, fine.lat.values, fine.lon.values)
            i = int(np.nanargmin(dist))
            brg = bearing(fine.lat.values, fine.lon.values, s.lat, s.lon)
            quad = np.select([brg < 90, brg < 180, brg < 270], ["ne", "se", "sw"], "nw")
            qidx = (brg // 90).astype(int).clip(0, 3)
            cls = 0
            for k in (34, 50, 64):
                rq = np.column_stack([fine[f"r{k}_{q}"].values for q in ("ne", "se", "sw", "nw")])
                rk = rq[np.arange(len(quad)), qidx] * 1.852   # nautical miles to km
                if np.any(dist <= np.nan_to_num(rk, nan=-1)):
                    cls = k
            d0, v0 = dist[i], fine.vmax.values[i]
            out.append(dict(site=s.site, storm=name, storm_year=yr, closest_time=fine.index[i],
                            min_dist_km=round(d0, 2), vmax_closest_kt=v0,
                            wind_radius_class_kt=cls,
                            wind_proxy_kt=round(v0 * min(1.0, (40.0 / max(d0, 1e-6)) ** 0.5), 2),
                            within_radius=d0 <= STORM_RADIUS_KM))
    return pd.DataFrame(out), pd.DataFrame(slist)


def freeze(sites):
    stn_meta, rows = [], []
    cache = RAW_DIR / "ghcn" / "ghcnd-stations.txt"   # downloaded in step 1
    coords = {}
    for ln in open(cache):
        if ln[:11] in GHCN_STATIONS:
            coords[ln[:11]] = (float(ln[12:20]), float(ln[21:30]))
    for sid, label in GHCN_STATIONS.items():
        with gzip.open(RAW_DIR / "ghcn" / f"{sid}.csv.gz", "rt") as f:
            g = pd.read_csv(f, header=None, usecols=[0, 1, 2, 3, 5], names=["id", "date", "el", "val", "qf"],
                            dtype={"date": str, "qf": str})
        g = g[(g.el == "TMIN") & g.qf.isna()]
        g["date"] = pd.to_datetime(g.date, format="%Y%m%d")
        w = g[(g.date >= FREEZE_START) & (g.date <= FREEZE_END)]
        tmin = w.val / 10.0
        # winter context: lowest TMIN in every Dec-Feb winter 2006-2026
        g["winter"] = g.date.dt.year.where(g.date.dt.month <= 2, g.date.dt.year + 1)
        wmins = g[g.date.dt.month.isin([12, 1, 2]) & g.winter.between(2006, 2026)].groupby("winter").val.min() / 10
        stn_meta.append(dict(ghcn_id=sid, station=label, lat=coords[sid][0], lon=coords[sid][1],
                             days_with_data=len(w), tmin_min_c=tmin.min(), days_le_0c=int((tmin <= 0).sum()),
                             days_le_m4c=int((tmin <= -4).sum()),
                             winters_as_cold_or_colder_2006_2026=int((wmins <= tmin.min()).sum()),
                             n_winters=len(wmins)))
    st = pd.DataFrame(stn_meta)
    st.to_csv(TAB_DIR / "freeze_stations_feb2021.csv", index=False, float_format="%.2f")
    for s in sites.itertuples():
        dd = hav(s.lat, s.lon, st.lat.values, st.lon.values)
        j = int(np.argmin(dd))
        rows.append(dict(site=s.site, ghcn_id=st.ghcn_id.iat[j], ghcn_dist_km=round(dd[j], 1),
                         freeze_tmin_c=st.tmin_min_c.iat[j], freeze_days_le_0c=st.days_le_0c.iat[j]))
    return pd.DataFrame(rows), st


def climdiv(sites):
    def read(var, code):
        f = sorted((RAW_DIR / "climdiv").glob(f"climdiv-{var}-v1.0.0-*"))[-1]
        recs = []
        for ln in open(f):
            p = ln.split()
            c = p[0]
            if c[:2] == CLIMDIV_STATE and c[2:4] in ("06", "07", "08", "09") and c[4:6] == code:
                yr = int(c[6:10])
                for m, v in enumerate(p[1:13], 1):
                    v = float(v)
                    recs.append(dict(div=int(c[2:4]), year=yr, month=m, var=var, value=np.nan if v <= -99 else v))
        return recs
    recs = read("pdsidv", "05") + read("sp06dv", "74") + read("pcpndv", "01")
    cd = pd.DataFrame(recs).pivot_table(index=["div", "year", "month"], columns="var", values="value").reset_index()
    cd = cd.rename(columns={"pdsidv": "pdsi", "sp06dv": "spi6", "pcpndv": "precip_in"})
    shp = sorted((RAW_DIR / "climdiv").rglob("*.shp"))[0]
    dv = gpd.read_file(shp)
    dv = dv[dv.ST_ABBRV == "LA"].to_crs(26915)
    dv["div"] = dv.CD_2DIG.astype(int)
    pts = gpd.GeoDataFrame(sites[["site"]], geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=4326).to_crs(26915)
    j = gpd.sjoin_nearest(pts, dv[["div", "geometry"]], how="left", distance_col="dist_m")
    j = j.drop_duplicates("site")[["site", "div", "dist_m"]]
    return cd, j


def main():
    sy = pd.read_csv(WORK_DIR / "veg_site_year.csv")
    sites = sy.groupby("site", as_index=False)[["lat", "lon"]].mean()
    hd = read_hurdat(sorted((RAW_DIR / "hurdat").glob("hurdat2-1851-*.txt"))[-1])
    se, sl = storm_exposure(sites, hd)
    se.to_csv(WORK_DIR / "storm_exposure.csv", index=False)
    sl.to_csv(TAB_DIR / "storm_list.csv", index=False)
    print(sl.to_string(index=False))
    print(se.groupby("storm").agg(sites_le_50km=("min_dist_km", lambda x: (x <= 50).sum()),
                                  sites_64kt=("wind_radius_class_kt", lambda x: (x == 64).sum()),
                                  sites_34kt=("wind_radius_class_kt", lambda x: (x >= 34).sum())))
    fz, st = freeze(sites)
    fz.to_csv(WORK_DIR / "freeze_exposure.csv", index=False)
    print(st.to_string(index=False))
    cd, j = climdiv(sites)
    cd.to_csv(WORK_DIR / "climdiv_monthly.csv", index=False)
    j.to_csv(WORK_DIR / "site_divisions.csv", index=False)
    print("sites by climate division:", j["div"].value_counts().to_dict(), "max join distance m:", j.dist_m.max())


if __name__ == "__main__":
    main()
