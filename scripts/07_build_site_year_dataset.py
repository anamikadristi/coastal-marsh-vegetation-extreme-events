"""Step 7: build the main derived dataset, one row per CRMS site and survey year.

Joins vegetation (step 3), growing-season hydrology (step 4), event exposure and drought
indices (step 5) and MODIS NDVI (step 6). Also derives:
  marsh_type_base   modal CRMS community in the site's first three survey years
  salinity_index    community-weighted salinity optimum (ppt). Each species' optimum is its
                    cover-weighted mean of growing-season site salinity over all plot visits
                    (weighted averaging, species with >= 20 plot visits); the plot index is the
                    cover-weighted mean of the optima of species present, averaged by site-year.
                    It is a composition index on the salinity scale, not a measured salinity.
  ndvi_gs_mean      mean April-September MOD13Q1 NDVI of good or marginal pixels within 1 km,
                    after dropping pixels whose all-years median June-September NDVI is < 0.2
                    (mostly open water). A year needs 9 or more valid April-September composites.
Outputs: data/crms_marsh_site_year.csv, data/crms_veg_sites.csv, data/crms_veg_sites.geojson,
         WORK_DIR/ndvi_site_composite.csv.gz, outputs/tables/species_salinity_optima.csv
Run: python -I scripts/07_build_site_year_dataset.py
"""
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import WORK_DIR, DATA_DIR, TAB_DIR, GROW_MONTHS  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
MIN_SPECIES_VISITS = 20
WATER_NDVI = 0.2


def salinity_index(hy):
    sp = pd.read_csv(WORK_DIR / "veg_species_cover.csv.gz", parse_dates=["date"])
    sp = sp.merge(hy[["site", "year", "sal_gs_mean"]], on=["site", "year"], how="left")
    cal = sp.dropna(subset=["sal_gs_mean"])
    cal = cal.assign(w=cal.cover_num * cal.sal_gs_mean)
    opt = cal.groupby("species").agg(n_visits=("station", "size"), cover_sum=("cover_num", "sum"),
                                     wsum=("w", "sum"), halophyte=("halophyte", "first"))
    opt["sal_optimum_ppt"] = opt.wsum / opt.cover_sum
    opt = opt[opt.n_visits >= MIN_SPECIES_VISITS]
    opt[["n_visits", "sal_optimum_ppt", "halophyte"]].sort_values("sal_optimum_ppt") \
        .to_csv(TAB_DIR / "species_salinity_optima.csv", float_format="%.3f")
    s = sp.merge(opt[["sal_optimum_ppt"]], left_on="species", right_index=True, how="inner")
    s["w"] = s.cover_num * s.sal_optimum_ppt
    pv = s.groupby(["site", "station", "year"]).agg(w=("w", "sum"), c=("cover_num", "sum"))
    pv["idx"] = pv.w / pv.c
    print(f"species with salinity optima: {len(opt)} of {sp.species.nunique()}")
    return pv.groupby(["site", "year"]).idx.mean().rename("salinity_index").reset_index()


def ndvi(sites_order):
    z = np.load(WORK_DIR / "modis" / "ndvi_pixels.npz")
    px = pd.read_csv(WORK_DIR / "modis" / "pixels.csv")
    dates = pd.to_datetime([f"{d[:4]}-01-01" for d in z["dates"]]) + \
        pd.to_timedelta([int(d[4:]) - 1 for d in z["dates"]], unit="D")
    v = z["ndvi"].astype(float) / 10000.0
    v[(z["ndvi"] == -3000) | ~np.isin(z["rel"], [0, 1])] = np.nan
    jun_sep = dates.month.isin([6, 7, 8, 9])
    med = np.nanmedian(v[jun_sep], axis=0)
    key = pd.Series(np.arange(len(z["grow"])), index=pd.MultiIndex.from_arrays([z["grow"], z["gcol"]]))
    px["j"] = key.loc[list(zip(px.grow, px.gcol))].values
    px["land"] = med[px.j.values] >= WATER_NDVI
    rows = []
    for site, g in px[px.land].groupby("site"):
        m = np.nanmean(v[:, g.j.values], axis=1)
        n = np.sum(~np.isnan(v[:, g.j.values]), axis=1)
        rows.append(pd.DataFrame(dict(site=site, date=dates, ndvi=m, n_pix=n)))
    nc = pd.concat(rows)
    nc = nc[nc.n_pix >= 3]
    nc.to_csv(WORK_DIR / "ndvi_site_composite.csv.gz", index=False, compression="gzip", float_format="%.4f")
    gs = nc[nc.date.dt.month.between(*GROW_MONTHS)]
    y = gs.groupby(["site", gs.date.dt.year.rename("year")]).agg(
        ndvi_gs_mean=("ndvi", "mean"), ndvi_gs_n=("ndvi", "count")).reset_index()
    y.loc[y.ndvi_gs_n < 9, "ndvi_gs_mean"] = np.nan   # needs 9 of about 11 Apr-Sep composites
    cnt = px.groupby("site").agg(ndvi_pixels_total=("j", "size"), ndvi_pixels_land=("land", "sum")).reset_index()
    print(f"NDVI: {nc.site.nunique()} sites with land pixels, {len(nc)} site-composites, "
          f"land pixels {int(px.land.sum())} of {len(px)}")
    return y.drop(columns="ndvi_gs_n"), cnt


def main():
    sy = pd.read_csv(WORK_DIR / "veg_site_year.csv")
    hy = pd.read_csv(WORK_DIR / "hydro_site_year.csv")
    first = sy.sort_values("year").groupby("site").head(3)
    base = first.groupby("site").community.agg(lambda x: x.mode().iat[0]).rename("marsh_type_base")
    df = sy.merge(base, on="site").rename(columns={"community": "marsh_type_year"})
    df = df.merge(hy[["site", "year", "sal_gs_mean", "sal_gs_p95", "flood_gs_pct", "wlm_gs_mean"]],
                  on=["site", "year"], how="left")
    df = df.merge(salinity_index(hy), on=["site", "year"], how="left")
    nd, cnt = ndvi(df.site.unique())
    df = df.merge(nd, on=["site", "year"], how="left").merge(cnt, on="site", how="left")

    # drought indices for the site's climate division
    dv = pd.read_csv(WORK_DIR / "site_divisions.csv")[["site", "div"]].rename(columns={"div": "climdiv"})
    cd = pd.read_csv(WORK_DIR / "climdiv_monthly.csv")
    g = cd[cd.month.between(*GROW_MONTHS)].groupby(["div", "year"]).agg(
        pdsi_gs_mean=("pdsi", "mean"), precip_gs_in=("precip_in", "sum")).reset_index()
    s9 = cd[cd.month == GROW_MONTHS[1]][["div", "year", "spi6"]].rename(columns={"spi6": "spi6_sep"})
    g = g.merge(s9, on=["div", "year"]).rename(columns={"div": "climdiv"})
    df = df.merge(dv, on="site", how="left").merge(g, on=["climdiv", "year"], how="left")

    # storms: strongest exposure in each year, and whether it came before the survey
    se = pd.read_csv(WORK_DIR / "storm_exposure.csv", parse_dates=["closest_time"])
    se = se.sort_values(["site", "storm_year", "wind_proxy_kt"], ascending=[True, True, False])
    top = se.groupby(["site", "storm_year"]).first().reset_index()
    top = top.rename(columns={"storm_year": "year", "storm": "storm_name", "min_dist_km": "storm_min_dist_km",
                              "wind_radius_class_kt": "storm_wind_class_kt", "wind_proxy_kt": "storm_wind_proxy_kt"})
    df = df.merge(top[["site", "year", "storm_name", "storm_min_dist_km", "storm_wind_class_kt",
                       "storm_wind_proxy_kt", "closest_time"]], on=["site", "year"], how="left")
    df["storm_before_survey"] = (df.closest_time < pd.to_datetime(df.survey_date)).where(df.storm_name.notna())
    df = df.drop(columns="closest_time")
    fz = pd.read_csv(WORK_DIR / "freeze_exposure.csv")
    df = df.merge(fz[["site", "ghcn_id", "ghcn_dist_km", "freeze_tmin_c"]].rename(
        columns={"ghcn_id": "freeze_ghcn_id", "ghcn_dist_km": "freeze_ghcn_dist_km",
                 "freeze_tmin_c": "freeze2021_tmin_c"}), on="site", how="left")

    cols = ["site", "year", "survey_date", "lat", "lon", "marsh_type_base", "marsh_type_year", "n_plots",
            "cover_total", "richness_plot", "site_richness", "halo_share", "salinity_index",
            "avicennia_present", "avicennia_cover", "sal_gs_mean", "sal_gs_p95", "flood_gs_pct", "wlm_gs_mean",
            "ndvi_gs_mean", "ndvi_pixels_total", "ndvi_pixels_land", "climdiv", "pdsi_gs_mean", "precip_gs_in",
            "spi6_sep", "storm_name", "storm_min_dist_km", "storm_wind_class_kt", "storm_wind_proxy_kt",
            "storm_before_survey", "freeze_ghcn_id", "freeze_ghcn_dist_km", "freeze2021_tmin_c"]
    df = df[cols].sort_values(["site", "year"])
    df.to_csv(DATA_DIR / "crms_marsh_site_year.csv", index=False, float_format="%.4f")
    sites = df.groupby("site").agg(lat=("lat", "mean"), lon=("lon", "mean"), marsh_type_base=("marsh_type_base", "first"),
                                   first_year=("year", "min"), last_year=("year", "max"),
                                   n_years=("year", "nunique"), climdiv=("climdiv", "first")).reset_index()
    sites.to_csv(DATA_DIR / "crms_veg_sites.csv", index=False, float_format="%.6f")
    gpd.GeoDataFrame(sites, geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=4326) \
        .to_file(DATA_DIR / "crms_veg_sites.geojson", driver="GeoJSON")
    print(f"site-year rows {len(df)}, sites {df.site.nunique()}, years {df.year.min()}-{df.year.max()}")
    print("non-missing:", df[["sal_gs_mean", "flood_gs_pct", "salinity_index", "ndvi_gs_mean", "pdsi_gs_mean"]]
          .notna().sum().to_dict())
    print(sites.marsh_type_base.value_counts().to_dict())


if __name__ == "__main__":
    main()
