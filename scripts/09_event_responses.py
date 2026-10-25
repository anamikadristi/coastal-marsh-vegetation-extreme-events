"""Step 9: vegetation and greenness responses to hurricanes, the February 2021 freeze and the
2022-2023 drought, plus resistance and recovery metrics.

Hurricane seasons are treated as events: 2008 (Gustav, Ike), 2012 (Isaac), 2020 (Laura,
Delta, Zeta), 2021 (Ida), 2024 (Francine). For each site and event:
  exposure   highest HURDAT2 wind-radius class reached at the site (0 = outside 34 kt,
             34, 50, 64 kt) and the largest wind proxy among the season's storms
  pre        the last vegetation survey before the first storm's closest approach (<= 400 days)
  post       the first survey after the last storm's closest approach (<= 400 days)
  delta      post minus pre (total cover, halophyte share, richness)
Before-after-control-impact (BACI) style model, pooled over events:
  delta ~ C(event) + C(exposure class) * marsh type (sum contrasts not used; reference =
          outside 34 kt) + (1 | site)
so the exposure coefficients are differences in change relative to unexposed sites in the
same event. A second model uses log(closest-approach distance). Event-by-event OLS models of
change on the wind proxy (per 10 kt) plus marsh type are also reported.
Resistance = post / pre and recovery = value two surveys after the event / pre (ratios of
total cover and of growing-season NDVI), summarised for sites inside the 64 kt radius.
NDVI anomalies: each 16-day composite minus the site's 2006-2025 median for that composite
day of year, in NDVI units.
Freeze: black mangrove (Avicennia germinans) presence before (2018-2020) and after
(2021-2025); spring 2021 NDVI anomaly versus the nearest station minimum temperature, with and
without adjustment for the site's 2020 hurricane wind-radius class.
Drought: monthly salinity anomalies (2022-2023 versus the site's 2006-2021 monthly means)
and change in cover from 2021 to 2023 versus the 2023 growing-season salinity anomaly.
All results are associations.
Outputs (outputs/tables): hurricane_before_after.csv, hurricane_baci_models.csv,
  hurricane_event_models.csv,
  resistance_recovery.csv, ndvi_event_anomalies.csv, freeze_mangrove.csv,
  freeze_ndvi_spring2021.csv, drought_salinity_anomaly.csv, drought_vegetation_models.csv,
  and data/hurricane_site_events.csv
Run: python -I scripts/09_event_responses.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_DIR, TAB_DIR, WORK_DIR, MARSH_TYPES  # noqa: E402

warnings.filterwarnings("ignore")
EVENTS = {2008: ["GUSTAV", "IKE"], 2012: ["ISAAC"], 2020: ["LAURA", "DELTA", "ZETA"], 2021: ["IDA"],
          2024: ["FRANCINE"]}
VARS = ["cover_total", "halo_share", "richness_plot"]


def mixed(formula, data):
    return smf.mixedlm(formula, data, groups=data["site"]).fit(reml=True, method=["lbfgs", "powell"])


def coef_rows(m, label, n):
    out = []
    for k in m.params.index:
        if k in ("Group Var",) or k.startswith("C(event)") or k == "Intercept":
            continue
        b, se = m.params[k], m.bse[k]
        out.append(dict(model=label, term=k, estimate=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se,
                        p=m.pvalues[k], n_obs=n))
    return out


def ndvi_anomalies():
    nc = pd.read_csv(WORK_DIR / "ndvi_site_composite.csv.gz", parse_dates=["date"])
    nc["doy"] = nc.date.dt.day_of_year
    clim = nc.groupby(["site", "doy"]).ndvi.median().rename("clim")
    nc = nc.join(clim, on=["site", "doy"])
    nc["anom"] = nc.ndvi - nc.clim
    return nc


def hurricanes(df, nc):
    se = pd.read_csv(WORK_DIR / "storm_exposure.csv", parse_dates=["closest_time"])
    veg = df[["site", "year", "survey_date", "marsh_type_base", "ndvi_gs_mean"] + VARS].copy()
    veg["survey_date"] = pd.to_datetime(veg.survey_date)
    veg = veg.sort_values(["site", "survey_date"])
    rows = []
    for ev, names in EVENTS.items():
        s = se[se.storm.isin(names) & (se.storm_year == ev)]
        g = s.groupby("site").agg(t0=("closest_time", "min"), t1=("closest_time", "max"),
                                  wind_class=("wind_radius_class_kt", "max"),
                                  wind_proxy=("wind_proxy_kt", "max"), min_dist_km=("min_dist_km", "min"))
        for site, r in g.iterrows():
            v = veg[veg.site == site]
            pre = v[(v.survey_date < r.t0) & (v.survey_date >= r.t0 - pd.Timedelta(days=400))].tail(1)
            post_all = v[v.survey_date > r.t1]
            post = post_all[post_all.survey_date <= r.t1 + pd.Timedelta(days=400)].head(1)
            if pre.empty or post.empty:
                continue
            later = post_all.iloc[2:3]          # two surveys after the first post-event survey
            rec = dict(site=site, event=ev, storms="/".join(names), marsh_type=pre.marsh_type_base.iat[0],
                       wind_class=int(r.wind_class), wind_proxy=r.wind_proxy, min_dist_km=r.min_dist_km,
                       pre_date=pre.survey_date.iat[0].date(), post_date=post.survey_date.iat[0].date())
            for x in VARS + ["ndvi_gs_mean"]:
                rec[f"{x}_pre"] = pre[x].iat[0]
                rec[f"{x}_post"] = post[x].iat[0]
                rec[f"{x}_rec"] = later[x].iat[0] if len(later) else np.nan
            # NDVI anomaly in the three composites after the last storm
            a = nc[(nc.site == site) & (nc.date > r.t1) & (nc.date <= r.t1 + pd.Timedelta(days=48))]
            rec["ndvi_anom_post48d"] = a.anom.mean() if len(a) else np.nan
            rows.append(rec)
    h = pd.DataFrame(rows)
    for x in VARS + ["ndvi_gs_mean"]:
        h[f"d_{x}"] = h[f"{x}_post"] - h[f"{x}_pre"]
    h["exposure"] = pd.Categorical(h.wind_class.map({0: "lt34", 34: "kt34", 50: "kt50", 64: "kt64"}),
                                   ["lt34", "kt34", "kt50", "kt64"])
    h["log_dist"] = np.log(h.min_dist_km.clip(lower=1))
    h.to_csv(DATA_DIR / "hurricane_site_events.csv", index=False, float_format="%.4f")

    summ = h.groupby(["event", "exposure"], observed=True).agg(
        n_sites=("site", "size"), d_cover_mean=("d_cover_total", "mean"), d_cover_median=("d_cover_total", "median"),
        d_halo_share_mean=("d_halo_share", "mean"), d_richness_mean=("d_richness_plot", "mean"),
        ndvi_anom_post48d=("ndvi_anom_post48d", "mean")).reset_index()
    summ.to_csv(TAB_DIR / "hurricane_before_after.csv", index=False, float_format="%.4f")

    out = []
    for y in ("d_cover_total", "d_halo_share", "d_richness_plot", "ndvi_anom_post48d"):
        d = h.dropna(subset=[y])
        m = mixed(f"{y} ~ C(event) + C(exposure)", d)
        out += coef_rows(m, f"{y} ~ event + exposure class", len(d))
        form = "log_dist * C(marsh_type)" if y in ("d_cover_total", "ndvi_anom_post48d") else "log_dist + C(marsh_type)"
        m = mixed(f"{y} ~ C(event) + {form}", d)
        out += coef_rows(m, f"{y} ~ event + {form.replace('C(marsh_type)', 'marsh type').replace('log_dist', 'log(distance)')}", len(d))
    # same exposure-class model without 2021, whose unexposed sites were recovering from the 2020 storms
    d = h[(h.event != 2021)].dropna(subset=["d_cover_total"])
    m = mixed("d_cover_total ~ C(event) + C(exposure)", d)
    out += coef_rows(m, "d_cover_total ~ event + exposure class, 2021 excluded", len(d))
    pd.DataFrame(out).to_csv(TAB_DIR / "hurricane_baci_models.csv", index=False, float_format="%.4g")

    # event-by-event: one observation per site, OLS on the wind proxy (per 10 kt) and marsh type
    ev_rows = []
    for ev, g in h.groupby("event"):
        for y in ("d_cover_total", "ndvi_anom_post48d"):
            d = g.dropna(subset=[y]).assign(wp10=lambda x: x.wind_proxy / 10)
            m = smf.ols(f"{y} ~ wp10 + C(marsh_type)", d).fit()
            ev_rows.append(dict(event=ev, storms=g.storms.iat[0], response=y, n_sites=len(d),
                                slope_per_10kt=m.params["wp10"], ci_low=m.conf_int().loc["wp10", 0],
                                ci_high=m.conf_int().loc["wp10", 1], p=m.pvalues["wp10"],
                                wind_proxy_range=f"{d.wind_proxy.min():.0f}-{d.wind_proxy.max():.0f}"))
    evm = pd.DataFrame(ev_rows)
    evm.to_csv(TAB_DIR / "hurricane_event_models.csv", index=False, float_format="%.4g")
    print(evm.round(4).to_string(index=False))

    # resistance and recovery
    rr = []
    for (mt, ex), g in h.groupby(["marsh_type", "exposure"], observed=True):
        for x in ("cover_total", "ndvi_gs_mean"):
            ok = g[g[f"{x}_pre"] > 0]
            res = ok[f"{x}_post"] / ok[f"{x}_pre"]
            rec = (ok[f"{x}_rec"] / ok[f"{x}_pre"]).dropna()
            rr.append(dict(marsh_type=mt, exposure=ex, variable=x, n_site_events=len(res),
                           resistance_median=res.median(), resistance_q25=res.quantile(.25),
                           resistance_q75=res.quantile(.75), n_recovery=len(rec),
                           recovery_median=rec.median(), recovery_q25=rec.quantile(.25),
                           recovery_q75=rec.quantile(.75)))
    rr = pd.DataFrame(rr)
    rr.to_csv(TAB_DIR / "resistance_recovery.csv", index=False, float_format="%.4f")
    return h, summ, pd.DataFrame(out), rr


def freeze(df, nc):
    av = df.groupby(["year"]).agg(sites_surveyed=("site", "nunique"),
                                  sites_with_avicennia=("avicennia_present", "sum"),
                                  mean_inplot_cover_at_avicennia_sites=("avicennia_cover",
                                                                        lambda x: x[df.loc[x.index, "avicennia_present"]].mean()))
    av = av.reset_index()
    pre_sites = set(df[(df.year.between(2018, 2020)) & df.avicennia_present].site)
    post = df[df.site.isin(pre_sites) & df.year.between(2021, 2025)]
    keep = post.groupby("year").agg(pre_sites_surveyed=("site", "nunique"),
                                    pre_sites_still_present=("avicennia_present", "sum")).reset_index()
    av = av.merge(keep, on="year", how="left")
    av.to_csv(TAB_DIR / "freeze_mangrove.csv", index=False, float_format="%.3f")

    sp = nc[(nc.date >= "2021-02-18") & (nc.date <= "2021-05-31")].groupby("site").anom.mean().rename("ndvi_anom_spring2021")
    ref = nc[(nc.date.dt.month.between(3, 5)) & (nc.date.dt.year.between(2006, 2020))]
    ref = ref.groupby([ref.site, ref.date.dt.year]).anom.mean().groupby("site").std().rename("sd_spring_anom")
    s = df.groupby("site").agg(marsh_type=("marsh_type_base", "first"), tmin=("freeze2021_tmin_c", "first"),
                               ghcn=("freeze_ghcn_id", "first")).join(sp).join(ref).dropna(subset=["ndvi_anom_spring2021"])
    rows = []
    for mt in MARSH_TYPES + ["All"]:
        x = s if mt == "All" else s[s.marsh_type == mt]
        r, p = stats.spearmanr(x.tmin, x.ndvi_anom_spring2021)
        rows.append(dict(marsh_type=mt, n_sites=len(x), mean_ndvi_anom=x.ndvi_anom_spring2021.mean(),
                         sd_ndvi_anom=x.ndvi_anom_spring2021.std(),
                         spearman_r_vs_station_tmin=r, p=p))
    by_stn = s.groupby("ghcn").agg(tmin=("tmin", "first"), n_sites=("tmin", "size"),
                                   mean_ndvi_anom=("ndvi_anom_spring2021", "mean")).reset_index()
    # the coldest stations (Lake Charles, Lafayette) are in the area hit by Laura and Delta in 2020,
    # so adjust for each site's 2020 wind-radius class
    se = pd.read_csv(WORK_DIR / "storm_exposure.csv")
    w20 = se[se.storm_year == 2020].groupby("site").wind_radius_class_kt.max().rename("wind2020")
    s = s.join(w20)
    m = smf.ols("ndvi_anom_spring2021 ~ tmin + C(wind2020) + C(marsh_type)", s.dropna(subset=["wind2020"])).fit()
    m0 = smf.ols("ndvi_anom_spring2021 ~ tmin + C(marsh_type)", s).fit()
    for lab, mm in (("marsh type", m0), ("marsh type and 2020 wind class", m)):
        rows.append(dict(marsh_type=f"All, OLS adjusted for {lab}", n_sites=int(mm.nobs),
                         ols_slope_per_1C_warmer=mm.params["tmin"], p=mm.pvalues["tmin"],
                         ci_low=mm.conf_int().loc["tmin", 0], ci_high=mm.conf_int().loc["tmin", 1]))
    t = pd.DataFrame(rows)
    t.to_csv(TAB_DIR / "freeze_ndvi_spring2021.csv", index=False, float_format="%.4f")
    return av, t, by_stn


def drought(df):
    hm = pd.read_csv(WORK_DIR / "hydro_site_month.csv.gz")
    hm = hm[hm.sal_days >= 15]
    base = hm[hm.year.between(2006, 2021)].groupby(["site", "month"]).sal_mean.mean().rename("sal_clim")
    hm = hm.join(base, on=["site", "month"])
    hm["sal_anom"] = hm.sal_mean - hm.sal_clim
    mt = df.groupby("site").marsh_type_base.first()
    hm["marsh_type"] = hm.site.map(mt)
    d = hm[hm.year.between(2021, 2024) & hm.marsh_type.notna()]
    tab = d.groupby(["marsh_type", "year", "month"]).agg(n_sites=("site", "nunique"),
                                                        sal_anom_mean=("sal_anom", "mean"),
                                                        sal_mean=("sal_mean", "mean")).reset_index()
    tab.to_csv(TAB_DIR / "drought_salinity_anomaly.csv", index=False, float_format="%.3f")
    peak = d[d.year.between(2022, 2023)].groupby("marsh_type").apply(
        lambda g: g.groupby(["year", "month"]).sal_anom.mean().agg(["max", "idxmax"])).reset_index()

    # vegetation change 2021 -> 2023 versus 2023 growing-season salinity anomaly
    gs = hm[hm.month.between(4, 9)]
    a23 = gs[gs.year == 2023].groupby("site").sal_anom.mean().rename("sal_anom_gs2023")
    w = df.pivot_table(index="site", columns="year", values=["cover_total", "halo_share", "ndvi_gs_mean"])
    rows, models = [], []
    v = pd.DataFrame({"marsh_type": mt}).join(a23)
    for x in ("cover_total", "halo_share", "ndvi_gs_mean"):
        for y0, y1 in ((2021, 2023), (2023, 2025)):
            if (x, y0) in w and (x, y1) in w:
                v[f"d_{x}_{y0}_{y1}"] = w[(x, y1)] - w[(x, y0)]
    v = v.dropna(subset=["sal_anom_gs2023"])
    for col in [c for c in v.columns if c.startswith("d_")]:
        for m in MARSH_TYPES + ["All"]:
            x = v.dropna(subset=[col]) if m == "All" else v[v.marsh_type == m].dropna(subset=[col])
            if len(x) < 8:
                continue
            sl = stats.linregress(x.sal_anom_gs2023, x[col])
            rows.append(dict(change=col, marsh_type=m, n_sites=len(x), mean_change=x[col].mean(),
                             ci95_half=1.96 * x[col].std() / np.sqrt(len(x)),
                             mean_sal_anom_gs2023=x.sal_anom_gs2023.mean(),
                             slope_vs_sal_anom=sl.slope, slope_p=sl.pvalue, r=sl.rvalue))
    t = pd.DataFrame(rows)
    t.to_csv(TAB_DIR / "drought_vegetation_models.csv", index=False, float_format="%.4g")
    return tab, peak, t


def main():
    pd.set_option("display.width", 220)
    df = pd.read_csv(DATA_DIR / "crms_marsh_site_year.csv")
    nc = ndvi_anomalies()
    h, summ, mods, rr = hurricanes(df, nc)
    print(f"hurricane site-events: {len(h)}")
    print(summ.round(3).to_string(index=False))
    print(mods.round(4).to_string(index=False))
    print(rr[rr.exposure.isin(["kt64", "lt34"])].round(3).to_string(index=False))
    av, ft, by_stn = freeze(df, nc)
    print(av.to_string(index=False))
    print(ft.round(4).to_string(index=False))
    print(by_stn.round(4).to_string(index=False))
    tab, peak, dt = drought(df)
    print(peak)
    print(dt.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
