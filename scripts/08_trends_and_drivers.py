"""Step 8: vegetation change by marsh type, marsh-type transitions, hydrologic drivers, and
satellite greenness versus field cover.

Models are linear mixed-effects models (statsmodels MixedLM, REML) with a random intercept
for CRMS site. All results are associations; nothing here identifies causes.
  Trends:   y ~ year (per decade) + (1 | site), fitted separately for each baseline marsh type
  Drivers:  y ~ salinity_within + salinity_between + flooding_within + flooding_between
                + marsh type + (1 | site)
            "within" = site-year value minus the site mean (year-to-year association);
            "between" = site mean (among-site association). Predictors are in SD units.
  NDVI vs cover: Pearson r across site-years and for within-site deviations, by marsh type.
Outputs (outputs/tables): trends_by_marsh_type.csv, marsh_type_transitions.csv,
  driver_models.csv, ndvi_vs_cover.csv
Run: python -I scripts/08_trends_and_drivers.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_DIR, TAB_DIR, MARSH_TYPES  # noqa: E402

warnings.filterwarnings("ignore")
RESP = {"cover_total": "Total cover (%)", "richness_plot": "Species per 4 m2 plot",
        "site_richness": "Species per site", "halo_share": "Halophyte share of cover",
        "salinity_index": "Salinity index (ppt)", "ndvi_gs_mean": "Growing-season NDVI"}


def mixed(formula, data, group="site"):
    m = smf.mixedlm(formula, data, groups=data[group]).fit(reml=True, method=["lbfgs", "powell"])
    return m


def trends(df):
    rows = []
    for mt in MARSH_TYPES:
        d0 = df[df.marsh_type_base == mt]
        for y in RESP:
            d = d0.dropna(subset=[y]).copy()
            d["dec"] = (d.year - 2006) / 10.0
            m = mixed(f"{y} ~ dec", d)
            b, se = m.params["dec"], m.bse["dec"]
            first = d[d.year <= 2008][y].mean()
            last = d[d.year >= 2023][y].mean()
            rows.append(dict(marsh_type=mt, response=y, n_sites=d.site.nunique(), n_obs=len(d),
                             slope_per_decade=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se,
                             p=m.pvalues["dec"], mean_2006_2008=first, mean_2023_2025=last))
    t = pd.DataFrame(rows)
    t.to_csv(TAB_DIR / "trends_by_marsh_type.csv", index=False, float_format="%.4g")
    return t


def transitions(df):
    order = {m: i for i, m in enumerate(MARSH_TYPES)}
    last = df.sort_values("year").groupby("site").tail(3)
    recent = last.groupby("site").marsh_type_year.agg(lambda x: x.mode().iat[0]).rename("marsh_type_recent")
    s = df.groupby("site").agg(marsh_type_base=("marsh_type_base", "first"), last_year=("year", "max")).join(recent)
    s = s[s.last_year >= 2023]
    s["shift"] = s.marsh_type_recent.map(order) - s.marsh_type_base.map(order)
    tab = pd.crosstab(s.marsh_type_base, s.marsh_type_recent).reindex(index=MARSH_TYPES, columns=MARSH_TYPES, fill_value=0)
    tab["n_sites"] = tab.sum(axis=1)
    tab["pct_saltier"] = [100 * (s[s.marsh_type_base == m]["shift"] > 0).mean() for m in tab.index]
    tab["pct_fresher"] = [100 * (s[s.marsh_type_base == m]["shift"] < 0).mean() for m in tab.index]
    tab.to_csv(TAB_DIR / "marsh_type_transitions.csv", float_format="%.1f")
    return tab, s


def drivers(df):
    rows = []
    d = df.dropna(subset=["sal_gs_mean", "flood_gs_pct"]).copy()
    for v in ("sal_gs_mean", "flood_gs_pct"):
        mu = d.groupby("site")[v].transform("mean")
        d[f"{v}_w"] = d[v] - mu
        d[f"{v}_b"] = mu
    for v in ("sal_gs_mean_w", "sal_gs_mean_b", "flood_gs_pct_w", "flood_gs_pct_b"):
        d[v + "z"] = (d[v] - d[v].mean()) / d[v].std()
    sd = {v: d[v].std() for v in ("sal_gs_mean_w", "sal_gs_mean_b", "flood_gs_pct_w", "flood_gs_pct_b")}
    for y in ("cover_total", "richness_plot", "halo_share", "ndvi_gs_mean"):
        dd = d.dropna(subset=[y])
        m = mixed(f"{y} ~ sal_gs_mean_wz + sal_gs_mean_bz + flood_gs_pct_wz + flood_gs_pct_bz + C(marsh_type_base)", dd)
        for term in ("sal_gs_mean_wz", "sal_gs_mean_bz", "flood_gs_pct_wz", "flood_gs_pct_bz"):
            b, se = m.params[term], m.bse[term]
            rows.append(dict(response=y, term=term.replace("z", ""), sd_of_predictor=sd[term[:-1]],
                             estimate_per_sd=b, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se,
                             p=m.pvalues[term], n_sites=dd.site.nunique(), n_obs=len(dd)))
    # marsh-type-specific within-site salinity association with total cover
    for mt in MARSH_TYPES:
        dd = d[d.marsh_type_base == mt].dropna(subset=["cover_total"])
        m = mixed("cover_total ~ sal_gs_mean_w + flood_gs_pct_w", dd)
        for term in ("sal_gs_mean_w", "flood_gs_pct_w"):
            b, se = m.params[term], m.bse[term]
            rows.append(dict(response=f"cover_total [{mt}]", term=term + " (per unit: ppt or % time)",
                             sd_of_predictor=np.nan, estimate_per_sd=b, ci_low=b - 1.96 * se,
                             ci_high=b + 1.96 * se, p=m.pvalues[term], n_sites=dd.site.nunique(), n_obs=len(dd)))
    t = pd.DataFrame(rows)
    t.to_csv(TAB_DIR / "driver_models.csv", index=False, float_format="%.4g")
    return t


def ndvi_cover(df):
    rows = []
    d = df.dropna(subset=["ndvi_gs_mean", "cover_total"]).copy()
    d["ndvi_w"] = d.ndvi_gs_mean - d.groupby("site").ndvi_gs_mean.transform("mean")
    d["cov_w"] = d.cover_total - d.groupby("site").cover_total.transform("mean")
    for mt in MARSH_TYPES + ["All"]:
        x = d if mt == "All" else d[d.marsh_type_base == mt]
        r, p = stats.pearsonr(x.ndvi_gs_mean, x.cover_total)
        rw, pw = stats.pearsonr(x.ndvi_w, x.cov_w)
        sm = x.groupby("site")[["ndvi_gs_mean", "cover_total"]].mean()
        rb, pb = stats.pearsonr(sm.ndvi_gs_mean, sm.cover_total)
        rows.append(dict(marsh_type=mt, n_sites=x.site.nunique(), n_site_years=len(x), r_all=r, p_all=p,
                         r_within_site=rw, p_within=pw, r_between_site_means=rb, p_between=pb,
                         mean_ndvi=x.ndvi_gs_mean.mean()))
    t = pd.DataFrame(rows)
    t.to_csv(TAB_DIR / "ndvi_vs_cover.csv", index=False, float_format="%.4g")
    return t


def main():
    df = pd.read_csv(DATA_DIR / "crms_marsh_site_year.csv")
    pd.set_option("display.width", 200)
    t = trends(df)
    print(t[["marsh_type", "response", "n_sites", "slope_per_decade", "ci_low", "ci_high", "p",
             "mean_2006_2008", "mean_2023_2025"]].round(4).to_string(index=False))
    tab, s = transitions(df)
    print(tab.round(1).to_string())
    print(drivers(df).round(4).to_string(index=False))
    print(ndvi_cover(df).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
