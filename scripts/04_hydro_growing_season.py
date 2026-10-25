"""Step 4: CRMS hydrology summaries by site and growing season (April to September).

From the station-day file of step 2, per site and year (stations at a site are averaged by day):
  sal_gs_mean   mean daily salinity (ppt)
  sal_gs_p95    95th percentile of daily salinity (ppt)
  flood_gs_pct  percent of valid hours with water above the marsh surface
  wlm_gs_mean   mean water elevation relative to the marsh surface (ft)
A season needs MIN_GROW_DAYS valid days for each variable.
Also writes monthly site salinity (used for the 2022-2023 drought analysis).
Outputs: WORK_DIR/hydro_site_year.csv, WORK_DIR/hydro_site_month.csv.gz
Run: python -I scripts/04_hydro_growing_season.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import WORK_DIR, GROW_MONTHS, MIN_GROW_DAYS  # noqa: E402


def main():
    d = pd.read_csv(WORK_DIR / "crms_hydro_daily.csv.gz", parse_dates=["date"])
    s = d.groupby(["site", "date"], as_index=False)[["sal", "wl_marsh_ft", "flood_frac"]].mean()
    s["year"], s["month"] = s.date.dt.year, s.date.dt.month
    m = s.groupby(["site", "year", "month"]).agg(sal_mean=("sal", "mean"), sal_days=("sal", "count"),
                                                 flood_pct=("flood_frac", "mean"),
                                                 flood_days=("flood_frac", "count")).reset_index()
    m["flood_pct"] *= 100
    m.to_csv(WORK_DIR / "hydro_site_month.csv.gz", index=False, compression="gzip", float_format="%.4f")

    g = s[s.month.between(*GROW_MONTHS)]
    y = g.groupby(["site", "year"]).agg(
        sal_gs_mean=("sal", "mean"), sal_gs_p95=("sal", lambda x: x.quantile(0.95)), sal_days=("sal", "count"),
        flood_gs_pct=("flood_frac", "mean"), wlm_gs_mean=("wl_marsh_ft", "mean"),
        flood_days=("flood_frac", "count")).reset_index()
    y["flood_gs_pct"] *= 100
    for v, n in (("sal_gs_mean", "sal_days"), ("sal_gs_p95", "sal_days"),
                 ("flood_gs_pct", "flood_days"), ("wlm_gs_mean", "flood_days")):
        y.loc[y[n] < MIN_GROW_DAYS, v] = float("nan")
    y = y.dropna(subset=["sal_gs_mean", "flood_gs_pct"], how="all")
    y.to_csv(WORK_DIR / "hydro_site_year.csv", index=False, float_format="%.4f")
    print(f"site-years with growing-season hydrology: {len(y)} ({y.site.nunique()} sites, "
          f"{y.year.min()}-{y.year.max()}); salinity {y.sal_gs_mean.notna().sum()}, "
          f"flooding {y.flood_gs_pct.notna().sum()}")


if __name__ == "__main__":
    main()
