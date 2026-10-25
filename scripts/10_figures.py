"""Step 10: maps and figures (PNG, 200 dpi) in outputs/figures.

Marsh types use a fixed categorical order: Freshwater, Intermediate, Brackish, Saline.
Run: python -I scripts/10_figures.py
"""
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_DIR, FIG_DIR, TAB_DIR, WORK_DIR, RAW_DIR, MARSH_TYPES, STORMS  # noqa: E402

COL = dict(zip(MARSH_TYPES, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
EVENT_YEARS = {2008: "Gustav/Ike", 2012: "Isaac", 2020: "Laura/Delta/Zeta", 2021: "Ida + freeze",
               2022: "drought", 2023: "drought", 2024: "Francine"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.facecolor": "white"})


def save(fig, name):
    fig.savefig(FIG_DIR / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def shade_events(ax):
    for y in EVENT_YEARS:
        ax.axvspan(y - 0.5, y + 0.5, color="#f1efe8", zorder=0, lw=0)


def fig_map(df):
    sites = pd.read_csv(DATA_DIR / "crms_veg_sites.csv")
    basins = gpd.read_file(next((RAW_DIR / "gis" / "Basins").rglob("*.shp")))
    basins = basins.to_crs(4326)
    from importlib import import_module
    ev = import_module("05_event_exposure")
    hd = ev.read_hurdat(sorted((RAW_DIR / "hurdat").glob("hurdat2-1851-*.txt"))[-1])
    stn = pd.read_csv(TAB_DIR / "freeze_stations_feb2021.csv")
    fig, ax = plt.subplots(figsize=(9, 5.2))
    basins.plot(ax=ax, color="#f4f3ef", edgecolor="#bdbcb6", lw=0.5)
    for name, yr in STORMS:
        t = hd[(hd.name == name) & (hd.time.dt.year == yr)]
        ax.plot(t.lon, t.lat, color=INK2, lw=1, alpha=0.7)
        inb = t[t.lon.between(-93.9, -89.0) & t.lat.between(29.0, 30.9)]
        if len(inb):
            p = inb.iloc[len(inb) // 2]
            ax.text(p.lon + 0.05, p.lat, f"{name.title()} {yr}", fontsize=6.5, color=INK2)
    for mt in MARSH_TYPES:
        s = sites[sites.marsh_type_base == mt]
        ax.scatter(s.lon, s.lat, s=14, color=COL[mt], edgecolor="white", lw=0.4, label=f"{mt} ({len(s)})", zorder=3)
    ax.scatter(stn.lon, stn.lat, marker="^", s=40, color=INK, label="GHCN station (freeze)", zorder=4)
    ax.set_xlim(-94.0, -88.8)
    ax.set_ylim(28.85, 30.85)
    ax.set_aspect(1 / np.cos(np.radians(29.8)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.legend(loc="lower left", fontsize=7, title="Baseline marsh type (sites)", title_fontsize=7)
    ax.set_title("CRMS herbaceous vegetation sites, basins and HURDAT2 tracks of the study storms", loc="left")
    save(fig, "fig01_site_map.png")


def fig_trends(df):
    vars_ = [("cover_total", "Total cover (%)"), ("richness_plot", "Species per 4 m2 plot"),
             ("halo_share", "Halophyte share of cover"), ("salinity_index", "Salinity index (ppt)"),
             ("sal_gs_mean", "Growing-season salinity (ppt)"), ("ndvi_gs_mean", "Growing-season NDVI")]
    fig, axs = plt.subplots(2, 3, figsize=(11, 6), sharex=True)
    for ax, (v, lab) in zip(axs.flat, vars_):
        shade_events(ax)
        for mt in MARSH_TYPES:
            g = df[df.marsh_type_base == mt].groupby("year")[v]
            m, se = g.mean(), g.sem()
            ax.fill_between(m.index, m - 1.96 * se, m + 1.96 * se, color=COL[mt], alpha=0.15, lw=0)
            ax.plot(m.index, m.values, color=COL[mt], lw=2, label=mt)
        ax.set_title(lab, loc="left")
    axs[0, 0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Site-year means by baseline marsh type (95% CI band); shaded years had a hurricane, the freeze or drought",
                 x=0.01, ha="left", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig02_trends_by_marsh_type.png")


def fig_transitions():
    t = pd.read_csv(TAB_DIR / "marsh_type_transitions.csv", index_col=0)
    m = t[MARSH_TYPES].values
    share = m / m.sum(axis=1, keepdims=True)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.imshow(share, cmap="Blues", vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{m[i, j]}\n{100 * share[i, j]:.0f}%", ha="center", va="center", fontsize=8,
                    color="white" if share[i, j] > 0.6 else INK)
    ax.set_xticks(range(4), MARSH_TYPES, rotation=20)
    ax.set_yticks(range(4), MARSH_TYPES)
    ax.set_xlabel("Recent marsh type (modal, last 3 surveys)")
    ax.set_ylabel("Baseline marsh type (first 3 surveys)")
    ax.grid(False)
    ax.set_title("CRMS marsh-type transitions (sites surveyed through 2023 or later)", loc="left", fontsize=9)
    save(fig, "fig03_marsh_type_transitions.png")


def fig_hurricane():
    h = pd.read_csv(DATA_DIR / "hurricane_site_events.csv")
    lab = {"lt34": "<34 kt", "kt34": "34-49 kt", "kt50": "50-63 kt", "kt64": ">=64 kt"}
    evs = sorted(h.event.unique())
    fig, axs = plt.subplots(1, len(evs), figsize=(12, 3.6), sharey=True)
    for ax, ev in zip(axs, evs):
        g = h[h.event == ev]
        for k, ex in enumerate(lab):
            x = g[g.exposure == ex].d_cover_total.dropna()
            if len(x) == 0:
                continue
            jit = np.random.default_rng(1).uniform(-0.18, 0.18, len(x))
            ax.scatter(k + jit, x, s=8, color="#9fc1ea", lw=0)
            ax.plot([k - 0.28, k + 0.28], [x.median()] * 2, color=INK, lw=2)
            ax.text(k, ax.get_ylim()[0] if False else -98, f"n={len(x)}", ha="center", fontsize=6.5, color=INK2)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xticks(range(4), list(lab.values()), rotation=30, fontsize=7)
        ax.set_title(f"{ev}: {g.storms.iat[0].title()}", loc="left", fontsize=8.5)
    axs[0].set_ylabel("Change in total cover (percentage points)\npost minus pre-storm survey")
    axs[0].set_ylim(-100, 80)
    fig.suptitle("Plot cover change across hurricane seasons by HURDAT2 wind-radius class at the site (line = median)",
                 x=0.01, ha="left", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig04_hurricane_cover_change.png")


def fig_resist():
    rr = pd.read_csv(TAB_DIR / "resistance_recovery.csv")
    rr = rr[rr.exposure.isin(["lt34", "kt64"])]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, var, lab in zip(axs, ["cover_total", "ndvi_gs_mean"], ["Total cover", "Growing-season NDVI"]):
        for k, mt in enumerate(MARSH_TYPES):
            for j, (ex, mk) in enumerate([("lt34", "o"), ("kt64", "s")]):
                r = rr[(rr.marsh_type == mt) & (rr.exposure == ex) & (rr.variable == var)]
                if r.empty:
                    continue
                r = r.iloc[0]
                x0 = k + (j - 0.5) * 0.3
                ax.errorbar(x0 - 0.05, r.resistance_median, yerr=[[r.resistance_median - r.resistance_q25],
                            [r.resistance_q75 - r.resistance_median]], fmt=mk, color=COL[mt], ms=6,
                            mfc="white" if ex == "lt34" else COL[mt], lw=1.2)
                ax.errorbar(x0 + 0.05, r.recovery_median, yerr=[[r.recovery_median - r.recovery_q25],
                            [r.recovery_q75 - r.recovery_median]], fmt=mk, color=COL[mt], ms=6, alpha=0.55,
                            mfc="white" if ex == "lt34" else COL[mt], lw=1.2)
        ax.axhline(1, color=INK2, lw=0.8)
        ax.set_xticks(range(4), MARSH_TYPES)
        ax.set_title(f"{lab}: resistance (post/pre, solid) and recovery (two surveys later/pre, faded)",
                     loc="left", fontsize=8)
    axs[0].set_ylabel("Ratio to pre-storm value (median, IQR)")
    from matplotlib.lines import Line2D
    axs[1].legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, label="outside 34 kt radius"),
                           Line2D([], [], marker="s", ls="", color=INK2, label="inside 64 kt radius")],
                  fontsize=7, loc="lower right")
    fig.tight_layout()
    save(fig, "fig05_resistance_recovery.png")


def fig_drought():
    t = pd.read_csv(TAB_DIR / "drought_salinity_anomaly.csv")
    t["date"] = pd.to_datetime(dict(year=t.year, month=t.month, day=15))
    cd = pd.read_csv(WORK_DIR / "climdiv_monthly.csv")
    cd = cd[cd.year.between(2021, 2024)]
    cd["date"] = pd.to_datetime(dict(year=cd.year, month=cd.month, day=15))
    fig, axs = plt.subplots(2, 1, figsize=(9, 5.5), sharex=True, height_ratios=[1.6, 1])
    for mt in MARSH_TYPES:
        g = t[t.marsh_type == mt]
        axs[0].plot(g.date, g.sal_anom_mean, color=COL[mt], lw=2, label=mt)
    axs[0].axhline(0, color=INK2, lw=0.8)
    axs[0].set_ylabel("Salinity anomaly (ppt)\nvs 2006-2021 monthly mean")
    axs[0].legend(fontsize=7, ncol=4, loc="upper left")
    axs[0].set_title("CRMS monthly salinity anomaly by baseline marsh type, 2021-2024", loc="left")
    names = {6: "East Central", 7: "Southwest", 8: "South Central", 9: "Southeast"}
    for k, (dv, g) in enumerate(cd.groupby("div")):
        axs[1].plot(g.date, g.pdsi, lw=1.6, color=["#4a3aa7", "#008300", "#e87ba4", "#52514e"][k],
                    label=f"Div {dv} {names[dv]}")
    axs[1].axhline(-3, color=INK2, lw=0.8, ls="--")
    axs[1].text(cd.date.min(), -3.3, "severe drought (-3)", fontsize=6.5, color=INK2, va="top")
    axs[1].set_ylabel("Palmer Drought\nSeverity Index")
    axs[1].legend(fontsize=7, ncol=4, loc="upper left")
    axs[1].set_title("nClimDiv PDSI, Louisiana climate divisions", loc="left")
    fig.tight_layout()
    save(fig, "fig06_drought_salinity_pdsi.png")


def fig_freeze(df):
    av = pd.read_csv(TAB_DIR / "freeze_mangrove.csv")
    nc = pd.read_csv(WORK_DIR / "ndvi_site_composite.csv.gz", parse_dates=["date"])
    nc["doy"] = nc.date.dt.day_of_year
    nc = nc.join(nc.groupby(["site", "doy"]).ndvi.median().rename("clim"), on=["site", "doy"])
    sp = nc[(nc.date >= "2021-02-18") & (nc.date <= "2021-05-31")].assign(a=lambda x: x.ndvi - x.clim)
    s = sp.groupby("site").a.mean().rename("anom").to_frame().join(
        df.groupby("site").agg(tmin=("freeze2021_tmin_c", "first"), mt=("marsh_type_base", "first")))
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    ax.bar(av.year, av.sites_with_avicennia, color="#1baf7a", width=0.7)
    ax.axvline(2020.5, color=INK2, lw=0.8, ls="--")
    ax.text(2020.6, av.sites_with_avicennia.max() * 0.95, "Feb 2021\nfreeze", fontsize=7, color=INK2)
    ax.set_ylabel("Sites with black mangrove recorded")
    ax.set_xticks(range(2006, 2026, 3))
    ax.set_title("Avicennia germinans records at CRMS vegetation sites", loc="left", fontsize=9)
    ax = axs[1]
    rng = np.random.default_rng(2)
    for mt in MARSH_TYPES:
        g = s[s.mt == mt]
        ax.scatter(g.tmin + rng.uniform(-0.12, 0.12, len(g)), g.anom, s=10, color=COL[mt], label=mt, alpha=0.8, lw=0)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("Minimum TMIN, 14-20 Feb 2021, nearest GHCN station (C)")
    ax.set_ylabel("NDVI anomaly, late Feb-May 2021")
    ax.legend(fontsize=7)
    ax.set_title("Spring 2021 greenness anomaly vs freeze severity", loc="left", fontsize=9)
    fig.tight_layout()
    save(fig, "fig07_freeze_mangrove_ndvi.png")


def fig_ndvi_cover(df):
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.3), sharex=True, sharey=True)
    t = pd.read_csv(TAB_DIR / "ndvi_vs_cover.csv").set_index("marsh_type")
    for ax, mt in zip(axs, MARSH_TYPES):
        g = df[df.marsh_type_base == mt].dropna(subset=["ndvi_gs_mean", "cover_total"])
        ax.scatter(g.cover_total, g.ndvi_gs_mean, s=4, color=COL[mt], alpha=0.35, lw=0)
        ax.set_title(f"{mt}: r = {t.loc[mt, 'r_all']:.2f} (within-site {t.loc[mt, 'r_within_site']:.2f})",
                     loc="left", fontsize=8)
        ax.set_xlabel("CRMS plot total cover (%)")
    axs[0].set_ylabel("MODIS growing-season NDVI (1 km)")
    fig.tight_layout()
    save(fig, "fig08_ndvi_vs_field_cover.png")


def fig_drivers():
    d = pd.read_csv(TAB_DIR / "driver_models.csv")
    d = d[~d.response.str.contains(r"\[")]
    labels = {"sal_gs_mean_w": "Salinity, within site", "sal_gs_mean_b": "Salinity, between sites",
              "flood_gs_pct_w": "Flooding, within site", "flood_gs_pct_b": "Flooding, between sites"}
    resp = {"cover_total": "Total cover (%)", "richness_plot": "Species per plot",
            "halo_share": "Halophyte share", "ndvi_gs_mean": "NDVI"}
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
    for ax, (r, lab) in zip(axs, resp.items()):
        g = d[d.response == r].reset_index(drop=True)
        y = np.arange(len(g))[::-1]
        ax.errorbar(g.estimate_per_sd, y, xerr=[g.estimate_per_sd - g.ci_low, g.ci_high - g.estimate_per_sd],
                    fmt="o", color="#2a78d6", ms=6, lw=1.5)
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_yticks(y, [labels[t] for t in g.term])
        ax.set_title(f"{lab}, per 1 SD", loc="left", fontsize=8.5)
    fig.suptitle("Mixed-model associations of growing-season salinity and flooding with vegetation (95% CI)",
                 x=0.01, ha="left", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig09_driver_associations.png")


def fig_ndvi_anom(df):
    nc = pd.read_csv(WORK_DIR / "ndvi_site_composite.csv.gz", parse_dates=["date"])
    nc["doy"] = nc.date.dt.day_of_year
    nc = nc.join(nc.groupby(["site", "doy"]).ndvi.median().rename("clim"), on=["site", "doy"])
    nc["anom"] = nc.ndvi - nc.clim
    nc["mt"] = nc.site.map(df.groupby("site").marsh_type_base.first())
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for mt in MARSH_TYPES:
        g = nc[nc.mt == mt].groupby("date").anom.mean().rolling(3, center=True).mean()
        ax.plot(g.index, g.values, color=COL[mt], lw=1.4, label=mt)
    st = pd.read_csv(TAB_DIR / "storm_list.csv", parse_dates=["la_landfall_time"])
    for r in st.dropna(subset=["la_landfall_time"]).itertuples():
        ax.axvline(r.la_landfall_time, color=INK2, lw=0.6, ls=":")
        if r.storm in ("DELTA", "ZETA"):
            continue
        lab = "Laura, Delta, Zeta" if r.storm == "LAURA" else r.storm.title()
        ax.text(r.la_landfall_time, 0.105, lab, rotation=90, fontsize=6, color=INK2, va="top", ha="right")
    ax.axvline(pd.Timestamp("2021-02-15"), color="#4a3aa7", lw=0.8, ls="--")
    ax.text(pd.Timestamp("2021-02-15"), 0.105, "Feb 2021 freeze", fontsize=6, color="#4a3aa7", rotation=90, va="top", ha="right")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_ylim(-0.12, 0.11)
    ax.set_ylabel("NDVI anomaly (3-composite mean)")
    ax.legend(fontsize=7, ncol=4, loc="upper left", bbox_to_anchor=(0, -0.08))
    ax.set_title("MODIS NDVI anomaly at CRMS site buffers (vs each site's composite-date median, 2006-2025)", loc="left")
    save(fig, "fig10_ndvi_anomaly_timeseries.png")


def main():
    df = pd.read_csv(DATA_DIR / "crms_marsh_site_year.csv")
    fig_map(df)
    fig_trends(df)
    fig_transitions()
    fig_hurricane()
    fig_resist()
    fig_drought()
    fig_freeze(df)
    fig_ndvi_cover(df)
    fig_drivers()
    fig_ndvi_anom(df)


if __name__ == "__main__":
    main()
