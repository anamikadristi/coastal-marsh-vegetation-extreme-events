"""Step 3: CRMS herbaceous marsh vegetation, plot visits and site-year summaries.

Input: CIMS CRMS_Marsh_Vegetation.csv (4 m2 plots, Braun-Blanquet cover by species).
Kept: herbaceous marsh plots (station IDs with -V), communities Freshwater, Intermediate,
Brackish and Saline (Swamp and forested -F stations are excluded), survey years 2006-2025.
Sites with no coordinates in the export are dropped (reported when run).
Per plot visit (station and date): total cover, species richness (species recorded inside the
plot), and halophyte cover share (see config.STRICT_HALOPHYTES).
Cover codes: "<1" is coded 0.5 percent and "Solitary" 0.1 percent.
Black mangrove (Avicennia germinans) presence counts any record, including species seen
just outside the plot ("Out"), because mangroves are often near but not inside a plot.
Outputs: WORK_DIR/veg_species_cover.csv.gz (long table) and WORK_DIR/veg_plot_visits.csv.gz,
WORK_DIR/veg_site_year.csv
Run: python -I scripts/03_vegetation_site_year.py
"""
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RAW_DIR, WORK_DIR, VEG_YEARS, MARSH_TYPES, STRICT_HALOPHYTES  # noqa: E402

ZIP = RAW_DIR / "crms" / "CRMS_Marsh_Vegetation.zip"
NOT_SPECIES = {"Bare Ground", "Water", "Dead Vegetation", "Wrack"}


def cover_value(s):
    s = str(s).strip()
    if s in ("nan", ""):
        return np.nan
    if s.startswith("<"):
        return 0.5
    if s.lower() == "solitary":
        return 0.1
    try:
        return float(s)
    except ValueError:
        return np.nan


def is_halophyte(name):
    return any(str(name).startswith(h) for h in STRICT_HALOPHYTES)


def main():
    with zipfile.ZipFile(ZIP) as zf:
        with zf.open("CRMS_Marsh_Vegetation.csv") as f:
            d = pd.read_csv(f, encoding="latin-1", low_memory=False)
    n_raw = len(d)
    d = d.rename(columns={"Station ID": "station", "Community": "community",
                          "Scientific Name As Currently Recognized": "species", "In/Out": "inout",
                          "% Cover Total": "cover_total", "% Cover": "cover",
                          "Latitude": "lat", "Longitude": "lon"})
    d["date"] = pd.to_datetime(d["Collection Date (mm/dd/yyyy)"], format="%m/%d/%Y", errors="coerce")
    d["year"] = d.date.dt.year
    d["site"] = d.station.str[:8]
    d = d[d.station.str.contains("-V", na=False) & d.community.isin(MARSH_TYPES)
          & d.year.between(*VEG_YEARS)].copy()
    print(f"raw rows {n_raw}; herbaceous marsh rows {len(d)}")

    # mangrove presence per plot visit (in plot or just outside)
    d["avicennia_any"] = d.species.str.startswith("Avicennia", na=False)
    d["cover_num"] = d.cover.map(cover_value)
    d["cover_total_num"] = pd.to_numeric(d.cover_total, errors="coerce")
    inplot = d[d.inout.isin(["Both", "In"]) & ~d.species.isin(NOT_SPECIES) & d.species.notna()].copy()
    inplot["cover_num"] = inplot.cover_num.fillna(0.1)
    inplot["halophyte"] = inplot.species.map(is_halophyte)

    inplot["halo_c"] = inplot.cover_num * inplot.halophyte
    sp = inplot[["site", "station", "date", "year", "community", "species", "cover_num", "halophyte", "halo_c"]]
    sp.to_csv(WORK_DIR / "veg_species_cover.csv.gz", index=False, compression="gzip")

    keys = ["site", "station", "date", "year"]
    visits = d.groupby(keys, as_index=False).agg(
        community=("community", "first"), cover_total=("cover_total_num", "max"),
        lat=("lat", "first"), lon=("lon", "first"), avicennia_any=("avicennia_any", "any"))
    s = sp.groupby(keys).agg(richness=("species", "nunique"), sum_cover=("cover_num", "sum"),
                             halo_cover=("halo_c", "sum"))
    av = sp[sp.species.str.startswith("Avicennia")].groupby(keys).cover_num.sum().rename("avicennia_cover")
    visits = visits.merge(s, on=keys, how="left").merge(av, on=keys, how="left")
    visits[["richness", "sum_cover", "halo_cover", "avicennia_cover"]] = \
        visits[["richness", "sum_cover", "halo_cover", "avicennia_cover"]].fillna(0)
    visits["halo_share"] = (visits.halo_cover / visits.sum_cover).where(visits.sum_cover > 0)
    visits.to_csv(WORK_DIR / "veg_plot_visits.csv.gz", index=False, compression="gzip")

    # site-year summary
    site_rich = sp.groupby(["site", "year"]).species.nunique().rename("site_richness")
    sy = visits.groupby(["site", "year"]).agg(
        n_plots=("station", "nunique"), survey_date=("date", "median"),
        community=("community", lambda x: x.mode().iat[0]),
        cover_total=("cover_total", "mean"), richness_plot=("richness", "mean"),
        halo_share=("halo_share", "mean"), avicennia_present=("avicennia_any", "any"),
        avicennia_cover=("avicennia_cover", "mean"),
        lat=("lat", "mean"), lon=("lon", "mean")).join(site_rich).reset_index()
    sy["survey_date"] = sy.survey_date.dt.date
    nocoord = sy[sy.lat.isna() | sy.lon.isna()]
    if len(nocoord):
        print(f"dropping {nocoord.site.nunique()} site(s) without coordinates "
              f"({', '.join(nocoord.site.unique())}; {len(nocoord)} site-years)")
        sy = sy.dropna(subset=["lat", "lon"])
    sy.to_csv(WORK_DIR / "veg_site_year.csv", index=False, float_format="%.4f")
    print(f"plot visits {len(visits)}, plots {visits.station.nunique()}, sites {sy.site.nunique()}, "
          f"site-years {len(sy)}, years {sy.year.min()}-{sy.year.max()}, species {sp.species.nunique()}")
    print(sy.groupby("community").site.nunique().rename("sites by modal community (any year)"))


if __name__ == "__main__":
    main()
