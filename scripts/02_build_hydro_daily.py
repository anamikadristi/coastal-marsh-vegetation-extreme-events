"""Step 2: reduce the 15 GB CRMS hourly hydrographic CSV to station-day values.

Streams the CSV out of the zip with `unzip -p` into awk (scripts/02_hydro_daily.awk), which is
far faster than reading it in Python. Requires unzip and awk (macOS, Linux, or WSL).
Daily values: mean salinity (ppt), mean water elevation relative to the marsh surface (ft),
and the fraction of valid hours with water above the marsh surface (flooded fraction).
A daily value needs at least 12 valid hourly readings.
Output: RAW_DIR/work/crms_hydro_daily.csv.gz
Run: python -I scripts/02_build_hydro_daily.py
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RAW_DIR, WORK_DIR  # noqa: E402

ZIP = RAW_DIR / "crms" / "CRMS_Continuous_Hydrographic.zip"
AWK = Path(__file__).resolve().parent / "02_hydro_daily.awk"
RAW_OUT = WORK_DIR / "crms_daily_flood.csv"
OUT = WORK_DIR / "crms_hydro_daily.csv.gz"


def main():
    if not RAW_OUT.exists():
        cmd = (f'unzip -p "{ZIP}" CRMS_Continuous_Hydrographic.csv | '
               f'LC_ALL=C awk -f "{AWK}" > "{RAW_OUT}.part" && mv "{RAW_OUT}.part" "{RAW_OUT}"')
        print("streaming hourly CSV (several minutes)...")
        subprocess.run(cmd, shell=True, check=True)
    d = pd.read_csv(RAW_OUT)
    d["date"] = pd.to_datetime(d["date"], format="%m/%d/%Y")
    g = d.groupby(["station", "date"], as_index=False)[["sal_sum", "sal_n", "wlm_sum", "wlm_n", "flood_h"]].sum()
    g["sal"] = (g.sal_sum / g.sal_n).where(g.sal_n >= 12)
    g["wl_marsh_ft"] = (g.wlm_sum / g.wlm_n).where(g.wlm_n >= 12)
    g["flood_frac"] = (g.flood_h / g.wlm_n).where(g.wlm_n >= 12)
    g["site"] = g.station.str[:8]
    out = g[["site", "station", "date", "sal", "wl_marsh_ft", "flood_frac"]]
    out = out.dropna(subset=["sal", "wl_marsh_ft"], how="all")
    out.to_csv(OUT, index=False, compression="gzip", float_format="%.4f")
    print(f"hydro stations: {out.station.nunique()}, sites: {out.site.nunique()}, station-days: {len(out)}, "
          f"dates {out.date.min().date()} to {out.date.max().date()}")
    print(f"station-days with salinity: {out.sal.notna().sum()}, with flooding: {out.flood_frac.notna().sum()}")


if __name__ == "__main__":
    main()
