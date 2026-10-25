#!/usr/bin/env bash
# Run the full pipeline. Large downloads go to $CMV_RAW_DIR (default: ./raw_data).
# The CRMS hourly hydrographic zip is about 1.3 GB; step 2 streams it (needs unzip and awk).
set -euo pipefail
cd "$(dirname "$0")/.."
for s in 01_download_data 02_build_hydro_daily 03_vegetation_site_year 04_hydro_growing_season \
         05_event_exposure 06_modis_ndvi 07_build_site_year_dataset 08_trends_and_drivers \
         09_event_responses 10_figures 11_fgdc_metadata; do
  echo "=== $s"
  python -I "scripts/$s.py"
done
