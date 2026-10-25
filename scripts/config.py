"""Shared paths, sources and study settings for the pipeline.

Raw downloads go to RAW_DIR (default: <repo>/raw_data, ignored by git).
Set the environment variable CMV_RAW_DIR to keep large downloads elsewhere.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RAW_DIR = Path(os.environ.get("CMV_RAW_DIR", REPO / "raw_data"))
WORK_DIR = RAW_DIR / "work"            # large intermediate files (not committed)
DATA_DIR = REPO / "data"               # small cleaned files (committed)
FIG_DIR = REPO / "outputs" / "figures"
TAB_DIR = REPO / "outputs" / "tables"
META_DIR = REPO / "metadata"

for d in (RAW_DIR, WORK_DIR, DATA_DIR, FIG_DIR, TAB_DIR, META_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---------------- public, keyless sources ----------------
CIMS = "https://cims.coastal.la.gov/RequestedDownloads/ZippedFiles"
CIMS_GIS = "https://cims.coastal.louisiana.gov/viewer/metadata/zips"
CLIMDIV = "https://www.ncei.noaa.gov/monitoring-content/data/us/climdiv/monthly/current/"
GHCN = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station"
HURDAT_DIR = "https://www.nhc.noaa.gov/data/hurdat/"
FGDC_XSD = "https://www.fgdc.gov/schemas/metadata/"

SOURCES = {
    "crms_veg": (f"{CIMS}/CRMS_Marsh_Vegetation.zip", "crms/CRMS_Marsh_Vegetation.zip"),
    # about 1.3 GB zipped, 15 GB CSV; streamed by step 2, never extracted
    "crms_hydro": (f"{CIMS}/CRMS_Continuous_Hydrographic.zip", "crms/CRMS_Continuous_Hydrographic.zip"),
    "basins": (f"{CIMS_GIS}/Basins.zip", "gis/Basins.zip"),
    "climdiv_shp": (f"{CLIMDIV}CONUS_CLIMATE_DIVISIONS.shp.zip", "climdiv/CONUS_CLIMATE_DIVISIONS.shp.zip"),
}

# GHCN-Daily first-order stations used for the February 2021 freeze
GHCN_STATIONS = {
    "USW00003937": "Lake Charles Regional AP",
    "USW00013976": "Lafayette Regional AP",
    "USC00163433": "Galliano (cooperative station)",
    "USW00012916": "New Orleans Intl AP (MSY)",
    "USW00053865": "Slidell AP",
}

# nClimDiv Louisiana coastal climate divisions (state code 16)
CLIMDIV_STATE = "16"
CLIMDIV_DIVS = {"16 06": "East Central", "16 07": "Southwest", "16 08": "South Central", "16 09": "Southeast"}

# ---------------- study design ----------------
VEG_YEARS = (2006, 2025)          # 2026 surveys were incomplete at download time
MARSH_TYPES = ["Freshwater", "Intermediate", "Brackish", "Saline"]
GROW_MONTHS = (4, 9)              # April to September growing season for hydrology and NDVI
MIN_GROW_DAYS = 90                # minimum valid days per growing season for a hydro summary
BUFFER_M = 1000                   # station buffer radius for MODIS NDVI

# Named storms (HURDAT2 names and years); landfall dates are taken from the track file
STORMS = [("GUSTAV", 2008), ("IKE", 2008), ("ISAAC", 2012), ("LAURA", 2020), ("DELTA", 2020),
          ("ZETA", 2020), ("IDA", 2021), ("FRANCINE", 2024)]
STORM_RADIUS_KM = 400             # track points within this distance are used for exposure
FREEZE_START, FREEZE_END = "2021-02-14", "2021-02-20"
DROUGHT_YEARS = (2022, 2023)

# Halophytes typical of Louisiana saline marsh (genus or species prefix match on the
# 'Scientific Name As Currently Recognized' field). Used for the salt-tolerant cover share.
STRICT_HALOPHYTES = ["Spartina alterniflora", "Sporobolus alterniflorus", "Distichlis spicata",
                     "Juncus roemerianus", "Avicennia germinans", "Batis maritima", "Salicornia",
                     "Sarcocornia", "Lycium carolinianum", "Borrichia frutescens", "Sesuvium",
                     "Limonium", "Suaeda"]
