"""Step 11: FGDC CSDGM (FGDC-STD-001-1998) XML metadata for data/crms_marsh_site_year.csv,
validated against the official CSDGM XML schema from fgdc.gov (downloaded in step 1).

Sections written: idinfo, dataqual (with lineage), spdoinfo, spref, eainfo, distinfo, metainfo.
Attribute ranges, counts and the bounding box are computed from the CSV at run time.
Validation uses two independent XSD validators (lxml/libxml2 and xmlschema); both results are
written to metadata/fgdc_validation_report.txt.
Run: python -I scripts/11_fgdc_metadata.py
"""
import datetime as dt
import sys
from pathlib import Path

import pandas as pd
import xmlschema
from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_DIR, META_DIR, RAW_DIR, VEG_YEARS  # noqa: E402

CSV = DATA_DIR / "crms_marsh_site_year.csv"
XML = META_DIR / "crms_marsh_site_year_fgdc.xml"
REPORT = META_DIR / "fgdc_validation_report.txt"
XSD = RAW_DIR / "schemas" / "fgdc-std-001-1998.xsd"
AUTHOR = "Anamika Dristi"

CIMS = "https://cims.coastal.la.gov"
ATTR = {
    "site": ("CRMS site identifier (first eight characters of the CRMS station ID).", "CPRA CIMS", None),
    "year": ("Vegetation survey year.", "CPRA CIMS", "year"),
    "survey_date": ("Median collection date of the site's vegetation plots in that year (YYYY-MM-DD).", "CPRA CIMS", None),
    "lat": ("Mean latitude of the site's vegetation plots, decimal degrees.", "CPRA CIMS", "decimal degrees"),
    "lon": ("Mean longitude of the site's vegetation plots, decimal degrees.", "CPRA CIMS", "decimal degrees"),
    "marsh_type_base": ("Baseline marsh type: modal CRMS community in the site's first three survey years.", "Derived from CPRA CIMS", "enum"),
    "marsh_type_year": ("Modal CRMS community recorded at the site's plots in this year.", "CPRA CIMS", "enum"),
    "n_plots": ("Number of 4 m2 herbaceous vegetation plots surveyed at the site in the year.", "Derived", "count"),
    "cover_total": ("Mean total vegetation cover of the plots.", "CPRA CIMS (% Cover Total)", "percent"),
    "richness_plot": ("Mean number of species recorded inside a plot.", "Derived", "species per plot"),
    "site_richness": ("Number of species recorded inside any plot at the site in the year.", "Derived", "species"),
    "halo_share": ("Mean plot share of summed species cover from halophytes (list in docs/METADATA_NOTES.md), 0 to 1.", "Derived", "proportion"),
    "salinity_index": ("Community-weighted mean of species salinity optima (weighted averaging on CRMS growing-season salinity).", "Derived", "ppt"),
    "avicennia_present": ("True if Avicennia germinans was recorded in or next to any plot at the site.", "Derived", "boolean"),
    "avicennia_cover": ("Mean in-plot cover of Avicennia germinans across the site's plots.", "Derived", "percent"),
    "sal_gs_mean": ("Mean daily surface-water salinity, April to September, CRMS continuous recorders at the site.", "Derived from CPRA CIMS hydrographic data", "ppt"),
    "sal_gs_p95": ("95th percentile of daily salinity, April to September.", "Derived from CPRA CIMS hydrographic data", "ppt"),
    "flood_gs_pct": ("Percent of valid hours, April to September, with water above the marsh surface.", "Derived from CPRA CIMS hydrographic data", "percent"),
    "wlm_gs_mean": ("Mean water elevation relative to the marsh surface, April to September.", "Derived from CPRA CIMS hydrographic data", "feet"),
    "ndvi_gs_mean": ("Mean April-September MODIS MOD13Q1 NDVI of good or marginal land pixels within 1 km.", "Derived from NASA LP DAAC MOD13Q1 v061 via Microsoft Planetary Computer", "NDVI"),
    "ndvi_pixels_total": ("Number of 250 m MODIS pixels whose centers fall within 1 km of the site.", "Derived", "count"),
    "ndvi_pixels_land": ("Pixels kept after removing those with all-years median June-September NDVI below 0.2.", "Derived", "count"),
    "climdiv": ("NOAA nClimDiv Louisiana climate division number joined to the site.", "NOAA NCEI nClimDiv", "division"),
    "pdsi_gs_mean": ("Mean April-September Palmer Drought Severity Index for the climate division.", "NOAA NCEI nClimDiv", "index"),
    "precip_gs_in": ("Total April-September precipitation for the climate division.", "NOAA NCEI nClimDiv", "inches"),
    "spi6_sep": ("Six-month Standardized Precipitation Index ending in September.", "NOAA NCEI nClimDiv", "index"),
    "storm_name": ("Study hurricane of that year with the largest wind proxy at the site (blank if none).", "NOAA NHC HURDAT2", None),
    "storm_min_dist_km": ("Closest approach of that storm's track to the site.", "Derived from HURDAT2", "kilometers"),
    "storm_wind_class_kt": ("Highest HURDAT2 wind-radius class (0, 34, 50, 64 kt) whose quadrant radius reached the site.", "Derived from HURDAT2", "knots"),
    "storm_wind_proxy_kt": ("Rough wind exposure index: Vmax times min(1, (40 km / distance)^0.5).", "Derived from HURDAT2", "knots"),
    "storm_before_survey": ("True if the storm's closest approach came before the year's survey date.", "Derived", "boolean"),
    "freeze_ghcn_id": ("GHCN-Daily station nearest to the site used for the February 2021 freeze.", "NOAA NCEI GHCN-Daily", None),
    "freeze_ghcn_dist_km": ("Distance from the site to that station.", "Derived", "kilometers"),
    "freeze2021_tmin_c": ("Lowest daily minimum temperature at that station, 14-20 February 2021.", "NOAA NCEI GHCN-Daily", "degrees Celsius"),
}
SOURCES = [
    ("Coastal Protection and Restoration Authority of Louisiana (CPRA)", "CRMS herbaceous marsh vegetation (CRMS_Marsh_Vegetation.zip)",
     f"{CIMS}/RequestedDownloads/ZippedFiles/CRMS_Marsh_Vegetation.zip", "CRMS_VEG", "2006", "2025",
     "Species percent cover, community and plot coordinates."),
    ("Coastal Protection and Restoration Authority of Louisiana (CPRA)", "CRMS continuous hourly hydrographic data (CRMS_Continuous_Hydrographic.zip)",
     f"{CIMS}/RequestedDownloads/ZippedFiles/CRMS_Continuous_Hydrographic.zip", "CRMS_HYDRO", "2006", "2025",
     "Hourly salinity and water elevation relative to the marsh surface."),
    ("NOAA National Hurricane Center", "Atlantic hurricane database (HURDAT2)",
     "https://www.nhc.noaa.gov/data/hurdat/", "HURDAT2", "2008", "2024", "Best-track positions, intensity and wind radii."),
    ("NOAA National Centers for Environmental Information", "Global Historical Climatology Network - Daily (GHCN-Daily)",
     "https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/", "GHCND", "2021", "2021", "Daily minimum temperature."),
    ("NOAA National Centers for Environmental Information", "nClimDiv climate division data (PDSI, SPI-6, precipitation)",
     "https://www.ncei.noaa.gov/monitoring-content/data/us/climdiv/monthly/current/", "NCLIMDIV", "2006", "2025",
     "Monthly divisional drought indices and precipitation."),
    ("NASA LP DAAC (served by Microsoft Planetary Computer)", "MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 (MOD13Q1)",
     "https://planetarycomputer.microsoft.com/dataset/modis-13Q1-061", "MOD13Q1", "2006", "2025", "16-day NDVI and pixel reliability."),
]
STEPS = [
    "Downloaded the public inputs without accounts or keys (scripts/01_download_data.py).",
    "Reduced the hourly CRMS hydrographic export to station-day salinity, water level relative to the marsh surface and flooded fraction, requiring 12 or more valid hours per day (scripts/02_build_hydro_daily.py).",
    "Selected herbaceous marsh plots (Freshwater, Intermediate, Brackish, Saline; years 2006-2025), coded cover classes (<1 as 0.5, Solitary as 0.1) and summarised plots to site-years (scripts/03_vegetation_site_year.py).",
    "Summarised hydrology by site for April-September, requiring 90 or more valid days (scripts/04_hydro_growing_season.py).",
    "Computed hurricane, freeze and drought exposure for each site from HURDAT2, GHCN-Daily and nClimDiv (scripts/05_event_exposure.py).",
    "Read MOD13Q1 NDVI windows at 1 km site buffers from Microsoft Planetary Computer and kept good or marginal pixels (scripts/06_modis_ndvi.py).",
    "Joined all inputs, derived the baseline marsh type, salinity index and growing-season NDVI, and wrote the site-year CSV (scripts/07_build_site_year_dataset.py).",
]


def sub(parent, tag, text=None):
    e = etree.SubElement(parent, tag)
    if text is not None:
        e.text = str(text)
    return e


def contact(parent):
    ci = sub(parent, "cntinfo")
    p = sub(ci, "cntperp")
    sub(p, "cntper", AUTHOR)
    sub(ci, "cntpos", "Independent researcher")
    a = sub(ci, "cntaddr")
    sub(a, "addrtype", "mailing")
    sub(a, "city", "Not provided")
    sub(a, "state", "Not provided")
    sub(a, "postal", "Not provided")
    sub(a, "country", "USA")
    sub(ci, "cntvoice", "Not provided; contact by email")
    sub(ci, "cntemail", "anamikadristi@gmail.com")


def citation(parent, origin, title, pubdate, online=None, geoform=None):
    c = sub(sub(parent, "citation") if parent.tag != "srccite" else parent, "citeinfo")
    sub(c, "origin", origin)
    sub(c, "pubdate", pubdate)
    sub(c, "title", title)
    if geoform:
        sub(c, "geoform", geoform)
    if online:
        sub(c, "onlink", online)
    return c


def build(df, today):
    root = etree.Element("metadata")
    # ---------------- idinfo ----------------
    idinfo = sub(root, "idinfo")
    citation(idinfo, AUTHOR, "CRMS herbaceous marsh vegetation, hydrology, extreme-event exposure and MODIS NDVI by site and year, coastal Louisiana, 2006-2025",
             today.strftime("%Y%m%d"), "https://github.com/anamikadristi/coastal-marsh-vegetation-extreme-events", "tabular digital data")
    d = sub(idinfo, "descript")
    sub(d, "abstract", f"One row per Coastwide Reference Monitoring System (CRMS) herbaceous marsh site and survey year "
        f"({len(df)} rows, {df.site.nunique()} sites, {df.year.min()}-{df.year.max()}). Each row combines plot vegetation "
        "summaries (total cover, species richness, halophyte cover share, a salinity index from species composition, "
        "black mangrove presence) with April-September salinity and flooding from CRMS continuous recorders, MODIS NDVI "
        "within 1 km, nClimDiv drought indices, HURDAT2 hurricane exposure and February 2021 freeze minimum temperature. "
        "It is an independent derived product built from public data.")
    sub(d, "purpose", "To study how Louisiana coastal marsh vegetation responds to hurricanes, the February 2021 freeze and the "
        "2022-2023 drought, and how vegetation relates to salinity, flooding and satellite greenness, by marsh type.")
    sub(d, "supplinf", "Associations only; the dataset supports observational analyses and does not establish causes. "
        "Code and documentation: scripts/ and docs/METADATA_NOTES.md in the repository.")
    tp = sub(sub(idinfo, "timeperd"), "timeinfo")
    rd = sub(tp, "rngdates")
    sub(rd, "begdate", str(df.survey_date.min()).replace("-", ""))
    sub(rd, "enddate", str(df.survey_date.max()).replace("-", ""))
    sub(idinfo.find("timeperd"), "current", "ground condition")
    st = sub(idinfo, "status")
    sub(st, "progress", "Complete")
    sub(st, "update", "As needed")
    b = sub(sub(idinfo, "spdom"), "bounding")
    sub(b, "westbc", f"{df.lon.min():.4f}")
    sub(b, "eastbc", f"{df.lon.max():.4f}")
    sub(b, "northbc", f"{df.lat.max():.4f}")
    sub(b, "southbc", f"{df.lat.min():.4f}")
    kw = sub(idinfo, "keywords")
    th = sub(kw, "theme")
    sub(th, "themekt", "ISO 19115 Topic Category")
    for k in ("biota", "environment", "inlandWaters", "climatologyMeteorologyAtmosphere"):
        sub(th, "themekey", k)
    th = sub(kw, "theme")
    sub(th, "themekt", "None")
    for k in ("coastal wetlands", "marsh vegetation", "CRMS", "hurricanes", "freeze", "drought", "salinity",
              "flooding duration", "NDVI", "MODIS", "resistance", "recovery"):
        sub(th, "themekey", k)
    pl = sub(kw, "place")
    sub(pl, "placekt", "None")
    for k in ("Louisiana", "northern Gulf coast", "Mississippi River Delta", "Chenier Plain"):
        sub(pl, "placekey", k)
    tm = sub(kw, "temporal")
    sub(tm, "tempkt", "None")
    for k in ("2006-2025", "2008", "2012", "2020", "2021", "2022", "2023", "2024"):
        sub(tm, "tempkey", k)
    sub(idinfo, "accconst", "None")
    sub(idinfo, "useconst", "None for this derived product. Source data keep their terms: CPRA CIMS data are provided "
        "for general informational purposes and CPRA asks that users acquire the data from CPRA; NOAA and NASA data are public domain. "
        "Cite the original sources listed in the lineage.")
    contact(sub(idinfo, "ptcontac"))
    sub(idinfo, "native", "Python 3 (pandas, statsmodels, rasterio, geopandas); scripts in the repository")
    # ---------------- dataqual ----------------
    dq = sub(root, "dataqual")
    aa = sub(dq, "attracc")
    sub(aa, "attraccr", "Vegetation values are as reported by CRMS field crews (visual cover estimates on 4 m2 plots). "
        "Hydrologic summaries depend on recorder completeness; site-years with fewer than 90 valid growing-season days are "
        "left blank. NDVI is a 250 m, 16-day composite and mixes marsh, water and other cover within 1 km. Storm exposure "
        "is from best-track positions and wind radii, not modeled or measured site winds. Freeze temperature is from the "
        "nearest of five GHCN stations, up to "
        f"{df.freeze_ghcn_dist_km.max():.0f} km away.")
    sub(dq, "logic", "Each row is unique by site and year. Percent values are within 0-100 and proportions within 0-1. "
        "Missing values are empty fields. Checked by the pipeline before writing.")
    sub(dq, "complete", f"{df.site.nunique()} sites; vegetation for all rows; salinity for {df.sal_gs_mean.notna().sum()} rows, "
        f"flooding for {df.flood_gs_pct.notna().sum()}, NDVI for {df.ndvi_gs_mean.notna().sum()}, "
        f"salinity index for {df.salinity_index.notna().sum()}. Swamp (forested) sites are excluded. "
        "A CRMS Floristic Quality Index could not be computed because the species coefficient of conservatism list was "
        "not available from a keyless source at build time.")
    pa = sub(sub(dq, "posacc"), "horizpa")
    sub(pa, "horizpar", "Coordinates are the mean of the plot coordinates in the CPRA export; positional accuracy is not stated by the source.")
    lin = sub(dq, "lineage")
    for org, title, url, abbr, b0, b1, contrib in SOURCES:
        si = sub(lin, "srcinfo")
        sc = sub(si, "srccite")
        citation(sc, org, title, "Unknown", url)
        sub(si, "typesrc", "online")
        tim = sub(si, "srctime")
        rd = sub(sub(tim, "timeinfo"), "rngdates")
        sub(rd, "begdate", b0)
        sub(rd, "enddate", b1)
        sub(tim, "srccurr", "ground condition")
        sub(si, "srccitea", abbr)
        sub(si, "srccontr", contrib)
    for i, s in enumerate(STEPS):
        ps = sub(lin, "procstep")
        sub(ps, "procdesc", s)
        sub(ps, "procdate", today.strftime("%Y%m%d"))
    # ---------------- spdoinfo ----------------
    sp = sub(root, "spdoinfo")
    sub(sp, "direct", "Point")
    sd = sub(sub(sp, "ptvctinf"), "sdtsterm")
    sub(sd, "sdtstype", "Entity point")
    sub(sd, "ptvctcnt", str(len(df)))
    # ---------------- spref ----------------
    hs = sub(sub(root, "spref"), "horizsys")
    g = sub(hs, "geograph")
    sub(g, "latres", "0.0001")
    sub(g, "longres", "0.0001")
    sub(g, "geogunit", "Decimal degrees")
    # ---------------- eainfo ----------------
    ea = sub(sub(root, "eainfo"), "detailed")
    et = sub(ea, "enttyp")
    sub(et, "enttypl", CSV.name)
    sub(et, "enttypd", "Comma-separated table, one row per CRMS herbaceous marsh site and vegetation survey year.")
    sub(et, "enttypds", AUTHOR)
    for col, (defn, src, unit) in ATTR.items():
        a = sub(ea, "attr")
        sub(a, "attrlabl", col)
        sub(a, "attrdef", defn)
        sub(a, "attrdefs", src)
        dom = sub(a, "attrdomv")
        s = df[col]
        if unit == "enum":
            dom.getparent().remove(dom)
            for v in sorted(s.dropna().unique()):
                dom = sub(a, "attrdomv")
                ed = sub(dom, "edom")
                sub(ed, "edomv", v)
                sub(ed, "edomvd", f"CRMS community type: {v}")
                sub(ed, "edomvds", "CPRA CIMS")
        elif unit == "boolean":
            dom.getparent().remove(dom)
            for v in ("True", "False"):
                dom = sub(a, "attrdomv")
                ed = sub(dom, "edom")
                sub(ed, "edomv", v)
                sub(ed, "edomvd", "Condition met" if v == "True" else "Condition not met")
                sub(ed, "edomvds", AUTHOR)
        elif unit is not None and pd.api.types.is_numeric_dtype(s):
            r = sub(dom, "rdom")
            sub(r, "rdommin", f"{s.min():.4g}")
            sub(r, "rdommax", f"{s.max():.4g}")
            sub(r, "attrunit", unit)
        else:
            sub(dom, "udom", "Text identifier or date; see definition.")
    # ---------------- distinfo ----------------
    di = sub(root, "distinfo")
    contact(sub(di, "distrib"))
    sub(di, "resdesc", "coastal-marsh-vegetation-extreme-events/data/crms_marsh_site_year.csv")
    sub(di, "distliab", "The author provides this derived dataset as is, without warranty. Users should consult the original "
        "CPRA, NOAA and NASA sources before decisions.")
    so = sub(di, "stdorder")
    dfm = sub(so, "digform")
    ti = sub(dfm, "digtinfo")
    sub(ti, "formname", "CSV")
    sub(ti, "formcont", "UTF-8 comma-separated values with a header row")
    sub(ti, "transize", f"{CSV.stat().st_size / 1e6:.2f}")
    opt = sub(sub(dfm, "digtopt"), "onlinopt")
    ca = sub(opt, "computer")
    sub(sub(ca, "networka"), "networkr", "https://github.com/anamikadristi/coastal-marsh-vegetation-extreme-events")
    sub(so, "fees", "None")
    # ---------------- metainfo ----------------
    mi = sub(root, "metainfo")
    sub(mi, "metd", today.strftime("%Y%m%d"))
    contact(sub(mi, "metc"))
    sub(mi, "metstdn", "FGDC Content Standard for Digital Geospatial Metadata")
    sub(mi, "metstdv", "FGDC-STD-001-1998")
    return root


def main():
    df = pd.read_csv(CSV, dtype={"climdiv": "Int64"})
    today = dt.date.today()
    root = build(df, today)
    etree.ElementTree(root).write(str(XML), pretty_print=True, xml_declaration=True, encoding="UTF-8")
    lines = [f"FGDC CSDGM validation report, {today.isoformat()}",
             f"File: metadata/{XML.name}",
             "Schema: FGDC-STD-001-1998 XML schema (fgdc-std-001-1998.xsd and its ten section files)",
             "Schema source: https://www.fgdc.gov/schemas/metadata/ (downloaded keylessly by scripts/01_download_data.py)", ""]
    schema = etree.XMLSchema(etree.parse(str(XSD)))
    doc = etree.parse(str(XML))
    ok1 = schema.validate(doc)
    lines.append(f"Validator 1, lxml {etree.__version__} (libxml2): {'VALID' if ok1 else 'INVALID'}")
    lines += [f"  line {e.line}: {e.message}" for e in schema.error_log]
    xs = xmlschema.XMLSchema10(str(XSD))
    errs = list(xs.iter_errors(str(XML)))
    lines.append(f"Validator 2, xmlschema {xmlschema.__version__}: {'VALID' if not errs else 'INVALID'}")
    lines += [f"  {e.path}: {e.reason}" for e in errs[:50]]
    # negative controls: the schema must reject a record missing metainfo and a bad status value
    import copy
    bad1 = copy.deepcopy(doc)
    bad1.getroot().remove(bad1.getroot().find("metainfo"))
    bad2 = copy.deepcopy(doc)
    bad2.find(".//status/progress").text = "Done"
    lines.append(f"Negative control 1 (metainfo removed): {'rejected' if not schema.validate(bad1) else 'NOT rejected'}")
    lines.append(f"Negative control 2 (progress = 'Done', not an allowed value): "
                 f"{'rejected' if not schema.validate(bad2) else 'NOT rejected'}")
    n_attr = len(doc.findall(".//attr"))
    lines += ["", f"Sections present: " + ", ".join(t for t in ("idinfo", "dataqual", "spdoinfo", "spref", "eainfo",
                                                              "distinfo", "metainfo") if doc.find(t) is not None),
              f"Lineage sources: {len(doc.findall('.//srcinfo'))}; process steps: {len(doc.findall('.//procstep'))}; "
              f"attributes documented: {n_attr} of {df.shape[1]} columns",
              "Note: schema validity checks structure and element order only. It does not check that the content is",
              "complete or correct; a review against the CSDGM workbook or the USGS mp tool is still advisable."]
    REPORT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
