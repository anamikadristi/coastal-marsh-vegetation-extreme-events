# FGDC metadata notes

## What was produced

- `metadata/crms_marsh_site_year_fgdc.xml` describes `data/crms_marsh_site_year.csv`, the main derived dataset (one row per CRMS herbaceous marsh site and survey year).
- The record follows the FGDC Content Standard for Digital Geospatial Metadata, FGDC-STD-001-1998 (CSDGM).
- `scripts/11_fgdc_metadata.py` builds the record. Row counts, the date range, the bounding box, attribute value ranges and file size are computed from the CSV each time it runs, so the metadata stays in step with the data.

## Sections included

| CSDGM section | Content |
|---|---|
| `idinfo` | Citation, abstract, purpose, supplemental information, time period of content (survey dates), status, bounding coordinates, theme keywords (ISO 19115 topic categories and free terms), place and temporal keywords, access and use constraints, point of contact, native environment |
| `dataqual` | Attribute accuracy, logical consistency, completeness (including the reason FQI is absent), horizontal positional accuracy statement, **lineage**: 6 source citations (CRMS vegetation, CRMS hydrology, HURDAT2, GHCN-Daily, nClimDiv, MOD13Q1), each with time period and contribution, plus 7 process steps that map one-to-one to scripts 01-07 |
| `spdoinfo` | Direct spatial reference method `Point`, SDTS point count |
| `spref` | Geographic coordinates in decimal degrees, resolution 0.0001 |
| `eainfo` | Detailed entity and attribute information: all 34 columns, each with a definition and definition source. Numeric columns have range domains with units, marsh types and true/false fields have enumerated domains, and text identifiers and dates have unrepresentable domains |
| `distinfo` | Distributor, resource description, liability statement, standard order process (CSV, online at the repository URL, no fees) |
| `metainfo` | Metadata date, metadata contact, standard name and version |

## Validation

- The schema is the official CSDGM XML schema from the FGDC: `fgdc-std-001-1998.xsd` plus its ten section files (`-sect01` to `-sect10`). Step 1 downloads it without a key from <https://www.fgdc.gov/schemas/metadata/>.
- Two independent validators were run:
  - lxml 6.1.3 (libxml2): **VALID**;
  - xmlschema 4.3.2 (pure Python): **VALID**.
- As negative controls, the same schema rejected a copy with `metainfo` removed and a copy with an invalid `progress` value. This shows the validation actually checks the record.
- The full output is in `metadata/fgdc_validation_report.txt`.

Schema validity covers structure, required elements, element order and controlled values only. It does not show that the content is complete or accurate. Before formal submission, a review with the USGS `mp` metadata parser or the USGS Online Metadata Editor is advisable. Neither was run here.

## Choices and gaps to know about

- **Contact phone.** CSDGM requires `cntvoice`. To keep personal phone numbers out of a public repository, it is filled with "Not provided; contact by email". The mailing city, state and postal code are also "Not provided". An agency release would use an office address and phone.
- **Horizontal datum.** The CIMS vegetation export and its bundled metadata do not state the datum of the plot coordinates. The optional `geodetic` block (datum and ellipsoid) is therefore left out instead of guessed. If CPRA confirms the datum (CRMS work commonly uses NAD83), add it under `spref/horizsys/geodetic`.
- **Source publication dates.** The lineage citations use `pubdate` = "Unknown" because the online sources are continuously updated. The access date (24 October 2026) and the file names with version stamps are recorded in `raw_data/download_log.json` and the README.
- **FQI.** The Floristic Quality Index was not computed, because the coefficient of conservatism list could not be downloaded from a keyless source (lacoast.gov returned HTTP 502). The `complete` element says so.
- **Halophyte list** used for `halo_share`: Spartina (Sporobolus) alterniflora, Distichlis spicata, Juncus roemerianus, Avicennia germinans, Batis maritima, Salicornia, Sarcocornia, Lycium carolinianum, Borrichia frutescens, Sesuvium, Limonium, Suaeda (prefix match on "Scientific Name As Currently Recognized").
- **Source terms.** CPRA asks that CIMS data be obtained from CPRA directly. The metadata's use constraints point users to the original sources, and NOAA and NASA data are public domain.

## Reusing the template

`build()` in `scripts/11_fgdc_metadata.py` is a plain lxml tree builder. To document another table:
1. edit the `ATTR` dictionary (column, definition, source, unit or `enum`/`boolean`);
2. edit `SOURCES` and `STEPS`;
3. rerun the script. It re-validates against the same XSD and rewrites the report.
