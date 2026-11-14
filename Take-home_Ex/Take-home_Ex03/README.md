# Take-home Exercise 3: Geospatial Analytics for Urban Sustainability

A Huff model of weekday morning-peak boardings at Singapore's MRT and LRT stations,
calibrated on LTA passenger volumes, used to predict boardings at the Jurong Region Line
and Cross Island Line (phase 1) stations and their effect on the stations already open.

ISSS626 Geospatial Analytics and Applications, Take-home Exercise 3.

All commands below are run from **this directory**, `Take-home_Ex/Take-home_Ex03/`.

| Deliverable | Location |
|---|---|
| Technical report | `index.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex03/> |
| Executive summary | `slides.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex03/slides.html> |
| Analysis pipeline | `R/00_config.R` and `R/00_download.R` through `R/08_sensitivity.R` |
| Improvement audit | `IMPROVEMENT_AUDIT.md` (what was checked, corrected and verified on 14 Nov 2026) |

## Data

No data file is committed; `data/` is git-ignored across the repository.
`R/00_download.R` fetches everything it can by script. Three LTA files refuse scripted
requests and have to be saved by hand from the DataMall static data page, and the
passenger volumes need a free DataMall API key.

| File | Source | How it arrives |
|---|---|---|
| `transport_node_train_YYYYMM.zip` | LTA DataMall, Passenger Volume by Train Stations (API, latest three months) | script, needs `DATAMALL_KEY` in `~/.Renviron` |
| `TrainStation_Sep2026.zip` | LTA DataMall static data, Train Station | by hand |
| `BusStopLocation_Sep2026.zip` | LTA DataMall static data, Bus Stop Location | by hand (OpenStreetMap bus stops are used if absent) |
| `HDBExistingBuilding.geojson` | data.gov.sg, HDB Existing Building | script |
| `HDBPropertyInformation.csv` | data.gov.sg, HDB Property Information (dwelling units per block) | script |
| `MasterPlan2019SubzoneBoundaryNoSea.geojson` | data.gov.sg, URA Master Plan 2019 subzones | script |
| `operational_stations.csv`, `future_stations.csv` | github.com/elliotwutingfeng/singapore_train_station_coordinates, pinned commit `9f5972c` | script |
| `malaysia-singapore-brunei-latest.osm.pbf` | Geofabrik OpenStreetMap extract | script |
| `future_station_schedule.csv` (committed, not in `data/`) | Planned opening of every future station, compiled from LTA publications with a source on every row | in the repository |

## Reproducing it

```bash
Rscript R/00_download.R      # every input that can be fetched by script
Rscript R/01_street_codes.R  # link HDB dwelling units to building footprints (OneMap, about 15 min)
Rscript R/02_osm_poi.R       # points of interest and fallback bus stops from OpenStreetMap
Rscript R/03_prepare.R       # origins, stations, attraction variables, ridership, interaction frames
                             #   (station rules live in R/prep_functions.R, shared with 05b)
Rscript R/04_calibrate.R     # Huff calibration, cross-validation, bootstrap
Rscript R/05_scenario.R      # today + JRL + CRL phase 1: predicted boardings and market areas
Rscript R/05b_phased_scenarios.R  # the same model on the end-2030 and full planned networks
Rscript R/04b_diagnostics.R  # data and model checks (needs 04 and 05)
Rscript R/06_value_added.R   # regional distance decay and exponential decay
Rscript R/07_scenario_r5.R   # the expanded network re-run with regional decay
Rscript R/08_sensitivity.R   # all-metric model comparison; exclusion-rule and radius sensitivity
quarto render index.qmd && quarto render slides.qmd
```

Rendering `index.qmd` on a fresh checkout runs the whole chain itself if its outputs are
missing (the setup chunk), so the report builds end to end from the raw data.

Run times on a 2020s laptop: `01` about 15 minutes (OneMap requests), `04` about 8 minutes
(200 bootstrap refits), `06` about 25 minutes (100 regional refits); everything else in
seconds to two minutes.

R 4.x with sf, tidyverse, MCI, spdep, httr2, patchwork, ggrepel and ragg (exact versions are printed at the end of the report). The report
and slides read the `.rds` files the scripts write to `data/derived/`, so they render
in seconds once the pipeline has run.
