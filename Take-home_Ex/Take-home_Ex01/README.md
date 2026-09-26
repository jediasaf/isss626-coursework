# Take-home Exercise 1: Geospatial Analytics for Public Safety

First- and second-order point-pattern analysis of VIIRS active-fire detections in
**Pulang Pisau Regency, Central Kalimantan**, 1 July to 25 September 2026.

ISSS626 Geospatial Analytics and Applications, Take-home Exercise 1.

This directory sits inside the ISSS626 coursework site. All commands below are run
from **this directory**, `Take-home_Ex/Take-home_Ex01/`.

| Deliverable | Location |
|---|---|
| Technical report | `index.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex01/> |
| Executive summary, 10 slides | `slides.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex01/slides.html> |
| Analysis pipeline | `R/00_config.R` through `R/08_extended_range.R` |

## The argument in one paragraph

A VIIRS hotspot is a 375 m pixel observed during one satellite overpass, not a
fire. One large peat fire fills many pixels and is re-detected roughly four times
a day, so the raw detection set overstates the number of distinct burning
locations by a factor of about 3.4 and any test against complete spatial
randomness rejects before the analysis starts. This project changes the unit of
analysis to the 375 m cell and the null model to an inhomogeneous Poisson process
whose intensity is re-estimated for every simulated pattern. Under that harder
null the clustering survives, out to roughly 2 km; a fitted Thomas process puts
the operational cluster radius near 1 km by an independent route; and a Knox
permutation test under random labelling finds space-time interaction concentrated
within a few hundred metres and a few days.

## Reproducing it

```bash
# 1. FIRMS needs a free key: https://firms.modaps.eosdis.nasa.gov/api/map_key/
cp .Renviron.example .Renviron && $EDITOR .Renviron

# 2. Covariates come from a Geofabrik extract (~390 MB, downloaded once)
curl -L -o data/raw/kalimantan-latest-free.gpkg.zip \
  https://download.geofabrik.de/asia/indonesia/kalimantan-latest-free.gpkg.zip

# 3. Run the pipeline in order
Rscript R/01_acquire.R          # download, retention probe, regency ranking
Rscript R/02_prepare.R          # screening, audit trail, the three patterns
Rscript R/03_firstorder.R       # quadrats, bandwidths, KDE, FRP, temporal
Rscript R/04_secondorder.R      # G, F, L, pcf, inhomogeneous versions, global tests
Rscript R/05_spacetime.R        # Knox test, relative risk, trajectory, ST kernel
Rscript R/06_covariates_models.R  # OSM covariates, ppm, AIC, residual K, scan test
Rscript R/07_cluster_scale.R      # Thomas process by explicit minimum contrast
Rscript R/08_extended_range.R     # second-order statistics to a 20 km lag

# 4. Render, from the repository root
quarto render Take-home_Ex/Take-home_Ex01/index.qmd
quarto render Take-home_Ex/Take-home_Ex01/slides.qmd
```

`data/` is git-ignored, matching the convention used by the hands-on exercises in
this repository. Nothing in it is required to read the published pages; it is
required only to re-run the analysis.

Each script caches to `data/derived/*.rds`; the two documents read the caches and
compute no analysis inline, so rendering cannot silently disagree with what was
analysed. The seed is fixed at 626 in `R/00_config.R`, so all Monte Carlo results
reproduce exactly.

## Requirements

R 4.5 with `sf`, `spatstat` (`.geom`, `.explore`, `.model`, `.random`),
`tidyverse`, `ggplot2`, `patchwork`, `scales`, `knitr`, and Quarto 1.9.

`sparr` and `stpp` are deliberately **not** used. Both load `tcltk` at startup and
fail without XQuartz on macOS. The space-time kernel estimate and the space-time
interaction test are implemented directly on `spatstat` primitives instead, which
also makes the separability assumption and the permutation null explicit rather
than hidden inside a package default.

## Data sources

| Source | Detail | Licence |
|---|---|---|
| NASA FIRMS | VIIRS 375 m C2 NRT, S-NPP and NOAA-20, `area` API | Open, attribution requested |
| geoBoundaries | gbOpen Indonesia ADM1/ADM2, commit `9469f09` | CC BY 3.0 IGO |
| OpenStreetMap | Geofabrik Kalimantan extract, roads and waterways | ODbL |

`data/raw/` is git-ignored: it is large and fully re-creatable by `R/01_acquire.R`
plus the Geofabrik download above.

## Two things this analysis cannot do

The FIRMS `area` REST endpoint retains near-real-time granules for roughly 90 days
and the standard-processing product is not built for 2026 yet, so the scripted
record begins on 1 July 2026. There is no wet-season baseline and no inter-annual
comparison, and nothing here speaks to whether 2026 was unusual.

This is a limit of the access route, not of the data. FIRMS publishes the full
archive (VIIRS S-NPP from January 2012) through its download service at
<https://firms.modaps.eosdis.nasa.gov/download/>, which authenticates by Earthdata
login or emailed code and delivers by email. It cannot be driven from an API key,
which is why this pipeline does not use it, but it would lift the limitation
entirely.

OpenStreetMap does not map the minor canal network of the former Mega Rice
Project. The drainage covariate therefore measures mapping coverage, and its flat
relationship with intensity must not be read as evidence that drainage does not
matter.
