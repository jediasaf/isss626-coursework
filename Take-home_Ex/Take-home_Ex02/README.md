# Take-home Exercise 2: Geospatial Analytics for Peace

Local spatial association and emerging hot spots of armed conflict across
Myanmar's 330 townships, 1 January 2021 to 30 September 2025.

ISSS626 Geospatial Analytics and Applications, Take-home Exercise 2.

All commands below are run from **this directory**, `Take-home_Ex/Take-home_Ex02/`.

| Deliverable | Location |
|---|---|
| Technical report | `index.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex02/> |
| Executive summary | `slides.qmd` → <https://isss626-coursework.vercel.app/Take-home_Ex/Take-home_Ex02/slides.html> |
| Analysis pipeline | `R/00_config.R` through `R/04_diffusion.R` (plus `R/03b_ehsa_gi.R`) |

## Data

Neither input is committed. ACLED's terms do not permit redistribution of the
event file, and `data/` is git-ignored across the repository.

| File | Source | Place it at |
|---|---|---|
| `ACLED_Data_Myanmar_Jan2021-Sep2025.csv` | ACLED extract supplied on eLearn (87,109 rows, 31 columns) | `data/raw/aspatial/` |
| `mmr_polbnda_adm3_250k_mimu_1.*` | MIMU township boundaries, PCode 9.4, MIMU GeoNode WFS | `data/raw/geospatial/` (downloaded automatically by `R/01_prepare.R` if absent) |

## Reproducing it

```bash
Rscript R/01_prepare.R     # audit, clean, point-in-polygon join, quarterly and monthly cubes
Rscript R/02_lmsa.R        # weights, global Moran, Local Moran's I, Gi*, sensitivity
Rscript R/03_ehsa.R        # Gi* per quarter, Mann-Kendall, EHSA, location check, sensitivity (~15 min)
Rscript R/03b_ehsa_gi.R    # keeps EHSA's space-time Gi* series for the trajectory figure (~1 min)
Rscript R/04_diffusion.R   # value-added: spatial Markov diffusion test
quarto render index.qmd
quarto render slides.qmd
```

Each script caches its results to `data/derived/*.rds`; the report and slides
only read those files, so rendering takes seconds once the pipeline has run.

Tested with R 4.5.3, sf 1.1.1, sfdep 0.2.5, spdep 1.4-2, Kendall 2.2.2.

## A note on sfdep 0.2.5

`emerging_hotspot_analysis()` computes its results in alphabetical order of the
location column but labels them in the geometry's order. If the two differ, the
classifications land on the wrong townships with no error. `R/01_prepare.R` sorts
the townships by PCode so the orders coincide, and `R/03_ehsa.R` asserts that
they do and checks the labels against an independently computed Gi* series.
