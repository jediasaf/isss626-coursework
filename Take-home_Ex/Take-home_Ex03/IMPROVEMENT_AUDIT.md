# Improvement audit: Take-home Exercise 3

Audit of the version published on 14 November 2026 (commit `b220ad7`), carried out the
same evening against the official brief (course site, last modified 30 October 2026).
Every number below was recomputed from the pipeline; none was copied from the old slides.

## 1. Does the work do what the brief asks?

| Task in the brief | Before | After |
|---|---|---|
| Spatial interaction data frame, existing network | `si_existing` (127,174 block-station pairs) | unchanged, now shown in the report (section 4.5) |
| Spatial interaction data frame, with proposed stations | `si_2030`, but it included JS2A, which is not one of the brief's JRL stations | JS2A removed; phased variants added (`R/05b_phased_scenarios.R`) |
| Origin table with local market potential | `origins`, 10,817 HDB blocks, 1,179,156 dwelling units | unchanged, checks added |
| Destination table with observed ridership and attractions | `stations` | unchanged for existing stations (verified identical) |
| Calibrate a probabilistic spatial interaction model | Huff model with MCI-style attractiveness, four nested models | unchanged; error measures widened |
| Estimate and map ridership | yes | yes, plus an observed-against-predicted chart with factor-of-two lines |
| Predict existing and proposed stations | yes | yes, on the corrected network |
| Show the impact on existing stations' market areas | yes | yes, plus a ranked change chart and a clearer catchment map |

## 2. The research design, stated plainly

- **Unit:** a station complex (an interchange counts once); 184 operating today.
- **Time:** weekday tap-ins 06:00 to 09:00, mean of August, September and October 2026, per month.
- **Who is counted:** residents of HDB flats, through the dwelling units of each block. Private housing is not.
- **Trip assumption:** morning boardings start at home.
- **Observed ridership:** tap-ins from LTA's Passenger Volume by Train Stations.
- **Attractiveness:** lines served, bus stops and points of interest within 400 m, LRT-only flag.
- **Distance:** straight-line, metres, SVY21; power decay; floor at 100 m.
- **Calibration:** least squares on log boardings; trip rate *k* solved exactly.
- **Network scenario:** today's network against today's plus the JRL (24 stations) and CRL phase 1 (12).
- **Validation:** 10-fold cross-validation over stations, all models on the same folds.

## 3. Issues found, and what was done

### 3.1 The "2030 network" included a station that opens in the mid-2030s

- **Before:** the expanded network took every station coded JE, JS or JW, which includes JS2A.
- **Problem:** JS2A is an infill station LTA announced on 4 March 2026 and targets for the mid-2030s, "in tandem with" surrounding development. It is not one of the JRL's 24 stations that the brief refers to, and LTA has not named it ("Forest Hill" comes from the station list, not LTA).
- **Evidence:** LTA factsheet, 4 March 2026; LTA JRL page ("JS1 to JS8 (excluding JS2A)").
- **Correction:** `JRL_LATER <- "JS2A"` in `R/00_config.R`; `R/03_prepare.R` assigns it to the later network. New stations in the brief's network: 30 to 29.
- **Files:** `R/00_config.R`, `R/03_prepare.R`; reruns of `05`, `05b`, `07`, `04b`, `08`.
- **Verified:** yes. Today's stations are identical before and after (`all.equal` on every existing row), so the calibration and the regional model did not need refitting.

### 3.2 Opening schedules were asserted, not checked

- **Before:** the deck said the new stations would exist "by 2030" without a source.
- **Problem:** the brief's "completed by 2030" needed checking against LTA's current statements.
- **Evidence:** every planned station's official schedule, with sources, now sits in `future_station_schedule.csv`. Key facts were spot-checked against the saved LTA pages: JRL stage 1 mid-2028 (moved from end-2027); stages 2 and 3 in 2028 and 2029, which LTA says it will update; CRL phase 1 by 2030; CRL phase 2 and Punggol extension 2032; TEL stage 5 and DTL3 extension on 26 December 2026; Circle Line stage 6 opened 12 July 2026.
- **Correction:** three networks are built from the schedule table: the brief's; end-2030 (adds Bedok South, Sungei Bedok, Xilin and Founders' Memorial); and the full planned network, labelled as not a forecast of any year. The brief's comparison stays the main one.
- **Result:** adding the four end-2030 stations leaves the JRL and CRL predictions unchanged (rank correlation 1.00; total change under 0.1%). Their main effect is on Bayshore (-33%), far from both lines.
- **Files:** `future_station_schedule.csv` (new), `R/05b_phased_scenarios.R` (new), `R/00_config.R`, `R/prep_functions.R` (new, shared station rules), report section 7.1.
- **Verified:** yes, with these caveats. The Founders' Memorial year is inferred from the memorial's own end-2028 date. No official date was found for Bukit Brown (CC18) or Tagore (TE4A), so they appear only in the full network.

### 3.3 Accuracy was reported with log-scale measures only

- **Before:** R² on logs, RMSLE, MdAPE and Spearman.
- **Problem:** these say nothing about errors in boardings, or bias in the total.
- **Correction:** `R/08_sensitivity.R` adds MAE, RMSE, WAPE and bias for every model (M0 to M3, R5, EXP), all on the same folds.
- **Result:** the full model's WAPE is 33% and it over-predicts the fitted stations' total by 9%. The regional model brings WAPE to 30% and bias to 2%. The textbook model under-predicts the total by 24%.
- **Verified:** yes. The old headline numbers reproduce exactly: R² 0.41 (M0), 0.75 (M3) and 0.77 (R5); MdAPE 34%; Spearman 0.90; exponent -2.48.

### 3.4 The exclusion of 36 stations was asserted, not tested

- **Before:** stations with fewer than 1,000 HDB dwelling units within 1 km were left out of the fit, with no sensitivity check.
- **Problem:** this could be selection bias, and an arbitrary threshold.
- **Evidence and correction:** refits with thresholds of 0, 500 and 2,000 (`R/08_sensitivity.R`).
  - 500 or 2,000 give almost the same model (exponent -2.2 or -2.7; R² on the common 148 stations 0.75 or 0.76; new-station totals within 5%; rank correlation 0.99 or better).
  - No threshold at all means fitting 180 stations. The exponent flattens to -1.17, R² on the residential stations falls from 0.74 to 0.64, and new-station totals rise 16%.
  - Under the main model, the excluded stations get a median of about 6% of their observed boardings, and the four Tuas stations get none, since no HDB block can reach them. 27 of the 36 receive more people than they send out in the morning.
- **Conclusion recorded in the report:** the rule keeps the model to stations its data can explain. The price is that it says nothing reliable about employment-led stations.
- **Verified:** yes.

### 3.5 The choice-set radius was never varied

- **Correction:** refits at 2 km and 5 km.
- **Result:** 5 km fits better (R² 0.79 against 0.74), with almost the same new-station predictions (totals -2%, rank 0.99). 3 km is kept because it was chosen before seeing results, and the report says so. Switching on test scores would overstate accuracy.
- **Verified:** yes.

### 3.6 Validation design: leakage, and what it can and cannot test

- **Checked:** within each fold, the exponents and *k* are fitted on training stations only. The exclusion rule and the attraction variables do not use ridership. No leakage was found.
- **Problem:** the old text implied the cross-validation tests new stations.
- **Correction:** the report and the slides now say it hides mature stations, not new ones. A true before-and-after test is impossible because DataMall serves only three months. The Circle Line stage 6 stations, still ramping up, show why: Prince Edward Road's morning boardings rose 45% from August to October.
- **Verified:** yes (`R/04b_diagnostics.R`).

### 3.7 Attractiveness effects were worded as causes

- **Before:** "people walk further to an interchange because it saves a transfer later".
- **Correction:** these are now reported as associations, holding the other variables fixed, with saving a transfer given as a plausible but unidentified reason.
- **New finding:** the attraction variables are not collinear (largest correlation 0.34, all VIFs below 1.4). But in the bootstrap the distance exponent and the LRT effect move together (correlation 0.86), so only their combined effect is firm. The report and the slide say this.
- **Multipliers reproduced:** a second line ×2.6; ten times the points of interest ×1.8; LRT-only 6%; twice the distance ÷5.6.
- **Verified:** yes.

### 3.8 Residual analysis: the test and the explanations

- **Before:** Moran's I 0.16 from an analytical test with 6 neighbours; explanations stated as fact.
- **Correction:** a 999-run permutation test (p = 0.001), significant for 4, 6, 8 and 10 neighbours (I from 0.13 to 0.22). The explanations are now "consistent with where the errors fall, not proven causes". The ramp-up evidence for Keppel is added.
- **Verified:** yes.

### 3.9 Uncertainty intervals were open to misreading

- **Before:** "95% bootstrap intervals" on predicted boardings, without saying what they cover.
- **Correction:** labelled everywhere as parameter uncertainty only. The report lists what they do not cover: model form, station scatter, spatial error, future land use, and which stations open and when.
- **Verified:** yes.

### 3.10 Redistribution was described as "relief"

- **Before:** "for the operator this is relief"; feeder buses "will need to be redrawn".
- **Problem:** the model holds homes and trips fixed, so it can only redistribute boardings. Entries are not train loads.
- **Evidence:** total predicted boardings are identical in both networks (difference about 1e-14%).
- **Correction:** changes are now "reallocation of modelled HDB-origin morning boardings". The report separates station access from crowding. The feeder-bus recommendation is replaced by "needs a bus-network study".
- **Verified:** yes.

### 3.11 Low-demand new stations were called the "ghost-station risk"

- **Before:** "In a home-based model they are the ghost-station risk"; "their predictions are a floor".
- **Problem:** few HDB homes nearby is not the same as low use. Employment, campus and airport demand is outside the model, and today's surroundings could understate or overstate 2030's.
- **Correction:** the text now says they have few HDB residents nearby, not that they will be little used, and today's figure is "neither a floor nor a ceiling". The opening recommendation is replaced by the evidence still needed.
- **Verified:** yes.

### 3.12 A factual error about Tengah

- **Before:** "Tengah's dwelling units currently have to reach stations two or more kilometres away".
- **Evidence:** the straight-line distance to the nearest station is 1.3 to 2.2 km, median 1.7 km.
- **Correction:** the report now computes and quotes this.
- **Verified:** yes.

### 3.13 Reproducibility

- **Before:** the report read saved results and showed little of the code that made them.
- **Correction:**
  - every script now appears in the report as a folded, full code listing beside the section it supports;
  - the report runs the whole pipeline itself if the saved results are missing;
  - software versions are printed;
  - the station rules live in one shared file, `R/prep_functions.R`, and the refactor was verified to reproduce all outputs exactly;
  - the README lists the new scripts, run times and the schedule file.
- **Verified:** yes for the refactor and the rendering. A from-scratch run was not repeated, because it needs about 50 minutes and a DataMall key.

### 3.14 Slides

- **Before:** section-heading titles; a closing slide with unsupported recommendations.
- **After:**
  - ten substantive slides with one claim per title, one main visual, a short plain reading and a marked qualification;
  - a worked example of the model on one Tengah block;
  - a ranked change chart for existing stations;
  - a closing slide that separates what is supported, what is suggested and what is still needed, across four themes.
- **Fixed during inspection:**
  - the dashed-line explanation on slide 2 was reversed;
  - the closing-slide grid broke, and its headings were being turned into extra slides.

## 4. Updated results, before and after

| Quantity | Before | After |
|---|---|---|
| New stations in the brief's network | 30 | 29 (JS2A removed) |
| Cross-validated R² (log): M0 / M3 / R5 | 0.41 / 0.75 / 0.77 | same (reproduced) |
| Rank correlation (M3) | 0.90 | 0.90 |
| Typical miss, one station (MdAPE, M3) | 34% | 34% |
| WAPE / bias in total (M3) | not reported | 33% / +9% |
| Distance exponent (M3) | -2.48 | -2.48 (interval -3.24 to -1.82) |
| Moran's I of errors | 0.16 (analytical) | 0.16 (permutation p = 0.001) |
| Tengah Plantation, modelled boardings per month | 178,800 | 179,900 |
| Choa Chu Kang West | 111,900 | 113,600 |
| Largest losses (Chinese Garden / Lakeside / Pioneer) | -65% / -61% / -56% | -65% / -61% / -56% |
| Blocks changing main station / to a new one | 1,342 / 1,042 | 1,342 / 1,042 |

## 5. Verification summary

- **Rendering:** `index.qmd` and `slides.qmd` render without errors.
- **Slide count:** 12 in total, the cover, contents and 10 substantive slides.
- **Visual inspection:** every slide was viewed in a browser and fixed where needed, including text overflow, a broken grid and a reversed sentence. The report's figures, tables and code listings were checked in the browser: 11 figures, none broken, and no missing values or leftover placeholders.
- **Maths:** the Huff equations display in the report and on slide 3.
- **Consistency:** every number on the slides is computed from the same saved results as the report, with no hand-typed results. The abstract was corrected to name the same top stations as the slides.
- **Not done:**
  - a from-scratch rerun of the full pipeline;
  - a check on a second browser or device;
  - a check of the published site, since nothing has been pushed.

## 6. Unresolved limitations

- The origin-destination passenger data could not be downloaded (DataMall quota). Station-level tap-ins are used instead, which is all this model needs, but it rules out a destination-side model.
- LTA's bus stop file could not be downloaded by script. OpenStreetMap bus stops are used.
- Private housing is not represented. This is the most likely cause of the clustered errors.
- Straight-line distances stand in for walking-network distances.
- No before-and-after test at a real opening is possible with three months of data.
- JRL stages 2 and 3 dates are provisional, the JS2A name is unofficial, the Founders' Memorial year is inferred, and CC18 and TE4A have no official date.

## 7. Remaining manual actions before submission

1. Review the report and slides locally (the browser preview on port 8788 serves `_site`).
2. Approve commit, push and deploy. Nothing has been pushed or deployed in this audit.
3. Submit the three links on eLearn by 15 November 2026, 11:59pm.
