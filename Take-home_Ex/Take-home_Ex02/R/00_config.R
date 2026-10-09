# ---------------------------------------------------------------------------
# ISSS626 Take-home Exercise 2 - shared configuration
# Every parameter that governs the study design lives here, so the report, the
# slides and the scripts cannot drift out of agreement.
# ---------------------------------------------------------------------------

## --- Study window -----------------------------------------------------------
# Fixed by the brief. The window ends on the last day of a calendar quarter, so
# a quarterly cube has 19 complete quarters and no partial final period.
WINDOW_START <- as.Date("2021-01-01")
WINDOW_END   <- as.Date("2025-09-30")

# The coup. Used only for annotation, never to cut the data.
COUP_DATE    <- as.Date("2021-02-01")
# Operation 1027, the Three Brotherhood Alliance offensive in northern Shan.
# Used to annotate the series and to test whether the diffusion model is stable.
OP1027_DATE  <- as.Date("2023-10-27")

## --- Inputs -----------------------------------------------------------------
# ACLED Myanmar extract supplied on eLearn. ACLED's terms do not allow the raw
# file to be redistributed, so data/ is git-ignored; README.md says how to
# obtain it.
ACLED_FILE <- "data/raw/aspatial/ACLED_Data_Myanmar_Jan2021-Sep2025.csv"
# MIMU township boundaries, PCode release 9.4, fetched from the MIMU GeoNode
# WFS endpoint recorded in R/01_prepare.R.
MIMU_FILE  <- "data/raw/geospatial/mmr_polbnda_adm3_250k_mimu_1.shp"
MIMU_WFS   <- paste0(
  "http://geonode.themimu.info/geoserver/geonode/ows?format_options=charset%3AUTF-8",
  "&outputFormat=SHAPE-ZIP&version=1.0.0&service=WFS&request=GetFeature",
  "&typename=geonode:mmr_polbnda_adm3_250k_mimu_1")

## --- What counts as armed conflict -----------------------------------------
# ACLED's three political-violence event types. Protests and riots are excluded
# because they are not armed conflict, and strategic developments because they
# are not violent events (arrests, base changes, territory transfers).
PV_TYPES <- c("Battles", "Explosions/Remote violence", "Violence against civilians")
VAR_LABELS <- c(
  pv         = "All political violence",
  battles    = "Battles",
  erv        = "Explosions and remote violence",
  vac        = "Violence against civilians",
  fatalities = "Reported fatalities")

## --- Inclusion rules for an event ------------------------------------------
# geo_precision 3 means ACLED could only place the event in a larger region and
# coded it to that region's centroid. Assigning such an event to a township
# would be arbitrary, so it is dropped from township counts.
GEO_PRECISION_KEEP <- c(1L, 2L)
# Events whose coordinates fall just outside every polygon (coastline
# generalisation at 1:250k) are snapped to the nearest township if within this
# distance. Beyond it they are treated as offshore and dropped.
SNAP_TOLERANCE_M <- 5000

## --- Projection -------------------------------------------------------------
# WGS 84 / UTM zone 47N. Myanmar's centroid (~96.5 E) lies in zone 47
# (96-102 E). A metric CRS is needed for the snapping distance, the k-nearest
# neighbour fallback and areas.
STUDY_CRS <- 32647L

## --- Spatial weights --------------------------------------------------------
# Queen contiguity, row-standardised. Townships with no contiguous neighbour
# (islands) are linked to their K_ISLAND nearest neighbours by centroid
# distance, so no township is dropped and no zero-weight row enters the tests.
K_ISLAND <- 1L
# Sensitivity check: adaptive k-nearest neighbours with this k, about the
# median number of queen neighbours.
K_SENS   <- 6L

## --- Inference --------------------------------------------------------------
ALPHA      <- 0.05      # fixed by the brief
NSIM_LMSA  <- 999L      # permutations for local Moran and Gi*
NSIM_EHSA  <- 199L      # permutations per space-time Gi*; floor 1/200 = 0.005
SEED       <- 626L

## --- Diffusion analysis (value-added) ---------------------------------------
# A township-quarter is "active" if it records at least ACTIVE_MIN political
# violence events. Chosen as one event a month: it separates sustained activity
# from a single isolated incident. Sensitivity to 1 and 6 is reported.
ACTIVE_MIN <- 3L
ACTIVE_SENS <- c(1L, 6L)

## --- Paths ------------------------------------------------------------------
DIR_DERIVED <- "data/derived"
