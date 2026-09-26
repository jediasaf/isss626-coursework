# ---------------------------------------------------------------------------
# ISSS626 Take-home Exercise 1 - shared configuration
# Every parameter that governs the study design lives here, so that the
# report, the slides and the scripts cannot drift out of agreement.
# ---------------------------------------------------------------------------

## --- Temporal study window -------------------------------------------------
# The brief requires 2026 observations. Two constraints set the window.
#
# Lower bound, imposed by the ACCESS ROUTE rather than by the data, and verified
# empirically. The FIRMS `area` REST API serves near-real-time granules on a
# rolling ~90-day retention, and the standard-processing (SP) product has not
# been built for 2026 because NRT is only replaced by standard science quality
# after a ~5 month lag.
#
# Note carefully: FIRMS DOES publish the full archive (VIIRS S-NPP from
# 2012-01-20, NOAA-20 from 2018-04-01) through its archive download service at
# https://firms.modaps.eosdis.nasa.gov/download/ . That service authenticates by
# Earthdata login or emailed code and delivers by email, so it cannot be driven
# from an API key inside this pipeline. The 87-day window is therefore a
# consequence of choosing a scriptable route, not a hard limit on the record.
# Queries
# for a Kalimantan-wide box return a header with zero data rows for every
# sampled date from January to 30 June 2026 and non-zero counts from 1 July
# onward. Six consecutive months of exactly zero detections across a
# 10.6 deg x 9.7 deg box is not a physically credible result, so the zeros are
# read as absence of coverage rather than absence of fire. 1 July 2026 is the
# earliest retrievable day. See R/01_acquire.R, which re-runs the probe.
#
# Upper bound, a design choice: the last complete UTC day before access, so no
# partially-observed day enters the daily count series.
#
# The resulting 87-day window is not merely what was available. It spans the
# onset, escalation and peak of the 2026 Kalimantan dry-season fire episode,
# which is the period over which a fire-management agency would act.
WINDOW_START <- as.Date("2026-07-01")
WINDOW_END   <- as.Date("2026-09-25")
ACCESS_DATE  <- as.Date("2026-09-26")
# Earliest date the API returned data for, recorded by the retention probe.
ARCHIVE_FLOOR <- as.Date("2026-07-01")

## --- Sensor / product choice ----------------------------------------------
# VIIRS 375 m (I-band) Collection 2 from both operational platforms. MODIS is
# excluded: its 1 km footprint would confound the second-order analysis, whose
# distances of interest start below 1 km.
# The standard-processing (SP) archive is empty for 2026, verified by probe in
# R/01_acquire.R, so the NRT products are the only option.
FIRMS_SENSORS_NRT <- c(SNPP = "VIIRS_SNPP_NRT", NOAA20 = "VIIRS_NOAA20_NRT")
# The area API rejects day_range outside [1..5]; documented limit, tested.
FIRMS_MAX_DAY_RANGE <- 5L

# Kalimantan bounding box (W,S,E,N) used for the FIRMS area query.
KALIMANTAN_BBOX <- c(west = 108.5, south = -5.2, east = 119.1, north = 4.5)

## --- Inclusion rules for a point event ------------------------------------
# See report section "Point event definition" for the justification of each.
# Confidence. VIIRS confidence is categorical, not a percentage. The area API
# returns single letters (l/n/h); the keyless regional files spell the words.
# Both encodings are normalised in R/02_prepare.R. Low-confidence detections are
# dropped: they are dominated by small or cool anomalies and by edge-of-scan
# artefacts, and retaining them would inflate apparent intensity in a way that
# varies with scan geometry rather than with fire activity.
CONF_KEEP      <- c("nominal", "high")

# Source type. The `type` field (0 vegetation fire, 1 volcano, 2 other static
# land source, 3 offshore) exists ONLY in the standard-processing product. The
# NRT granules used here do not carry it, verified against the API response
# header. Three substitutes are used instead, documented in R/02_prepare.R:
#   type 1  no Holocene volcano exists on Borneo, so the risk is nil
#   type 3  removed by the clip to the terrestrial regency polygon
#   type 2  addressed by the persistence screen below
TYPE_KEEP      <- 0L          # applied only if the column is present
TYPE_AVAILABLE_IN_NRT <- FALSE

# Persistence screen, standing in for the missing `type` field. A 375 m cell
# detected on more than this fraction of the days in the window behaves like a
# fixed industrial heat source (flare stack, mill, coal face), not a vegetation
# fire, which burns out or moves. Kalimantan has both palm-oil mills and coal
# operations, so this is a live contamination risk rather than a hypothetical.
# Sensitivity to the threshold is reported.
PERSIST_DAY_FRAC <- 0.60
DEDUP_GRID_M   <- 375                   # VIIRS I-band nominal pixel size
# Study regency, chosen by a stated rule applied to the ranking computed in
# R/01_acquire.R over the whole 87-day window (not a convenience sample).
#
# Rule: the highest detection density among regencies large enough for
# second-order estimation to a 5 km lag, operationalised as area > 5,000 km2.
# Pulang Pisau wins this outright at 3,656 detections per 1,000 km2.
#
# Two alternatives were rejected, and the reasons are recorded because the
# choice is only reviewable alongside them:
#   Ketapang leads on ABSOLUTE count (65,901 against 35,679) but is 29,722 km2,
#     three times the area, and ranks second on density. A sparser pattern in a
#     much larger window is a weaker basis for second-order work.
#   Kota Palangka Raya has the highest density of any unit (8,222) but covers
#     only 2,606 km2, so a 5 km lag approaches the window's own dimension and
#     edge effects would dominate every estimate.
# Pulang Pisau is also rank 2 of 56 on raw count, so it is not a fringe pick.
#
# It is the degraded peatland of the former Mega Rice Project, so the
# public-safety framing (peat smoke, respiratory burden, suppression access) is
# substantive rather than incidental.
STUDY_REGENCY  <- "Pulang Pisau"
STUDY_PROVINCE <- "Central Kalimantan"
# WGS 84 / UTM zone 49S. The regency centroid is 113.96 E, inside zone 49
# (108-114 E). A projected, equal-unit CRS is required because every
# second-order statistic is a function of Euclidean distance in metres.
STUDY_CRS      <- 32749L

## --- Analysis parameters --------------------------------------------------
KDE_SIGMAS     <- c("bw.diggle", "bw.ppl", "bw.CvL", "bw.scott")
NSIM_ENVELOPE  <- 199L   # Monte Carlo simulations -> 0.01 two-sided alpha
RMAX_FRACTION  <- 0.25   # r range capped at a quarter of the window's shortest side
SEED           <- 626

## --- Paths ----------------------------------------------------------------
DIR_RAW     <- "data/raw"
DIR_DERIVED <- "data/derived"
