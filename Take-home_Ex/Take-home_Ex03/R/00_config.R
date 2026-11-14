# ---------------------------------------------------------------------------
# ISSS626 Take-home Exercise 3 - shared configuration
# Every parameter that governs the study design lives here, so the report, the
# slides and the scripts cannot drift out of agreement.
# ---------------------------------------------------------------------------

## --- Projection -------------------------------------------------------------
# SVY21 / Singapore TM. Every distance, buffer and area is computed in metres
# in this CRS; inputs arrive in WGS 84 and are transformed on read.
STUDY_CRS <- 3414L

## --- Inputs -----------------------------------------------------------------
RAW <- "data/raw"
DER <- "data/derived"
F_BLDG     <- file.path(RAW, "HDBExistingBuilding.geojson")
F_PROP     <- file.path(RAW, "HDBPropertyInformation.csv")
F_SUBZONE  <- file.path(RAW, "MasterPlan2019SubzoneBoundaryNoSea.geojson")
F_STN_OP   <- file.path(RAW, "operational_stations.csv")
F_STN_FUT  <- file.path(RAW, "future_stations.csv")
F_LTA_STN  <- file.path(RAW, "TrainStation_Sep2026.zip")
F_LTA_EXIT <- file.path(RAW, "TrainStationExit.zip")
F_BUSSTOP  <- file.path(RAW, "BusStopLocation_Sep2026.zip")
F_POI      <- file.path(DER, "osm_poi_sg.gpkg")
PV_GLOB    <- "transport_node_train_*.zip"
OD_GLOB    <- "origin_destination_train_*.zip"

## --- What ridership is modelled ----------------------------------------------
# Weekday tap-ins in the morning peak. Between 6 and 9 am most trips start at
# home, so a station's boardings are drawn from the homes around it: the
# assumption a residential Huff model needs. Over the whole day a CBD station's
# boardings are mostly people going home from work, whom no residential
# catchment can explain. The all-day series is kept as a sensitivity check.
DAY_TYPE <- "WEEKDAY"
AM_HOURS <- 6:8          # TIME_PER_HOUR 6, 7 and 8, i.e. 06:00 to 08:59

## --- Origins ------------------------------------------------------------------
# HDB residential blocks, weighted by dwelling units (the local market
# potential). Private housing is not in either HDB file; see the limitations.
MIN_DU <- 1L

## --- Choice set ---------------------------------------------------------------
# A block chooses among the stations within CHOICE_RADIUS_M of it, and always
# among at least its K_MIN nearest, so blocks far from the network (Tengah)
# keep a choice. 3 km covers a walk plus a short feeder-bus ride.
CHOICE_RADIUS_M <- 3000
K_MIN <- 3L
# Distances below this are floored, so a block beside a station does not get
# an effectively infinite pull under a power decay.
MIN_DIST_M <- 100

## --- Which stations the model is judged on ----------------------------------------
# The market potential is HDB dwelling units, so the model can only explain
# boardings at stations with HDB homes around them. Stations with fewer than
# this many HDB dwelling units within 1 km (private estates, industry, the
# airport, the central area) stay in every block's choice set, so they still
# compete for HDB riders, but their observed boardings do not enter the fit.
CAL_MIN_DU_1KM <- 1000

## --- Attractiveness -----------------------------------------------------------
# Counted around each station: about a five-minute walk.
CATCH_M <- 400
POI_GROUPS <- c("retail", "food", "office", "education", "health", "civic")

## --- Planned lines --------------------------------------------------------------
# The brief's 2030 scenario: the Jurong Region Line (JE, JS and JW codes) and
# the Cross Island Line phase 1 (CR2 Aviation Park to CR13 Bright Hill).
JRL_PREFIX <- c("JE", "JS", "JW")
CRL1_CODES <- paste0("CR", 2:13)

## --- Line colours (official LTA line colours, used only to mark lines) ----------
LINE_COLS <- c(NS = "#D42E12", EW = "#009645", CG = "#009645", NE = "#9900AA",
               CC = "#FA9E0D", CE = "#FA9E0D", DT = "#005EC4", TE = "#9D5B25",
               JRL = "#0099AA", CRL = "#97C616", LRT = "#748477")

SEED <- 626L
