# ---------------------------------------------------------------------------
# 00_download.R - fetch every input dataset into data/raw/
#
# Run from Take-home_Ex/Take-home_Ex03/:  Rscript R/00_download.R
# Each file is skipped if it already exists, so the script can be re-run.
#
# Three LTA DataMall static files cannot be fetched by script (the server
# refuses non-browser requests). Download them by hand from
#   https://datamall.lta.gov.sg/content/datamall/en/static-data.html
# and save them in data/raw/ under these names:
#   TrainStation_Sep2026.zip, BusStopLocation_Sep2026.zip, TrainStationExit.zip
#
# Passenger volume by origin-destination train station comes from the DataMall
# API and needs a free account key (request it at
# https://datamall.lta.gov.sg/content/datamall/en/request-for-api.html).
# Put it in ~/.Renviron as DATAMALL_KEY=... before running.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(httr2); library(jsonlite) })

raw <- "data/raw"
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
have <- function(f) file.exists(file.path(raw, f)) && file.size(file.path(raw, f)) > 0
fetch <- function(url, f, timeout = 600) {
  if (have(f)) return(invisible(message("have ", f)))
  request(url) |> req_timeout(timeout) |> req_perform(path = file.path(raw, f))
  message("got  ", f, sprintf(" (%.1f MB)", file.size(file.path(raw, f)) / 1e6))
}

# 1. data.gov.sg ---------------------------------------------------------------
# poll-download returns a short-lived signed URL. The open API allows a few
# requests a minute without a key, so retry politely on 429.
datagov <- function(id, f) {
  if (have(f)) return(invisible(message("have ", f)))
  url <- paste0("https://api-open.data.gov.sg/v1/public/api/datasets/", id, "/poll-download")
  for (attempt in 1:6) {
    j <- request(url) |> req_error(is_error = \(r) FALSE) |> req_perform() |> resp_body_json()
    if (!is.null(j$data$url)) return(fetch(j$data$url, f))
    Sys.sleep(15)
  }
  stop("data.gov.sg did not return a download link for ", id)
}
datagov("d_16b157c52ed637edd6ba1232e026258d", "HDBExistingBuilding.geojson")
datagov("d_17f5382f26140b1fdae0ba2ef6239d2f", "HDBPropertyInformation.csv")
datagov("d_8594ae9ff96d0c708bc2af633048edfb", "MasterPlan2019SubzoneBoundaryNoSea.geojson")

# 2. Station coordinates, operational and planned -------------------------------
# Pinned to a commit so the station list cannot change under the analysis.
sha <- "9f5972c71db8a529801eee9c527946a45a78e409"
for (f in c("operational_stations.csv", "future_stations.csv")) {
  fetch(paste0("https://raw.githubusercontent.com/elliotwutingfeng/",
               "singapore_train_station_coordinates/", sha, "/", f), f)
}

# 3. OpenStreetMap (Geofabrik extract, about 250 MB) -----------------------------
fetch("https://download.geofabrik.de/asia/malaysia-singapore-brunei-latest.osm.pbf",
      "malaysia-singapore-brunei-latest.osm.pbf", timeout = 1800)

# 4. LTA DataMall static files (manual, see header) -------------------------------
for (f in c("TrainStation_Sep2026.zip", "BusStopLocation_Sep2026.zip", "TrainStationExit.zip")) {
  if (!have(f)) message("MISSING ", f, ": download it by hand from the DataMall static data page")
}

# 5. Passenger volume by origin-destination train station ------------------------
# DataMall serves only the three most recent months. Each call returns a link
# to a zipped CSV for one month.
key <- Sys.getenv("DATAMALL_KEY")
if (!nzchar(key)) {
  message("DATAMALL_KEY not set: skipping origin-destination passenger volume")
} else {
  months <- format(seq(as.Date(cut(Sys.Date(), "month")), by = "-1 month", length.out = 4), "%Y%m")
  # Passenger Volume by Train Stations: tap-ins and tap-outs per station and
  # hour. This is what the model is calibrated on. The origin-destination
  # files below are optional; on 14 Nov 2026 that endpoint answered this key
  # with a quota violation while this one did not.
  for (m in months) {
    f <- paste0("transport_node_train_", m, ".zip")
    if (have(f)) { message("have ", f); next }
    Sys.sleep(5)
    r <- request("https://datamall2.mytransport.sg/ltaodataservice/PV/Train") |>
      req_url_query(Date = m) |> req_headers(AccountKey = key, accept = "application/json") |>
      req_error(is_error = \(r) FALSE) |> req_perform()
    link <- if (resp_status(r) == 200) resp_body_json(r)$value[[1]]$Link else NULL
    if (is.null(link)) { message("no station volume file for ", m, " (HTTP ", resp_status(r), ")"); next }
    fetch(link, f)
  }
  for (m in months) {
    f <- paste0("origin_destination_train_", m, ".zip")
    if (have(f)) { message("have ", f); next }
    # DataMall answers a burst of calls with a 500 "rate limit quota violation";
    # space the calls and retry.
    Sys.sleep(5)
    r <- request("https://datamall2.mytransport.sg/ltaodataservice/PV/ODTrain") |>
      req_url_query(Date = m) |> req_headers(AccountKey = key, accept = "application/json") |>
      req_retry(max_tries = 2, backoff = \(i) 60,
                is_transient = \(r) resp_status(r) %in% c(429, 503)) |>
      req_error(is_error = \(r) FALSE) |> req_perform()
    # The current month is not published until the next one starts (404).
    link <- if (resp_status(r) == 200) resp_body_json(r)$value[[1]]$Link else NULL
    if (is.null(link)) { message("no OD file for ", m); next }
    fetch(link, f)
  }
}
message("done")
