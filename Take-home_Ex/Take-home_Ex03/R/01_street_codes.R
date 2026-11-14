# ---------------------------------------------------------------------------
# 01_street_codes.R - link HDB Property Information to HDB Existing Building
#
# The two data.gov.sg files describe the same blocks but share no key:
#   Property Information: blk_no, street (abbreviated name), dwelling units
#   Existing Building:    BLK_NO, ST_COD (a street code), POSTAL_COD, footprint
# ST_COD is one code per street (657 street names, 660 codes), so one geocoded
# block per street is enough to learn the street -> ST_COD crosswalk. Each
# street's block is looked up in the OneMap search API, its postal code is
# matched to a footprint, and that footprint's ST_COD becomes the street's code.
#
# About 660 API calls, spaced 1.2 s apart; results are cached in
# data/derived/street_codes.csv, so a re-run only fetches what is missing.
# Run from Take-home_Ex/Take-home_Ex03/:  Rscript R/01_street_codes.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(tidyverse); library(httr2)
})

bldg <- st_read("data/raw/HDBExistingBuilding.geojson", quiet = TRUE) |>
  st_drop_geometry() |> select(BLK_NO, ST_COD, POSTAL_COD)
prop <- read_csv("data/raw/HDBPropertyInformation.csv", show_col_types = FALSE)

cache_f <- "data/derived/street_codes.csv"
dir.create("data/derived", showWarnings = FALSE)
cache <- if (file.exists(cache_f)) read_csv(cache_f, show_col_types = FALSE,
                                            col_types = cols(.default = "c")) else
  tibble(street = character(), blk_no = character(), postal = character(),
         road_name = character(), ST_COD = character())

onemap <- function(q) {
  r <- request("https://www.onemap.gov.sg/api/common/elastic/search") |>
    req_url_query(searchVal = q, returnGeom = "Y", getAddrDetails = "Y", pageNum = 1) |>
    req_retry(max_tries = 4, backoff = \(i) 5 * i,
              is_transient = \(resp) resp_status(resp) %in% c(429, 500, 502, 503)) |>
    req_error(is_error = \(resp) FALSE) |>
    req_perform()
  if (!grepl("json", resp_content_type(r))) return(NULL)
  resp_body_json(r)$results
}

todo <- setdiff(unique(prop$street), cache$street)
message(length(todo), " streets to look up")

for (s in todo) {
  # Try up to three of the street's blocks until one resolves to a footprint.
  blks <- head(prop$blk_no[prop$street == s], 3)
  hit <- NULL
  for (b in blks) {
    Sys.sleep(1.2)
    res <- onemap(paste(b, s))
    if (is.null(res)) next
    for (x in res) {
      if (!identical(x$BLK_NO, b) || x$POSTAL == "NIL") next
      m <- bldg[bldg$POSTAL_COD == x$POSTAL & bldg$BLK_NO == b, ]
      if (nrow(m) == 1) { hit <- tibble(street = s, blk_no = b, postal = x$POSTAL,
                                        road_name = x$ROAD_NAME, ST_COD = m$ST_COD); break }
    }
    if (!is.null(hit)) break
  }
  if (is.null(hit)) hit <- tibble(street = s, blk_no = NA, postal = NA, road_name = NA, ST_COD = NA)
  cache <- bind_rows(cache, hit)
  write_csv(cache, cache_f)
}

message(sum(!is.na(cache$ST_COD)), " of ", nrow(cache), " streets matched to a street code")

# ---------------------------------------------------------------------------
# Second pass: residential blocks the crosswalk could not place on a footprint
# (two streets OneMap could not resolve, and blocks newer than the footprint
# file, mostly in Tengah). Each is geocoded on its own and kept as a point.
# Cached in data/derived/block_points.csv.
# ---------------------------------------------------------------------------
res <- prop |> filter(residential == "Y") |>
  left_join(select(cache, street, ST_COD), by = "street") |>
  left_join(bldg, by = c("blk_no" = "BLK_NO", "ST_COD"))
miss <- res |> filter(is.na(POSTAL_COD)) |> distinct(blk_no, street)

pts_f <- "data/derived/block_points.csv"
pts <- if (file.exists(pts_f)) read_csv(pts_f, show_col_types = FALSE, col_types = cols(.default = "c")) else
  tibble(blk_no = character(), street = character(), postal = character(), x = character(), y = character())
todo2 <- anti_join(miss, pts, by = c("blk_no", "street"))
message(nrow(todo2), " unplaced blocks to geocode")
for (k in seq_len(nrow(todo2))) {
  b <- todo2$blk_no[k]; s <- todo2$street[k]
  Sys.sleep(1.2)
  res_k <- onemap(paste(b, s))
  hit <- Filter(\(x) identical(x$BLK_NO, b), res_k %||% list())
  row <- if (length(hit)) tibble(blk_no = b, street = s, postal = hit[[1]]$POSTAL,
                                 x = hit[[1]]$X, y = hit[[1]]$Y) else
    tibble(blk_no = b, street = s, postal = NA, x = NA, y = NA)
  pts <- bind_rows(pts, row)
  write_csv(pts, pts_f)
}
message(sum(!is.na(pts$x)), " of ", nrow(pts), " unplaced blocks geocoded to a point")
