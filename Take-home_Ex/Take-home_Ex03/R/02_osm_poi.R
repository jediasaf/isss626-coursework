# ---------------------------------------------------------------------------
# 02_osm_poi.R - points of interest for Singapore from the Geofabrik extract
#
# The extract covers Malaysia, Singapore and Brunei (about 250 MB). GDAL's OSM
# driver reads it directly; a bounding-box filter keeps Singapore and a strip
# of Johor, which is clipped away against the Master Plan 2019 subzones.
#
# Each POI is the point itself or the centroid of a mapped area (many shops,
# schools and malls are drawn as building outlines, not points), classified
# into a few broad groups from its OSM tags.
# Run from Take-home_Ex/Take-home_Ex03/:  Rscript R/02_osm_poi.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
sf_use_s2(FALSE)

pbf <- "data/raw/malaysia-singapore-brunei-latest.osm.pbf"
out <- "data/derived/osm_poi_sg.gpkg"
if (!file.exists(out)) {
  bb <- c("-spat", "103.59", "1.15", "104.10", "1.48")
  tmp <- tempfile(fileext = ".gpkg")
  sf::gdal_utils("vectortranslate", pbf, tmp,
                 options = c(bb, "-nln", "points", "points"),
                 config_options = c(OSM_MAX_TMPFILE_SIZE = "1024"))
  sf::gdal_utils("vectortranslate", pbf, tmp,
                 options = c(bb, "-update", "-nln", "areas",
                             "-where", "amenity IS NOT NULL OR shop IS NOT NULL OR office IS NOT NULL OR tourism IS NOT NULL OR leisure IS NOT NULL",
                             "multipolygons"),
                 config_options = c(OSM_MAX_TMPFILE_SIZE = "1024"))

  tag <- function(other, key) str_match(other, paste0('"', key, '"=>"([^"]*)"'))[, 2]
  pts <- st_read(tmp, "points", quiet = TRUE) |>
    mutate(amenity = tag(other_tags, "amenity"), shop = tag(other_tags, "shop"),
           office = tag(other_tags, "office"), tourism = tag(other_tags, "tourism"),
           leisure = tag(other_tags, "leisure"), healthcare = tag(other_tags, "healthcare"),
           source = "point") |>
    filter(!is.na(amenity) | !is.na(shop) | !is.na(office) | !is.na(tourism) |
             !is.na(leisure) | !is.na(healthcare))
  areas <- st_read(tmp, "areas", quiet = TRUE) |>
    mutate(healthcare = tag(other_tags, "healthcare"), source = "area") |>
    st_make_valid() |> st_transform(3414) |> st_point_on_surface() |> st_transform(4326)
  keep <- c("osm_id", "name", "amenity", "shop", "office", "tourism", "leisure", "healthcare", "source")
  areas$osm_id <- coalesce(areas$osm_id, areas$osm_way_id)
  poi <- bind_rows(select(pts, any_of(keep)), select(areas, any_of(keep))) |>
    st_transform(3414)

  sg <- st_read("data/raw/MasterPlan2019SubzoneBoundaryNoSea.geojson", quiet = TRUE) |>
    st_transform(3414) |> st_union() |> st_make_valid()
  poi <- poi[lengths(st_intersects(poi, sg)) > 0, ]

  food <- c("restaurant", "fast_food", "cafe", "food_court", "bar", "pub", "ice_cream", "bbq")
  edu <- c("school", "college", "university", "kindergarten", "childcare", "library", "prep_school")
  health <- c("clinic", "hospital", "doctors", "dentist", "pharmacy")
  civic <- c("place_of_worship", "community_centre", "police", "post_office", "townhall", "bank", "marketplace")
  poi <- poi |> mutate(group = case_when(
    !is.na(shop) ~ "retail",
    amenity %in% food ~ "food",
    !is.na(office) ~ "office",
    amenity %in% edu ~ "education",
    amenity %in% health | !is.na(healthcare) ~ "health",
    amenity %in% civic ~ "civic",
    !is.na(leisure) | !is.na(tourism) ~ "leisure",
    TRUE ~ "other"))
  write_sf(poi, out, delete_dsn = TRUE)
}
poi <- st_read(out, quiet = TRUE)
print(count(st_drop_geometry(poi), group, sort = TRUE))

# Bus stops from OpenStreetMap: a fallback only, used if LTA's
# BusStopLocation file is not in data/raw.
bus_out <- "data/derived/osm_busstops_sg.gpkg"
if (!file.exists(bus_out)) {
  tmp <- tempfile(fileext = ".gpkg")
  sf::gdal_utils("vectortranslate", pbf, tmp,
                 options = c("-spat", "103.59", "1.15", "104.10", "1.48",
                             "-where", "highway = 'bus_stop'", "-nln", "bus", "points"))
  sg <- st_read("data/raw/MasterPlan2019SubzoneBoundaryNoSea.geojson", quiet = TRUE) |>
    st_transform(3414) |> st_union() |> st_make_valid()
  bus <- st_read(tmp, "bus", quiet = TRUE) |> st_transform(3414) |> select(osm_id, name)
  write_sf(bus[lengths(st_intersects(bus, sg)) > 0, ], bus_out, delete_dsn = TRUE)
}
message(nrow(st_read(bus_out, quiet = TRUE)), " OpenStreetMap bus stops")
