# ---------------------------------------------------------------------------
# Shared plotting helpers. One theme and one palette for every figure, so the
# maps and the statistical graphics read as one system. The structure follows
# Take-home Exercises 1 and 2.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(sf); library(scales)
})

# Page tokens follow the coursework site. Ridership is one blue-green ramp;
# losses to the new lines are a warm ramp, so gain and loss never share a hue.
# The official line colours (R/00_config.R) mark lines and nothing else.
PAL <- list(
  ink        = "#16181D",
  ground     = "#F7F8FA",
  rule       = "#E5E7EB",
  muted      = "#667085",
  land       = "#ECEEF1",
  land_line  = "#FFFFFF",
  null       = "#D5D9DF",
  accent     = "#137C8B",   # site accent
  accent_lt  = "#9CCFD3",
  loss       = "#B4432F",
  loss_lt    = "#EDB9A6",
  seq        = c("#E6F2F3", "#B9DCDF", "#86C0C6", "#4FA0A9", "#137C8B", "#0B5560", "#06343B"),
  div        = c("#9A2C1E", "#C9644C", "#EDB9A6", "#F2F2F2", "#9CCFD3", "#4FA0A9", "#0B5560")
)

theme_ex03 <- function(base_size = 11, grid = "y") {
  th <- theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(family = "Georgia", size = rel(1.3),
                                   colour = PAL$ink, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = PAL$muted, size = rel(0.92),
                                   margin = margin(b = 10), lineheight = 1.25),
      plot.caption  = element_text(colour = PAL$muted, size = rel(0.78), hjust = 0,
                                   margin = margin(t = 8)),
      plot.caption.position = "plot",
      plot.title.position = "plot",
      axis.title    = element_text(colour = PAL$muted, size = rel(0.88)),
      axis.text     = element_text(colour = PAL$muted, size = rel(0.82)),
      legend.title  = element_text(colour = PAL$muted, size = rel(0.82)),
      legend.text   = element_text(colour = PAL$ink, size = rel(0.78)),
      legend.position = "top", legend.justification = "left",
      legend.key.height = unit(8, "pt"), legend.key.width = unit(22, "pt"),
      strip.text    = element_text(colour = PAL$ink, hjust = 0, size = rel(0.9)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = PAL$rule, linewidth = 0.25),
      plot.background  = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA))
  if (grid == "y") th <- th + theme(panel.grid.major.x = element_blank())
  if (grid == "x") th <- th + theme(panel.grid.major.y = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major = element_blank())
  th
}

theme_map <- function(base_size = 11) {
  theme_ex03(base_size, grid = "none") +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          panel.grid = element_blank())
}

# Land underlay for every map: subzones in a flat tone, white boundaries.
base_land <- function(subzone, lwd = 0.12) {
  geom_sf(data = subzone, fill = PAL$land, colour = PAL$land_line, linewidth = lwd)
}

# Crop every island-wide map to the area with homes and stations, so the
# offshore islands do not shrink the part that carries the information.
sg_crop <- function(x, pad = 1200) {
  b <- sf::st_bbox(x)
  coord_sf(xlim = c(b[["xmin"]] - pad, b[["xmax"]] + pad),
           ylim = c(b[["ymin"]] - pad, b[["ymax"]] + pad), expand = FALSE, datum = NA)
}

fmt_n <- function(x, acc = 1) label_comma(accuracy = acc)(x)
fmt_pct <- function(x, acc = 1) label_percent(accuracy = acc)(x)
