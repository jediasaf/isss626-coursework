# ---------------------------------------------------------------------------
# Shared plotting helpers. One theme, one palette, applied to every figure so
# the statistical graphics and the maps read as a single system rather than as
# a collection of library defaults.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(sf); library(scales)
})

# Palette. The page chrome follows the coursework site's own tokens so the report
# does not read as a foreign document dropped into the site. The sequential ramp
# stays warm, deliberately: it encodes fire intensity, and a teal heat map would
# fight its own subject. The site's teal is kept as the secondary reference colour
# for envelopes, reference lines and anything that is not the fire itself.
PAL <- list(
  ink        = "#16181D",   # site --text-primary
  ground     = "#F7F8FA",   # site --bg
  ground_alt = "#EEF1F3",   # site --code-bg
  rule       = "#E5E7EB",   # site --border
  muted      = "#667085",   # site --text-secondary
  accent     = "#C2410C",   # warm, semantic: fire
  accent_dk  = "#7C2408",
  accent_lt  = "#F3A26A",
  cool       = "#137C8B",   # site --accent, used for reference and contrast
  seq        = c("#F7F8FA", "#F6DFC8", "#F0B97F", "#E2803C", "#C2410C", "#8A2A06", "#4A1403")
)

theme_ex01 <- function(base_size = 11, grid = "y") {
  th <- theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(family = "Georgia", size = rel(1.35),
                                   colour = PAL$ink, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = PAL$muted, size = rel(0.95),
                                   margin = margin(b = 12), lineheight = 1.25),
      plot.caption  = element_text(colour = PAL$muted, size = rel(0.8), hjust = 0,
                                   margin = margin(t = 10)),
      plot.caption.position = "plot",
      plot.title.position = "plot",
      axis.title    = element_text(colour = PAL$muted, size = rel(0.9)),
      axis.text     = element_text(colour = PAL$muted, size = rel(0.85)),
      legend.title  = element_text(colour = PAL$muted, size = rel(0.85)),
      legend.text   = element_text(colour = PAL$muted, size = rel(0.8)),
      legend.position = "top", legend.justification = "left",
      legend.key.height = unit(8, "pt"), legend.key.width = unit(28, "pt"),
      strip.text    = element_text(colour = PAL$ink, face = "plain", hjust = 0,
                                   size = rel(0.9)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = PAL$rule, linewidth = 0.25),
      plot.background  = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA)
    )
  if (grid == "y") th <- th + theme(panel.grid.major.x = element_blank())
  if (grid == "x") th <- th + theme(panel.grid.major.y = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major = element_blank())
  th
}

theme_map <- function(base_size = 11) {
  theme_ex01(base_size, grid = "none") +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank())
}

# spatstat im -> tidy data frame for geom_raster, dropping cells outside the
# window so the map does not show a rectangle of interpolated nothing.
im_df <- function(im, value = "intensity") {
  d <- as.data.frame(im)
  names(d) <- c("x", "y", value)
  d[is.finite(d[[value]]), , drop = FALSE]
}

# Intensity images are stored per square metre. Everything reported to a reader
# is per 1,000 km2, which is the unit a regency-scale map can actually be read in.
per_1000km2 <- function(v) v * 1e9

# A scale bar drawn as data, so it inherits the palette and survives resizing.
scalebar_layer <- function(bbox, km = 20, pad = 0.06, colour = PAL$ink) {
  w <- as.numeric(bbox["xmax"] - bbox["xmin"]); h <- as.numeric(bbox["ymax"] - bbox["ymin"])
  x0 <- as.numeric(bbox["xmin"]) + pad * w; y0 <- as.numeric(bbox["ymin"]) + pad * h
  list(
    annotate("segment", x = x0, xend = x0 + km * 1000, y = y0, yend = y0,
             colour = colour, linewidth = 0.7),
    annotate("text", x = x0 + km * 500, y = y0 + 0.018 * h,
             label = paste0(km, " km"), colour = colour, size = 2.9,
             family = "Helvetica", vjust = 0)
  )
}

# North arrow, deliberately minimal: a tick and a letter, not a compass rose.
north_layer <- function(bbox, pad = 0.06, colour = PAL$ink) {
  w <- as.numeric(bbox["xmax"] - bbox["xmin"]); h <- as.numeric(bbox["ymax"] - bbox["ymin"])
  x0 <- as.numeric(bbox["xmax"]) - pad * w; y0 <- as.numeric(bbox["ymin"]) + pad * h
  list(
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.035 * h,
             colour = colour, linewidth = 0.7,
             arrow = grid::arrow(length = unit(4, "pt"), type = "closed")),
    annotate("text", x = x0, y = y0 + 0.045 * h, label = "N",
             colour = colour, size = 2.9, family = "Helvetica", vjust = 0)
  )
}

# Square-root colour scales compress the upper range, so evenly spaced breaks
# collide at the top of the legend. These breaks are chosen on the transformed
# scale so the labels are readable and still land on round numbers.
sqrt_breaks <- function(v, n = 5) {
  hi <- max(v, na.rm = TRUE)
  raw <- (seq(0, 1, length.out = n)^2) * hi
  br <- unique(signif(raw, 1))
  br[br <= hi]
}

scale_fill_intensity <- function(name = "Burned cells per 1,000 km\u00b2", v = NULL, ...) {
  ggplot2::scale_fill_gradientn(
    colours = PAL$seq, trans = "sqrt",
    breaks = if (is.null(v)) waiver() else sqrt_breaks(v),
    labels = scales::label_number(big.mark = ","), name = name, ...)
}

fmt_p <- function(p, nperm = NULL) {
  if (is.na(p)) return("--")
  floor_p <- if (!is.null(nperm)) 1 / (nperm + 1) else NULL
  if (!is.null(floor_p) && p <= floor_p) paste0("< ", signif(floor_p, 2))
  else if (p < 0.001) "< 0.001" else format(round(p, 3), nsmall = 3)
}
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
