# ---------------------------------------------------------------------------
# Shared plotting helpers. One theme and one palette for every figure, so the
# maps and the statistical graphics read as one system rather than as a set of
# library defaults. The structure follows Take-home Exercise 1.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(sf); library(scales)
})

# Page tokens follow the coursework site. Violence is encoded in a single
# crimson ramp; the site's teal is reserved for the opposite pole (cold spots,
# Low-Low clusters) and for reference marks. Red and teal also stay
# distinguishable for the common forms of colour-vision deficiency.
PAL <- list(
  ink        = "#16181D",
  ground     = "#F7F8FA",
  ground_alt = "#EEF1F3",
  rule       = "#E5E7EB",
  muted      = "#667085",
  null       = "#E4E6EA",   # not significant: present, but recessive
  border     = "#FFFFFF",
  accent     = "#A3212A",   # violence
  accent_dk  = "#5E1016",
  accent_lt  = "#E6A19A",
  cool       = "#137C8B",   # site accent: the opposite pole
  cool_lt    = "#9CCFD3",
  seq        = c("#F7F8FA", "#F4D9D3", "#EBAEA4", "#DA7A6F", "#C04A43", "#9A2228", "#5E1016")
)

LISA_COLS <- function() c(
  "High-High" = PAL$accent, "Low-High" = PAL$accent_lt,
  "High-Low" = PAL$cool_lt, "Low-Low" = PAL$cool, "Not significant" = PAL$null)

GI_COLS <- function() c("Hot spot" = PAL$accent, "Cold spot" = PAL$cool,
                        "Not significant" = PAL$null)

# EHSA classes, ordered from the most to the least operationally urgent.
EHSA_LEVELS <- c("intensifying hotspot", "persistent hotspot", "consecutive hotspot",
                 "new hotspot", "oscilating hotspot", "sporadic hotspot",
                 "diminishing hotspot", "historical hotspot",
                 "new coldspot", "consecutive coldspot", "intensifying coldspot",
                 "persistent coldspot", "diminishing coldspot", "sporadic coldspot",
                 "oscilating coldspot", "historical coldspot", "no pattern detected")
EHSA_COLS <- function() c(
  "intensifying hotspot" = PAL$accent_dk, "persistent hotspot" = PAL$accent,
  "consecutive hotspot" = "#D0543F", "new hotspot" = "#E98B3A",
  "oscilating hotspot" = "#DE8A92", "sporadic hotspot" = PAL$accent_lt,
  "diminishing hotspot" = "#B98A86", "historical hotspot" = "#CDB9B6",
  "new coldspot" = "#5FB3BC", "consecutive coldspot" = "#3E9AA5",
  "intensifying coldspot" = "#0B4F59", "persistent coldspot" = PAL$cool,
  "diminishing coldspot" = "#8AB2B6", "sporadic coldspot" = PAL$cool_lt,
  "oscilating coldspot" = "#B7D8DB", "historical coldspot" = "#C9D9DB",
  "no pattern detected" = PAL$null)
# sfdep spells "oscillating" with one l; labels are corrected for display only.
ehsa_label <- function(x) sub("oscilating", "oscillating", x)

theme_ex02 <- function(base_size = 11, grid = "y") {
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
  theme_ex02(base_size, grid = "none") +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank())
}

# Myanmar is tall and narrow: 1,900 km north to south, 900 km across. Maps are
# drawn portrait and small multiples are laid out in a row, so each township
# stays large enough to read.
map_base <- function(states, bbox = NULL) {
  list(
    geom_sf(data = states, fill = NA, colour = PAL$ink, linewidth = 0.28),
    coord_sf(expand = FALSE, datum = NA),
    theme_map()
  )
}

# A scale bar drawn as data, so it inherits the palette and survives resizing.
scalebar_layer <- function(bbox, km = 200, pad = 0.05, colour = PAL$ink) {
  w <- as.numeric(bbox["xmax"] - bbox["xmin"]); h <- as.numeric(bbox["ymax"] - bbox["ymin"])
  x0 <- as.numeric(bbox["xmin"]) + pad * w; y0 <- as.numeric(bbox["ymin"]) + pad * h
  list(
    annotate("segment", x = x0, xend = x0 + km * 1000, y = y0, yend = y0,
             colour = colour, linewidth = 0.6),
    annotate("text", x = x0 + km * 500, y = y0 + 0.012 * h,
             label = paste0(km, " km"), colour = colour, size = 2.6,
             family = "Helvetica", vjust = 0)
  )
}

north_layer <- function(bbox, pad = 0.05, colour = PAL$ink) {
  w <- as.numeric(bbox["xmax"] - bbox["xmin"]); h <- as.numeric(bbox["ymax"] - bbox["ymin"])
  x0 <- as.numeric(bbox["xmax"]) - pad * w; y0 <- as.numeric(bbox["ymax"]) - 0.08 * h
  list(
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.03 * h,
             colour = colour, linewidth = 0.6,
             arrow = grid::arrow(length = unit(4, "pt"), type = "closed")),
    annotate("text", x = x0, y = y0 + 0.037 * h, label = "N",
             colour = colour, size = 2.6, family = "Helvetica", vjust = 0)
  )
}

sqrt_breaks <- function(v, n = 5) {
  hi <- max(v, na.rm = TRUE)
  raw <- (seq(0, 1, length.out = n)^2) * hi
  br <- unique(signif(raw, 1))
  br[br <= hi]
}

scale_fill_count <- function(name = "Events", v = NULL, ...) {
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
fmt_pct <- function(x, d = 0) paste0(format(round(100 * x, d), nsmall = d), "%")
