#!/usr/bin/env Rscript
# =============================================================================
# NewFig. 8 | Cxcl10 is confined to the vascular-contact microglial compartment
#             and decays along the perivascular gradient
#
# Core conclusion
#   In one EAE brain section the vascular-contact (EAE_POS / MgIV) and
#   parenchymal (EAE_NEG / MgIVC) programmes partition in opposite directions;
#   vascular-contact microglia sit closest to the Vwf+Cdh5+ large-vessel mask;
#   Cxcl10 detection is confined to that compartment and decays with distance
#   from large vessels; and in an independent pMCAO single-cell data set the
#   same compartment carries a glycolytic shift.
#
# Archetype  asymmetric mixed-modality figure
# Backend    R only (ggplot2 + patchwork + scattermore); no cross-language render
# Size       183 mm wide, <= 247 mm tall; 5 pt rendered-glyph floor
#
# Panel map
#   a  Tcell EAE_POS / EAE_NEG prediction score over the whole section
#   b  Mg    EAE_POS / EAE_NEG prediction score over the whole section
#   c  H&E reference plus the same four scores inside the perivascular window
#   d  within-spot EAE_POS minus EAE_NEG difference, lesion vs extra-lesional
#   e  spatial Cxcl10 in EAE_POS and EAE_NEG microglial spots
#   f  proportion of Cxcl10-positive microglial spots with Wilson intervals
#   g  Cxcl10+ spots on the large-vessel mask, perivascular enrichment and the
#      Cxcl10+ fraction along the distance-to-vessel gradient
#   h  microglia on the vessel mask, the distance-to-vessel map and the
#      cell-to-vessel distance of MgIV against MgIVC
#   i  glycolysis / OXPHOS / balance score by vascular-contact status and the
#      mean glycolysis score per cluster x contact (independent pMCAO scRNA-seq)
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(grid)
  library(scattermore)
  library(RANN)
})
# --- path resolver (GitHub repo IVstaining_scRNA) ---
.iv <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)[1])
  if (!length(f) || is.na(f) || !nzchar(f)) {
    fig <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    list(script = NA_character_, fig = fig, root = NA_character_,
         fig_root = dirname(fig))
  } else {
    sp <- normalizePath(f, winslash = "/", mustWork = TRUE)
    d <- dirname(sp)
    if (basename(d) %in% c("upstream", "qa", "original_scripts")) d <- dirname(d)
    fig <- if (identical(basename(d), "scripts")) dirname(d) else d
    root <- fig
    repeat {
      if (dir.exists(file.path(root, "figures")) &&
          file.exists(file.path(root, "README.md"))) break
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    list(script = sp, fig = fig, root = root,
         fig_root = file.path(root, "figures"))
  }
})
if (!exists("OUTDIR") || is.null(OUTDIR)) OUTDIR <- .iv$fig
IV_ROOT     <- .iv$root
IV_FIG_ROOT <- .iv$fig_root
IV_FIG      <- .iv$fig
# --- end path resolver ---

set.seed(42)

SRC_PRED <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig5_Nature_vector/source_data")
SRC_CX   <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig6_Nature_vector2/source_data")
SRC_TAB  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/analysis_figs_run01")
OUTDIR   <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig8_Nature_v3")
OUT_SD   <- file.path(OUTDIR, "source_data")
QA       <- file.path(OUTDIR, "qa")
for (d in c(OUT_SD, QA, file.path(OUTDIR, "scripts"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

FIG_W_MM <- 183
FAM      <- "Arial"
LABEL_MM <- 5.2 / ggplot2::.pt   # geom_text size that renders at 5.2 pt
SMALL_MM <- 5.0 / ggplot2::.pt   # n-per-bin annotations, at the 5 pt floor

# Row heights in mm. The c row is derived from the true crop aspect so the H&E
# tile and the four score maps stay pixel-registered; every other row is fixed.
# The sum is the page height, so a millimetre here is a millimetre on the page.
H_ROW <- c(ab = 36.5, c = NA, d = 31, ef = 36.5, g = 42, h = 41, i = 29)
C_TITLE_MM <- 4.5   # title strip above a c tile
C_W_MM     <- 35    # width of one c tile
MAP_MM     <- 32    # square section map in rows a, b, e
LEG_MM     <- 9     # colour-bar column
D_MM       <- 140   # panel d is centred rather than stretched to full width
SPAT_MM    <- 34    # square section map in rows g, h

# ---------------------------------------------------------------------------
# One palette for the whole figure.
#   neutral  light-grey tissue silhouette and grey scope reference
#   signal   sequential warm ramp for prediction-score / count magnitude
#   identity red = vascular contact (EAE_POS / MgIV / exposed)
#            blue = parenchyma      (EAE_NEG / MgIVC / parenchymal)
#   accent   orange for Cxcl10, purple for the perivascular lesion windows
# The identity pair is used in every categorical panel (f, h, i) and the ROI
# outlines are purple everywhere they appear (a, b, c, e), so the reader learns
# the vocabulary once.
# ---------------------------------------------------------------------------
seq_cols    <- c("#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#7F0000")
TISSUE_GREY <- "#DCDCDC"
BG_GREY     <- "#E2E5E8"
C_CONTACT   <- "#C0392B"   # vascular contact / EAE_POS / MgIV / exposed
C_PAREN     <- "#3E6D9C"   # parenchyma / EAE_NEG / MgIVC / parenchymal
C_PAREN_LT  <- "#8FA9C4"   # distal distance class in g
C_CXCL10    <- "#E8720C"   # Cxcl10 accent
C_VESSEL    <- "#111111"   # Vwf+ Cdh5+ large-vessel mask
ROI_COL     <- "#6C3483"   # perivascular lesion window
prog_cols   <- c(EAE_POS = C_CONTACT, EAE_NEG = C_PAREN)
scope_cols  <- c(`Extra-lesional` = "#9E9E9E", `Perivascular lesion` = ROI_COL)

fm <- systemfonts::match_fonts(FAM)
cat("font:", FAM, "->", fm$path, "\n")
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM; do.call(reg, args)
  }
}

theme_set(
  theme_classic(base_size = 6.5, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(0.8, "pt"),
      axis.title = element_text(size = 6.5, colour = "black"),
      axis.text = element_text(size = 5.6, colour = "black"),
      legend.title = element_text(size = 5.6),
      legend.text = element_text(size = 5.2),
      legend.key.size = unit(2.2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 5.6, face = "bold", colour = "black")
    )
)

tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.004, 0.99)
)

# --- formatting -------------------------------------------------------------
# Exact P/q values are printed everywhere; the figure carries no asterisk code.
fmt_q <- function(p) {
  if (!is.finite(p)) return("q = NA")
  if (p < 1e-300) return("q < 1e-300")
  sprintf("q = %.0e", p)
}
fmt_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 1e-300) return("P < 1e-300")
  if (p < 1e-4)  return(sprintf("P = %.0e", p))
  sprintf("P = %.3g", p)
}
signed_txt <- function(v, nd = 3) sprintf("%s%.*f", ifelse(v < 0, "\u2212", "+"),
                                          nd, abs(v))

# --- effect sizes -----------------------------------------------------------
# Matched-pairs rank-biserial correlation: the effect size that belongs with a
# Wilcoxon signed-rank test. Positive means EAE_POS exceeds EAE_NEG.
rb_paired <- function(d) {
  d <- d[is.finite(d)]
  nz <- d[d != 0]
  if (!length(nz)) return(0)
  r <- rank(abs(nz))
  (sum(r[nz > 0]) - sum(r[nz < 0])) / sum(r)
}
# Rank-biserial correlation for a rank-sum test; positive means lv1 is higher.
rb_rs <- function(x, g, lv1) {
  r <- rank(x); n1 <- sum(g == lv1); n2 <- sum(g != lv1)
  U <- sum(r[g == lv1]) - n1 * (n1 + 1) / 2
  2 * (U / (n1 * n2)) - 1
}
# Wilson score interval; the normal approximation is unusable at 1% positivity.
wilson <- function(k, n) {
  if (n == 0) return(c(NA_real_, NA_real_, NA_real_))
  p <- k / n; z <- 1.959964; d <- 1 + z^2 / n
  c(p,
    max(0, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d),
    min(1, (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d))
}
wilson_df <- function(k, n) {
  m <- t(mapply(wilson, k, n))
  data.frame(pct = 100 * m[, 1], lo = 100 * m[, 2], hi = 100 * m[, 3])
}
# Two-sided Fisher exact test written so that the odds ratio is the odds of
# detection in the FIRST group relative to the second. table() on a factor plus
# a logical returns the reciprocal, which would invert the reported direction.
fisher2 <- function(k1, n1, k2, n2, nm1, nm2) {
  tb <- matrix(c(k1, n1 - k1, k2, n2 - k2), nrow = 2, byrow = TRUE,
               dimnames = list(group = c(nm1, nm2),
                               outcome = c("positive", "negative")))
  ft <- fisher.test(tb)
  list(tab = tb, or = as.numeric(ft$estimate), lo = ft$conf.int[1],
       hi = ft$conf.int[2], p = ft$p.value)
}

# =============================================================================
# Load and cross-check the two parent data sets
# =============================================================================
cat("Loading source data ...\n")
gj    <- readRDS(file.path(SRC_PRED, "Fig5GJ_spatial_scores.rds"))
les   <- readRDS(file.path(SRC_PRED, "Fig5HJ_lesion_scores.rds"))
he    <- readRDS(file.path(SRC_PRED, "HE_Brain_lowres.rds"))
rects <- readRDS(file.path(SRC_PRED, "Fig5E_roi_rects.rds"))
cx    <- as.data.frame(readRDS(file.path(SRC_CX, "Fig6H_spatial_cxcl10.rds")))
he_cx <- readRDS(file.path(SRC_CX, "Fig6H_HE_image.rds"))
sp    <- read.csv(file.path(SRC_TAB, "spatial_eae_export.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
mc    <- read.csv(file.path(SRC_TAB, "mcao_microglia_export.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(he, he_cx))

score_cols <- c("Tcell EAE_POS", "Tcell EAE_NEG", "Mg EAE_POS", "Mg EAE_NEG")
stopifnot(all(score_cols %in% names(gj)), all(score_cols %in% names(les)))
stopifnot(all(is.finite(as.matrix(gj[, score_cols]))))

# The two halves of this figure come from different exports. They may only be
# combined if they describe the same spots in the same image frame, so that is
# asserted here rather than assumed.
gj$cell <- as.character(gj$cell)
sp$cell_id <- as.character(sp$cell_id)
stopifnot(identical(sort(gj$cell), sort(sp$cell_id)))
sp <- sp[match(gj$cell, sp$cell_id), ]
stopifnot(identical(sp$cell_id, gj$cell))
stopifnot(max(abs(sp$x - gj$x_img)) < 1e-6, max(abs(sp$y - gj$y_img)) < 1e-6)
cat("  frame check: same", nrow(gj), "spots, identical x/y image coordinates\n")

gj$in_lesion <- gj$cell %in% les$cell
cx$cell <- as.character(cx$cell)
cx$in_lesion <- cx$cell %in% les$cell
cx$group_label <- factor(cx$group_label,
                         levels = c("Mg_EAE_POS", "Mg_EAE_NEG"),
                         labels = c("EAE_POS", "EAE_NEG"))
n_spot <- nrow(gj); n_les <- sum(gj$in_lesion); n_out <- n_spot - n_les
cat(sprintf("  prediction spots=%d (lesion=%d, extra-lesional=%d); none excluded\n",
            n_spot, n_les, n_out))

# The microglial EAE_POS / EAE_NEG split and the MgIV / MgIVC cluster split are
# asserted to be the same partition. Without this the figure would use two
# different definitions of "vascular contact" under one colour.
part <- table(cx$group_label,
              sp$cluster[match(cx$cell, sp$cell_id)])
cat("  microglial programme x cluster cross-check:\n"); print(part)
stopifnot(part["EAE_POS", "MgIV"] == sum(cx$group_label == "EAE_POS"),
          part["EAE_NEG", "MgIVC"] == sum(cx$group_label == "EAE_NEG"),
          sum(part) == nrow(cx))
cat("  EAE_POS == MgIV == vascular contact; EAE_NEG == MgIVC == parenchymal\n")
cx_csv <- sp$Cxcl10[match(cx$cell, sp$cell_id)]
stopifnot(identical(as.numeric(cx$Cxcl10), as.numeric(cx_csv)))
cat("  Cxcl10 counts identical in both exports\n")

# Shared colour cap per cell type so the two maps of a pair can be compared
# directly. Capping at the 99th percentile of the pooled pair keeps the ramp off
# a handful of extreme spots without hiding any spot.
cap_tcell <- as.numeric(quantile(c(gj[["Tcell EAE_POS"]], gj[["Tcell EAE_NEG"]]),
                                 0.99, names = FALSE))
cap_mg    <- as.numeric(quantile(c(gj[["Mg EAE_POS"]], gj[["Mg EAE_NEG"]]),
                                 0.99, names = FALSE))
cat(sprintf("  shared colour cap: Tcell %.4f, Mg %.4f\n", cap_tcell, cap_mg))

# =============================================================================
# Map helper. Every spot is drawn first in light grey so the tissue silhouette
# is always visible - without it the panels whose scores sit near zero render as
# an almost blank page and the reader cannot tell a weak programme from missing
# tissue. Spots with a positive score are then drawn on top in ascending order,
# so a sparse hotspot is never overdrawn by the surrounding low values.
# =============================================================================
spat_score <- function(df, score_col, title, cap, show_legend = FALSE,
                       rects_df = NULL, raster = TRUE, pt = 2.0,
                       pixels = c(820, 820), xlim = NULL, ylim = NULL,
                       bar_title = NULL) {
  d <- df
  d$val <- pmin(d[[score_col]], cap)
  d_pos <- d[d$val > 0, , drop = FALSE]
  d_pos <- d_pos[order(d_pos$val), , drop = FALSE]
  cat(sprintf("  %-14s all=%5d positive=%5d cap=%.4f\n",
              score_col, nrow(d), nrow(d_pos), cap))

  p <- ggplot()
  if (raster) {
    p <- p +
      scattermore::geom_scattermore(
        data = d, aes(x_img, -y_img), colour = TISSUE_GREY,
        pointsize = pt, pixels = pixels, interpolate = FALSE) +
      scattermore::geom_scattermore(
        data = d_pos, aes(x_img, -y_img, colour = val),
        pointsize = pt, pixels = pixels, interpolate = FALSE)
  } else {
    p <- p +
      geom_point(data = d, aes(x_img, -y_img), colour = TISSUE_GREY,
                 size = pt, stroke = 0) +
      geom_point(data = d_pos, aes(x_img, -y_img, colour = val),
                 size = pt, stroke = 0)
  }
  p <- p +
    scale_colour_gradientn(colours = seq_cols, limits = c(0, cap),
                           name = bar_title, oob = squish,
                           breaks = pretty_breaks(n = 3)) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = is.null(xlim)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 6.5, base_family = FAM) +
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                family = FAM, lineheight = 1.05),
      legend.title = element_text(size = 5.6, family = FAM),
      legend.text = element_text(size = 5.2, family = FAM),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.6, "mm"),
      plot.margin = margin(0.8, 0.8, 0.8, 0.8)
    )
  if (show_legend) {
    p <- p + guides(colour = guide_colourbar(theme = theme(
      legend.key.width = unit(1.4, "mm"),
      legend.key.height = unit(11, "mm"),
      legend.frame = element_blank(),
      legend.ticks = element_line(linewidth = 0.3, colour = "white")
    ))) + theme(legend.position = "right")
  } else {
    p <- p + guides(colour = "none")
  }
  if (!is.null(rects_df) && nrow(rects_df) > 0) {
    p <- p + geom_rect(
      data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
      inherit.aes = FALSE, fill = NA, colour = ROI_COL,
      linewidth = 0.3, linetype = "dashed")
  }
  p
}

# Histological reference tile. The stored low-resolution H&E is strongly
# saturated, so it is shown once, unmodified, next to the score maps rather than
# placed underneath them where it would swamp the colour ramp. No contrast,
# gamma or channel adjustment is applied.
he_tile <- function(he_crop, nr, nc, title, rects_df = NULL) {
  p <- ggplot() +
    annotation_raster(he_crop, xmin = 0, xmax = nc, ymin = -nr, ymax = 0)
  if (!is.null(rects_df) && nrow(rects_df) > 0) {
    p <- p + geom_rect(
      data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
      inherit.aes = FALSE, fill = NA, colour = ROI_COL,
      linewidth = 0.3, linetype = "dashed")
  }
  p +
    coord_fixed(xlim = c(0, nc), ylim = c(-nr, 0), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 6.5, base_family = FAM) +
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                family = FAM),
      plot.margin = margin(0.8, 0.8, 0.8, 0.8)
    )
}

# Categorical section map used by g and h. The dense all-spot layer is a raster
# mark so every one of the 32,811 spots is kept without emitting 32,811 vector
# circles per panel; the sparse overlays stay vector. The "All spots" key is
# carried by the discrete scale limits rather than by a dummy point.
cat_map <- function(bg, layers, keys, cols, title, bar_px = 100) {
  p <- ggplot() +
    scattermore::geom_scattermore(
      data = bg, aes(x, -y), colour = BG_GREY, pointsize = 1.6,
      pixels = c(820, 820), interpolate = FALSE)
  for (L in layers) {
    p <- p + geom_point(data = L$d, aes(x, -y, colour = L$key),
                        size = L$size, stroke = 0, alpha = L$alpha)
  }
  # Coordinate reference bar in the empty bottom-right corner of the section.
  # It is labelled in low-resolution image pixels: the archived image carries no
  # calibrated pixel size, so a micrometre scale bar cannot be drawn honestly.
  x1 <- max(bg$x) - 0.03 * diff(range(bg$x))
  y0 <- -max(bg$y) + 0.03 * diff(range(bg$y))
  p +
    annotate("segment", x = x1 - bar_px, xend = x1, y = y0, yend = y0,
             linewidth = 0.7, colour = "black", lineend = "butt") +
    annotate("text", x = x1 - bar_px / 2, y = y0, label = "100 px",
             vjust = -0.45, size = SMALL_MM, family = FAM) +
    scale_colour_manual(name = NULL, values = cols, limits = keys,
                        drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 0.9, alpha = 1))) +
    coord_fixed(expand = TRUE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 6.5, base_family = FAM) +
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                family = FAM, lineheight = 1.05),
      legend.text = element_text(size = 5.2, family = FAM),
      # The section is centred in the frame, so the top-left corner is empty.
      # The inside anchor needs legend.justification.inside; plain
      # legend.justification is ignored for inside legends.
      # The section sits centre-high in the frame, so the bottom-left corner is
      # the only region with no spots; the coordinate bar takes bottom-right.
      legend.position = "inside",
      legend.position.inside = c(0.0, 0.0),
      legend.justification.inside = c(0, 0),
      legend.key.spacing.y = unit(0, "mm"),
      legend.key.size = unit(1.6, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0.8, 0.8, 0.8, 0.8)
    )
}

# =============================================================================
# a-b  whole-section prediction score maps
# =============================================================================
cat("Panels a-b ...\n")
mA1 <- spat_score(gj, "Tcell EAE_POS", "Tcell EAE_POS", cap_tcell, rects_df = rects)
mA2 <- spat_score(gj, "Tcell EAE_NEG", "Tcell EAE_NEG", cap_tcell,
                  show_legend = TRUE, rects_df = rects)
mB1 <- spat_score(gj, "Mg EAE_POS", "Mg EAE_POS", cap_mg, rects_df = rects)
mB2 <- spat_score(gj, "Mg EAE_NEG", "Mg EAE_NEG", cap_mg,
                  show_legend = TRUE, rects_df = rects)

pA <- wrap_elements(full = (mA1 | mA2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "a") + tag_theme
pB <- wrap_elements(full = (mB1 | mB2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "b") + tag_theme

# =============================================================================
# c  the same four scores inside the perivascular lesion window
# =============================================================================
cat("Panel c ...\n")
nr_all <- dim(he)[1]; nc_all <- dim(he)[2]
pad <- 15
xmin <- max(0, min(les$x_img) - pad); xmax <- min(nc_all, max(les$x_img) + pad)
ymin <- max(0, min(les$y_img) - pad); ymax <- min(nr_all, max(les$y_img) + pad)
r1 <- max(1, floor(ymin) + 1); r2 <- min(nr_all, ceiling(ymax))
c1 <- max(1, floor(xmin) + 1); c2 <- min(nc_all, ceiling(xmax))
he_crop <- he[r1:r2, c1:c2, , drop = FALSE]
nrc <- dim(he_crop)[1]; ncc <- dim(he_crop)[2]
crop_ar <- ncc / nrc

# Every spot inside the window is shown, not only the lesion spots: with the
# lesion set alone the tile collapses to three disconnected patches floating on
# white, and the reader cannot see that the lesions sit inside continuous
# tissue. The dashed boxes mark the lesion outlines within that tissue.
win_c <- gj %>%
  filter(x_img >= c1 - 1, x_img <= c2, y_img >= r1 - 1, y_img <= r2) %>%
  mutate(x_img = x_img - (c1 - 1), y_img = y_img - (r1 - 1))
rects_c <- transform(rects, xmin = xmin - (c1 - 1), xmax = xmax - (c1 - 1),
                     ymin = ymin - (r1 - 1), ymax = ymax - (r1 - 1))
cat(sprintf("  window %d x %d px (aspect %.3f); %d spots, of which %d lesion\n",
            ncc, nrc, crop_ar, nrow(win_c), sum(win_c$in_lesion)))

# All five tiles share the crop extent, so the H&E reference registers
# pixel-for-pixel with the four score maps beside it. The row height follows
# from the tile width and the crop aspect, so the row can never overflow.
H_ROW[["c"]] <- C_W_MM / crop_ar + C_TITLE_MM
c_pad <- (FIG_W_MM - 5 * C_W_MM) / 2
FIG_H_MM <- round(sum(H_ROW), 1)
stopifnot(FIG_H_MM <= 247)   # Nature full-page height limit
cat(sprintf("  c tile %.1f mm; row %.1f mm; page %.1f x %.1f mm\n",
            C_W_MM, H_ROW[["c"]], FIG_W_MM, FIG_H_MM))

c_xlim <- c(0, ncc); c_ylim <- c(-nrc, 0)
mk_c <- function(col, cap) {
  spat_score(win_c, col, col, cap, rects_df = rects_c, raster = FALSE,
             pt = 0.42, xlim = c_xlim, ylim = c_ylim)
}
cC <- list(
  he_tile(he_crop, nrc, ncc, "H&E", rects_c),
  mk_c("Tcell EAE_POS", cap_tcell),
  mk_c("Tcell EAE_NEG", cap_tcell),
  mk_c("Mg EAE_POS", cap_mg),
  mk_c("Mg EAE_NEG", cap_mg)
)
pC <- wrap_elements(panel = (
  plot_spacer() | cC[[1]] | cC[[2]] | cC[[3]] | cC[[4]] | cC[[5]] | plot_spacer()
) + plot_layout(widths = c(c_pad, rep(C_W_MM, 5), c_pad))) +
  labs(tag = "c") + tag_theme

# =============================================================================
# d  within-spot EAE_POS - EAE_NEG difference
#
# The two scores are measured at the same spot, so the scientific question is
# paired. An unpaired side-by-side comparison of the two marginal distributions
# would also hide that these scores are compositional across the reference
# label set, which sets their absolute magnitudes.
# =============================================================================
cat("Panel d ...\n")
gj$scope <- factor(ifelse(gj$in_lesion, "Perivascular lesion", "Extra-lesional"),
                   levels = c("Extra-lesional", "Perivascular lesion"))
diff_df <- gj %>%
  transmute(scope,
            `T cell` = `Tcell EAE_POS` - `Tcell EAE_NEG`,
            Microglia = `Mg EAE_POS` - `Mg EAE_NEG`) %>%
  pivot_longer(-scope, names_to = "celltype", values_to = "diff")
n_long <- nrow(diff_df)
diff_df <- diff_df %>% filter(is.finite(diff))
cat(sprintf("  paired differences: %d of %d finite\n", nrow(diff_df), n_long))
diff_df$celltype <- factor(diff_df$celltype, levels = c("T cell", "Microglia"))

paired_tab <- diff_df %>%
  group_by(celltype, scope) %>%
  summarise(n = n(), median_diff = median(diff), mean_diff = mean(diff),
            pct_pos_gt_neg = 100 * mean(diff > 0),
            p = wilcox.test(diff, mu = 0, exact = FALSE)$p.value,
            r_rb = rb_paired(diff), .groups = "drop")
paired_tab$q <- p.adjust(paired_tab$p, method = "BH")
cat("  paired signed-rank (EAE_POS - EAE_NEG vs 0):\n")
print(as.data.frame(paired_tab), digits = 4)

shift_tab <- diff_df %>%
  group_by(celltype) %>%
  summarise(median_lesion = median(diff[scope == "Perivascular lesion"]),
            median_outside = median(diff[scope == "Extra-lesional"]),
            p = wilcox.test(diff ~ scope, exact = FALSE)$p.value,
            r_rb = rb_rs(diff, as.character(scope), "Perivascular lesion"),
            .groups = "drop")
shift_tab$q <- p.adjust(shift_tab$p, method = "BH")
cat("  lesion vs extra-lesional shift of the difference:\n")
print(as.data.frame(shift_tab), digits = 4)

mean_d <- diff_df %>% group_by(celltype, scope) %>%
  summarise(m = mean(diff), .groups = "drop")

# Annotation is anchored on the rendered data envelope of each facet, not on a
# fixed y, so the labels cannot land on a violin.
env_d <- diff_df %>% group_by(celltype) %>%
  summarise(lo = min(diff), hi = max(diff), .groups = "drop")
ann_d <- paired_tab %>% left_join(env_d, by = "celltype") %>%
  mutate(label = sprintf("%s\nr_rb = %+.2f", vapply(q, fmt_q, character(1)), r_rb),
         y = hi + 0.10 * (hi - lo))
ann_shift <- shift_tab %>% left_join(env_d, by = "celltype") %>%
  mutate(label = sprintf("lesion vs extra-lesional: %s, r_rb = %+.2f",
                         vapply(q, fmt_q, character(1)), r_rb),
         y = lo - 0.06 * (hi - lo))

pD <- ggplot(diff_df, aes(scope, diff, fill = scope)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
             colour = "grey35") +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.22,
              colour = "black", alpha = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.25,
               fill = "white", colour = "black") +
  geom_point(data = mean_d, aes(scope, m), inherit.aes = FALSE,
             shape = 23, size = 0.9, stroke = 0.25,
             fill = "white", colour = "black") +
  geom_text(data = ann_d, aes(x = scope, y = y, label = label),
            inherit.aes = FALSE, size = LABEL_MM, family = FAM,
            vjust = 0, lineheight = 1.05) +
  geom_text(data = ann_shift, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = LABEL_MM, family = FAM,
            vjust = 1, lineheight = 1.05) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = scope_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.24))) +
  # Two lines: on one line the axis title is taller than the panel and cairo
  # clips both ends of it.
  labs(x = NULL, y = "EAE_POS \u2212 EAE_NEG\nprediction score") +
  theme(
    axis.title.y = element_text(size = 6.5, colour = "black", lineheight = 1.05),
    axis.text.x = element_text(size = 5.6, colour = "black"),
    strip.text = element_text(size = 6.5, face = "bold"),
    panel.spacing = unit(3.5, "mm")
  )

# =============================================================================
# e  spatial Cxcl10 among microglial spots
# =============================================================================
cat("Panel e ...\n")
cap_cx <- max(cx$Cxcl10)
cat(sprintf("  Cxcl10 counts 0-%d; no cap applied (integer UMI counts)\n", cap_cx))
# The two groups differ four-fold in spot number, so n is placed in the title:
# without it the sparser EAE_NEG map could be read as a weaker signal rather
# than as fewer assigned spots.
e_xlim <- range(gj$x_img); e_ylim <- rev(-range(gj$y_img))
mE1 <- spat_score(cx %>% filter(group_label == "EAE_POS"), "Cxcl10",
                  sprintf("EAE_POS (MgIV)\nn = %s spots",
                          format(sum(cx$group_label == "EAE_POS"), big.mark = ",")),
                  cap_cx, rects_df = rects, raster = FALSE, pt = 0.45,
                  xlim = e_xlim, ylim = e_ylim, bar_title = "Cxcl10")
mE2 <- spat_score(cx %>% filter(group_label == "EAE_NEG"), "Cxcl10",
                  sprintf("EAE_NEG (MgIVC)\nn = %s spots",
                          format(sum(cx$group_label == "EAE_NEG"), big.mark = ",")),
                  cap_cx, show_legend = TRUE, rects_df = rects,
                  raster = FALSE, pt = 0.45, xlim = e_xlim, ylim = e_ylim,
                  bar_title = "Cxcl10")
pE <- wrap_elements(full = (mE1 | mE2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "e") + tag_theme

# =============================================================================
# f  Cxcl10 positivity among microglial spots
#
# 95% of microglial spots have zero Cxcl10 counts, so a mean or a crossbar would
# be carried by a handful of spots. The detectable quantity is the proportion of
# positive spots, reported with a Wilson interval and a Fisher odds ratio.
# =============================================================================
cat("Panel f ...\n")
pos_tab <- cx %>% group_by(group_label) %>%
  summarise(n = n(), k = sum(Cxcl10 > 0), .groups = "drop")
pos_tab <- bind_cols(pos_tab, wilson_df(pos_tab$k, pos_tab$n))
cat("  section-wide Cxcl10 positivity:\n"); print(as.data.frame(pos_tab), digits = 4)

fF <- fisher2(pos_tab$k[pos_tab$group_label == "EAE_POS"],
              pos_tab$n[pos_tab$group_label == "EAE_POS"],
              pos_tab$k[pos_tab$group_label == "EAE_NEG"],
              pos_tab$n[pos_tab$group_label == "EAE_NEG"],
              "EAE_POS", "EAE_NEG")
cat("  Fisher 2 x 2 (odds of detection, EAE_POS relative to EAE_NEG):\n")
print(fF$tab)
stopifnot(fF$or > 1)
wtF <- wilcox.test(Cxcl10 ~ group_label, data = cx, exact = FALSE)
rb_cx <- rb_rs(cx$Cxcl10, as.character(cx$group_label), "EAE_POS")
cat(sprintf("  Fisher OR=%.2f (%.2f-%.2f) P=%.3g; rank-sum P=%.3g r_rb=%+.3f\n",
            fF$or, fF$lo, fF$hi, fF$p, wtF$p.value, rb_cx))

# The perivascular window holds only 3 EAE_NEG microglial spots, so it is
# reported as a number rather than drawn as an underpowered second comparison.
les_cx <- cx %>% filter(in_lesion) %>% group_by(group_label) %>%
  summarise(n = n(), k = sum(Cxcl10 > 0), .groups = "drop")
cat("  inside the perivascular window (reported in text, not plotted):\n")
print(as.data.frame(les_cx))

f_lab <- sprintf("Fisher OR = %.1f\n95%% CI %.1f\u2013%.1f\n%s",
                 fF$or, fF$lo, fF$hi, fmt_p(fF$p))
pF <- ggplot(pos_tab, aes(group_label, pct, fill = group_label)) +
  geom_col(width = 0.5, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.3) +
  geom_text(aes(y = hi, label = sprintf("%.1f%%", pct)),
            vjust = -0.55, size = LABEL_MM, family = FAM) +
  annotate("text", x = 1.5, y = 9.1, label = f_lab,
           size = LABEL_MM, family = FAM, lineheight = 1.1, vjust = 1) +
  scale_fill_manual(values = prog_cols, guide = "none") +
  # The positive/total counts ride in the tick labels rather than as text at
  # y = 0, where they collided with the axis line at this panel height.
  scale_x_discrete(labels = setNames(
    sprintf("%s\n%s / %s spots",
            c(EAE_POS = "EAE_POS (MgIV)", EAE_NEG = "EAE_NEG (MgIVC)")[
              as.character(pos_tab$group_label)],
            format(pos_tab$k, big.mark = ",", trim = TRUE),
            format(pos_tab$n, big.mark = ",", trim = TRUE)),
    as.character(pos_tab$group_label))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.04)),
                     limits = c(0, 9.2), breaks = seq(0, 8, 2)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Cxcl10-positive\nmicroglial spots (%)",
       title = "Cxcl10 detection\n(whole section)") +
  theme(
    axis.text.x = element_text(size = 5.6, colour = "black"),
    axis.title.y = element_text(size = 6.5, lineheight = 1.05),
    plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                              lineheight = 1.05),
    plot.margin = margin(1.2, 2.5, 4.5, 2.5)
  ) +
  labs(tag = "f") + tag_theme

# =============================================================================
# Distance to the large-vessel mask
#
# The mask is the Vwf+ Cdh5+ double-positive spot set. It is deliberately
# sparse: 123 of 32,811 spots. Every distance below is the Euclidean distance
# from a spot to its nearest mask spot in low-resolution image pixels. The
# archived image carries no calibrated pixel size, so distances stay in pixels
# and no micrometre scale bar is drawn anywhere in the figure.
# =============================================================================
cat("Distance to vessel mask ...\n")
NEAR_PX <- 25
sp$is_vessel <- (sp$Vwf > 0) & (sp$Cdh5 > 0)
ves_xy <- as.matrix(sp[sp$is_vessel, c("x", "y")])
stopifnot(nrow(ves_xy) > 0)
sp$dist <- as.numeric(RANN::nn2(ves_xy, as.matrix(sp[, c("x", "y")]), k = 1)$nn.dists)
stopifnot(all(is.finite(sp$dist)), min(sp$dist) >= 0)
sp$cx_pos <- sp$Cxcl10 > 0
sp$near <- sp$dist <= NEAR_PX
cat(sprintf("  vessel mask %d of %d spots (%.2f%%)\n",
            sum(sp$is_vessel), nrow(sp), 100 * mean(sp$is_vessel)))
cat(sprintf("  Cxcl10+ spots %d (%.2f%%); perivascular (<= %d px) %d, distal %d\n",
            sum(sp$cx_pos), 100 * mean(sp$cx_pos), NEAR_PX,
            sum(sp$near), sum(!sp$near)))

# =============================================================================
# g  Cxcl10 co-localises with the large-vessel mask and decays with distance
# =============================================================================
cat("Panel g ...\n")
g_keys <- c("All spots", "Cxcl10+ spot", "Vessel mask")
pG1 <- cat_map(
  bg = sp,
  layers = list(
    list(d = sp[sp$cx_pos, ], key = "Cxcl10+ spot", size = 0.30, alpha = 0.85),
    list(d = sp[sp$is_vessel, ], key = "Vessel mask", size = 0.42, alpha = 0.95)
  ),
  keys = g_keys,
  cols = c("All spots" = BG_GREY, "Cxcl10+ spot" = C_CXCL10,
           "Vessel mask" = C_VESSEL),
  title = "Cxcl10+ on vessel mask"
)

# g2  perivascular enrichment. Two-sided Fisher with an odds ratio and a 95%
# CI, and Wilson intervals on the proportions, so the uncertainty definition is
# the same as in panel f.
near_tab <- sp %>%
  mutate(cls = factor(ifelse(near, "perivascular", "distal"),
                      levels = c("perivascular", "distal"))) %>%
  group_by(cls) %>% summarise(n = n(), k = sum(cx_pos), .groups = "drop")
near_tab <- bind_cols(near_tab, wilson_df(near_tab$k, near_tab$n))
cat("  perivascular vs distal Cxcl10 positivity:\n")
print(as.data.frame(near_tab), digits = 4)
fG <- fisher2(near_tab$k[near_tab$cls == "perivascular"],
              near_tab$n[near_tab$cls == "perivascular"],
              near_tab$k[near_tab$cls == "distal"],
              near_tab$n[near_tab$cls == "distal"],
              "perivascular", "distal")
print(fG$tab)
cat(sprintf("  Fisher OR=%.2f (%.2f-%.2f) P=%.3g\n", fG$or, fG$lo, fG$hi, fG$p))
stopifnot(fG$or > 1)

g2_top <- max(near_tab$hi) * 1.62
pG2 <- ggplot(near_tab, aes(cls, pct, fill = cls)) +
  geom_col(width = 0.5, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.3) +
  geom_text(aes(y = hi, label = sprintf("%.2f%%", pct)), vjust = -0.55,
            size = LABEL_MM, family = FAM) +
  annotate("text", x = 1.5, y = g2_top,
           label = sprintf("Fisher OR = %.2f\n95%% CI %.2f\u2013%.2f\n%s",
                           fG$or, fG$lo, fG$hi, fmt_p(fG$p)),
           size = LABEL_MM, family = FAM, lineheight = 1.1, vjust = 1) +
  scale_fill_manual(values = c(perivascular = C_CXCL10, distal = C_PAREN_LT),
                    guide = "none") +
  scale_x_discrete(labels = c(
    perivascular = sprintf("Perivascular\n\u2264 %d px\nn = %s", NEAR_PX,
                           format(near_tab$n[near_tab$cls == "perivascular"],
                                  big.mark = ",")),
    distal = sprintf("Distal\n> %d px\nn = %s", NEAR_PX,
                     format(near_tab$n[near_tab$cls == "distal"],
                            big.mark = ",")))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.04)),
                     limits = c(0, g2_top)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Cxcl10+ spots (%)",
       title = "Perivascular\nenrichment") +
  theme(axis.text.x = element_text(size = 5.4, colour = "black",
                                   lineheight = 1.15),
        plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                  lineheight = 1.05),
        plot.margin = margin(1.2, 2, 1.2, 1.2))

# g3  Cxcl10+ fraction along the distance gradient. Eight EQUAL-COUNT distance
# bins, not equal-width bins: the distance distribution is strongly
# right-skewed, so equal-width bins put 18,726 spots in the first bin and 31 in
# the last, and the Wilson interval of that last bin alone is wider than the
# whole biological effect. Equal-count bins give every point comparable
# precision and spend the resolution where the spots are, which is the
# perivascular range. Every spot enters exactly one bin and each bin carries a
# Wilson interval, so the uncertainty definition matches f and g2.
NBIN <- 8
q_brk <- as.numeric(quantile(sp$dist, probs = seq(0, 1, length.out = NBIN + 1),
                             names = FALSE))
stopifnot(!any(duplicated(q_brk)))
sp$dbin <- cut(sp$dist, breaks = q_brk, include.lowest = TRUE)
binned <- sp %>% group_by(dbin) %>%
  summarise(n = n(), k = sum(cx_pos), dmid = mean(dist),
            dlo = min(dist), dhi = max(dist), .groups = "drop")
binned <- bind_cols(binned, wilson_df(binned$k, binned$n))
stopifnot(sum(binned$n) == nrow(sp), !any(is.na(sp$dbin)))
cat(sprintf("  %d equal-count distance bins over %.2f-%.2f px:\n",
            NBIN, min(sp$dist), max(sp$dist)))
print(as.data.frame(binned), digits = 4)

sr <- suppressWarnings(cor.test(sp$dist, sp$Cxcl10, method = "spearman",
                                exact = FALSE))
cat(sprintf("  Spearman rho=%.4f P=%.3g (two-sided, n=%d)\n",
            sr$estimate, sr$p.value, nrow(sp)))

g3_top <- max(binned$hi) * 1.34
# Equal-count binning makes the per-bin n identical to within one spot, so the
# eight n labels of the parent panel are replaced by one caption; the freed
# space carries horizontal bars showing the distance span each bin covers,
# which is what the reader otherwise cannot recover from a quantile binning.
pG3 <- ggplot(binned, aes(dmid, pct)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_line(linewidth = 0.5, colour = C_CXCL10) +
  geom_segment(aes(x = dlo, xend = dhi, y = pct, yend = pct),
               linewidth = 0.3, colour = C_CXCL10, alpha = 0.55) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = diff(range(binned$dmid)) / 55, linewidth = 0.3,
                colour = C_CXCL10) +
  geom_point(size = 0.9, colour = C_CXCL10, fill = "white", shape = 21,
             stroke = 0.35) +
  annotate("text", x = max(binned$dhi), y = g3_top,
           label = sprintf("Spearman \u03c1 = %s\n%s\nn = %s spots, %s per bin",
                           signed_txt(sr$estimate), fmt_p(sr$p.value),
                           format(nrow(sp), big.mark = ","),
                           format(min(binned$n), big.mark = ",")),
           hjust = 1, vjust = 1, size = LABEL_MM, family = FAM,
           lineheight = 1.1) +
  scale_y_continuous(limits = c(0, g3_top), breaks = pretty_breaks(n = 4),
                     expand = expansion(mult = c(0, 0.02))) +
  # Square-root x axis, as in h. Equal-count bins concentrate seven of the
  # eight points below 40 px while the last bin span reaches 138 px; on a linear
  # axis the informative perivascular range is squeezed into the left third.
  scale_x_sqrt(breaks = c(0, 5, 10, 20, 40, 80, 140),
               limits = c(0, max(binned$dhi)),
               expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Distance to nearest vessel-mask spot (px, square-root scale)",
       y = "Cxcl10+ spots (%)",
       title = "Cxcl10 decays along the perivascular gradient") +
  theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
        plot.margin = margin(1.2, 2, 1.2, 1.2))

pGrow <- wrap_elements(panel = (pG1 | pG2 | pG3) +
  plot_layout(widths = c(38, 42, 103))) + labs(tag = "g") + tag_theme

# =============================================================================
# h  vascular-contact microglia lie closest to the large-vessel mask
# =============================================================================
cat("Panel h ...\n")
mg <- sp[sp$cluster %in% c("MgIV", "MgIVC"), ]
mg$grp <- factor(mg$cluster, levels = c("MgIV", "MgIVC"))
cat(sprintf("  microglial spots %d (MgIV %d, MgIVC %d)\n", nrow(mg),
            sum(mg$grp == "MgIV"), sum(mg$grp == "MgIVC")))

h_keys <- c("All spots", "MgIVC", "MgIV", "Vessel mask")
pH1 <- cat_map(
  bg = sp,
  layers = list(
    list(d = mg[mg$grp == "MgIVC", ], key = "MgIVC", size = 0.30, alpha = 0.85),
    list(d = mg[mg$grp == "MgIV", ], key = "MgIV", size = 0.30, alpha = 0.75),
    list(d = sp[sp$is_vessel, ], key = "Vessel mask", size = 0.42, alpha = 0.95)
  ),
  keys = h_keys,
  cols = c("All spots" = BG_GREY, "MgIVC" = C_PAREN, "MgIV" = C_CONTACT,
           "Vessel mask" = C_VESSEL),
  title = "Microglia on vessel mask"
)

# h2  distance-to-vessel map. A sequential viridis ramp keeps this quantity
# visually distinct from the warm prediction-score ramp of a-c and e, and it
# stays interpretable in greyscale. The mask spots are the darkest class by
# construction (distance 0), so no separate mask overlay is drawn.
d_cap <- as.numeric(quantile(sp$dist, 0.98, names = FALSE))
cat(sprintf("  distance 0-%.2f px; colour capped at the 98th percentile %.2f px\n",
            max(sp$dist), d_cap))
pH2 <- ggplot() +
  scattermore::geom_scattermore(
    data = sp[order(-sp$dist), ], aes(x, -y, colour = pmin(dist, d_cap)),
    pointsize = 1.8, pixels = c(820, 820), interpolate = FALSE) +
  # Viridis truncated at 0.75: the full ramp ends in yellow, which is almost
  # invisible on the white page, and the far-from-vessel spots are the majority
  # of the section. Dark violet (at the mask) to green (distal) keeps both ends
  # legible, stays sequential in greyscale, and is clearly a different quantity
  # from the warm prediction-score ramp of a-c and e.
  scale_colour_viridis_c(name = "Distance\nto vessel (px)", limits = c(0, d_cap),
                         oob = squish, breaks = pretty_breaks(n = 3),
                         begin = 0, end = 0.75, direction = 1) +
  guides(colour = guide_colourbar(theme = theme(
    legend.key.width = unit(1.4, "mm"), legend.key.height = unit(11, "mm"),
    legend.frame = element_blank(),
    legend.ticks = element_line(linewidth = 0.3, colour = "white")))) +
  coord_fixed(expand = TRUE) +
  labs(title = "Distance-to-vessel map", x = NULL, y = NULL) +
  theme_void(base_size = 6.5, base_family = FAM) +
  theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                  family = FAM, lineheight = 1.05),
        legend.title = element_text(size = 5.6, family = FAM),
        legend.text = element_text(size = 5.2, family = FAM),
        legend.box.spacing = unit(0.6, "mm"),
        legend.position = "right",
        plot.margin = margin(0.8, 0.8, 0.8, 0.8))

# h3  cell-to-vessel distance. Two-sided rank-sum with a rank-biserial effect
# size; the parent used a one-sided test, but the direction is the result being
# claimed and should not be assumed by the test.
wtH <- wilcox.test(dist ~ grp, data = mg, exact = FALSE)
rb_H <- rb_rs(mg$dist, as.character(mg$grp), "MgIV")
med_H <- mg %>% group_by(grp) %>%
  summarise(n = n(), median = median(dist), mean = mean(dist), .groups = "drop")
cat("  cell-to-vessel distance by microglial cluster:\n")
print(as.data.frame(med_H), digits = 4)
cat(sprintf("  rank-sum P=%.3g, r_rb=%+.3f (positive = MgIV further)\n",
            wtH$p.value, rb_H))

# Distances are strongly right-skewed and non-negative, so a square-root axis
# keeps every observation on screen while leaving the informative 0-20 px range
# legible. The transform is valid because min(dist) >= 0 was asserted above.
h3_max <- max(mg$dist)
h3_brk <- c(0, 5, 10, 20, 40, 80, 160, 320)
h3_brk <- h3_brk[h3_brk <= h3_max * 1.02]
pH3 <- ggplot(mg, aes(grp, dist, fill = grp)) +
  geom_point(position = position_jitter(width = 0.17, height = 0, seed = 42),
             size = 0.12, stroke = 0, alpha = 0.12,
             aes(colour = grp), show.legend = FALSE) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.22,
              colour = "black", alpha = 0.85) +
  geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
               fill = "white", colour = "black") +
  geom_point(data = med_H, aes(grp, mean), inherit.aes = FALSE, shape = 23,
             size = 0.9, stroke = 0.25, fill = "white", colour = "black") +
  annotate("text", x = 1.5, y = h3_max * 1.5,
           label = sprintf("Wilcoxon rank-sum, %s\nr_rb = %+.2f (MgIV closer)",
                           fmt_p(wtH$p.value), rb_H),
           size = LABEL_MM, family = FAM, lineheight = 1.1, vjust = 1) +
  scale_fill_manual(values = c(MgIV = C_CONTACT, MgIVC = C_PAREN),
                    guide = "none") +
  scale_colour_manual(values = c(MgIV = C_CONTACT, MgIVC = C_PAREN),
                      guide = "none") +
  scale_x_discrete(labels = c(
    MgIV = sprintf("MgIV (vascular contact)\nn = %s, median %.2f px",
                   format(med_H$n[med_H$grp == "MgIV"], big.mark = ","),
                   med_H$median[med_H$grp == "MgIV"]),
    MgIVC = sprintf("MgIVC (parenchymal)\nn = %s, median %.2f px",
                    format(med_H$n[med_H$grp == "MgIVC"], big.mark = ","),
                    med_H$median[med_H$grp == "MgIVC"]))) +
  scale_y_sqrt(breaks = h3_brk, limits = c(0, h3_max * 1.55),
               expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Distance to vessel (px, square-root scale)",
       title = "Cell-to-vessel distance") +
  theme(axis.text.x = element_text(size = 5.4, colour = "black",
                                   lineheight = 1.15),
        axis.title.y = element_text(size = 6.2, colour = "black"),
        plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
        plot.margin = margin(1.2, 2, 1.2, 1.2))

pHrow <- wrap_elements(panel = (pH1 | pH2 | pH3) +
  plot_layout(widths = c(40, 46, 97))) + labs(tag = "h") + tag_theme

# =============================================================================
# i  glycolytic reprogramming in the vascular-contact compartment
#    (independent pMCAO single-cell data set - a different injury model and a
#     different assay, so this is corroboration, not a replicate)
# =============================================================================
cat("Panel i ...\n")
GLYCOLYTIC <- c("Pkm", "Hk2", "Ldha", "Pgk1", "Eno1", "Aldoa", "Tpi1",
                "Slc2a1", "Gapdh", "Slc2a3", "Hk1", "Pfkl", "Pfkp", "Pfkm")
OXPHOS <- c("Cox5a", "Atp5a1", "Sdha", "Ndufa1", "Ndufb8", "Uqcrc1")
gly <- intersect(GLYCOLYTIC, names(mc))
oxp <- intersect(OXPHOS, names(mc))
gly_missing <- setdiff(GLYCOLYTIC, names(mc))
oxp_missing <- setdiff(OXPHOS, names(mc))
cat(sprintf("  glycolytic genes used %d/%d (absent: %s)\n", length(gly),
            length(GLYCOLYTIC),
            ifelse(length(gly_missing), paste(gly_missing, collapse = ", "), "none")))
cat(sprintf("  OXPHOS genes used %d/%d (absent: %s)\n", length(oxp),
            length(OXPHOS),
            ifelse(length(oxp_missing), paste(oxp_missing, collapse = ", "), "none")))

# Mean z score across the gene set. Guard against a constant gene, which would
# divide by zero and silently produce NaN for every cell.
z_mean <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X)
  sdv <- sqrt(colMeans(sweep(X, 2, mu)^2))
  stopifnot(all(is.finite(mu)))
  if (any(sdv <= 0)) {
    cat("  constant gene dropped from the score:",
        paste(colnames(X)[sdv <= 0], collapse = ", "), "\n")
    X <- X[, sdv > 0, drop = FALSE]; mu <- mu[sdv > 0]; sdv <- sdv[sdv > 0]
  }
  rowMeans(sweep(sweep(X, 2, mu), 2, sdv, "/"))
}
mc$glyc <- z_mean(mc[, gly])
mc$oxph <- z_mean(mc[, oxp])
mc$gmo  <- mc$glyc - mc$oxph
stopifnot(all(is.finite(c(mc$glyc, mc$oxph, mc$gmo))))
mc$contact <- factor(mc$vascular_contact_type,
                     levels = c("parenchymal", "exposed"),
                     labels = c("Parenchymal", "Vascular contact"))
stopifnot(!any(is.na(mc$contact)))
mc$cluster <- factor(mc$cluster, levels = sort(unique(mc$cluster)))
cat(sprintf("  pMCAO microglia %d cells (parenchymal %d, vascular contact %d)\n",
            nrow(mc), sum(mc$contact == "Parenchymal"),
            sum(mc$contact == "Vascular contact")))

i_specs <- data.frame(
  col = c("glyc", "oxph", "gmo"),
  lab = c("Glycolysis score", "OXPHOS score", "Glycolysis \u2212 OXPHOS"),
  stringsAsFactors = FALSE
)
i_tab <- do.call(rbind, lapply(seq_len(nrow(i_specs)), function(i) {
  cl <- i_specs$col[i]
  v <- mc[[cl]]
  wt <- wilcox.test(v ~ mc$contact, exact = FALSE)
  data.frame(score = cl, label = i_specs$lab[i],
             n_par = sum(mc$contact == "Parenchymal"),
             n_vc = sum(mc$contact == "Vascular contact"),
             mean_par = mean(v[mc$contact == "Parenchymal"]),
             mean_vc = mean(v[mc$contact == "Vascular contact"]),
             median_par = median(v[mc$contact == "Parenchymal"]),
             median_vc = median(v[mc$contact == "Vascular contact"]),
             p = wt$p.value,
             r_rb = rb_rs(v, as.character(mc$contact), "Vascular contact"))
}))
i_tab$q <- p.adjust(i_tab$p, method = "BH")
i_tab$delta <- i_tab$mean_vc - i_tab$mean_par
cat("  vascular contact vs parenchymal (two-sided rank-sum, BH over 3 tests):\n")
print(i_tab, digits = 4)

mk_i <- function(k) {
  cl <- i_tab$score[k]
  d <- data.frame(contact = mc$contact, v = mc[[cl]])
  mn <- d %>% group_by(contact) %>% summarise(m = mean(v), .groups = "drop")
  lo <- min(quantile(d$v, 0.002)); hi <- max(quantile(d$v, 0.998))
  ggplot(d, aes(contact, v, fill = contact)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.22,
                colour = "black", alpha = 0.85) +
    geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
                 fill = "white", colour = "black") +
    geom_point(data = mn, aes(contact, m), inherit.aes = FALSE, shape = 23,
               size = 0.9, stroke = 0.25, fill = "white", colour = "black") +
    annotate("text", x = 1.5, y = hi + (hi - lo) * 0.30,
             label = sprintf("%s\n\u0394 = %s",
                             fmt_q(i_tab$q[k]), signed_txt(i_tab$delta[k])),
             size = LABEL_MM, family = FAM, lineheight = 1.1, vjust = 1) +
    scale_fill_manual(values = c(Parenchymal = C_PAREN,
                                 `Vascular contact` = C_CONTACT),
                      guide = "none") +
    scale_x_discrete(labels = c(
      Parenchymal = sprintf("Parenchymal\nn = %s\nmean %s",
                            format(i_tab$n_par[k], big.mark = ","),
                            signed_txt(i_tab$mean_par[k])),
      `Vascular contact` = sprintf("Vascular contact\nn = %s\nmean %s",
                                   format(i_tab$n_vc[k], big.mark = ","),
                                   signed_txt(i_tab$mean_vc[k])))) +
    # coord_cartesian, not scale limits: a scale limit would drop the cells
    # outside the view BEFORE the violin density is computed, which silently
    # changes the distribution being drawn. Here every cell enters the stat and
    # only the view is trimmed to the 0.2-99.8 percentile envelope.
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    coord_cartesian(ylim = c(lo - (hi - lo) * 0.06, hi + (hi - lo) * 0.32),
                    expand = FALSE) +
    labs(x = NULL, y = i_tab$label[k]) +
    theme(axis.text.x = element_text(size = 5.2, colour = "black",
                                     lineheight = 1.2),
          plot.margin = margin(1.2, 1.6, 1.2, 1.2))
}
pI1 <- mk_i(1); pI2 <- mk_i(2); pI3 <- mk_i(3)

# i4  mean glycolysis score per cluster x contact, with the per-cluster
# vascular-contact test. A diverging scale centred on zero with in-cell values
# so the tile survives greyscale print.
piv <- mc %>% group_by(cluster, contact) %>%
  summarise(m = mean(glyc), n = n(), .groups = "drop")
clu_p <- mc %>% group_by(cluster) %>%
  summarise(p = wilcox.test(glyc ~ contact, exact = FALSE)$p.value,
            r_rb = rb_rs(glyc, as.character(contact), "Vascular contact"),
            .groups = "drop")
clu_p$q <- p.adjust(clu_p$p, method = "BH")
cat("  mean glycolysis score per cluster x contact:\n")
print(as.data.frame(piv), digits = 4)
cat("  per-cluster vascular-contact effect (two-sided rank-sum, BH):\n")
print(as.data.frame(clu_p), digits = 4)

i4_lim <- max(abs(piv$m)) * 1.05
pI4 <- ggplot(piv, aes(contact, cluster, fill = m)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", m),
                colour = abs(m) > i4_lim * 0.55),
            size = LABEL_MM, family = FAM, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-i4_lim, i4_lim),
                       name = "Mean\nscore", breaks = pretty_breaks(n = 3)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1A1A1A"),
                      guide = "none") +
  scale_x_discrete(labels = c(Parenchymal = "Parenchymal",
                              `Vascular contact` = "Vascular\ncontact"),
                   expand = expansion(0)) +
  scale_y_discrete(limits = rev(levels(piv$cluster)), expand = expansion(0)) +
  guides(fill = guide_colourbar(theme = theme(
    legend.key.width = unit(1.4, "mm"), legend.key.height = unit(9, "mm"),
    legend.frame = element_blank(),
    legend.ticks = element_line(linewidth = 0.3, colour = "white")))) +
  labs(x = NULL, y = "Microglia cluster",
       title = "Glycolysis score:\ncluster \u00d7 contact") +
  theme(axis.text.x = element_text(size = 5.4, colour = "black",
                                   lineheight = 1.15),
        axis.line = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5,
                                  lineheight = 1.05),
        legend.box.spacing = unit(0.8, "mm"),
        plot.margin = margin(1.2, 1.2, 1.2, 1.2))

pIrow <- wrap_elements(panel = (pI1 | pI2 | pI3 | pI4) +
  plot_layout(widths = c(41, 41, 41, 60))) + labs(tag = "i") + tag_theme

# =============================================================================
# Assemble. Row widths are padded so each map cell keeps its true aspect ratio.
# =============================================================================
cat("Assemble ...\n")
ab_cell <- 2 * MAP_MM + LEG_MM
ab_pad  <- max((FIG_W_MM - 2 * ab_cell) / 2, 0.1)
row_ab <- wrap_elements(panel = (
  plot_spacer() | pA | pB | plot_spacer()
) + plot_layout(widths = c(ab_pad, ab_cell, ab_cell, ab_pad)))

# The tag rides on the wrapper, not on pD: inside the padded row it would
# otherwise land on top of the rotated y-axis title.
d_pad <- max((FIG_W_MM - D_MM) / 2, 0.1)
row_d <- wrap_elements(panel = (
  plot_spacer() | pD | plot_spacer()
) + plot_layout(widths = c(d_pad, D_MM, d_pad))) +
  labs(tag = "d") + tag_theme

row_ef <- (pE | pF | plot_spacer()) +
  plot_layout(widths = c(ab_cell, 68, FIG_W_MM - ab_cell - 68))

fig <- row_ab / pC / row_d / row_ef / pGrow / pHrow / pIrow +
  plot_layout(heights = as.numeric(H_ROW))

# =============================================================================
# Source-data tables
# =============================================================================
cat("Writing source data ...\n")
write.csv(gj %>% select(cell, x_img, y_img, in_lesion, all_of(score_cols)),
          file.path(OUT_SD, "NewFig8abcd_spatial_prediction_scores.csv"),
          row.names = FALSE)
write.csv(as.data.frame(paired_tab),
          file.path(OUT_SD, "NewFig8d_paired_signed_rank.csv"), row.names = FALSE)
write.csv(as.data.frame(shift_tab),
          file.path(OUT_SD, "NewFig8d_lesion_vs_extralesional.csv"),
          row.names = FALSE)
write.csv(cx %>% select(cell, x_img, y_img, cluster, group_label, in_lesion,
                        Cxcl10),
          file.path(OUT_SD, "NewFig8ef_spatial_cxcl10.csv"), row.names = FALSE)
write.csv(as.data.frame(pos_tab) %>%
            mutate(fisher_or = fF$or, fisher_ci_lo = fF$lo,
                   fisher_ci_hi = fF$hi, fisher_p = fF$p,
                   ranksum_p = wtF$p.value, ranksum_r_rb = rb_cx),
          file.path(OUT_SD, "NewFig8f_cxcl10_positivity.csv"), row.names = FALSE)
write.csv(as.data.frame(les_cx),
          file.path(OUT_SD, "NewFig8f_cxcl10_positivity_lesion_window.csv"),
          row.names = FALSE)
write.csv(sp %>% select(cell = cell_id, x, y, cluster, vascular_contact_type,
                        Vwf, Cdh5, is_vessel, dist, Cxcl10, cx_pos, near),
          file.path(OUT_SD, "NewFig8gh_vessel_distance_spots.csv"),
          row.names = FALSE)
write.csv(as.data.frame(near_tab) %>%
            mutate(fisher_or = fG$or, fisher_ci_lo = fG$lo,
                   fisher_ci_hi = fG$hi, fisher_p = fG$p),
          file.path(OUT_SD, "NewFig8g_perivascular_enrichment.csv"),
          row.names = FALSE)
write.csv(as.data.frame(binned) %>%
            mutate(spearman_rho = as.numeric(sr$estimate),
                   spearman_p = sr$p.value),
          file.path(OUT_SD, "NewFig8g_distance_gradient_bins.csv"),
          row.names = FALSE)
write.csv(as.data.frame(med_H) %>%
            mutate(ranksum_p = wtH$p.value, ranksum_r_rb = rb_H),
          file.path(OUT_SD, "NewFig8h_cell_to_vessel_distance.csv"),
          row.names = FALSE)
write.csv(mc %>% select(cell_id, cluster, vascular_contact_type, glyc, oxph, gmo),
          file.path(OUT_SD, "NewFig8i_mcao_metabolic_scores.csv"),
          row.names = FALSE)
write.csv(i_tab, file.path(OUT_SD, "NewFig8i_contact_tests.csv"),
          row.names = FALSE)
write.csv(as.data.frame(piv) %>%
            left_join(as.data.frame(clu_p), by = "cluster"),
          file.path(OUT_SD, "NewFig8i_cluster_by_contact.csv"), row.names = FALSE)
write.csv(rects, file.path(OUT_SD, "NewFig8_lesion_roi_rects.csv"),
          row.names = FALSE)
write.csv(data.frame(gene_set = c(rep("glycolysis", length(gly)),
                                  rep("oxphos", length(oxp))),
                     gene = c(gly, oxp)),
          file.path(OUT_SD, "NewFig8i_gene_sets.csv"), row.names = FALSE)

# =============================================================================
# Export
# =============================================================================
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base_path <- file.path(OUTDIR, "NewFigure8")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base_path, ".svg"), width = w_in, height = h_in,
                 bg = "white")
print(fig); dev.off()

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base_path, ".pdf"), width = w_in, height = h_in,
                     family = FAM)
print(fig); dev.off()

cat("Writing PNG ...\n")
ragg::agg_png(paste0(base_path, ".png"), width = w_in, height = h_in,
              units = "in", res = 600, background = "white")
print(fig); dev.off()

# Per-row QA renders at true final size, so each panel can be audited without
# squinting at the whole page. Same backend, same devices, same millimetres.
cat("Writing QA row renders ...\n")
qa_rows <- list(row_ab = list(p = row_ab, h = H_ROW[["ab"]]),
                row_c  = list(p = pC,     h = H_ROW[["c"]]),
                row_d  = list(p = row_d,  h = H_ROW[["d"]]),
                row_ef = list(p = row_ef, h = H_ROW[["ef"]]),
                row_g  = list(p = pGrow,  h = H_ROW[["g"]]),
                row_h  = list(p = pHrow,  h = H_ROW[["h"]]),
                row_i  = list(p = pIrow,  h = H_ROW[["i"]]))
for (nm in names(qa_rows)) {
  ragg::agg_png(file.path(QA, paste0(nm, ".png")), width = FIG_W_MM / 25.4,
                height = qa_rows[[nm]]$h / 25.4, units = "in", res = 600,
                background = "white")
  print(qa_rows[[nm]]$p); dev.off()
}

saveRDS(list(FIG_W_MM = FIG_W_MM, FIG_H_MM = FIG_H_MM, H_ROW = H_ROW,
             crop_ar = crop_ar, cap_tcell = cap_tcell, cap_mg = cap_mg,
             cap_cx = cap_cx, d_cap = d_cap, NEAR_PX = NEAR_PX, NBIN = NBIN,
             paired_tab = paired_tab, shift_tab = shift_tab, pos_tab = pos_tab,
             fF = fF, wtF_p = wtF$p.value, rb_cx = rb_cx, les_cx = les_cx,
             near_tab = near_tab, fG = fG, binned = binned,
             sr = list(rho = as.numeric(sr$estimate), p = sr$p.value),
             med_H = med_H, wtH_p = wtH$p.value, rb_H = rb_H,
             i_tab = i_tab, piv = piv, clu_p = clu_p,
             gly = gly, oxp = oxp, gly_missing = gly_missing,
             n_spot = n_spot, n_les = n_les, n_out = n_out,
             n_win = nrow(win_c), n_win_les = sum(win_c$in_lesion),
             ncc = ncc, nrc = nrc, n_cx = nrow(cx), n_mc = nrow(mc),
             n_vessel = sum(sp$is_vessel), n_mg = nrow(mg)),
        file.path(QA, "figure_state.rds"))

# =============================================================================
# Stats / legend / provenance / QA / README
# =============================================================================
gp <- function(ct, sc, field) {
  paired_tab[[field]][paired_tab$celltype == ct & paired_tab$scope == sc]
}
gs <- function(ct, field) shift_tab[[field]][shift_tab$celltype == ct]
gi <- function(sc, field) i_tab[[field]][i_tab$score == sc]
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

k_pos <- pos_tab$k[pos_tab$group_label == "EAE_POS"]
n_pos <- pos_tab$n[pos_tab$group_label == "EAE_POS"]
k_neg <- pos_tab$k[pos_tab$group_label == "EAE_NEG"]
n_neg <- pos_tab$n[pos_tab$group_label == "EAE_NEG"]
k_near <- near_tab$k[near_tab$cls == "perivascular"]
n_near <- near_tab$n[near_tab$cls == "perivascular"]
k_far  <- near_tab$k[near_tab$cls == "distal"]
n_far  <- near_tab$n[near_tab$cls == "distal"]
pct_pos <- pos_tab$pct[pos_tab$group_label == "EAE_POS"]
pct_neg <- pos_tab$pct[pos_tab$group_label == "EAE_NEG"]
pct_near <- near_tab$pct[near_tab$cls == "perivascular"]
pct_far  <- near_tab$pct[near_tab$cls == "distal"]
n_vessel <- sum(sp$is_vessel)
n_mgiv  <- med_H$n[med_H$grp == "MgIV"]
n_mgivc <- med_H$n[med_H$grp == "MgIVC"]
med_iv  <- med_H$median[med_H$grp == "MgIV"]
med_ivc <- med_H$median[med_H$grp == "MgIVC"]
les_pos_k <- les_cx$k[les_cx$group_label == "EAE_POS"]
les_pos_n <- les_cx$n[les_cx$group_label == "EAE_POS"]
les_neg_k <- les_cx$k[les_cx$group_label == "EAE_NEG"]
les_neg_n <- les_cx$n[les_cx$group_label == "EAE_NEG"]

writeLines(c(
  "NewFig. 8 statistics / integrity",
  "",
  sprintf("Spatial spots: %s in total; %s inside the three perivascular lesion",
          fmt_n(n_spot), fmt_n(n_les)),
  sprintf("  windows and %s extra-lesional. No spot was excluded from any panel.",
          fmt_n(n_out)),
  "Panels a-h come from ONE EAE brain section: n is spatial spots, not",
  "  biological replicates.",
  sprintf("Panel i: %s pMCAO single-cell microglia (parenchymal %s, vascular contact %s).",
          fmt_n(nrow(mc)),
          fmt_n(sum(mc$contact == "Parenchymal")),
          fmt_n(sum(mc$contact == "Vascular contact"))),
  "",
  "Panels a-c: label-transfer prediction scores. Within each cell type the two",
  sprintf("  maps share one colour scale (T cell 0-%.3f, microglia 0-%.3f, the 99th",
          cap_tcell, cap_mg),
  "  percentile of the pooled pair) so EAE_POS and EAE_NEG can be compared by eye.",
  "  Values above the cap are squished onto the top colour. Every spot is first",
  "  drawn in light grey to give the tissue silhouette, then positive spots are",
  "  drawn over it in ascending value order.",
  sprintf("  Panel c uses the same two scales as a and b and shows all %s spots of the",
          fmt_n(nrow(win_c))),
  sprintf("  perivascular window (%d x %d px), of which %s lie inside a lesion outline.",
          ncc, nrc, fmt_n(sum(win_c$in_lesion))),
  "  Its first tile is the unmodified stored low-resolution H&E of the same window,",
  "  shown as a histological reference; it is not used as a backdrop under the",
  "  score maps.",
  "",
  "Panel d: within-spot difference EAE_POS - EAE_NEG. Two-sided Wilcoxon",
  "  signed-rank test against zero, Benjamini-Hochberg across the four tests;",
  "  r_rb is the matched-pairs rank-biserial correlation (positive = EAE_POS",
  "  exceeds EAE_NEG). Box, median and interquartile range; white diamond, mean;",
  "  violin, full distribution at equal width scaling (width is not proportional",
  "  to n). The lesion-versus-extra-lesional line reports a two-sided Wilcoxon",
  "  rank-sum test on the same difference with its rank-sum r_rb.",
  sprintf("  T cell    extra-lesional  median %+.4f  r_rb %+.3f  %s",
          gp("T cell", "Extra-lesional", "median_diff"),
          gp("T cell", "Extra-lesional", "r_rb"),
          fmt_q(gp("T cell", "Extra-lesional", "q"))),
  sprintf("  T cell    lesion          median %+.4f  r_rb %+.3f  %s",
          gp("T cell", "Perivascular lesion", "median_diff"),
          gp("T cell", "Perivascular lesion", "r_rb"),
          fmt_q(gp("T cell", "Perivascular lesion", "q"))),
  sprintf("  Microglia extra-lesional  median %+.4f  r_rb %+.3f  %s",
          gp("Microglia", "Extra-lesional", "median_diff"),
          gp("Microglia", "Extra-lesional", "r_rb"),
          fmt_q(gp("Microglia", "Extra-lesional", "q"))),
  sprintf("  Microglia lesion          median %+.4f  r_rb %+.3f  %s",
          gp("Microglia", "Perivascular lesion", "median_diff"),
          gp("Microglia", "Perivascular lesion", "r_rb"),
          fmt_q(gp("Microglia", "Perivascular lesion", "q"))),
  sprintf("  Lesion shift: T cell r_rb %+.3f (%s); microglia r_rb %+.3f (%s).",
          gs("T cell", "r_rb"), fmt_q(gs("T cell", "q")),
          gs("Microglia", "r_rb"), fmt_q(gs("Microglia", "q"))),
  "  Both differences move DOWNWARD inside lesions, i.e. towards EAE_NEG.",
  "",
  sprintf("Panels e-f: %s microglial spots with an EAE_POS or EAE_NEG assignment",
          fmt_n(nrow(cx))),
  sprintf("  (EAE_POS %s, EAE_NEG %s). Cxcl10 is an integer UMI count, 0-%d, and is",
          fmt_n(n_pos), fmt_n(n_neg), cap_cx),
  sprintf("  exactly zero in %s of %s spots (%.1f%%).",
          fmt_n(sum(cx$Cxcl10 == 0)), fmt_n(nrow(cx)), 100 * mean(cx$Cxcl10 == 0)),
  sprintf("  Positivity: EAE_POS %d/%d = %.2f%% (Wilson 95%% CI %.2f-%.2f);",
          k_pos, n_pos, pct_pos,
          pos_tab$lo[pos_tab$group_label == "EAE_POS"],
          pos_tab$hi[pos_tab$group_label == "EAE_POS"]),
  sprintf("              EAE_NEG %d/%d = %.2f%% (Wilson 95%% CI %.2f-%.2f).",
          k_neg, n_neg, pct_neg,
          pos_tab$lo[pos_tab$group_label == "EAE_NEG"],
          pos_tab$hi[pos_tab$group_label == "EAE_NEG"]),
  sprintf("  Two-sided Fisher exact test: OR %.2f (95%% CI %.2f-%.2f), P = %.2e.",
          fF$or, fF$lo, fF$hi, fF$p),
  sprintf("  Wilcoxon rank-sum on the raw counts: P = %.2e, r_rb %+.3f.",
          wtF$p.value, rb_cx),
  sprintf("  Inside the perivascular window: EAE_POS %d/%d positive, EAE_NEG %d/%d.",
          les_pos_k, les_pos_n, les_neg_k, les_neg_n),
  "  With only 3 EAE_NEG spots in that window the contrast is uninformative and",
  "  is deliberately not drawn as a panel.",
  "",
  sprintf("Panels g-h: Vwf+ Cdh5+ large-vessel mask = %d of %s spots (%.2f%%).",
          n_vessel, fmt_n(n_spot), 100 * n_vessel / n_spot),
  "  Distance is Euclidean nearest-neighbour distance in low-resolution image",
  "  pixels. The archived image has no calibrated pixel size.",
  sprintf("  Perivascular (<= %d px): %s spots, Cxcl10+ = %.2f%%.",
          NEAR_PX, fmt_n(n_near), pct_near),
  sprintf("  Distal (> %d px): %s spots, Cxcl10+ = %.2f%%.",
          NEAR_PX, fmt_n(n_far), pct_far),
  sprintf("  Two-sided Fisher exact (perivascular vs distal): OR %.2f (95%% CI %.2f-%.2f), P = %.2e.",
          fG$or, fG$lo, fG$hi, fG$p),
  sprintf("  Spearman (distance vs Cxcl10 count, all spots): rho = %+.4f, P = %.2e.",
          as.numeric(sr$estimate), sr$p.value),
  sprintf("  Distance gradient: %d equal-count bins (n = %s to %s per bin).",
          NBIN, fmt_n(min(binned$n)), fmt_n(max(binned$n))),
  sprintf("  Cell-to-vessel distance: MgIV n = %s, median %.2f px; MgIVC n = %s, median %.2f px.",
          fmt_n(n_mgiv), med_iv, fmt_n(n_mgivc), med_ivc),
  sprintf("  Two-sided Wilcoxon rank-sum: P = %.2e, r_rb %+.3f (positive = MgIV further).",
          wtH$p.value, rb_H),
  "",
  sprintf("Panel i: independent pMCAO scRNA-seq microglia, %s cells.", fmt_n(nrow(mc))),
  sprintf("  Glycolytic genes used (%d): %s.", length(gly), paste(gly, collapse = ", ")),
  sprintf("  Glycolytic genes absent: %s.",
          ifelse(length(gly_missing), paste(gly_missing, collapse = ", "), "none")),
  sprintf("  OXPHOS genes used (%d): %s.", length(oxp), paste(oxp, collapse = ", ")),
  "  Score = mean z-score across the gene set. Tests are two-sided Wilcoxon",
  "  rank-sum with Benjamini-Hochberg across the three scores.",
  sprintf("  Glycolysis: mean vascular-contact %+.3f vs parenchymal %+.3f, Delta = %s, %s, r_rb %+.3f.",
          gi("glyc", "mean_vc"), gi("glyc", "mean_par"),
          signed_txt(gi("glyc", "delta")), fmt_q(gi("glyc", "q")), gi("glyc", "r_rb")),
  sprintf("  OXPHOS:     mean vascular-contact %+.3f vs parenchymal %+.3f, Delta = %s, %s, r_rb %+.3f.",
          gi("oxph", "mean_vc"), gi("oxph", "mean_par"),
          signed_txt(gi("oxph", "delta")), fmt_q(gi("oxph", "q")), gi("oxph", "r_rb")),
  sprintf("  Glyc-OXPHOS: mean vascular-contact %+.3f vs parenchymal %+.3f, Delta = %s, %s, r_rb %+.3f.",
          gi("gmo", "mean_vc"), gi("gmo", "mean_par"),
          signed_txt(gi("gmo", "delta")), fmt_q(gi("gmo", "q")), gi("gmo", "r_rb")),
  "",
  sprintf("Figure size: %.1f x %.1f mm", FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "NewFigure8_stats.txt"))

writeLines(c(
  paste0("Fig. 8 | Vascular-contact and parenchymal programmes partition in opposite ",
         "directions across one EAE brain section; Cxcl10 detection is restricted to ",
         "vascular-contact microglia and decays with distance from large vessels."),
  "",
  sprintf(paste0("a, Label-transfer prediction scores for the Tcell EAE_POS (left) and ",
                 "Tcell EAE_NEG (right) programmes across all %s spatial transcriptomic ",
                 "spots of one EAE brain section. Both maps share one colour scale spanning ",
                 "0 to %.2f, the 99th percentile of the pooled pair; values above the cap ",
                 "take the top colour. All spots are drawn in light grey to outline the ",
                 "tissue and spots with a positive score are drawn over that outline. ",
                 "Dashed boxes mark the three perivascular lesion windows enlarged in c."),
          fmt_n(n_spot), cap_tcell),
  sprintf(paste0("b, Mg EAE_POS and Mg EAE_NEG prediction scores for the same spots, ",
                 "displayed as in a on a shared scale spanning 0 to %.2f."), cap_mg),
  sprintf(paste0("c, The perivascular region boxed in a and b at higher magnification ",
                 "(%s spots, of which %s lie inside a lesion outline). Left, the unmodified ",
                 "low-resolution H&E histology of that window; right, the four scores of a ",
                 "and b for the same spots on the same two colour scales. Dashed boxes are ",
                 "the lesion outlines of a and b."),
          fmt_n(nrow(win_c)), fmt_n(sum(win_c$in_lesion))),
  sprintf(paste0("d, Within-spot difference between the EAE_POS and EAE_NEG prediction ",
                 "scores for each cell type, in extra-lesional tissue (n = %s spots) and ",
                 "inside the perivascular lesions (n = %s spots). The dashed line marks ",
                 "equal scores. Violins show the full distribution at equal width scaling, ",
                 "boxes the median and interquartile range, and white diamonds the mean. ",
                 "Upper labels, two-sided Wilcoxon signed-rank test of the difference ",
                 "against zero with Benjamini-Hochberg correction across the four tests, ",
                 "and the matched-pairs rank-biserial correlation r_rb (positive = EAE_POS ",
                 "above EAE_NEG). Lower labels, two-sided Wilcoxon rank-sum comparison of ",
                 "the same difference between lesion and extra-lesional tissue. Microglial ",
                 "spots favour EAE_POS and T cell spots favour EAE_NEG, and both differences ",
                 "move towards EAE_NEG inside lesions."),
          fmt_n(n_out), fmt_n(n_les)),
  sprintf(paste0("e, Spatial Cxcl10 counts in the %s EAE_POS (MgIV) and %s EAE_NEG (MgIVC) ",
                 "microglial spots of the same section, on a shared scale spanning 0 to %d ",
                 "counts. Grey marks every assigned spot and colour marks the spots with at ",
                 "least one count; the two maps carry different numbers of spots, so the ",
                 "sparser EAE_NEG map reflects fewer assigned spots as well as less signal. ",
                 "Dashed boxes are the lesion windows of a and b."),
          fmt_n(n_pos), fmt_n(n_neg), cap_cx),
  sprintf(paste0("f, Percentage of microglial spots with at least one Cxcl10 count. Bars ",
                 "show the observed percentage, error bars the Wilson 95%% score interval, ",
                 "and the counts below each bar the positive and total spots. Two-sided ",
                 "Fisher exact test on the 2 x 2 table of programme against Cxcl10 ",
                 "detection: odds ratio %.1f (95%% CI %.1f-%.1f), %s. Only 3 EAE_NEG ",
                 "microglial spots fall inside the perivascular windows, so the ",
                 "lesion-restricted contrast is reported in the statistics file rather than ",
                 "plotted."),
          fF$or, fF$lo, fF$hi, fmt_p(fF$p)),
  sprintf(paste0("g, Cxcl10-positive spots overlaid on the sparse Vwf+ Cdh5+ large-vessel ",
                 "mask (left; %d of %s spots). Middle, Cxcl10-positive fraction in ",
                 "perivascular tissue (distance to the nearest mask spot \u2264 %d px; n = %s) ",
                 "versus distal tissue (n = %s); bars show the observed percentage and error ",
                 "bars the Wilson 95%% score interval. Two-sided Fisher exact test: odds ",
                 "ratio %.2f (95%% CI %.2f-%.2f), %s. Right, Cxcl10-positive fraction in %d ",
                 "equal-count distance bins (n = %s spots per bin); points are the bin mean ",
                 "distance, horizontal bars the bin span, and vertical bars the Wilson 95%% ",
                 "interval. Two-sided Spearman correlation of distance against Cxcl10 count: ",
                 "\u03c1 = %s, %s. Distances are low-resolution image pixels; the archived ",
                 "image has no calibrated micrometre scale."),
          n_vessel, fmt_n(n_spot), NEAR_PX, fmt_n(n_near), fmt_n(n_far),
          fG$or, fG$lo, fG$hi, fmt_p(fG$p), NBIN, fmt_n(min(binned$n)),
          signed_txt(as.numeric(sr$estimate)), fmt_p(sr$p.value)),
  sprintf(paste0("h, Microglial spots on the same vessel mask (left) and the distance-to-vessel ",
                 "map of every spot (middle; colour capped at the 98th percentile, %.1f px). ",
                 "Right, cell-to-vessel distance for MgIV (vascular contact; n = %s; median ",
                 "%.2f px) versus MgIVC (parenchymal; n = %s; median %.2f px). Violins show ",
                 "the full distribution, boxes the median and interquartile range, and white ",
                 "diamonds the mean; the y axis is square-root scaled so the informative ",
                 "0-20 px range stays readable. Two-sided Wilcoxon rank-sum test, %s, rank-biserial ",
                 "r_rb = %+.2f (negative = MgIV closer)."),
          d_cap, fmt_n(n_mgiv), med_iv, fmt_n(n_mgivc), med_ivc,
          fmt_p(wtH$p.value), rb_H),
  sprintf(paste0("i, Independent pMCAO single-cell microglia (n = %s cells). Glycolysis, ",
                 "OXPHOS and glycolysis-minus-OXPHOS scores (mean z-score across %d and %d ",
                 "genes respectively) in parenchymal (n = %s) versus vascular-contact ",
                 "(n = %s) cells. Violins, boxes and white diamonds as in d. Upper labels, ",
                 "two-sided Wilcoxon rank-sum test with Benjamini-Hochberg correction across ",
                 "the three scores, and the mean difference \u0394 (vascular contact minus ",
                 "parenchymal). Right, mean glycolysis score per cluster \u00d7 contact; tile ",
                 "values are the cell-group means. This panel is a different injury model and ",
                 "a different assay from a-h and is corroboration, not a replicate."),
          fmt_n(nrow(mc)), length(gly), length(oxp),
          fmt_n(sum(mc$contact == "Parenchymal")),
          fmt_n(sum(mc$contact == "Vascular contact"))),
  "",
  paste0("Panels a-h derive from a single EAE brain section, so n refers to spatial ",
         "spots rather than biological replicates. Source data are provided as a ",
         "Source Data file."),
  "",
  paste0("Abbreviations: EAE_POS / MgIV, vascular-contact (intravascular-stained) ",
         "programme; EAE_NEG / MgIVC, parenchyma-associated programme; r_rb, ",
         "rank-biserial correlation.")
), file.path(OUTDIR, "NewFigure8_legend.txt"))

writeLines(c(
  "NewFig. 8 | provenance",
  "",
  "New folder. No earlier figure, script or source file is overwritten.",
  "",
  "Backend: R only (ggplot2 + patchwork + scattermore; svglite / cairo_pdf / ragg).",
  "Reuse level: structural adaptation. Every panel is rebuilt from the stored",
  "  spot-level source tables of the two parent figures; nothing is cropped out",
  "  of the parent SVG files.",
  "",
  "Panel mapping from the parent figures",
  "  first parent (spatial-prediction composite), panels a-f  -> a-f",
  "  second parent (vascular-contact composite), panel b      -> g",
  "  second parent (vascular-contact composite), panel c      -> h",
  "  second parent (vascular-contact composite), panel d      -> i",
  "",
  "Both spatial parents describe the same EAE section in the same low-resolution",
  "  image frame; the stored H&E arrays were compared and are identical, and the",
  "  32,811 spot barcodes and x/y coordinates match between the two exports.",
  "  EAE_POS microglia == MgIV and EAE_NEG microglia == MgIVC, so one colour",
  "  vocabulary is used for vascular contact versus parenchyma throughout.",
  "",
  "Changes relative to the parent panels",
  "  a-f  inherited from the first parent with the same shared colour caps,",
  "       paired within-spot difference, Wilson/Fisher positivity, and unmodified",
  "       H&E reference tile. Programme colours are remapped onto the",
  "       vascular-contact (red) / parenchyma (blue) pair used in g-i so the",
  "       nine-panel page has one identity vocabulary.",
  "  g    the parent used equal-width distance bins (18,726 spots in the first",
  "       bin and 31 in the last) and a one-sided Fisher test with binomial SE.",
  "       This panel uses eight equal-count bins, Wilson intervals on every",
  "       proportion, and a two-sided Fisher test with a 95% CI, matching panel f.",
  "  h    the parent used a one-sided Mann-Whitney test. The direction is the",
  "       result being claimed, so the test is two-sided and carries a rank-biserial",
  "       effect size. The distance map uses a truncated viridis ramp so it cannot",
  "       be confused with the warm prediction-score ramp of a-c and e.",
  "  i    the parent used asterisk codes and a one-line P summary. This panel",
  "       prints exact q values, mean differences and a BH correction across the",
  "       three scores. Cluster tiles keep in-cell values so the heat map survives",
  "       greyscale print.",
  "",
  "Image integrity: the H&E tile in c is the stored low-resolution image of the",
  "  section, cropped to the lesion bounding box only. No contrast, gamma, channel",
  "  or local adjustment was applied. The stored image carries no calibrated pixel",
  "  size, so distances stay in pixels and no micrometre scale bar is drawn."
), file.path(OUTDIR, "NewFigure8_provenance.txt"))

writeLines(c(
  "NewFig. 8 | QA record",
  "",
  sprintf("Backend            R only (ggplot2 %s); drawing, export and QA crops all in R",
          as.character(utils::packageVersion("ggplot2"))),
  sprintf("Final size         %.1f x %.1f mm (double column; Nature full-page limit 247 mm)",
          FIG_W_MM, FIG_H_MM),
  "Exports            NewFigure8.pdf (cairo_pdf), NewFigure8.svg (svglite), NewFigure8.png (600 dpi)",
  "Editable text      yes; live text in both SVG and PDF",
  "Glyph floor        5.2 pt minimum (axis, legend, annotation); panel tags 8 pt.",
  "                   cairo_pdf emits '1 Tf' and scales through the text matrix, so",
  "                   Tf-based PDF scanners report a false failure; sizes verified in SVG.",
  "Vector/raster      a, b, g, h rasterized point marks (scattermore), all spots retained;",
  "                   c, e vector points; the first tile of c is an unmodified H&E",
  "                   raster; d, f, i violins/bars/tiles are vector",
  "Colour             one sequential warm ramp over a light-grey tissue silhouette",
  "                   (a, b, c, e); one red/blue identity pair for vascular contact",
  "                   versus parenchyma (f, h, i); orange Cxcl10 accent (g); truncated",
  "                   viridis distance ramp (h). No rainbow or jet.",
  "Shared scales      EAE_POS and EAE_NEG maps of one cell type share a colour cap,",
  "                   and c reuses the caps of a and b. Proportions in f and g share",
  "                   Wilson intervals and two-sided Fisher tests.",
  "",
  "Panel audit (unique claim | summary | spread | replicate unit | risks)",
  sprintf("a  Which T cell programme occupies the section? | %s spots | none (census) |",
          fmt_n(n_spot)),
  "   one section | EAE_POS is genuinely faint on the shared scale; that is the result",
  sprintf("b  Which microglial programme occupies the section? | %s spots | none |",
          fmt_n(n_spot)),
  "   one section | shared cap as in a",
  sprintf("c  What does the lesion window look like? | %s spots | none | one section |",
          fmt_n(nrow(win_c))),
  "   the H&E reference tile is unmodified and therefore looks strongly saturated;",
  "   no calibrated scale bar exists for the stored image, so none is drawn",
  "d  Which programme wins per spot and does the lesion change it? | median + IQR + mean |",
  "   signed-rank q with matched-pairs r_rb, plus a rank-sum lesion contrast | spot |",
  sprintf("   r_rb  T cell %+.3f / %+.3f, microglia %+.3f / %+.3f (extra-lesional / lesion)",
          gp("T cell", "Extra-lesional", "r_rb"), gp("T cell", "Perivascular lesion", "r_rb"),
          gp("Microglia", "Extra-lesional", "r_rb"), gp("Microglia", "Perivascular lesion", "r_rb")),
  sprintf("e  Where is Cxcl10 detected? | %d of %d microglial spots positive | none |",
          sum(cx$Cxcl10 > 0), nrow(cx)),
  "   one section | very sparse; grey = every assigned spot, colour = at least one count",
  sprintf("f  Is detection restricted to EAE_POS? | %.2f%% vs %.2f%% | Wilson 95%% CI |",
          pct_pos, pct_neg),
  sprintf("   spot | Fisher OR %.2f (%.2f-%.2f); group sizes unbalanced (%s vs %s)",
          fF$or, fF$lo, fF$hi, fmt_n(n_pos), fmt_n(n_neg)),
  sprintf("g  Does Cxcl10 sit at large vessels? | %.2f%% vs %.2f%% near/far | Wilson + Fisher |",
          pct_near, pct_far),
  sprintf("   spot | mask is only %d spots; equal-count bins keep per-bin n comparable",
          n_vessel),
  sprintf("h  Do vascular-contact microglia sit closer to vessels? | median %.2f vs %.2f px |",
          med_iv, med_ivc),
  sprintf("   rank-sum r_rb %+.3f | spot | square-root axis; sparse mask as in g", rb_H),
  sprintf("i  Is the same compartment glycolytic in pMCAO? | Delta glyc %s | rank-sum q + r_rb |",
          signed_txt(gi("glyc", "delta"))),
  "   cell | different injury model and assay from a-h; corroboration, not a replicate",
  "",
  "OUTSTANDING BEFORE SUBMISSION",
  "- DIRECTION: the EAE_POS minus EAE_NEG difference moves DOWNWARD inside lesions",
  "  for both cell types. Do not write that the vascular-contact programme is",
  "  enriched in perivascular lesions relative to the parenchymal one.",
  "- T cell EAE_NEG exceeds EAE_POS at almost every spot. Label-transfer scores are",
  "  compositional across the reference label set, so this partly reflects reference",
  "  composition. Frame the T cell result as a balance, not as an absolute abundance.",
  sprintf("- The Cxcl10 result is the clean one: a %.1f-fold odds ratio for detection in",
          fF$or),
  sprintf("  EAE_POS microglial spots, and a %.2f-fold odds ratio for perivascular versus",
          fG$or),
  "  distal tissue. It is still one section.",
  "- n in a-h is spatial spots from a single EAE section, not biological replicates.",
  "  State this explicitly; it is the main reviewer risk in this figure.",
  "- q values near zero are driven by 32,811 spots or 33,867 cells. Always quote the",
  "  effect size alongside them.",
  "- Panel i is a different injury model. Do not call it a replicate of a-h.",
  "- Nine panels is at the upper limit of one Nature display item. a-f and g-i can",
  "  be split if the editor asks."
), file.path(OUTDIR, "NewFigure8_QA.txt"))

writeLines(c(
  "NewFig. 8 | Nature-format multipanel figure",
  "==========================================",
  "",
  "New folder. Nothing outside this folder was written to, and no earlier figure,",
  "script or source file was overwritten.",
  "",
  "Deliverables",
  "------------",
  "  NewFigure8.pdf              vector, cairo_pdf, editable Arial text",
  "  NewFigure8.svg              vector, svglite, editable Arial text",
  "  NewFigure8.png              600 dpi raster preview",
  "  NewFigure8_legend.txt       submission-ready figure legend",
  "  NewFigure8_stats.txt        statistics, sample sizes and integrity notes",
  "  NewFigure8_QA.txt           QA record, panel audit and pre-submission warnings",
  "  NewFigure8_provenance.txt   panel mapping and every change from the parent panels",
  "  NewFigure8_contract.txt     the figure contract written before any plotting code",
  "  source_data/                per-panel CSV source data and test tables",
  "  scripts/                    the R script that produces everything here",
  "  qa/                         per-row crops used for visual QA",
  "",
  sprintf("Final size %.1f x %.1f mm (double column). Backend: R only.",
          FIG_W_MM, FIG_H_MM),
  "",
  "Panels",
  "------",
  "  a  Tcell EAE_POS / EAE_NEG prediction score, whole section",
  "  b  Mg EAE_POS / EAE_NEG prediction score, whole section",
  "  c  H&E reference plus the four scores inside the perivascular lesion window",
  "  d  within-spot EAE_POS minus EAE_NEG difference, lesion vs extra-lesional",
  "  e  spatial Cxcl10 in EAE_POS and EAE_NEG microglial spots",
  "  f  proportion of Cxcl10-positive microglial spots with Wilson intervals",
  "  g  Cxcl10+ spots on the large-vessel mask; perivascular enrichment; distance gradient",
  "  h  microglia on the vessel mask; distance-to-vessel map; MgIV vs MgIVC distance",
  "  i  glycolysis / OXPHOS / balance by vascular-contact status (pMCAO scRNA-seq)",
  "",
  "Reproducing",
  "-----------",
  "  Rscript scripts/make_newfigure8.R",
  "",
  "  Requires R with ggplot2, dplyr, tidyr, patchwork, scales, svglite, ragg,",
  "  systemfonts, grid, scattermore and RANN.",
  "",
  "Read NewFigure8_QA.txt before writing the manuscript text."
), file.path(OUTDIR, "README.txt"))

cat("DONE\n")
