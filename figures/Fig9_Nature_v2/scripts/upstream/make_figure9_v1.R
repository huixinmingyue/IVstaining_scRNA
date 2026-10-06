#!/usr/bin/env Rscript
# =============================================================================
# Fig. 9 | L1-L3 and L4-L6 windows, drawn separately
#          Cxcl10 stays in the vascular-contact compartment
#
# Core conclusion
#   In one EAE brain section the vascular-contact (EAE_POS / MgIV) and
#   parenchymal (EAE_NEG / MgIVC) programmes partition in opposite directions;
#   vascular-contact microglia sit closest to the Vwf+Cdh5+ large-vessel mask;
#   Cxcl10 detection is confined to that compartment and decays with distance
#   from large vessels. L4-L6 is a second window set, never pooled with L1-L3.
#   The independent pMCAO panel is not part of this figure.
#
# Archetype  asymmetric mixed-modality figure
# Backend    R only (ggplot2 + patchwork + scattermore); no cross-language render
# Size       183 mm wide, <= 247 mm tall; 5 pt rendered-glyph floor
#
# Panel map
#   a  Tcell EAE_POS / EAE_NEG prediction score over the whole section
#   b  Mg    EAE_POS / EAE_NEG prediction score over the whole section
#   c  H&E plus the four scores inside L1-L3
#   d  the same five-tile layout inside L4-L6
#   e  within-spot EAE_POS minus EAE_NEG, Outside vs L1-L3 vs L4-L6
#   f  spatial Cxcl10 in EAE_POS and EAE_NEG microglial spots
#   g  Cxcl10-positive fraction, whole section and inside L4-L6
#   h  Cxcl10+ spots on the large-vessel mask, enrichment and distance gradient
#   i  microglia on the vessel mask and cell-to-vessel distance,
#      whole section and inside L4-L6
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

OUTDIR   <- file.path(IV_FIG_ROOT, "Fig9_Nature_v1")
SRC_PRED <- file.path(OUTDIR, "source_data")
SRC_CX   <- SRC_PRED
SRC_TAB  <- SRC_PRED
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
H_ROW <- c(ab = 36.5, c = NA, d = NA, e = 34, fg = 36, h = 40, i = 40)
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
scope_cols  <- c(Outside = "#9E9E9E", `L1-L3` = ROI_COL, `L4-L6` = "#AF7AC5")

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
rects_add <- read.csv(file.path(SRC_TAB, "lesion_rects_expanded.csv"),
                      stringsAsFactors = FALSE)
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

gj$in_l13 <- gj$cell %in% les$cell
cx$cell <- as.character(cx$cell)
cx$in_l13 <- cx$cell %in% les$cell
in_box <- function(x, y, r) {
  x >= r$xmin & x <= r$xmax & y >= r$ymin & y <= r$ymax
}
rects_l46 <- rects_add[rects_add$set == "added", c("lesion", "xmin", "xmax", "ymin", "ymax")]
gj$in_l46 <- FALSE
for (i in seq_len(nrow(rects_l46))) {
  gj$in_l46 <- gj$in_l46 | in_box(gj$x_img, gj$y_img, rects_l46[i, ])
}
stopifnot(!any(gj$in_l13 & gj$in_l46))
cx$in_l46 <- cx$cell %in% gj$cell[gj$in_l46]
gj$win <- factor(ifelse(gj$in_l13, "L1-L3", ifelse(gj$in_l46, "L4-L6", "Outside")),
                 levels = c("Outside", "L1-L3", "L4-L6"))
cx$win <- gj$win[match(cx$cell, gj$cell)]
sp$win <- gj$win[match(sp$cell_id, gj$cell)]
rects_all <- rbind(
  rects[, c("lesion", "xmin", "xmax", "ymin", "ymax")],
  rects_l46
)
n_spot <- nrow(gj)
n_l13 <- sum(gj$in_l13)
n_l46 <- sum(gj$in_l46)
n_out <- sum(gj$win == "Outside")
cat(sprintf("  spots=%d  L1-L3=%d  L4-L6=%d  outside=%d\n",
            n_spot, n_l13, n_l46, n_out))
stopifnot(n_l13 + n_l46 + n_out == n_spot)
cx$group_label <- factor(cx$group_label,
                         levels = c("Mg_EAE_POS", "Mg_EAE_NEG"),
                         labels = c("EAE_POS", "EAE_NEG"))

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
cat_map <- function(bg, layers, keys, cols, title, bar_px = 100, rects_df = NULL) {
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
  p <- p +
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
      legend.position = "inside",
      legend.position.inside = c(0.0, 0.0),
      legend.justification.inside = c(0, 0),
      legend.key.spacing.y = unit(0, "mm"),
      legend.key.size = unit(1.6, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0.8, 0.8, 0.8, 0.8)
    )
  if (!is.null(rects_df) && nrow(rects_df) > 0) {
    p <- p + geom_rect(
      data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
      inherit.aes = FALSE, fill = NA, colour = ROI_COL,
      linewidth = 0.3, linetype = "dashed")
  }
  p
}

# =============================================================================
# a-b  whole-section prediction score maps
# =============================================================================
cat("Panels a-b ...\n")
mA1 <- spat_score(gj, "Tcell EAE_POS", "Tcell EAE_POS", cap_tcell, rects_df = rects_all)
mA2 <- spat_score(gj, "Tcell EAE_NEG", "Tcell EAE_NEG", cap_tcell,
                  show_legend = TRUE, rects_df = rects_all)
mB1 <- spat_score(gj, "Mg EAE_POS", "Mg EAE_POS", cap_mg, rects_df = rects_all)
mB2 <- spat_score(gj, "Mg EAE_NEG", "Mg EAE_NEG", cap_mg,
                  show_legend = TRUE, rects_df = rects_all)

pA <- wrap_elements(full = (mA1 | mA2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "a") + tag_theme
pB <- wrap_elements(full = (mB1 | mB2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "b") + tag_theme

# =============================================================================
# c  L1-L3 zoom    d  L4-L6 zoom, same five-tile layout
# =============================================================================
cat("Panels c-d ...\n")
nr_all <- dim(he)[1]; nc_all <- dim(he)[2]
pad <- 15
c_pad <- (FIG_W_MM - 5 * C_W_MM) / 2

make_zoom <- function(anchor, rects_sub, he_title, which) {
  xmin <- max(0, min(anchor$x_img) - pad)
  xmax <- min(nc_all, max(anchor$x_img) + pad)
  ymin <- max(0, min(anchor$y_img) - pad)
  ymax <- min(nr_all, max(anchor$y_img) + pad)
  r1 <- max(1, floor(ymin) + 1); r2 <- min(nr_all, ceiling(ymax))
  c1 <- max(1, floor(xmin) + 1); c2 <- min(nc_all, ceiling(xmax))
  he_crop <- he[r1:r2, c1:c2, , drop = FALSE]
  nrc <- dim(he_crop)[1]; ncc <- dim(he_crop)[2]
  win <- gj %>%
    filter(x_img >= c1 - 1, x_img <= c2, y_img >= r1 - 1, y_img <= r2) %>%
    mutate(x_img = x_img - (c1 - 1), y_img = y_img - (r1 - 1))
  rects_z <- transform(rects_sub,
                       xmin = xmin - (c1 - 1), xmax = xmax - (c1 - 1),
                       ymin = ymin - (r1 - 1), ymax = ymax - (r1 - 1))
  cat(sprintf("  %s window %d x %d px; %d spots, of which %d in this set\n",
              which, ncc, nrc, nrow(win),
              if (which == "L1-L3") sum(win$in_l13) else sum(win$in_l46)))
  xlim <- c(0, ncc); ylim <- c(-nrc, 0)
  mk <- function(col, cap) {
    spat_score(win, col, col, cap, rects_df = rects_z, raster = FALSE,
               pt = 0.42, xlim = xlim, ylim = ylim)
  }
  tiles <- list(
    he_tile(he_crop, nrc, ncc, he_title, rects_z),
    mk("Tcell EAE_POS", cap_tcell),
    mk("Tcell EAE_NEG", cap_tcell),
    mk("Mg EAE_POS", cap_mg),
    mk("Mg EAE_NEG", cap_mg)
  )
  row <- wrap_elements(panel = (
    plot_spacer() | tiles[[1]] | tiles[[2]] | tiles[[3]] | tiles[[4]] | tiles[[5]] | plot_spacer()
  ) + plot_layout(widths = c(c_pad, rep(C_W_MM, 5), c_pad)))
  list(row = row, h = C_W_MM / (ncc / nrc) + C_TITLE_MM,
       n = nrow(win), n_set = if (which == "L1-L3") sum(win$in_l13) else sum(win$in_l46),
       ncc = ncc, nrc = nrc)
}

zoom_c <- make_zoom(les, rects, "H&E, L1\u2013L3", "L1-L3")
zoom_d <- make_zoom(gj[gj$in_l46, ], rects_l46, "H&E, L4\u2013L6", "L4-L6")
H_ROW[["c"]] <- zoom_c$h
H_ROW[["d"]] <- zoom_d$h
FIG_H_MM <- round(sum(H_ROW), 1)
cat(sprintf("  page %.1f x %.1f mm\n", FIG_W_MM, FIG_H_MM))
stopifnot(FIG_H_MM <= 247)
pC <- zoom_c$row + labs(tag = "c") + tag_theme
pDwin <- zoom_d$row + labs(tag = "d") + tag_theme
n_win_c <- zoom_c$n; n_set_c <- zoom_c$n_set
n_win_d <- zoom_d$n; n_set_d <- zoom_d$n_set

# =============================================================================
# e  within-spot EAE_POS - EAE_NEG difference
#    Outside, L1-L3 and L4-L6 stay separate.
# =============================================================================
cat("Panel e ...\n")
gj$scope <- gj$win
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

shift_one <- function(d, a) {
  x <- d$diff[d$scope == a]
  y <- d$diff[d$scope == "Outside"]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.frame(
    contrast = a, n = length(x),
    median_window = median(x), median_outside = median(y),
    p = wt$p.value,
    r_rb = rb_rs(c(x, y), c(rep(a, length(x)), rep("Outside", length(y))), a),
    stringsAsFactors = FALSE
  )
}
shift_tab <- bind_rows(lapply(c("T cell", "Microglia"), function(ct) {
  d <- diff_df[diff_df$celltype == ct, ]
  bind_rows(lapply(c("L1-L3", "L4-L6"), function(a) {
    cbind(celltype = ct, shift_one(d, a))
  }))
}))
shift_tab$q <- p.adjust(shift_tab$p, method = "BH")
cat("  each window set vs outside:\n")
print(as.data.frame(shift_tab), digits = 4)

mean_d <- diff_df %>% group_by(celltype, scope) %>%
  summarise(m = mean(diff), .groups = "drop")

env_d <- diff_df %>% group_by(celltype) %>%
  summarise(lo = min(diff), hi = max(diff), .groups = "drop")
ann_d <- paired_tab %>% left_join(env_d, by = "celltype") %>%
  mutate(label = sprintf("%s\nr_rb = %+.2f", vapply(q, fmt_q, character(1)), r_rb),
         y = hi + 0.08 * (hi - lo))
ann_shift <- shift_tab %>% left_join(env_d, by = "celltype") %>%
  mutate(ord = ifelse(contrast == "L1-L3", 1, 2),
         label = sprintf("%s vs outside: %s, r_rb = %+.2f",
                         contrast, vapply(q, fmt_q, character(1)), r_rb),
         y = lo - (0.02 + 0.10 * ord) * (hi - lo))

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
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = scope_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
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
                  cap_cx, rects_df = rects_all, raster = FALSE, pt = 0.45,
                  xlim = e_xlim, ylim = e_ylim, bar_title = "Cxcl10")
mE2 <- spat_score(cx %>% filter(group_label == "EAE_NEG"), "Cxcl10",
                  sprintf("EAE_NEG (MgIVC)\nn = %s spots",
                          format(sum(cx$group_label == "EAE_NEG"), big.mark = ",")),
                  cap_cx, show_legend = TRUE, rects_df = rects_all,
                  raster = FALSE, pt = 0.45, xlim = e_xlim, ylim = e_ylim,
                  bar_title = "Cxcl10")
pE <- wrap_elements(full = (mE1 | mE2) +
                      plot_layout(widths = c(MAP_MM, MAP_MM + LEG_MM))) +
  labs(tag = "f") + tag_theme

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

# L1-L3 holds only 3 EAE_NEG microglial spots, so that contrast is recorded
# and not drawn. L4-L6 has enough spots of both programmes to plot.
les_cx <- cx %>% filter(in_l13) %>% group_by(group_label) %>%
  summarise(n = n(), k = sum(Cxcl10 > 0), .groups = "drop")
cat("  L1-L3 microglia (reported, not plotted):\n")
print(as.data.frame(les_cx))
pos_l46 <- cx %>% filter(in_l46) %>% group_by(group_label) %>%
  summarise(n = n(), k = sum(Cxcl10 > 0), .groups = "drop")
pos_l46 <- bind_cols(pos_l46, wilson_df(pos_l46$k, pos_l46$n))
f46 <- fisher2(pos_l46$k[pos_l46$group_label == "EAE_POS"],
               pos_l46$n[pos_l46$group_label == "EAE_POS"],
               pos_l46$k[pos_l46$group_label == "EAE_NEG"],
               pos_l46$n[pos_l46$group_label == "EAE_NEG"],
               "EAE_POS", "EAE_NEG")
cat(sprintf("  L4-L6 Fisher OR=%.2f (%.2f-%.2f) P=%.3g\n",
            f46$or, f46$lo, f46$hi, f46$p))
pos_plot <- bind_rows(
  pos_tab %>% mutate(scope = "All spots"),
  pos_l46 %>% mutate(scope = "L4-L6")
)
pos_plot$scope <- factor(pos_plot$scope, levels = c("All spots", "L4-L6"))
pos_plot$group_label <- factor(pos_plot$group_label, levels = c("EAE_POS", "EAE_NEG"))
f_top <- max(pos_plot$hi) * 1.85
pF <- ggplot(pos_plot, aes(scope, pct, fill = group_label)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64,
           colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.72),
                width = 0.16, linewidth = 0.3) +
  geom_text(aes(y = hi, label = sprintf("%.1f%%", pct), group = group_label),
            position = position_dodge(width = 0.72),
            vjust = -0.4, size = LABEL_MM, family = FAM) +
  annotate("text", x = 1.5, y = f_top,
           label = sprintf("All spots  OR = %.1f (%.1f\u2013%.1f), %s\nL4\u2013L6  OR = %.1f (%.1f\u2013%.1f), %s",
                           fF$or, fF$lo, fF$hi, fmt_p(fF$p),
                           f46$or, f46$lo, f46$hi, fmt_p(f46$p)),
           size = LABEL_MM, family = FAM, lineheight = 1.05, vjust = 1) +
  scale_fill_manual(values = prog_cols, name = NULL,
                    labels = c(EAE_POS = "EAE_POS (MgIV)", EAE_NEG = "EAE_NEG (MgIVC)")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.02)),
                     limits = c(0, f_top)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Cxcl10-positive\nmicroglial spots (%)",
       title = "Cxcl10 detection") +
  theme(
    axis.text.x = element_text(size = 5.6, colour = "black"),
    axis.title.y = element_text(size = 6.5, lineheight = 1.05),
    plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 5.2),
    legend.key.size = unit(2.2, "mm"),
    plot.margin = margin(1.2, 1.5, 1.2, 1.5)
  ) +
  labs(tag = "g") + tag_theme

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
  title = "Cxcl10+ on vessel mask",
  rects_df = rects_all
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
  plot_layout(widths = c(38, 42, 103))) + labs(tag = "h") + tag_theme

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
  title = "Microglia on vessel mask",
  rects_df = rects_all
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
        plot.margin = margin(0.8, 0.8, 0.8, 0.8)) +
  geom_rect(data = rects_all,
            aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
            inherit.aes = FALSE, fill = NA, colour = ROI_COL,
            linewidth = 0.3, linetype = "dashed")

# h3  cell-to-vessel distance. Two-sided rank-sum with a rank-biserial effect
# size; the parent used a one-sided test, but the direction is the result being
# claimed and should not be assumed by the test.
dist_test <- function(d) {
  wt <- wilcox.test(dist ~ grp, data = d, exact = FALSE)
  data.frame(n_iv = sum(d$grp == "MgIV"), n_ivc = sum(d$grp == "MgIVC"),
             p = wt$p.value,
             r_rb = rb_rs(d$dist, as.character(d$grp), "MgIV"))
}
mg_h <- bind_rows(
  mg %>% mutate(scope = "All spots"),
  mg %>% filter(win == "L4-L6") %>% mutate(scope = "L4-L6")
)
mg_h$scope <- factor(mg_h$scope, levels = c("All spots", "L4-L6"))
dist_tab <- bind_rows(lapply(levels(mg_h$scope), function(s) {
  cbind(scope = s, dist_test(mg_h[mg_h$scope == s, ]))
}))
dist_tab$q <- p.adjust(dist_tab$p, method = "BH")
wtH <- list(p.value = dist_tab$p[dist_tab$scope == "All spots"])
rb_H <- dist_tab$r_rb[dist_tab$scope == "All spots"]
med_H <- mg_h %>% group_by(scope, grp) %>%
  summarise(n = n(), median = median(dist), mean = mean(dist), .groups = "drop")
cat("  cell-to-vessel distance:\n")
print(as.data.frame(med_H), digits = 4)
print(as.data.frame(dist_tab), digits = 4)
h3_max <- max(mg_h$dist)
h3_brk <- c(0, 5, 10, 20, 40, 80, 160, 320)
h3_brk <- h3_brk[h3_brk <= h3_max * 1.02]
ann_h <- dist_tab
ann_h$label <- sprintf("%s\nr_rb = %+.2f", vapply(dist_tab$p, fmt_p, ""), dist_tab$r_rb)
pH3 <- ggplot(mg_h, aes(grp, dist, fill = grp)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.22,
              colour = "black", alpha = 0.85) +
  geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
               fill = "white", colour = "black") +
  geom_point(data = med_H, aes(grp, mean), inherit.aes = FALSE, shape = 23,
             size = 0.9, stroke = 0.25, fill = "white", colour = "black") +
  geom_text(data = ann_h, aes(x = 1.5, y = h3_max * 1.35, label = label),
            inherit.aes = FALSE, size = LABEL_MM, family = FAM,
            lineheight = 1.05, vjust = 1) +
  facet_wrap(~ scope, nrow = 1) +
  scale_fill_manual(values = c(MgIV = C_CONTACT, MgIVC = C_PAREN),
                    guide = "none") +
  scale_x_discrete(labels = c(MgIV = "MgIV", MgIVC = "MgIVC")) +
  scale_y_sqrt(breaks = h3_brk, limits = c(0, h3_max * 1.7),
               expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Distance to vessel\n(px, square-root scale)",
       title = "Cell-to-vessel distance") +
  theme(axis.text.x = element_text(size = 5.2, colour = "black"),
        axis.title.y = element_text(size = 6.2, colour = "black", lineheight = 1.05),
        plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 5.6, face = "bold"),
        plot.margin = margin(1.2, 1.2, 1.2, 1.2))

pHrow <- wrap_elements(panel = (pH1 | pH2 | pH3) +
  plot_layout(widths = c(38, 44, 101))) + labs(tag = "i") + tag_theme


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
d_pad <- max((FIG_W_MM - 168) / 2, 0.1)
row_e <- wrap_elements(panel = (
  plot_spacer() | pD | plot_spacer()
) + plot_layout(widths = c(d_pad, 168, d_pad))) +
  labs(tag = "e") + tag_theme

row_fg <- (pE | pF | plot_spacer()) +
  plot_layout(widths = c(ab_cell, 78, FIG_W_MM - ab_cell - 78))

fig <- row_ab / pC / pDwin / row_e / row_fg / pGrow / pHrow +
  plot_layout(heights = as.numeric(H_ROW))

# =============================================================================
# Source-data tables, export, legend
# =============================================================================
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
cat("Writing source data ...\n")
write.csv(gj %>% select(cell, x_img, y_img, win, in_l13, in_l46, all_of(score_cols)),
          file.path(OUT_SD, "Fig9_spatial_prediction_scores.csv"), row.names = FALSE)
write.csv(as.data.frame(paired_tab),
          file.path(OUT_SD, "Fig9e_paired_signed_rank.csv"), row.names = FALSE)
write.csv(as.data.frame(shift_tab),
          file.path(OUT_SD, "Fig9e_window_vs_outside.csv"), row.names = FALSE)
write.csv(cx %>% select(cell, x_img, y_img, cluster, group_label, win, Cxcl10),
          file.path(OUT_SD, "Fig9f_spatial_cxcl10.csv"), row.names = FALSE)
write.csv(as.data.frame(pos_tab) %>%
            mutate(fisher_or = fF$or, fisher_ci_lo = fF$lo,
                   fisher_ci_hi = fF$hi, fisher_p = fF$p),
          file.path(OUT_SD, "Fig9g_cxcl10_positivity_all.csv"), row.names = FALSE)
write.csv(as.data.frame(pos_l46) %>%
            mutate(fisher_or = f46$or, fisher_ci_lo = f46$lo,
                   fisher_ci_hi = f46$hi, fisher_p = f46$p),
          file.path(OUT_SD, "Fig9g_cxcl10_positivity_L4L6.csv"), row.names = FALSE)
write.csv(as.data.frame(les_cx),
          file.path(OUT_SD, "Fig9g_cxcl10_positivity_L1L3_not_plotted.csv"),
          row.names = FALSE)
write.csv(sp %>% select(cell = cell_id, x, y, cluster, win, Vwf, Cdh5,
                        is_vessel, dist, Cxcl10, cx_pos, near),
          file.path(OUT_SD, "Fig9hi_vessel_distance_spots.csv"), row.names = FALSE)
write.csv(as.data.frame(near_tab) %>%
            mutate(fisher_or = fG$or, fisher_ci_lo = fG$lo,
                   fisher_ci_hi = fG$hi, fisher_p = fG$p),
          file.path(OUT_SD, "Fig9h_perivascular_enrichment.csv"), row.names = FALSE)
write.csv(as.data.frame(binned) %>%
            mutate(spearman_rho = as.numeric(sr$estimate),
                   spearman_p = sr$p.value),
          file.path(OUT_SD, "Fig9h_distance_gradient_bins.csv"), row.names = FALSE)
write.csv(as.data.frame(dist_tab),
          file.path(OUT_SD, "Fig9i_cell_to_vessel_distance.csv"), row.names = FALSE)
write.csv(as.data.frame(med_H),
          file.path(OUT_SD, "Fig9i_distance_medians.csv"), row.names = FALSE)
write.csv(rects_all, file.path(OUT_SD, "Fig9_window_rects.csv"), row.names = FALSE)

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base_path <- file.path(OUTDIR, "Figure9")
cat("Writing SVG ...\n")
svglite::svglite(paste0(base_path, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); dev.off()
cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base_path, ".pdf"), width = w_in, height = h_in, family = FAM)
print(fig); dev.off()
cat("Writing PNG ...\n")
ragg::agg_png(paste0(base_path, ".png"), width = w_in, height = h_in,
              units = "in", res = 600, background = "white")
print(fig); dev.off()
cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(base_path, ".tiff"), width = w_in, height = h_in,
               units = "in", res = 600, compression = "lzw", background = "white")
print(fig); dev.off()

cat("Writing QA row renders ...\n")
qa_rows <- list(
  row_ab = list(p = row_ab, h = H_ROW[["ab"]]),
  row_c  = list(p = pC,     h = H_ROW[["c"]]),
  row_d  = list(p = pDwin,  h = H_ROW[["d"]]),
  row_e  = list(p = row_e,  h = H_ROW[["e"]]),
  row_fg = list(p = row_fg, h = H_ROW[["fg"]]),
  row_h  = list(p = pGrow,  h = H_ROW[["h"]]),
  row_i  = list(p = pHrow,  h = H_ROW[["i"]])
)
for (nm in names(qa_rows)) {
  ragg::agg_png(file.path(QA, paste0(nm, ".png")), width = FIG_W_MM / 25.4,
                height = qa_rows[[nm]]$h / 25.4, units = "in", res = 300,
                background = "white")
  print(qa_rows[[nm]]$p)
  dev.off()
}

gp <- function(ct, sc, col) paired_tab[[col]][paired_tab$celltype == ct & paired_tab$scope == sc]
gs <- function(ct, contrast, col) shift_tab[[col]][shift_tab$celltype == ct & shift_tab$contrast == contrast]

writeLines(c(
  "Figure 9 statistics",
  sprintf("Spots: %s total; L1-L3 %s; L4-L6 %s; outside %s. None excluded.",
          fmt_n(n_spot), fmt_n(n_l13), fmt_n(n_l46), fmt_n(n_out)),
  "L1-L3 membership is the original lesion spot list. L4-L6 membership is",
  "spots inside the added rectangles. The two sets do not overlap and are not pooled.",
  sprintf("Panel c crop: %s spots in the L1-L3 frame, %s inside L1-L3.",
          fmt_n(n_win_c), fmt_n(n_set_c)),
  sprintf("Panel d crop: %s spots in the L4-L6 frame, %s inside L4-L6.",
          fmt_n(n_win_d), fmt_n(n_set_d)),
  "Panel e: two-sided Wilcoxon signed-rank of (EAE_POS - EAE_NEG) vs 0,",
  "Benjamini-Hochberg across the six cell-type x window tests.",
  sprintf("  T cell    Outside r_rb %+.3f %s", gp("T cell", "Outside", "r_rb"), fmt_q(gp("T cell", "Outside", "q"))),
  sprintf("  T cell    L1-L3   r_rb %+.3f %s", gp("T cell", "L1-L3", "r_rb"), fmt_q(gp("T cell", "L1-L3", "q"))),
  sprintf("  T cell    L4-L6   r_rb %+.3f %s", gp("T cell", "L4-L6", "r_rb"), fmt_q(gp("T cell", "L4-L6", "q"))),
  sprintf("  Microglia Outside r_rb %+.3f %s", gp("Microglia", "Outside", "r_rb"), fmt_q(gp("Microglia", "Outside", "q"))),
  sprintf("  Microglia L1-L3   r_rb %+.3f %s", gp("Microglia", "L1-L3", "r_rb"), fmt_q(gp("Microglia", "L1-L3", "q"))),
  sprintf("  Microglia L4-L6   r_rb %+.3f %s", gp("Microglia", "L4-L6", "r_rb"), fmt_q(gp("Microglia", "L4-L6", "q"))),
  "Panel e rank-sum, each window set vs outside, BH across the four tests.",
  "Positive r_rb means the window set has the higher (EAE_POS - EAE_NEG) difference.",
  sprintf("  T cell    L1-L3 vs outside r_rb %+.3f %s", gs("T cell", "L1-L3", "r_rb"), fmt_q(gs("T cell", "L1-L3", "q"))),
  sprintf("  T cell    L4-L6 vs outside r_rb %+.3f %s", gs("T cell", "L4-L6", "r_rb"), fmt_q(gs("T cell", "L4-L6", "q"))),
  sprintf("  Microglia L1-L3 vs outside r_rb %+.3f %s", gs("Microglia", "L1-L3", "r_rb"), fmt_q(gs("Microglia", "L1-L3", "q"))),
  sprintf("  Microglia L4-L6 vs outside r_rb %+.3f %s", gs("Microglia", "L4-L6", "r_rb"), fmt_q(gs("Microglia", "L4-L6", "q"))),
  sprintf("Panel g all-spots Fisher OR %.2f (%.2f-%.2f) %s", fF$or, fF$lo, fF$hi, fmt_p(fF$p)),
  sprintf("Panel g L4-L6 Fisher OR %.2f (%.2f-%.2f) %s", f46$or, f46$lo, f46$hi, fmt_p(f46$p)),
  "L1-L3 EAE_NEG microglia n = 3, so that contrast is not plotted.",
  sprintf("Panel h vessel-distance Fisher OR %.2f (%.2f-%.2f) %s", fG$or, fG$lo, fG$hi, fmt_p(fG$p)),
  sprintf("Panel h Spearman rho %+.3f %s", as.numeric(sr$estimate), fmt_p(sr$p.value)),
  sprintf("Panel i all-spots rank-sum %s r_rb %+.3f",
          fmt_p(dist_tab$p[dist_tab$scope == "All spots"]),
          dist_tab$r_rb[dist_tab$scope == "All spots"]),
  sprintf("Panel i L4-L6 rank-sum %s r_rb %+.3f (q after BH across the two tests %s)",
          fmt_p(dist_tab$p[dist_tab$scope == "L4-L6"]),
          dist_tab$r_rb[dist_tab$scope == "L4-L6"],
          fmt_q(dist_tab$q[dist_tab$scope == "L4-L6"])),
  "Positive r_rb in panel i means MgIV is farther from the vessel mask than MgIVC.",
  sprintf("Figure size: %.1f x %.1f mm", FIG_W_MM, FIG_H_MM),
  "One EAE section. n is spatial spots, not biological replicates.",
  "H&E crops are unmodified. Distances are low-resolution image pixels."
), file.path(OUTDIR, "Figure9_stats.txt"))

writeLines(c(
  "Fig. 9 | Vascular-contact and parenchymal programmes across L1-L3 and L4-L6 windows of one EAE brain section.",
  "",
  sprintf("a, T-cell EAE_POS and EAE_NEG prediction scores on all %s spots. Both maps share one colour scale (0 to %.2f). Dashed boxes mark L1-L3 and L4-L6.", fmt_n(n_spot), cap_tcell),
  sprintf("b, Microglial EAE_POS and EAE_NEG scores for the same spots, on a shared scale (0 to %.2f).", cap_mg),
  sprintf("c, L1-L3 at higher magnification (%s spots in the frame, %s inside L1-L3). Left, unmodified H&E. Right, the four scores in a and b.", fmt_n(n_win_c), fmt_n(n_set_c)),
  sprintf("d, L4-L6 in the same layout (%s spots in the frame, %s inside L4-L6).", fmt_n(n_win_d), fmt_n(n_set_d)),
  sprintf(paste0(
    "e, Within-spot EAE_POS minus EAE_NEG difference outside the windows (n = %s), ",
    "inside L1-L3 (n = %s) and inside L4-L6 (n = %s). Violins show the full distribution, ",
    "boxes the median and interquartile range, and white diamonds the mean. Labels give ",
    "the two-sided Wilcoxon signed-rank test against zero, Benjamini-Hochberg across six tests, ",
    "and the matched-pairs rank-biserial correlation. Two-sided rank-sum tests of each window ",
    "set against outside tissue, Benjamini-Hochberg across four tests: T cell L1-L3 r_rb %+.2f, ",
    "L4-L6 r_rb %+.2f; microglia L1-L3 r_rb %+.2f, L4-L6 r_rb %+.2f. The two window sets are not pooled."
  ), fmt_n(n_out), fmt_n(n_l13), fmt_n(n_l46),
     gs("T cell", "L1-L3", "r_rb"), gs("T cell", "L4-L6", "r_rb"),
     gs("Microglia", "L1-L3", "r_rb"), gs("Microglia", "L4-L6", "r_rb")),
  "f, Cxcl10 counts in EAE_POS (MgIV) and EAE_NEG (MgIVC) microglial spots. Grey marks every assigned spot. Dashed boxes are the six windows.",
  sprintf("g, Percentage of microglial spots with at least one Cxcl10 count, for all microglial spots and for L4-L6. Bars are observed percentages and error bars Wilson 95%% intervals. Two-sided Fisher exact tests: all spots, odds ratio %.1f (95%% CI %.1f-%.1f), %s; L4-L6, odds ratio %.1f (95%% CI %.1f-%.1f), %s. L1-L3 contains only 3 EAE_NEG microglial spots, so that contrast is not plotted.", fF$or, fF$lo, fF$hi, fmt_p(fF$p), f46$or, f46$lo, f46$hi, fmt_p(f46$p)),
  sprintf("h, Cxcl10-positive spots on the Vwf+ Cdh5+ vessel mask, the perivascular (<= %d px) versus distal contrast, and the Cxcl10-positive fraction in %d equal-count distance bins. Every spot, including L4-L6, enters this distance analysis. Fisher odds ratio %.2f (95%% CI %.2f-%.2f), %s. Spearman correlation of distance against Cxcl10 count: rho = %+.3f, %s. Distances are image pixels.", NEAR_PX, NBIN, fG$or, fG$lo, fG$hi, fmt_p(fG$p), as.numeric(sr$estimate), fmt_p(sr$p.value)),
  sprintf("i, Microglial spots on the vessel mask, the distance-to-vessel map, and cell-to-vessel distance for MgIV versus MgIVC among all microglial spots and inside L4-L6. Two-sided rank-sum tests, Benjamini-Hochberg across the two tests. Positive r_rb means MgIV is farther from the mask. L1-L3 is not drawn as its own violin because it contains only 3 MgIVC spots."),
  "",
  "Panels a-i are one EAE brain section, so n is spatial spots rather than biological replicates. Source data are provided as a Source Data file."
), file.path(OUTDIR, "Figure9_legend.txt"))

writeLines(c(
  "Figure 9 | new folder. Does not overwrite the source figure or any earlier folder.",
  "Structural adaptation of the previous Figure 8.",
  "Panel d repeats the panel-c five-tile layout on the L4-L6 bounding box.",
  "Panels e, g and i test L4-L6 separately from L1-L3.",
  "The previous pMCAO glycolysis panel is omitted.",
  "Backend: R exclusive."
), file.path(OUTDIR, "Figure9_provenance.txt"))

writeLines(c(
  "Figure 9 QA",
  sprintf("Size %.1f x %.1f mm. R only. 600 dpi PNG/TIFF, SVG, cairo_pdf.", FIG_W_MM, FIG_H_MM),
  "Glyph floor 5.2 pt. Panel tags 8 pt.",
  "cairo_pdf Tf = 1 is a CTM false positive; judge size from the SVG.",
  "L1-L3 and L4-L6 are not pooled.",
  "L1-L3 EAE_NEG microglia n = 3; that contrast is reported and not plotted.",
  "H&E tiles are unmodified crops. No micrometre scale bar.",
  "validate_figure.py is recorded after the render."
), file.path(OUTDIR, "Figure9_QA.txt"))

cat("DONE\n")
