#!/usr/bin/env Rscript
# =============================================================================
# Fig. 4 - Nature build (double column, 183 mm), R / ggplot2 exclusive
#
# Core conclusion
#   In the spleen Stereo-seq section, the follicular B-cell and plasma-cell
#   niches carry a higher CD45-IV-high reference signature than the matched
#   CD45-IV-low signature in the same spatial bins, whereas the CD4 T-cell
#   compartment shows the opposite, weaker shift.
#
# Evidence chain (one unique question per panel row)
#   a/b/c  B cells      : where is the niche | how does the score map | how big
#   d/e/f  Plasma cells : is the effect reproduced in a second lineage
#   g/h/i  CD4 T cells  : is the effect lineage-specific in direction
#
# Archetype : image plate + quant (3 lineages x [niche map | score maps | quant])
#
# Differences from the earlier build of this figure, and why
#   1. Paired statistics. Every spatial bin carries both the CD45-IV-low and the
#      CD45-IV-high reference score, so the two columns are matched, not
#      independent samples. The earlier build used an unpaired rank-sum test on
#      overlapping marginal distributions. This build tests the per-bin
#      difference with a two-sided Wilcoxon signed-rank test and reports the
#      median paired difference with its 95% confidence interval.
#   2. One signature per row. The B-cell row previously mapped one reference
#      build (spcd45) and quantified a different one (spcd452). Map and
#      quantification now use the same reference build within each row; the
#      unused columns are still written to source data.
#   3. One cluster palette for the whole figure. Cluster 5 was cyan in one niche
#      map and green in another; cluster identity now has a single colour.
#   4. Perceptually ordered score scale. The earlier rainbow ramp
#      (viridis segment spliced onto orange/red) is replaced by a single
#      sequential ramp in which high score is dark on a white page, and the
#      colour limits are taken from the plotted distribution instead of a
#      display cutoff that left most of the ramp unused.
#   5. Shared coordinate frame. All six spatial panels use the same x/y limits,
#      so positions are directly comparable, with the whole section drawn in
#      light grey behind every subset.
#   6. Group names are the biological contrast (CD45-IV low / high reference)
#      rather than POS/NEG codes, and n is printed on the panel.
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(grid)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_Nature_v2")
SRC    <- file.path(OUTDIR, "source_data")
SD     <- file.path(OUTDIR, "source_data_plotted")
QA     <- file.path(OUTDIR, "qa")
for (d in c(SD, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 197
FAM      <- "Arial"
BASE     <- 6.5

# Stereo-seq chip geometry: the bin grid steps by 50 chip units and the DNB
# pitch is 500 nm, i.e. bin 50 = 25 um. Confirm before submission; the scale bar
# below is the only place this constant is used.
UM_PER_UNIT <- 0.5
SCALEBAR_UM <- 1000

SZ_BG   <- 0.10
SZ_HI   <- 0.26
SZ_FEAT <- 0.34
SZ_JIT  <- 0.16

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args <- list(reg()$Helvetica); names(args) <- FAM
      try(do.call(reg, args), silent = TRUE)
    }
  }
}

# ---- shared vocabulary ------------------------------------------------------
# One colour per louvain cluster across every niche map (Okabe-Ito).
clus_cols <- c("1" = "#0072B2", "4" = "#E69F00", "5" = "#009E73",
               "8" = "#CC79A7", "9" = "#D55E00")
# Light enough that the palest scored bins still separate from the section.
BG_COL    <- "#EBEBEB"
# High score is dark on the white page and the ramp is monotone in lightness.
# The palest fifth is dropped: on white paper those bins read as absent rather
# than as low-scoring.
score_ramp <- grDevices::hcl.colors(256, "YlGnBu", rev = TRUE)[52:256]
# Neutral reference vs lineage accent, identical logic in all three quant panels.
LOW_FILL  <- "#BFBFBF"
lineage_col <- c(B = "#0072B2", Plasma = "#009E73", T4 = "#CC79A7")

LAB_LOW  <- "CD45-IV\nlow"
LAB_HIGH <- "CD45-IV\nhigh"
Y_LAB    <- "Prediction score"

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title   = element_text(size = BASE - 0.3, colour = "black"),
    axis.text    = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.5),
    legend.text  = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title    = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5, colour = "grey30"),
    plot.tag      = element_text(size = 8, face = "bold", family = FAM,
                                 hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    plot.margin  = margin(1, 2, 1, 1),
    panel.grid   = element_blank()
  )

theme_sp <- theme_void(base_family = FAM) +
  theme(
    plot.title   = element_text(size = BASE, face = "bold", hjust = 0.5,
                                family = FAM),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5,
                                 colour = "grey30", family = FAM),
    plot.tag     = element_text(size = 8, face = "bold", family = FAM,
                                hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    legend.title = element_text(size = BASE - 0.5, family = FAM),
    legend.text  = element_text(size = BASE - 0.8, family = FAM),
    plot.margin  = margin(1, 1, 1, 1)
  )

# Exact P wherever the value survives double precision; the floor is reported
# only when the computed value has underflowed to zero.
fmtp <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p <= 0) return("< 2.2e-16")
  sprintf("= %.1e", p)
}
stars_from_p <- function(p) {
  if (!is.finite(p)) "ns"
  else if (p < 1e-4) "****" else if (p < 1e-3) "***"
  else if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns"
}

# ---- load and define the shared coordinate frame ----------------------------
needed <- c("panel_a_slim.rds", "panel_bc_slim.rds", "panel_d_slim.rds",
            "panel_ef_plasma_slim.rds", "panel_g_slim.rds", "panel_hi_slim.rds")
miss <- needed[!file.exists(file.path(SRC, needed))]
if (length(miss)) stop("missing source data: ", paste(miss, collapse = ", "))

niche_a <- readRDS(file.path(SRC, "panel_a_slim.rds"))
niche_d <- readRDS(file.path(SRC, "panel_d_slim.rds"))
niche_g <- readRDS(file.path(SRC, "panel_g_slim.rds"))
bc      <- readRDS(file.path(SRC, "panel_bc_slim.rds"))
pl      <- readRDS(file.path(SRC, "panel_ef_plasma_slim.rds"))
tf      <- readRDS(file.path(SRC, "panel_hi_slim.rds"))

section <- niche_a$coords
section$xp <-  as.numeric(section$x)
section$yp <- -as.numeric(section$y)
PAD  <- 150
XLIM <- range(section$xp) + c(-PAD, PAD)
YLIM <- range(section$yp) + c(-PAD, PAD)
cat(sprintf("section: %d bins, %.2f x %.2f mm at %.2f um per unit\n",
            nrow(section), diff(range(section$xp)) * UM_PER_UNIT / 1000,
            diff(range(section$yp)) * UM_PER_UNIT / 1000, UM_PER_UNIT))

# The section outline is context, not evidence, and it is repeated behind all
# six spatial panels. Bins sit on a regular grid, so each row of bins is
# collapsed into its runs of consecutive occupied positions and drawn as
# rectangles. The covered area is identical to the bin cloud, but the file holds
# a few hundred rectangles instead of six copies of 30,624 points.
STEP <- 50
runs <- do.call(rbind, lapply(split(section, section$yp), function(g) {
  xs <- sort(unique(g$xp))
  brk <- c(0, which(diff(xs) != STEP), length(xs))
  do.call(rbind, lapply(seq_len(length(brk) - 1), function(i) {
    seg <- xs[(brk[i] + 1):brk[i + 1]]
    data.frame(xmin = min(seg) - STEP / 2, xmax = max(seg) + STEP / 2,
               ymin = g$yp[1] - STEP / 2, ymax = g$yp[1] + STEP / 2)
  }))
}))
cat(sprintf("section outline: %d rectangles for %d bins\n",
            nrow(runs), nrow(section)))

section_layer <- function() {
  geom_rect(data = runs,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = BG_COL, colour = NA, inherit.aes = FALSE)
}

scalebar_len <- SCALEBAR_UM / UM_PER_UNIT
sb_x1 <- XLIM[2] - 0.06 * diff(XLIM)
sb_x0 <- sb_x1 - scalebar_len
sb_y  <- YLIM[1] + 0.045 * diff(YLIM)

STATS <- list()

# -----------------------------------------------------------------------------
# Niche maps (a, d, g): the whole section in grey, the requested clusters in
# colour, drawn in one shared coordinate frame.
# -----------------------------------------------------------------------------
make_niche <- function(slim, title, tag, scalebar = FALSE) {
  df <- slim$coords
  df$xp <-  as.numeric(df$x)
  df$yp <- -as.numeric(df$y)
  keep <- slim$keep
  hi <- df[df$show %in% keep, , drop = FALSE]
  hi$show <- factor(hi$show, levels = keep)
  p <- ggplot() +
    section_layer() +
    geom_point(data = hi, aes(xp, yp, colour = show),
               size = SZ_HI, stroke = 0, shape = 16) +
    scale_colour_manual(values = clus_cols, name = "Louvain cluster",
                        drop = TRUE) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.1),
                                 title.position = "left")) +
    coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = title, tag = tag) +
    theme_sp +
    # The tag sits in the left margin: at this panel width a centred title long
    # enough to name the clusters otherwise runs into it.
    theme(plot.margin = margin(1, 1, 1, 4),
          plot.tag.position = c(0, 1),
          legend.position = "bottom",
          legend.key.size = unit(2.0, "mm"),
          legend.box.spacing = unit(0.8, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(colour = "grey80", fill = NA,
                                      linewidth = 0.25))
  if (scalebar) {
    p <- p +
      annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
               linewidth = 0.6, colour = "black") +
      annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y + 0.028 * diff(YLIM),
               label = paste0(SCALEBAR_UM / 1000, " mm"),
               size = (BASE - 1.0) / .pt, family = FAM, colour = "black")
  }
  p
}

# -----------------------------------------------------------------------------
# Score maps (b, e, h): one sequential ramp per row, shared limits within the
# row so that Low and High are directly comparable.
# -----------------------------------------------------------------------------
make_feat <- function(coords, vals, title, limits, tag = NULL) {
  df <- data.frame(xp =  as.numeric(coords$x),
                   yp = -as.numeric(coords$y),
                   val = pmin(pmax(as.numeric(vals), limits[1]), limits[2]))
  df <- df[order(df$val), , drop = FALSE]
  ggplot() +
    section_layer() +
    geom_point(data = df, aes(xp, yp, colour = val),
               size = SZ_FEAT, stroke = 0, shape = 16) +
    scale_colour_gradientn(colours = score_ramp, limits = limits, oob = squish,
                           name = Y_LAB,
                           breaks = pretty(limits, n = 4)) +
    guides(colour = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(22, "mm"), barheight = unit(1.6, "mm"),
      ticks.colour = "white", ticks.linewidth = 0.4,
      frame.colour = "black", frame.linewidth = 0.25)) +
    coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = title, tag = tag) +
    theme_sp +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "grey80", fill = NA,
                                      linewidth = 0.25))
}

# Colour limits come from the data actually plotted in the row, not from the
# display cutoff of the earlier build. Those cutoffs sat far above the bulk of
# the distribution, so almost the whole ramp went unused and the low and high
# maps looked alike. The 1st-99th percentile of the pooled low and high scores
# is used instead, with out-of-range bins squished to the end colours.
row_limits <- function(pred, low_col, high_col) {
  v <- c(as.numeric(pred[, low_col]), as.numeric(pred[, high_col]))
  q <- unname(quantile(v[is.finite(v)], c(0.01, 0.99)))
  lim <- c(floor(q[1] * 100) / 100, ceiling(q[2] * 100) / 100)
  cat(sprintf("  colour limits %.2f-%.2f (%.1f%% of bins squished)\n",
              lim[1], lim[2],
              100 * mean(v < lim[1] | v > lim[2], na.rm = TRUE)))
  lim
}

make_feat_pair <- function(coords, pred, low_col, high_col, limits, tag) {
  pL <- make_feat(coords, pred[, low_col],  "CD45-IV-low reference",  limits, tag)
  pH <- make_feat(coords, pred[, high_col], "CD45-IV-high reference", limits, NULL)
  (pL | pH) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box.spacing = unit(0.8, "mm"))
}

# -----------------------------------------------------------------------------
# Quantification (c, f, i)
# The two scores are measured in the same spatial bin, so the contrast is a
# paired one: the test is a two-sided Wilcoxon signed-rank test on the per-bin
# difference, and the reported effect is the median paired difference.
# -----------------------------------------------------------------------------
make_box <- function(pred, low_col, high_col, accent, tag, title) {
  low  <- as.numeric(pred[, low_col])
  high <- as.numeric(pred[, high_col])
  ok   <- is.finite(low) & is.finite(high)
  n_in <- length(low); low <- low[ok]; high <- high[ok]
  if (sum(!ok)) cat("  dropped", sum(!ok), "of", n_in, "bins with non-finite score\n")

  wt <- wilcox.test(high, low, paired = TRUE, conf.int = TRUE,
                    exact = FALSE, correct = TRUE)
  d  <- high - low

  df <- data.frame(
    score = c(low, high),
    group = factor(rep(c(LAB_LOW, LAB_HIGH), each = length(low)),
                   levels = c(LAB_LOW, LAB_HIGH))
  )
  y_max <- max(df$score)
  y_ann <- y_max * 1.05
  y_top <- y_max * 1.26

  ann <- sprintf("median \u0394 = %+.3f (%s)\n95%% CI %+.3f to %+.3f\nP %s, n = %s bins",
                 unname(wt$estimate), stars_from_p(wt$p.value),
                 wt$conf.int[1], wt$conf.int[2],
                 fmtp(wt$p.value), comma(length(low)))

  p <- ggplot(df, aes(group, score, fill = group)) +
    geom_jitter(width = 0.17, height = 0, size = SZ_JIT, alpha = 0.14,
                stroke = 0, colour = "grey30", shape = 16) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.92,
                 linewidth = 0.3, colour = "black") +
    annotate("text", x = 1.5, y = y_ann, label = ann,
             size = (BASE - 1.2) / .pt, family = FAM, vjust = 0,
             lineheight = 1.05, colour = "grey15") +
    scale_fill_manual(values = setNames(c(LOW_FILL, accent),
                                        c(LAB_LOW, LAB_HIGH)), guide = "none") +
    scale_y_continuous(limits = c(0, y_top), breaks = pretty(c(0, y_max), n = 5),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = NULL, y = Y_LAB, title = title, tag = tag) +
    theme_pub +
    theme(axis.text.x = element_text(size = BASE - 0.8, lineheight = 0.95))

  list(plot = p, wilcox = wt, n = length(low),
       med_low = median(low), med_high = median(high), delta = d,
       df = data.frame(score_low = low, score_high = high, delta = d))
}

# ---- row 1: follicular B cells ---------------------------------------------
cat("Row 1 (B cells) ...\n")
lim_b <- row_limits(bc$pred, "Bcell.spcd452.low", "Bcell.spcd452.high")
pA <- make_niche(niche_a, "B-cell niche (clusters 1, 5)", "a", scalebar = TRUE)
# The cluster numbers stay in the titles of d and g without the word "clusters",
# which the legend below each map already supplies.
pB <- make_feat_pair(bc$coords, bc$pred,
                     "Bcell.spcd452.low", "Bcell.spcd452.high", lim_b, "b")
qC <- make_box(bc$pred, "Bcell.spcd452.low", "Bcell.spcd452.high",
               lineage_col[["B"]], "c", "Follicular B cells")
pC <- qC$plot

# ---- row 2: plasma cells ----------------------------------------------------
cat("Row 2 (plasma cells) ...\n")
lim_p <- row_limits(pl$pred, "Plasma.spcd452.low", "Plasma.spcd452.high")
pD <- make_niche(niche_d, "Peri-follicular niche (4, 5, 8, 9)", "d")
pE <- make_feat_pair(pl$coords, pl$pred,
                     "Plasma.spcd452.low", "Plasma.spcd452.high", lim_p, "e")
qF <- make_box(pl$pred, "Plasma.spcd452.low", "Plasma.spcd452.high",
               lineage_col[["Plasma"]], "f", "Plasma cells")
pF <- qF$plot

# ---- row 3: CD4 T cells -----------------------------------------------------
cat("Row 3 (CD4 T cells) ...\n")
lim_t <- row_limits(tf$pred, "CD4T.spcd45.low", "CD4T.spcd45.high")
pG <- make_niche(niche_g, "B/T follicular niche (1, 4, 5, 8)", "g")
pH <- make_feat_pair(tf$coords, tf$pred,
                     "CD4T.spcd45.low", "CD4T.spcd45.high", lim_t, "h")
qI <- make_box(tf$pred, "CD4T.spcd45.low", "CD4T.spcd45.high",
               lineage_col[["T4"]], "i", "CD4 T cells")
pI <- qI$plot

# ---- statistics table -------------------------------------------------------
qlist <- list(c = qC, f = qF, i = qI)
stats <- do.call(rbind, lapply(names(qlist), function(k) {
  q <- qlist[[k]]
  data.frame(panel = k,
             lineage = c(c = "Follicular B cells", f = "Plasma cells",
                         i = "CD4 T cells")[[k]],
             n_bins = q$n,
             median_low = q$med_low, median_high = q$med_high,
             median_paired_delta = unname(q$wilcox$estimate),
             ci_low = q$wilcox$conf.int[1], ci_high = q$wilcox$conf.int[2],
             p_value = q$wilcox$p.value,
             test = "two-sided Wilcoxon signed-rank (paired by spatial bin)",
             stringsAsFactors = FALSE)
}))
print(stats, row.names = FALSE)
write.csv(stats, file.path(SD, "SourceData_Fig4_statistics.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Assemble. Rows are sealed so that no row distorts the column widths of another.
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
seal <- function(x) wrap_elements(full = x)
W <- c(1.02, 1.95, 0.95)
row1 <- seal((pA | pB | pC) + plot_layout(widths = W))
row2 <- seal((pD | pE | pF) + plot_layout(widths = W))
row3 <- seal((pG | pH | pI) + plot_layout(widths = W))

fig <- (row1 / row2 / row3) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(theme = theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", colour = NA)
  ))

# ---- export -----------------------------------------------------------------
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "Figure4")

cat("Writing SVG ...\n")
svglite::svglite(paste0(stem, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in,
                     family = FAM, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PNG ...\n")
ragg::agg_png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())

# ---- source data exactly as plotted -----------------------------------------
write.csv(cbind(bc$coords, as.data.frame(bc$pred)),
          file.path(SD, "SourceData_Fig4bc_Bcell_bins.csv"), row.names = FALSE)
write.csv(cbind(pl$coords, as.data.frame(pl$pred)),
          file.path(SD, "SourceData_Fig4ef_Plasma_bins.csv"), row.names = FALSE)
write.csv(cbind(tf$coords, as.data.frame(tf$pred)),
          file.path(SD, "SourceData_Fig4hi_CD4T_bins.csv"), row.names = FALSE)
write.csv(data.frame(barcode = niche_a$coords$barcode,
                     x = niche_a$coords$x, y = niche_a$coords$y,
                     cluster = niche_a$coords$cluster,
                     shown_in_a = niche_a$coords$show %in% niche_a$keep,
                     shown_in_d = niche_d$coords$show %in% niche_d$keep,
                     shown_in_g = niche_g$coords$show %in% niche_g$keep),
          file.path(SD, "SourceData_Fig4adg_section_clusters.csv"),
          row.names = FALSE)

# ---- legend -----------------------------------------------------------------
dir_word <- function(q) if (unname(q$wilcox$estimate) > 0) "higher" else "lower"
row_txt <- function(q) sprintf(
  "Boxes show the median and interquartile range, whiskers 1.5x the interquartile range, and each point one spatial bin. The CD45-IV-high reference score is %s than the CD45-IV-low score in the same bin (median paired difference %+.3f, 95%% CI %+.3f to %+.3f; two-sided Wilcoxon signed-rank P %s; n = %s bins).",
  dir_word(q), unname(q$wilcox$estimate), q$wilcox$conf.int[1],
  q$wilcox$conf.int[2], fmtp(q$wilcox$p.value), comma(q$n))

writeLines(c(
  "Fig. 4 | Spatial CD45-IV reference signatures in splenic B-cell, plasma-cell and CD4 T-cell compartments.",
  sprintf("a, Louvain clusters 1 and 5 (B-cell niche) on the spleen Stereo-seq section; the complete section is drawn in grey behind the highlighted bins. All six spatial panels share one coordinate frame. Scale bar, %s mm.", SCALEBAR_UM / 1000),
  sprintf("b, Per-bin prediction scores for the CD45-IV-low and CD45-IV-high follicular B-cell reference signatures (Bcell.spcd452). The two maps share one colour scale spanning %.2f to %.2f, the 1st to 99th percentile of the pooled scores; bins outside that range take the end colour.", lim_b[1], lim_b[2]),
  paste("c, Prediction scores in follicular B-cell bins.", row_txt(qC)),
  "d, Louvain clusters 4, 5, 8 and 9 (peri-follicular architecture containing the plasma-cell bins).",
  sprintf("e, As in b for the plasma-cell reference signatures (Plasma.spcd452); shared colour scale %.2f to %.2f.", lim_p[1], lim_p[2]),
  paste("f, Prediction scores in plasma-cell bins.", row_txt(qF)),
  "g, Louvain clusters 1, 4, 5 and 8 (B/T follicular architecture).",
  sprintf("h, As in b for the CD4 T-cell reference signatures (CD4T.spcd45); shared colour scale %.2f to %.2f.", lim_t[1], lim_t[2]),
  paste("i, Prediction scores in CD4 T-cell bins.", row_txt(qI)),
  "Statistics are paired because the two reference scores are measured in the same spatial bin; no bin was excluded. Source data are provided with this figure."
), file.path(OUTDIR, "Figure4_legend.txt"), useBytes = TRUE)

cat("DONE ->", OUTDIR, "\n")
