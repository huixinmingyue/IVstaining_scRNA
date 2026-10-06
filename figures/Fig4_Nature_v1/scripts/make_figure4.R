#!/usr/bin/env Rscript
# =============================================================================
# Fig. 4 v1 - Nature rebuild in final figures (R / ggplot2 exclusive)
#
# Core conclusion
#   In the spleen Stereo-seq section, follicular B-cell bins (cluster 1 and
#   cluster 5) and plasma-cell bins carry a higher CD45-IV-high than
#   CD45-IV-low reference score in the same spatial bin, whereas CD4 T-cell
#   bins (cluster 8 and cluster 4) shift in the opposite direction.
#
# Evidence chain (one lettered row = one unique subset)
#   a  Follicular B cells, clusters 1+5 niche, Bcell.spcd452 in cluster 1
#   b  Cluster 5 only, Bcell.spcd45 in cluster 5
#   c  Plasma cells, peri-follicular 4/5/8/9, Plasma.spcd452 in cluster 9
#   d  CD4 T cells, B/T follicular 1/4/5/8, CD4T.spcd45 in cluster 8
#   e  Clusters 4+8 niche, CD4T.spcd45 in cluster 4
#
# Archetype : image plate + quant
# Size      : 183 x 312 mm
# New folder — does not overwrite any prior Fig. 4 package
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "Fig4_Nature_v1")
}
SRC4 <- file.path(OUTDIR, "source_data", "fig4v2")
SRCS <- file.path(OUTDIR, "source_data", "s32")
SD   <- file.path(OUTDIR, "source_data_plotted")
QA   <- file.path(OUTDIR, "qa")
for (d in c(SD, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 312
FAM      <- "Arial"
BASE     <- 6.5

UM_PER_UNIT <- 0.5
SCALEBAR_UM <- 1000

SZ_HI   <- 0.26
SZ_FEAT <- 0.34
SZ_JIT  <- 0.16

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args_f <- list(reg()$Helvetica); names(args_f) <- FAM
      try(do.call(reg, args_f), silent = TRUE)
    }
  }
}

clus_cols <- c("1" = "#0072B2", "4" = "#E69F00", "5" = "#009E73",
               "8" = "#CC79A7", "9" = "#D55E00")
BG_COL     <- "#EBEBEB"
score_ramp <- grDevices::hcl.colors(256, "YlGnBu", rev = TRUE)[52:256]
LOW_FILL   <- "#BFBFBF"
lineage_col <- c(B = "#0072B2", B5 = "#0072B2", Plasma = "#009E73",
                 T8 = "#CC79A7", T4 = "#CC79A7")

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

# ---- load -------------------------------------------------------------------
needed4 <- c("panel_a_slim.rds", "panel_bc_slim.rds", "panel_d_slim.rds",
             "panel_ef_plasma_slim.rds", "panel_g_slim.rds", "panel_hi_slim.rds")
neededS <- c("panel_bc_slim.rds", "panel_ef_slim.rds")
miss <- c(needed4[!file.exists(file.path(SRC4, needed4))],
          neededS[!file.exists(file.path(SRCS, neededS))])
if (length(miss)) stop("missing source data: ", paste(miss, collapse = ", "))

niche_sec <- readRDS(file.path(SRC4, "panel_a_slim.rds"))
bc        <- readRDS(file.path(SRC4, "panel_bc_slim.rds"))
pl        <- readRDS(file.path(SRC4, "panel_ef_plasma_slim.rds"))
tf        <- readRDS(file.path(SRC4, "panel_hi_slim.rds"))
bc5       <- readRDS(file.path(SRCS, "panel_bc_slim.rds"))
tf4       <- readRDS(file.path(SRCS, "panel_ef_slim.rds"))

section <- niche_sec$coords
section$xp <-  as.numeric(section$x)
section$yp <- -as.numeric(section$y)
PAD  <- 150
XLIM <- range(section$xp) + c(-PAD, PAD)
YLIM <- range(section$yp) + c(-PAD, PAD)
cat(sprintf("section: %d bins, %.2f x %.2f mm at %.2f um per unit\n",
            nrow(section), diff(range(section$xp)) * UM_PER_UNIT / 1000,
            diff(range(section$yp)) * UM_PER_UNIT / 1000, UM_PER_UNIT))

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

# -----------------------------------------------------------------------------
make_niche <- function(coords, keep, title, tag, scalebar = FALSE) {
  df <- coords
  df$xp <-  as.numeric(df$x)
  df$yp <- -as.numeric(df$y)
  keep <- as.character(keep)
  hi <- df[as.character(df$cluster) %in% keep, , drop = FALSE]
  hi$show <- factor(as.character(hi$cluster), levels = keep)
  stopifnot(nrow(hi) > 0)
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

row_limits <- function(pred, low_col, high_col) {
  v <- c(as.numeric(pred[, low_col]), as.numeric(pred[, high_col]))
  q <- unname(quantile(v[is.finite(v)], c(0.01, 0.99)))
  lim <- c(floor(q[1] * 100) / 100, ceiling(q[2] * 100) / 100)
  cat(sprintf("  colour limits %.2f-%.2f (%.1f%% of bins squished)\n",
              lim[1], lim[2],
              100 * mean(v < lim[1] | v > lim[2], na.rm = TRUE)))
  lim
}

make_feat_pair <- function(coords, pred, low_col, high_col, limits) {
  pL <- make_feat(coords, pred[, low_col],  "CD45-IV-low reference",  limits, NULL)
  pH <- make_feat(coords, pred[, high_col], "CD45-IV-high reference", limits, NULL)
  (pL | pH) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box.spacing = unit(0.8, "mm"))
}

make_box <- function(pred, low_col, high_col, accent, title) {
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
    labs(x = NULL, y = Y_LAB, title = title) +
    theme_pub +
    theme(axis.text.x = element_text(size = BASE - 0.8, lineheight = 0.95))

  list(plot = p, wilcox = wt, n = length(low),
       med_low = median(low), med_high = median(high), delta = d,
       df = data.frame(score_low = low, score_high = high, delta = d))
}

# ---- rows -------------------------------------------------------------------
sec <- niche_sec$coords

cat("Row a (follicular B, clusters 1+5, Bcell.spcd452) ...\n")
lim_a <- row_limits(bc$pred, "Bcell.spcd452.low", "Bcell.spcd452.high")
pA_n <- make_niche(sec, c("1", "5"), "B-cell niche (clusters 1, 5)", "a",
                   scalebar = TRUE)
pA_f <- make_feat_pair(bc$coords, bc$pred,
                       "Bcell.spcd452.low", "Bcell.spcd452.high", lim_a)
qA   <- make_box(bc$pred, "Bcell.spcd452.low", "Bcell.spcd452.high",
                 lineage_col[["B"]], "Follicular B cells")

cat("Row b (cluster 5, Bcell.spcd45) ...\n")
lim_b <- row_limits(bc5$pred, "Bcell.spcd45.low", "Bcell.spcd45.high")
pB_n <- make_niche(sec, "5", "Cluster 5", "b")
pB_f <- make_feat_pair(bc5$coords, bc5$pred,
                       "Bcell.spcd45.low", "Bcell.spcd45.high", lim_b)
qB   <- make_box(bc5$pred, "Bcell.spcd45.low", "Bcell.spcd45.high",
                 lineage_col[["B5"]], "Cluster 5 B cells")

cat("Row c (plasma, clusters 4/5/8/9, Plasma.spcd452) ...\n")
lim_c <- row_limits(pl$pred, "Plasma.spcd452.low", "Plasma.spcd452.high")
pC_n <- make_niche(sec, c("4", "5", "8", "9"),
                   "Peri-follicular niche (4, 5, 8, 9)", "c")
pC_f <- make_feat_pair(pl$coords, pl$pred,
                       "Plasma.spcd452.low", "Plasma.spcd452.high", lim_c)
qC   <- make_box(pl$pred, "Plasma.spcd452.low", "Plasma.spcd452.high",
                 lineage_col[["Plasma"]], "Plasma cells")

cat("Row d (CD4 T, clusters 1/4/5/8, CD4T.spcd45 cluster 8) ...\n")
lim_d <- row_limits(tf$pred, "CD4T.spcd45.low", "CD4T.spcd45.high")
pD_n <- make_niche(sec, c("1", "4", "5", "8"),
                   "B/T follicular niche (1, 4, 5, 8)", "d")
pD_f <- make_feat_pair(tf$coords, tf$pred,
                       "CD4T.spcd45.low", "CD4T.spcd45.high", lim_d)
qD   <- make_box(tf$pred, "CD4T.spcd45.low", "CD4T.spcd45.high",
                 lineage_col[["T8"]], "CD4 T cells")

cat("Row e (clusters 4+8, CD4T.spcd45 cluster 4) ...\n")
lim_e <- row_limits(tf4$pred, "CD4T.spcd45.low", "CD4T.spcd45.high")
pE_n <- make_niche(sec, c("4", "8"), "Clusters 4, 8", "e")
pE_f <- make_feat_pair(tf4$coords, tf4$pred,
                       "CD4T.spcd45.low", "CD4T.spcd45.high", lim_e)
qE   <- make_box(tf4$pred, "CD4T.spcd45.low", "CD4T.spcd45.high",
                 lineage_col[["T4"]], "Cluster 4 CD4 T cells")

# ---- statistics -------------------------------------------------------------
qlist <- list(a = qA, b = qB, c = qC, d = qD, e = qE)
lin_lab <- c(a = "Follicular B cells (cluster 1, Bcell.spcd452)",
             b = "Cluster 5 B cells (Bcell.spcd45)",
             c = "Plasma cells (cluster 9, Plasma.spcd452)",
             d = "CD4 T cells (cluster 8, CD4T.spcd45)",
             e = "Cluster 4 CD4 T cells (CD4T.spcd45)")
stats <- do.call(rbind, lapply(names(qlist), function(k) {
  q <- qlist[[k]]
  data.frame(panel = k, lineage = lin_lab[[k]], n_bins = q$n,
             median_low = q$med_low, median_high = q$med_high,
             median_paired_delta = unname(q$wilcox$estimate),
             ci_low = q$wilcox$conf.int[1], ci_high = q$wilcox$conf.int[2],
             p_value = q$wilcox$p.value,
             test = "two-sided Wilcoxon signed-rank (paired by spatial bin)",
             stringsAsFactors = FALSE)
}))
print(stats, row.names = FALSE)
write.csv(stats, file.path(SD, "SourceData_Fig4_statistics.csv"), row.names = FALSE)

# ---- assemble ---------------------------------------------------------------
cat("Assemble ...\n")
seal <- function(x) wrap_elements(full = x)
W <- c(1.02, 1.95, 0.95)
rowA <- seal((pA_n | pA_f | qA$plot) + plot_layout(widths = W))
rowB <- seal((pB_n | pB_f | qB$plot) + plot_layout(widths = W))
rowC <- seal((pC_n | pC_f | qC$plot) + plot_layout(widths = W))
rowD <- seal((pD_n | pD_f | qD$plot) + plot_layout(widths = W))
rowE <- seal((pE_n | pE_f | qE$plot) + plot_layout(widths = W))

fig <- (rowA / rowB / rowC / rowD / rowE) +
  plot_layout(heights = c(1, 1, 1, 1, 1)) +
  plot_annotation(theme = theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", colour = NA)
  ))

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

cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(stem, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); invisible(dev.off())

write.csv(cbind(bc$coords, as.data.frame(bc$pred)),
          file.path(SD, "SourceData_Fig4a_Bcell_cluster1_bins.csv"), row.names = FALSE)
write.csv(cbind(bc5$coords, as.data.frame(bc5$pred)),
          file.path(SD, "SourceData_Fig4b_Bcell_cluster5_bins.csv"), row.names = FALSE)
write.csv(cbind(pl$coords, as.data.frame(pl$pred)),
          file.path(SD, "SourceData_Fig4c_Plasma_bins.csv"), row.names = FALSE)
write.csv(cbind(tf$coords, as.data.frame(tf$pred)),
          file.path(SD, "SourceData_Fig4d_CD4T_cluster8_bins.csv"), row.names = FALSE)
write.csv(cbind(tf4$coords, as.data.frame(tf4$pred)),
          file.path(SD, "SourceData_Fig4e_CD4T_cluster4_bins.csv"), row.names = FALSE)
write.csv(data.frame(barcode = sec$barcode, x = sec$x, y = sec$y,
                     cluster = sec$cluster,
                     shown_in_a = as.character(sec$cluster) %in% c("1", "5"),
                     shown_in_b = as.character(sec$cluster) %in% "5",
                     shown_in_c = as.character(sec$cluster) %in% c("4", "5", "8", "9"),
                     shown_in_d = as.character(sec$cluster) %in% c("1", "4", "5", "8"),
                     shown_in_e = as.character(sec$cluster) %in% c("4", "8")),
          file.path(SD, "SourceData_Fig4_section_clusters.csv"),
          row.names = FALSE)

cat("DONE ->", OUTDIR, "\n")
