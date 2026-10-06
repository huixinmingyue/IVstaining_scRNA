#!/usr/bin/env Rscript
# =============================================================================
# Fig. S3-2 — Nature SI true-vector (ggplot2 geom_point only)
# Updates vs prior truevector:
#   a: highlight cluster 5 only
#   d: highlight clusters 4 and 8 only
#   b/e titles: POS_B/NEG_B and POS_T/NEG_T (match c/f x-axis)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-2_NC_submission")
SRC <- file.path(OUTDIR, "source_data")
FIG_W_MM <- 183
FIG_H_MM <- 155
FAM <- "Arial"
BASE <- 6.5

SZ_BG <- 0.12
SZ_HI <- 0.28
SZ_FEAT <- 0.40
SZ_JIT <- 0.20

pal_feat <- c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725", "#FD8D3C", "#E31A1C")
pal_a <- c("5" = "#00BFC4")
pal_d <- c("4" = "#7CAE00", "8" = "#00BFC4")
fill_b <- c(POS_B = "#9B59B6", NEG_B = "#E67E22")
fill_t <- c(POS_T = "#3498DB", NEG_T = "#E74C3C")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = BASE),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.4),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.8, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid = element_blank()
  )

theme_sp <- theme_void(base_family = FAM) +
  theme(
    plot.title = element_text(size = BASE - 0.2, face = "bold", hjust = 0.5, family = FAM),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    legend.title = element_text(size = BASE - 0.5, family = FAM),
    legend.text = element_text(size = BASE - 0.8, family = FAM),
    legend.key.height = unit(2.2, "mm"),
    legend.key.width = unit(2.2, "mm"),
    plot.margin = margin(1, 1, 1, 1)
  )

fmtp <- function(p) if (!is.finite(p) || p < 2.2e-16) "<2.2e-16" else sprintf("%.3g", p)

prep_xy <- function(coords) {
  df <- coords
  df$y <- -as.numeric(df$y)
  df$x <- as.numeric(df$x)
  df$cluster <- as.character(df$cluster)
  df
}

# ---------- Panel a: cluster 5 only ----------
cat("Panel a (cluster 5 only) ...\n")
a <- readRDS(file.path(SRC, "panel_a_slim.rds"))
df_a <- prep_xy(a$coords)
df_a$hl <- ifelse(df_a$cluster == "5", "5", "bg")
df_bg <- df_a[df_a$hl == "bg", , drop = FALSE]
df_hi <- df_a[df_a$hl == "5", , drop = FALSE]
df_hi$hl <- factor(df_hi$hl, levels = "5")
pA <- ggplot() +
  geom_point(data = df_bg, aes(x, y), colour = "#B0B0B0", size = SZ_BG, stroke = 0, shape = 16) +
  geom_point(data = df_hi, aes(x, y, colour = hl), size = SZ_HI, stroke = 0, shape = 16) +
  scale_colour_manual(values = pal_a, name = "Cluster") +
  coord_fixed(expand = FALSE) +
  labs(tag = "a") + theme_sp +
  theme(legend.position = "right")

# ---------- Panel b: POS_B | NEG_B (match panel c) ----------
cat("Panel b (POS_B / NEG_B) ...\n")
bc <- readRDS(file.path(SRC, "panel_bc_slim.rds"))
lim_b <- c(0, 0.55)
make_feat <- function(coords, vals, title, letter = NULL, limits, show_leg = TRUE) {
  df <- data.frame(
    x = as.numeric(coords$x),
    y = -as.numeric(coords$y),
    val = pmin(pmax(as.numeric(vals), limits[1]), limits[2])
  )
  df <- df[order(df$val), , drop = FALSE]
  p <- ggplot(df, aes(x, y, colour = val)) +
    geom_point(size = SZ_FEAT, stroke = 0, shape = 16) +
    scale_colour_gradientn(colours = pal_feat, limits = limits, oob = squish, name = "Score") +
    guides(colour = guide_colourbar(
      display = "rectangles", nbin = 100,
      ticks.colour = "black", frame.colour = "black"
    )) +
    coord_fixed(expand = FALSE) +
    labs(title = title, tag = letter) + theme_sp +
    theme(
      legend.position = "right",
      legend.key.height = unit(6, "mm"),
      legend.key.width = unit(1.8, "mm")
    )
  if (!show_leg) p <- p + guides(colour = "none")
  p
}
# left POS_B = high; right NEG_B = low (same order as panel c)
pB <- make_feat(bc$coords, bc$pred[, "Bcell.spcd45.high"], "POS_B", "b", lim_b, FALSE) |
  make_feat(bc$coords, bc$pred[, "Bcell.spcd45.low"], "NEG_B", NULL, lim_b, TRUE)

# ---------- Panel c ----------
cat("Panel c ...\n")
df_c <- data.frame(
  score = c(bc$pred[, "Bcell.spcd45.high"], bc$pred[, "Bcell.spcd45.low"]),
  group = factor(
    rep(c("POS_B", "NEG_B"), each = nrow(bc$pred)),
    levels = c("POS_B", "NEG_B")
  )
)
df_c <- df_c[is.finite(df_c$score), , drop = FALSE]
wt_c <- wilcox.test(score ~ group, data = df_c)
stars_c <- if (wt_c$p.value < 1e-4) "****" else if (wt_c$p.value < 0.05) "*" else "ns"
y_c <- max(df_c$score, na.rm = TRUE)
pC <- ggplot(df_c, aes(group, score, fill = group)) +
  geom_jitter(width = 0.18, size = SZ_JIT, alpha = 0.18, stroke = 0, colour = "grey25") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
  scale_fill_manual(values = fill_b, guide = "none") +
  annotate("text", x = 1.5, y = min(0.78, y_c * 1.08), label = stars_c,
           size = BASE / ggplot2::.pt, fontface = "bold", family = FAM) +
  scale_y_continuous(limits = c(0, max(0.8, y_c * 1.12)), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Prediction score", tag = "c") +
  theme_pub

# ---------- Panel d: clusters 4 and 8 only ----------
cat("Panel d (clusters 4 + 8) ...\n")
d <- readRDS(file.path(SRC, "panel_d_slim.rds"))
df_d <- prep_xy(d$coords)
keep_d <- c("4", "8")
df_d$hl <- ifelse(df_d$cluster %in% keep_d, df_d$cluster, "bg")
df_dbg <- df_d[df_d$hl == "bg", , drop = FALSE]
df_dhi <- df_d[df_d$hl %in% keep_d, , drop = FALSE]
df_dhi$hl <- factor(df_dhi$hl, levels = keep_d)
pD <- ggplot() +
  geom_point(data = df_dbg, aes(x, y), colour = "#B0B0B0", size = SZ_BG, stroke = 0, shape = 16) +
  geom_point(data = df_dhi, aes(x, y, colour = hl), size = SZ_HI, stroke = 0, shape = 16) +
  scale_colour_manual(values = pal_d, name = "Cluster") +
  coord_fixed(expand = FALSE) +
  labs(tag = "d") + theme_sp +
  theme(legend.position = "right")

# ---------- Panel e: POS_T | NEG_T (match panel f) ----------
cat("Panel e (POS_T / NEG_T) ...\n")
tf <- readRDS(file.path(SRC, "panel_ef_slim.rds"))
lim_t <- c(0, 0.45)
pE <- make_feat(tf$coords, tf$pred[, "CD4T.spcd45.high"], "POS_T", "e", lim_t, FALSE) |
  make_feat(tf$coords, tf$pred[, "CD4T.spcd45.low"], "NEG_T", NULL, lim_t, TRUE)

# ---------- Panel f ----------
cat("Panel f ...\n")
df_f <- data.frame(
  score = c(tf$pred[, "CD4T.spcd45.high"], tf$pred[, "CD4T.spcd45.low"]),
  group = factor(
    rep(c("POS_T", "NEG_T"), each = nrow(tf$pred)),
    levels = c("POS_T", "NEG_T")
  )
)
df_f <- df_f[is.finite(df_f$score), , drop = FALSE]
wt_f <- wilcox.test(score ~ group, data = df_f)
stars_f <- if (wt_f$p.value < 1e-4) "****" else if (wt_f$p.value < 0.05) "*" else "ns"
pF <- ggplot(df_f, aes(group, score, fill = group)) +
  geom_jitter(width = 0.18, size = SZ_JIT, alpha = 0.18, stroke = 0, colour = "grey25") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
  scale_fill_manual(values = fill_t, guide = "none") +
  annotate("text", x = 1.5, y = 0.62, label = stars_f,
           size = BASE / ggplot2::.pt, fontface = "bold", family = FAM) +
  coord_cartesian(ylim = c(0, 0.65)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.65)) +
  labs(x = NULL, y = "Prediction score", tag = "f") +
  theme_pub

stats <- data.frame(
  panel = c("c", "f"),
  comparison = c("POS_B vs NEG_B (Bcell.spcd45 high vs low)",
                 "POS_T vs NEG_T (CD4T.spcd45 high vs low)"),
  n_POS = c(sum(df_c$group == "POS_B"), sum(df_f$group == "POS_T")),
  n_NEG = c(sum(df_c$group == "NEG_B"), sum(df_f$group == "NEG_T")),
  wilcox_p = c(wt_c$p.value, wt_f$p.value),
  test = "two-sided Wilcoxon rank-sum",
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(SRC, "FigS3-2_wilcoxon.csv"), row.names = FALSE)
cat("  a n_cluster5=", nrow(df_hi),
    " d n_4+8=", nrow(df_dhi),
    " c p=", fmtp(wt_c$p.value),
    " f p=", fmtp(wt_f$p.value), "\n")

row1 <- (pA | pB | pC) + plot_layout(widths = c(1.05, 1.45, 0.75))
row2 <- (pD | pE | pF) + plot_layout(widths = c(1.05, 1.45, 0.75))
fig <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(theme = theme(
    plot.margin = margin(3, 3, 3, 3),
    plot.background = element_rect(fill = "white", colour = NA)
  ))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS3-2")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); dev.off()

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM)
print(fig); grDevices::dev.off()

cat("Writing PNG 600 dpi ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in", res = 600, background = "white")
print(fig); dev.off()

writeLines(c(
  "Fig. S3-2 | Stereo-seq B and T niches with CD45 prediction scores.",
  "a, Spatial map of louvain cluster 5 on the spleen Stereo-seq section.",
  "b, Prediction-score maps for POS_B and NEG_B (Bcell.spcd45 high/low; max cutoff 0.55).",
  paste0(
    "c, Prediction scores for POS_B vs NEG_B. Boxes show median and IQR; whiskers 1.5×IQR. ",
    "Two-sided Wilcoxon; ****P < 1×10⁻⁴. Exact P = ", fmtp(wt_c$p.value),
    sprintf(". n = %d scores per group.", stats$n_POS[1])
  ),
  "d, Spatial map of louvain clusters 4 and 8.",
  "e, Prediction-score maps for POS_T and NEG_T (CD4T.spcd45 high/low; max cutoff 0.45).",
  paste0(
    "f, Prediction scores for POS_T vs NEG_T; y-axis limited to 0.65. ",
    "Two-sided Wilcoxon; ****P < 1×10⁻⁴. Exact P = ", fmtp(wt_f$p.value),
    sprintf(". n = %d scores per group.", stats$n_POS[2]),
    " Source data are provided with this figure."
  )
), file.path(OUTDIR, "FigureS3-2_legend.txt"), useBytes = TRUE)

writeLines(c(
  "Figure S3-2 | Nature SI submission (true vector ggplot)",
  "Core conclusion: Cluster 5 / 4+8 niches map POS vs NEG CD45 prediction programmes.",
  "Archetype: image plate + quant",
  "Backend: R 4.4.1 exclusive; geom_point only (no scattermore)",
  sprintf("Final size: %d x %d mm", FIG_W_MM, FIG_H_MM),
  "Panel a: cluster 5 only",
  "Panel d: clusters 4 and 8 only",
  "Panel b/e titles aligned to c/f: POS_B|NEG_B and POS_T|NEG_T",
  "Colourbars: guide_colourbar(display='rectangles')",
  "Exports: FigureS3-2.pdf / .svg / .png (600 dpi)",
  "New folder FigS3-2_NC_submission/; prior FigS3-2 folders untouched"
), file.path(OUTDIR, "FigureS3-2_QA.txt"), useBytes = TRUE)

writeLines(c(
  "Structural adaptation of FigS3-2_Nature_truevector with revised cluster/title rules.",
  "Outputs under FigS3-2_NC_submission/."
), file.path(OUTDIR, "FigureS3-2_provenance.txt"), useBytes = TRUE)

cat("DONE ->", OUTDIR, "\n")
