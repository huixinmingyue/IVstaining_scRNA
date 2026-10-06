#!/usr/bin/env Rscript
# =============================================================================
# Fig. S9 v1 — Nature multipanel (FigS9_Nature_v1)
# New folder; does not overwrite FigS5_NC_submission or earlier SI packages.
#
# Source figure: NC Fig. S5 submission (spatial markers + immune-lineage maps).
# Changed relative to that source
#   c-h SpatialFeaturePlot pt.size.factor 0.22 -> 1.00
#   h   microglia (Mg) IV | IVC prediction-score maps (new)
#   i   former panel h (NEG vs POS boxplots) with Mg added
#
# Core conclusion
#   EAE Visium resolves cluster marker programmes and maps IV (POS) versus
#   IVC (NEG) immune-lineage prediction scores, including microglia, onto
#   lesion territory with quantified group differences.
#
# Evidence chain
#   a  Which genes mark Visium clusters?                         heatmap
#   b  Where are proliferation / myeloid / Mg markers?           gene maps
#   c-h Where do DC, B, Neu1, PC, Mac and Mg scores sit?        spatial
#   i  Do POS and NEG scores differ for each lineage?            Wilcoxon
#
# Archetype : image plate + quant
# Backend   : R exclusive
# Size      : 183 x 235 mm
# Reuse     : structural adaptation of FigS5_NC_submission/make_figureS5.R
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS9_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
SRCFIG  <- file.path(OUTDIR, "source_figures")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 235
width_mm = 183
height_mm = 235
FAM  <- "Arial"
BASE <- 6.5
SMALL <- 5.4
PT_SIZE <- 1.00
has_sm <- requireNamespace("scattermore", quietly = TRUE)

pal_feat <- c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF8C00", "#FF0000")
clust_lvls <- as.character(0:8)
clust_pal <- setNames(
  c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
    "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"),
  clust_lvls
)
fill_box <- c(NEG = "#9B59B6", POS = "#E67E22")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = BASE),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.4),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(3.0, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid = element_blank()
  )
tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.01, 0.99)
)

fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}

fix_images <- function(sce) {
  for (nm in names(sce@images)) {
    img <- sce@images[[nm]]
    ok <- tryCatch({ slot(img, "misc"); TRUE }, error = function(e) FALSE)
    if (!ok) {
      slot(img, "misc", check = FALSE) <- list()
      sce@images[[nm]] <- img
    }
  }
  sce
}

fmtp <- function(p) {
  if (!is.finite(p) || p < 2.2e-16) "<2.2e-16" else sprintf("%.3g", p)
}

hm <- readRDS(file.path(SRC, "heatmap_slim.rds"))
full <- readRDS(file.path(SRC, "spatial_full_slim.rds"))
imm_slim <- readRDS(file.path(SRC, "spatial_immune_slim.rds"))
stopifnot(all(c("MgIV", "MgIVC") %in% rownames(imm_slim$pred)))

# -----------------------------------------------------------------------------
# a — cluster marker heatmap
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
set.seed(1)
mk <- hm$markers
mk$cluster <- as.character(mk$cluster)
gene_keep <- unique(unlist(lapply(split(mk, mk$cluster), function(d) {
  head(unique(as.character(d$gene)), 5)
}), use.names = FALSE))
gene_ord2 <- intersect(hm$gene_order, gene_keep)
mat <- hm$mat[gene_ord2, , drop = FALSE]
cl <- hm$cell_cluster
n_a_before <- length(cl)
keep <- unlist(lapply(split(seq_along(cl), cl), function(ix) {
  if (length(ix) <= 60) ix else sample(ix, 60)
}), use.names = FALSE)
mat <- mat[, keep, drop = FALSE]
cl <- cl[keep]
n_a_after <- length(cl)
mat_c <- pmax(pmin(mat, 2.5), -2.5)

write.csv(
  mk[, intersect(c("gene", "cluster", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"),
                 names(mk)), drop = FALSE],
  file.path(PLOTTED, "SourceData_FigS9a_markers.csv"),
  row.names = FALSE
)

df_a <- as.data.frame(as.table(mat_c), stringsAsFactors = FALSE)
colnames(df_a) <- c("gene", "cell", "z")
df_a$gene <- factor(df_a$gene, levels = rev(rownames(mat)))
df_a$cell <- factor(df_a$cell, levels = colnames(mat_c))
bar <- data.frame(
  cell = factor(colnames(mat_c), levels = colnames(mat_c)),
  cluster = factor(cl, levels = clust_lvls),
  y = 1
)
pA_bar <- ggplot(bar, aes(cell, y, fill = cluster)) +
  geom_tile(height = 1, width = 1, colour = NA) +
  scale_fill_manual(values = clust_pal, name = "Cluster", drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_family = FAM) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = BASE - 0.5),
    legend.text = element_text(size = BASE - 1),
    legend.key.size = unit(2.5, "mm"),
    plot.margin = margin(0, 0, 0, 0)
  )
pA_hm <- ggplot(df_a, aes(cell, gene, fill = z)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradientn(
    colours = c("#FF00FF", "#000000", "#FFFF00"),
    limits = c(-2.5, 2.5), oob = squish, name = "Scaled\nexpr."
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = "a") + theme_pub +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 5.2, face = "italic", family = FAM),
    legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.3)
  )
pA <- (pA_bar / pA_hm) + plot_layout(heights = c(0.06, 1), guides = "collect")

# -----------------------------------------------------------------------------
# b — gene spatial maps
# -----------------------------------------------------------------------------
spot_gene <- function(gene, letter = NULL, limits = NULL, pt = 0.12, show_legend = TRUE) {
  vals <- as.numeric(full$expr[gene, ]); names(vals) <- colnames(full$expr)
  if (is.null(limits)) {
    hi <- as.numeric(quantile(vals[vals > 0], 0.995, names = FALSE, na.rm = TRUE))
    if (!is.finite(hi) || hi <= 0) hi <- max(vals, na.rm = TRUE)
    limits <- c(0, max(hi, 0.1))
  }
  df <- data.frame(
    x = full$coords[, "x"], y = -full$coords[, "y"],
    val = pmin(pmax(vals[rownames(full$coords)], limits[1]), limits[2])
  )
  df <- df[order(df$val), , drop = FALSE]
  arr <- fade_he(full$he$image, 0.90)
  rast <- grDevices::as.raster(arr)
  ymin <- -(dim(arr)[1] + 0.5); ymax <- -0.5
  xmin <- 0.5; xmax <- dim(arr)[2] + 0.5
  p <- ggplot(df, aes(x, y, colour = val)) +
    annotation_raster(rast, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  if (has_sm) {
    p <- p + scattermore::geom_scattermore(pointsize = 1.05, pixels = c(800, 800),
                                           interpolate = TRUE)
  } else {
    p <- p + geom_point(size = pt, stroke = 0)
  }
  p <- p +
    scale_colour_gradientn(colours = pal_feat, limits = limits, oob = squish,
                           name = "Expr.") +
    coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    labs(title = gene, tag = letter) +
    theme_void(base_family = FAM) +
    theme(
      plot.title = element_text(size = BASE, face = "bold.italic", hjust = 0.5,
                                family = FAM),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.tag.position = c(0.01, 0.99),
      legend.title = element_text(size = BASE - 0.5),
      legend.text = element_text(size = BASE - 0.8),
      legend.key.height = unit(5, "mm"), legend.key.width = unit(2.0, "mm"),
      plot.margin = margin(1, 1, 1, 1)
    )
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

cat("Panel b ...\n")
lim_mki <- c(0, max(2, as.numeric(quantile(full$expr["Mki67", ], 0.995))))
lim_lga <- c(0, max(2, as.numeric(quantile(full$expr["Lgals3", ], 0.995))))
lim_ccl <- c(0, max(2, as.numeric(quantile(full$expr["Ccl2", ], 0.995))))
lim_hex <- c(0, max(2, as.numeric(quantile(full$expr["Hexb", ], 0.995))))
pB <- (spot_gene("Mki67", "b", lim_mki) | spot_gene("Lgals3", NULL, lim_lga)) /
  (spot_gene("Ccl2", NULL, lim_ccl) | spot_gene("Hexb", NULL, lim_hex))

# -----------------------------------------------------------------------------
# c–h — SpatialFeaturePlot (larger spots); h = Mg
# -----------------------------------------------------------------------------
cat("Loading predictions object for c-h ...\n")
spe_i <- fix_images(readRDS(file.path(SRC, "spe_immune_predictions.rds")))
DefaultAssay(spe_i) <- "predictions"
n_immune <- ncol(spe_i)
stopifnot(all(c("MgIV", "MgIVC") %in% rownames(spe_i)))
cat("spe_immune dim:", paste(dim(spe_i), collapse = "x"), "\n")

theme_sf <- theme_void(base_family = FAM) +
  theme(
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    legend.title = element_text(size = BASE - 0.5, family = FAM),
    legend.text = element_text(size = BASE - 0.8, family = FAM),
    legend.key.height = unit(5, "mm"),
    legend.key.width = unit(2, "mm"),
    plot.margin = margin(1, 1, 1, 1)
  )

sf_one <- function(feat, max.cutoff, title, letter = NULL, show_legend = TRUE) {
  args <- list(
    object = spe_i, features = feat, pt.size.factor = PT_SIZE,
    min.cutoff = 0, crop = TRUE, stroke = 0,
    image.alpha = 0.9, alpha = c(0.1, 1)
  )
  if (!is.null(max.cutoff)) args$max.cutoff <- max.cutoff
  p <- do.call(SpatialFeaturePlot, args)
  lim <- if (is.null(max.cutoff)) {
    v <- as.numeric(GetAssayData(spe_i, assay = "predictions", layer = "data")[feat, ])
    c(0, max(0.1, as.numeric(quantile(v, 0.995, names = FALSE, na.rm = TRUE))))
  } else c(0, max.cutoff)
  p <- p +
    scale_fill_gradientn(colours = pal_feat, limits = lim, name = "Score") +
    ggtitle(title) + labs(tag = letter) + theme_sf
  if (!show_legend) p <- p + guides(fill = "none")
  p
}

sf_pair <- function(f_pos, f_neg, max.cutoff, letter, lineage) {
  cat("  SpatialFeaturePlot", letter, lineage, "pt.size.factor", PT_SIZE, "\n")
  left <- sf_one(f_pos, max.cutoff, paste(lineage, "IV"), letter, FALSE)
  right <- sf_one(f_neg, max.cutoff, paste(lineage, "IVC"), NULL, TRUE)
  wrap_elements(full = (left | right))
}

q995 <- function(feat) {
  as.numeric(quantile(
    GetAssayData(spe_i, assay = "predictions", layer = "data")[feat, ],
    0.995, names = FALSE, na.rm = TRUE
  ))
}

cat("Building panels c-h from Seurat ...\n")
pC <- sf_pair("DCIV", "DCIVC", 0.65, "c", "DC")
pD <- sf_pair("BcellIV", "BcellIVC", 0.15, "d", "B cell")
pE <- sf_pair("Neu1IV", "Neu1IVC", 0.08, "e", "Neu1")
pc_hi <- max(q995("PCIV"), q995("PCIVC"), 0.1)
pF <- sf_pair("PCIV", "PCIVC", pc_hi, "f", "PC")
pG <- sf_pair("MacIV", "MacIVC", 0.45, "g", "Mac")
mg_hi <- max(q995("MgIV"), q995("MgIVC"), 0.1)
pH <- sf_pair("MgIV", "MgIVC", mg_hi, "h", "Mg")
rm(spe_i); gc(verbose = FALSE)

# -----------------------------------------------------------------------------
# i — boxplots (former h) including Mg
# -----------------------------------------------------------------------------
cat("Panel i ...\n")
make_box <- function(f_pos, f_neg, xlab, letter = NULL, show_y = TRUE) {
  df <- data.frame(
    score = c(as.numeric(imm_slim$pred[f_pos, ]), as.numeric(imm_slim$pred[f_neg, ])),
    group = factor(rep(c("POS", "NEG"), each = ncol(imm_slim$pred)),
                   levels = c("NEG", "POS"))
  )
  df <- df[is.finite(df$score), , drop = FALSE]
  n_neg <- sum(df$group == "NEG")
  n_pos <- sum(df$group == "POS")
  wt <- wilcox.test(score ~ group, data = df)
  stars <- if (wt$p.value < 1e-4) "****" else if (wt$p.value < 0.001) "***" else
    if (wt$p.value < 0.01) "**" else if (wt$p.value < 0.05) "*" else "ns"
  y_max <- max(df$score, na.rm = TRUE)
  med <- tapply(df$score, df$group, median, na.rm = TRUE)
  q25 <- tapply(df$score, df$group, quantile, probs = 0.25, na.rm = TRUE)
  q75 <- tapply(df$score, df$group, quantile, probs = 0.75, na.rm = TRUE)
  set.seed(1)
  df_j <- if (nrow(df) > 3000) df[sample.int(nrow(df), 3000), , drop = FALSE] else df
  p <- ggplot(df, aes(group, score, fill = group)) +
    geom_jitter(data = df_j, width = 0.18, size = 0.15, alpha = 0.15,
                stroke = 0, colour = "grey25") +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
    scale_fill_manual(values = fill_box, guide = "none") +
    annotate(
      "text", x = 1.5, y = min(1.15, y_max * 1.08 + 0.02), label = stars,
      size = BASE / ggplot2::.pt, fontface = "bold", family = FAM
    ) +
    labs(x = xlab, y = if (show_y) "Prediction score" else NULL, tag = letter) +
    coord_cartesian(ylim = c(0, max(1, y_max * 1.2))) +
    theme_pub
  list(
    plot = p, p = wt$p.value, n_neg = n_neg, n_pos = n_pos,
    med_NEG = unname(med["NEG"]), med_POS = unname(med["POS"]),
    q25_NEG = unname(q25["NEG"]), q75_NEG = unname(q75["NEG"]),
    q25_POS = unname(q25["POS"]), q75_POS = unname(q75["POS"])
  )
}

pairs <- list(
  DC = c("DCIV", "DCIVC"), B = c("BcellIV", "BcellIVC"),
  Neu1 = c("Neu1IV", "Neu1IVC"), PC = c("PCIV", "PCIVC"),
  Mac = c("MacIV", "MacIVC"), Mg = c("MgIV", "MgIVC")
)
labs_i <- c(DC = "DC", B = "B cell", Neu1 = "Neu1", PC = "PC", Mac = "Mac", Mg = "Mg")
box_out <- list(); pI_list <- list(); i <- 0
for (nm in names(pairs)) {
  i <- i + 1
  res <- make_box(
    pairs[[nm]][1], pairs[[nm]][2], labs_i[[nm]],
    letter = if (i == 1) "i" else NULL, show_y = i == 1
  )
  box_out[[nm]] <- res; pI_list[[nm]] <- res$plot
  cat(sprintf("  %s p=%s n_NEG=%d n_POS=%d\n", nm, fmtp(res$p), res$n_neg, res$n_pos))
}
pI <- wrap_plots(pI_list, nrow = 1)
stats_i <- data.frame(
  lineage = names(pairs),
  feature_POS = vapply(pairs, function(x) x[[1]], ""),
  feature_NEG = vapply(pairs, function(x) x[[2]], ""),
  n_NEG = vapply(box_out, function(x) x$n_neg, 0L),
  n_POS = vapply(box_out, function(x) x$n_pos, 0L),
  median_NEG = vapply(box_out, function(x) x$med_NEG, 0.0),
  median_POS = vapply(box_out, function(x) x$med_POS, 0.0),
  q25_NEG = vapply(box_out, function(x) x$q25_NEG, 0.0),
  q75_NEG = vapply(box_out, function(x) x$q75_NEG, 0.0),
  q25_POS = vapply(box_out, function(x) x$q25_POS, 0.0),
  q75_POS = vapply(box_out, function(x) x$q75_POS, 0.0),
  wilcox_p = vapply(box_out, function(x) x$p, 0.0),
  test = "two-sided Wilcoxon rank-sum",
  stringsAsFactors = FALSE
)
write.csv(stats_i, file.path(PLOTTED, "SourceData_FigS9i_wilcoxon.csv"),
          row.names = FALSE)
write.csv(stats_i, file.path(SRC, "FigS9i_wilcoxon.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- (pA | pB) + plot_layout(widths = c(1.0, 1.15))
row2 <- (pC | pD | pE)
row3 <- (pF | pG | pH)
fig <- (row1 / row2 / row3 / pI) +
  plot_layout(heights = c(1.28, 0.82, 0.82, 0.92)) +
  plot_annotation(theme = theme(
    plot.margin = margin(3, 3, 3, 3),
    plot.background = element_rect(fill = "white", colour = NA)
  ))

# -----------------------------------------------------------------------------
# Stats / legend / provenance
# -----------------------------------------------------------------------------
p_txt <- paste(sprintf("%s %s", stats_i$lineage, vapply(stats_i$wilcox_p, fmtp, "")),
               collapse = "; ")
stats <- c(
  "Figure S9 statistics / integrity (FigS9_Nature_v1)",
  sprintf("Immune Visium spots: %d (POS scores and NEG scores each n = %d).",
          n_immune, stats_i$n_NEG[1]),
  sprintf("Panel a heatmap display cells: %d -> %d (<=60 spots/cluster; markers use full table)",
          n_a_before, n_a_after),
  sprintf("Panels c-h: SpatialFeaturePlot pt.size.factor = %.2f (source Fig. S5 used 0.22)",
          PT_SIZE),
  sprintf("Panel h Mg colour scale max = %.3f (max of MgIV/MgIVC 99.5th percentiles)",
          mg_hi),
  sprintf("Panel i Wilcoxon (two-sided rank-sum): %s", p_txt),
  "Panel i jitter is display-subsampled to 3,000 points; Wilcoxon uses all finite scores.",
  "Image integrity: H&E global alpha 0.9; no contrast stretch.",
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
)
writeLines(stats, file.path(OUTDIR, "FigureS9_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. S9 | Spatial cluster markers and immune-lineage mapping, including microglia, in the EAE Visium section.",
  "",
  "a, Heatmap of top marker genes per Visium cluster (0-8); expression scaled and clipped to [-2.5, 2.5].",
  "b, Spatial maps of Mki67, Lgals3, Ccl2 and Hexb on Visium H&E (image alpha = 0.9).",
  "c-h, Prediction-score maps for DC (c), B cells (d), Neu1 (e), plasma cells (f), macrophages (g) and microglia (h); left, IV (POS); right, IVC (NEG); image alpha = 0.9.",
  paste0(
    "i, Prediction scores for the same lineages (NEG vs POS). Boxes show median and IQR; whiskers 1.5\u00d7IQR. ",
    "Two-sided Wilcoxon rank-sum test; ****P < 1\u00d710\u207b\u2074. Exact P: ",
    p_txt,
    sprintf(". n = %s spots per group (immune Visium subset).",
            format(stats_i$n_NEG[1], big.mark = ","))
  ),
  "",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS9_legend.txt"))

writeLines(c(
  "Figure S9 | new folder FigS9_Nature_v1 (does not overwrite FigS5_NC_submission).",
  "Reuse level: structural adaptation of NC Fig. S5 (make_figureS5.R).",
  "Panels a-g: same observations and cutoffs as source Fig. S5 a-g.",
  sprintf("Panels c-h: pt.size.factor increased 0.22 -> %.2f.", PT_SIZE),
  "Panel h: MgIV / MgIVC SpatialFeaturePlot (not in source Fig. S5).",
  "Panel i: former panel h plus Mg Wilcoxon on all 32,811 spots per group.",
  "Backend: R exclusive."
), file.path(OUTDIR, "FigureS9_provenance.txt"))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
base <- file.path(OUTDIR, "FigureS9")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig)
dev.off()

cat("Writing PDF ...\n")
pdf_tmp <- paste0(base, "_tmp.pdf")
tryCatch({
  grDevices::cairo_pdf(pdf_tmp, width = w_in, height = h_in, family = FAM)
  print(fig)
  dev.off()
  file.copy(pdf_tmp, paste0(base, ".pdf"), overwrite = TRUE)
  unlink(pdf_tmp)
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("cairo_pdf failed: ", conditionMessage(e))
  grDevices::pdf(paste0(base, ".pdf"), width = w_in, height = h_in,
                 useDingbats = FALSE)
  print(fig)
  dev.off()
})

cat("Writing PNG ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig)
dev.off()

cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(base, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig)
dev.off()

writeLines(c(
  "Build notes (full visual QA is in FigureS9_QA.txt)",
  "- Backend: R exclusive",
  "- a, b, i: native ggplot vector",
  sprintf("- c-h: Seurat SpatialFeaturePlot; H&E alpha 0.9; pt.size.factor %.2f", PT_SIZE),
  "- h: Mg added; i: former h plus Mg Wilcoxon",
  sprintf("- Exports: PDF, SVG, PNG, TIFF at %d x %d mm, raster 600 dpi",
          FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "FigureS9_build_notes.txt"))

cat("DONE\n")
