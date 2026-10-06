#!/usr/bin/env Rscript
# =============================================================================
# Fig. S11 — Nature Communications SI
# Archetype: image plate + quant
# Core conclusion: On one MCAO Visium section, inflammatory (cluster 5) and
# proliferative (cluster 3) microglia scores occupy restricted domains, and
# Lgals3, Mki67 and Cxcl10 are higher in the top POS spots than in the top
# NEG spots; Prdx1 is not.
# Reuse: exact geometry of the supplied Fig. S8 plate.
# width_mm = 183
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(grid)
})

set.seed(42)
options(warn = 1)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
ROOT <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else {
  normalizePath(file.path(getwd(), ".."))
}
SRC <- file.path(ROOT, "source_data")
PLOTTED <- file.path(ROOT, "source_data_plotted")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)

width_mm <- 183
height_mm <- 239
BASE <- 6.5
FAM <- "Arial"
fm <- systemfonts::match_fonts(FAM)
cat("font resolve:", FAM, "->", fm$path, "\n")
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args_font <- list(reg()$Helvetica)
    names(args_font) <- FAM
    do.call(reg, args_font)
  }
}
PT <- function(x) x / .pt

theme_set(
  theme_classic(base_size = 6.5, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid = element_blank(),
      text = element_text(family = FAM, size = 6.5),
      legend.position = "right"
    )
)

# Sequential viridis stops. Shared within each numeric range. Grayscale-safe.
scale_spot <- function(limits = NULL, name = NULL) {
  scale_colour_gradientn(
    colours = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
    limits = limits, oob = squish, name = name,
    guide = guide_colourbar(
      barwidth = unit(1.8, "mm"),
      barheight = unit(12, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  )
}

theme_spot <- function() {
  theme_void(base_size = 6.5, base_family = FAM) +
    theme(
      plot.title = element_text(
        size = 6.5, face = "bold", hjust = 0.5, family = FAM,
        margin = margin(b = 1, t = 0)
      ),
      plot.subtitle = element_text(
        size = 6, face = "plain", hjust = 0.5, family = FAM,
        colour = "grey25", margin = margin(b = 1, t = 0)
      ),
      legend.text = element_text(size = 5.5, family = FAM),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1),
      legend.key.height = unit(8, "mm"),
      legend.key.width = unit(1.8, "mm")
    )
}

hdr_plot <- function(txt) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = txt,
      size = PT(6.5), fontface = "bold", family = FAM
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

coord_df <- read.csv(file.path(SRC, "FigS11_coords.csv"), stringsAsFactors = FALSE)
pred_df <- read.csv(file.path(SRC, "FigS11_predictions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
expr_df <- read.csv(file.path(SRC, "FigS11_expression.csv"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(coord_df$barcode, pred_df$barcode))
stopifnot(identical(coord_df$barcode, expr_df$barcode))

coords <- data.frame(x = coord_df$x, y = coord_df$y, row.names = coord_df$barcode)
pred_mat <- as.matrix(pred_df[, setdiff(colnames(pred_df), "barcode"), drop = FALSE])
rownames(pred_mat) <- pred_df$barcode
expr_mat <- as.matrix(expr_df[, c("Lgals3", "Prdx1", "Mki67", "Cxcl10"), drop = FALSE])
rownames(expr_mat) <- expr_df$barcode

he5 <- readRDS(file.path(SRC, "he_immune5.rds"))
cat("spots", nrow(coords), " HE dim", paste(he5$dim, collapse = "x"), "\n")

# Global fade only. Recorded in the figure contract. Not a local contrast edit.
fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}

he_annotation <- function(he, alpha = 0.90) {
  arr <- fade_he(he$image, alpha = alpha)
  nr <- dim(arr)[1]
  nc <- dim(arr)[2]
  arr <- arr[rev(seq_len(nr)), , , drop = FALSE]
  list(
    raster = grDevices::as.raster(arr),
    xmin = 0.5, xmax = nc + 0.5,
    ymin = 0.5, ymax = nr + 0.5
  )
}

spot_plot <- function(values, cells = NULL, title = NULL, subtitle = NULL,
                      limits, pt = 0.55, show_legend = TRUE, he = NULL,
                      he_alpha = 0.90) {
  if (is.null(cells)) cells <- rownames(coords)
  cells <- intersect(cells, rownames(coords))
  df <- data.frame(
    x = coords[cells, "x"],
    y = coords[cells, "y"],
    val = as.numeric(values[cells])
  )
  df <- df[order(df$val), , drop = FALSE]
  p <- ggplot(df, aes(x, y, colour = val))
  if (!is.null(he)) {
    ha <- he_annotation(he, alpha = he_alpha)
    p <- p +
      annotation_raster(
        ha$raster, xmin = ha$xmin, xmax = ha$xmax,
        ymin = ha$ymin, ymax = ha$ymax
      ) +
      geom_point(size = pt, stroke = 0) +
      coord_fixed(
        xlim = c(ha$xmin, ha$xmax),
        ylim = c(ha$ymin, ha$ymax),
        expand = FALSE
      )
  } else {
    p <- p +
      geom_point(size = pt, stroke = 0) +
      coord_fixed()
  }
  p <- p +
    scale_spot(limits = limits) +
    labs(title = title, subtitle = subtitle) +
    theme_spot()
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

score_vec <- function(feat) {
  v <- pred_mat[, feat]
  names(v) <- rownames(pred_mat)
  v
}

# Top 100 by prediction score. Scientific contrast from the source plate, not a render sample.
top_cells <- function(feat, n = 100) {
  v <- score_vec(feat)
  names(sort(v, decreasing = TRUE))[seq_len(n)]
}

# =============================================================================
# Panel a — full-section prediction maps
# =============================================================================
pA5 <- (hdr_plot("Cluster 5 (inflammatory)") /
  (spot_plot(score_vec("n-5"), title = "n-5", limits = c(0, 2.0), pt = 0.42, show_legend = FALSE) |
   spot_plot(score_vec("p-5"), title = "p-5", limits = c(0, 2.0), pt = 0.42, show_legend = TRUE))) +
  plot_layout(heights = c(0.11, 1))

pA3 <- (hdr_plot("Cluster 3 (proliferative)") /
  (spot_plot(score_vec("n-3"), title = "n-3", limits = c(0, 0.5), pt = 0.42, show_legend = FALSE) |
   spot_plot(score_vec("p-3"), title = "p-3", limits = c(0, 0.5), pt = 0.42, show_legend = FALSE))) +
  plot_layout(heights = c(0.11, 1))

pA0 <- (hdr_plot("Cluster 0 (homeostatic)") /
  (spot_plot(score_vec("n-0"), title = "n-0", limits = c(0, 0.5), pt = 0.42, show_legend = FALSE) |
   spot_plot(score_vec("p-0"), title = "p-0", limits = c(0, 0.5), pt = 0.42, show_legend = TRUE))) +
  plot_layout(heights = c(0.11, 1))

pA <- wrap_plots(pA5, pA3, pA0, nrow = 1)

# =============================================================================
# Panel b — expression on the top-100 POS / NEG spots, over the section H&E
# =============================================================================
make_gene_col <- function(gene, n_feat, p_feat, cluster_lab,
                          limits = c(0, 2), show_legend = FALSE) {
  cells_pos <- top_cells(p_feat, 100)
  cells_neg <- top_cells(n_feat, 100)
  v <- expr_mat[, gene]
  names(v) <- rownames(expr_mat)
  p_pos <- spot_plot(
    v, cells_pos, title = gene,
    subtitle = paste0("MCAO_POS \u00b7 ", cluster_lab),
    limits = limits, pt = 1.35, show_legend = FALSE, he = he5
  ) +
    theme(plot.title = element_text(
      size = 7, face = "bold.italic", hjust = 0.5, family = FAM,
      margin = margin(b = 0, t = 0)
    ))
  p_neg <- spot_plot(
    v, cells_neg, title = "MCAO_NEG",
    limits = limits, pt = 1.35, show_legend = show_legend, he = he5
  )
  p_pos / p_neg
}

pB <- wrap_plots(
  make_gene_col("Lgals3", "n-5", "p-5", "cluster 5"),
  make_gene_col("Prdx1", "n-5", "p-5", "cluster 5"),
  make_gene_col("Mki67", "n-3", "p-3", "cluster 3"),
  make_gene_col("Cxcl10", "n-3", "p-3", "cluster 3", show_legend = TRUE),
  nrow = 1
)

# =============================================================================
# Panel c — same top-100 spot sets
# =============================================================================
wilcox_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 0.05) return("*")
  "ns"
}

make_box <- function(df, gene, levels_x) {
  df <- df %>% mutate(cluster = factor(cluster, levels = levels_x))
  wt <- wilcox.test(expression ~ cluster, data = df, exact = FALSE)
  star <- wilcox_stars(wt$p.value)
  y_max <- max(df$expression, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  y_star <- y_max * 1.10
  ggplot(df, aes(cluster, expression, fill = cluster)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(
      width = 0.18, height = 0, size = 0.45, stroke = 0,
      colour = "grey20", alpha = 0.55
    ) +
    scale_fill_manual(values = setNames(c("#D9D9D9", "#595959"), levels_x), guide = "none") +
    annotate(
      "text", x = 1.5, y = y_star, label = star,
      size = PT(7), fontface = "bold", family = FAM
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.28), clip = "off") +
    labs(x = NULL, y = "Expression", title = gene) +
    theme_classic(base_size = 6.5, base_family = FAM) +
    theme(
      plot.title = element_text(size = 7, face = "bold.italic", hjust = 0.5, family = FAM),
      axis.title.y = element_text(size = 6.5, family = FAM),
      axis.text = element_text(size = 6, family = FAM),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid = element_blank(),
      plot.margin = margin(6, 6, 2, 6)
    )
}

read_box <- function(stem) {
  read.csv(file.path(SRC, paste0("FigS8C_", stem, ".csv")), stringsAsFactors = FALSE)
}

lg <- read_box("Lgals3")
pr <- read_box("Prdx1")
mk <- read_box("Mki67")
cx <- read_box("Cxcl10")

pC <- (make_box(lg, "Lgals3", c("n-5", "p-5")) |
       make_box(pr, "Prdx1", c("n-5", "p-5")) |
       make_box(mk, "Mki67", c("n-3", "p-3")) |
       make_box(cx, "Cxcl10", c("n-3", "p-3")))

fig <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  wrap_elements(full = pC) +
  plot_layout(heights = c(1.05, 1.55, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99)
  )

pvals <- data.frame(
  gene = c("Lgals3", "Prdx1", "Mki67", "Cxcl10"),
  group_neg = c("n-5", "n-5", "n-3", "n-3"),
  group_pos = c("p-5", "p-5", "p-3", "p-3"),
  n_neg = c(nrow(lg[lg$cluster == "n-5", ]), nrow(pr[pr$cluster == "n-5", ]),
            nrow(mk[mk$cluster == "n-3", ]), nrow(cx[cx$cluster == "n-3", ])),
  n_pos = c(nrow(lg[lg$cluster == "p-5", ]), nrow(pr[pr$cluster == "p-5", ]),
            nrow(mk[mk$cluster == "p-3", ]), nrow(cx[cx$cluster == "p-3", ])),
  wilcox_p = c(
    wilcox.test(expression ~ cluster, data = lg, exact = FALSE)$p.value,
    wilcox.test(expression ~ cluster, data = pr, exact = FALSE)$p.value,
    wilcox.test(expression ~ cluster, data = mk, exact = FALSE)$p.value,
    wilcox.test(expression ~ cluster, data = cx, exact = FALSE)$p.value
  ),
  test = "two-sided Wilcoxon rank-sum, exact = FALSE, uncorrected",
  unit = "Visium spots on one section (top 100 by prediction score)",
  stringsAsFactors = FALSE
)
pvals$star <- vapply(pvals$wilcox_p, wilcox_stars, character(1))
write.csv(pvals, file.path(PLOTTED, "SourceData_FigS11c_wilcoxon.csv"), row.names = FALSE)
write.csv(lg, file.path(PLOTTED, "SourceData_FigS11c_Lgals3.csv"), row.names = FALSE)
write.csv(pr, file.path(PLOTTED, "SourceData_FigS11c_Prdx1.csv"), row.names = FALSE)
write.csv(mk, file.path(PLOTTED, "SourceData_FigS11c_Mki67.csv"), row.names = FALSE)
write.csv(cx, file.path(PLOTTED, "SourceData_FigS11c_Cxcl10.csv"), row.names = FALSE)

top_tbl <- lapply(c("n-5", "p-5", "n-3", "p-3", "n-0", "p-0"), function(feat) {
  cells <- top_cells(feat, 100)
  data.frame(
    feature = feat,
    rank = seq_along(cells),
    barcode = cells,
    score = score_vec(feat)[cells],
    stringsAsFactors = FALSE
  )
})
write.csv(bind_rows(top_tbl), file.path(PLOTTED, "SourceData_FigS11b_top100_spots.csv"), row.names = FALSE)

stats <- c(
  "Figure S11 | Nature Communications SI",
  "Archetype: image plate + quant",
  sprintf("spots in panel a: %d", nrow(coords)),
  sprintf(
    "Wilcoxon p: Lgals3=%.3g; Prdx1=%.3g; Mki67=%.3g; Cxcl10=%.3g",
    pvals$wilcox_p[1], pvals$wilcox_p[2], pvals$wilcox_p[3], pvals$wilcox_p[4]
  ),
  sprintf("n per group: %s", paste(pvals$n_neg, collapse = ",")),
  sprintf("Size: %d x %d mm; font: %s %.1f pt", width_mm, height_mm, FAM, BASE)
)
writeLines(stats, file.path(ROOT, "qa", "FigureS11_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
base <- file.path(ROOT, "FigureS11")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig)
dev.off()

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM, fallback_resolution = 600)
print(fig)
grDevices::dev.off()

cat("Writing PNG (600 dpi) ...\n")
ragg::agg_png(
  paste0(base, ".png"), width = w_in, height = h_in,
  units = "in", res = 600, background = "white"
)
print(fig)
dev.off()

cat("Writing TIFF (600 dpi) ...\n")
ragg::agg_tiff(
  paste0(base, ".tiff"), width = w_in, height = h_in,
  units = "in", res = 600, background = "white"
)
print(fig)
dev.off()

cat("DONE Figure S11 ->", ROOT, "\n")
