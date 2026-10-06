#!/usr/bin/env Rscript
# =============================================================================
# Fig. S7 — Nature multipanel (native ggplot2 vector; no scattermore)
#
# Core conclusion
#   Lineage markers and microglia gene-set scores map Stroke immune compartments
#   (a–b); n vs p DimPlots and cluster proportions show composition shifts (c–d);
#   GSEA of cluster 5 vs 0 recovers innate programmes in KEGG/Wiki (e–f).
#
# Panel map
#   a  Lineage FeaturePlots (7 genes)
#   b  Module-score FeaturePlots (blue scale)
#   c  DimPlot merged n (n1+n2) vs p (p1+p2)  [NEW]
#   d  Cluster proportions n vs p              [NEW]
#   e  KEGG signaling (was c)
#   f  WikiPathways (was d)
#
# Archetype : quantitative grid
# Backend   : R exclusive
# Folder    : FigS7_NC_np_blue (does not overwrite prior packages)
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS7_NC_np_blue")
SRC <- file.path(OUTDIR, "source_data")
dir.create(file.path(OUTDIR, "scripts"), showWarnings = FALSE, recursive = TRUE)
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 290
BASE <- 6.2
FAM <- "Arial"

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args <- list(reg()$Helvetica)
      names(args) <- FAM
      try(do.call(reg, args), silent = TRUE)
    }
  }
}

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = BASE),
      axis.text = element_text(size = BASE - 0.4, colour = "black"),
      legend.title = element_text(size = BASE - 0.2),
      legend.text = element_text(size = BASE - 0.6),
      legend.key.size = unit(2.8, "mm"),
      plot.title = element_text(size = BASE + 0.3, face = "bold", hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold", family = FAM)
    )
)

theme_feat <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.title = element_blank(),
  legend.key.height = unit(8, "mm"),
  legend.key.width = unit(1.8, "mm"),
  plot.title = element_text(size = BASE - 0.1, face = "bold", hjust = 0.5),
  plot.margin = margin(1, 1, 1, 1)
)

# Native ggplot FeaturePlot helper
feat_plot <- function(df, value_col, title, italic_title = FALSE,
                      low = "#a1a3a6", high = "#840228") {
  d <- df[order(df[[value_col]]), , drop = FALSE]
  ggplot(d, aes(umap_1, umap_2, colour = .data[[value_col]])) +
    geom_point(size = 0.12, stroke = 0, shape = 16) +
    scale_colour_gradient(
      low = low, high = high, name = NULL,
      guide = guide_colorbar(barwidth = unit(1.8, "mm"), barheight = unit(10, "mm"))
    ) +
    coord_equal(expand = TRUE) +
    labs(title = title) +
    theme_feat +
    theme(
      plot.title = element_text(
        size = BASE - 0.1,
        face = if (italic_title) "bold.italic" else "bold",
        hjust = 0.5, family = FAM
      )
    )
}

# Seurat-like categorical cluster palette (hue)
cluster_cols <- hue_pal()(10)
names(cluster_cols) <- as.character(0:9)

# -----------------------------------------------------------------------------
# Panel a — lineage markers (grey -> maroon)
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
df_a <- read.csv(file.path(SRC, "FigS7A_lineage_features.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
genes_a <- c("Tmem119", "Ms4a7", "Ccr2", "S100a9", "Cd3e", "Itgax", "Cd79a")
stopifnot(all(genes_a %in% names(df_a)))
plist_a <- lapply(genes_a, function(g) feat_plot(df_a, g, g, italic_title = TRUE))
pA <- wrap_plots(plist_a, ncol = 4) +
  plot_annotation(title = "Lineage markers") &
  theme(plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0.5, family = FAM))

# -----------------------------------------------------------------------------
# Panel b — module scores (BLUE scale)
# -----------------------------------------------------------------------------
cat("Panel b (blue module scores) ...\n")
df_b <- read.csv(file.path(SRC, "FigS7B_module_scores_StrokeUMAP.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
score_cols <- c("activation", "IAM", "MHC", "lipid")
display_titles <- c(
  "Macrophage / microglia activation",
  "Ischemic-associated microglia",
  "MHC class I genes",
  "Lipid metabolism"
)
stopifnot(all(score_cols %in% names(df_b)))
# blue projection: light grey -> deep blue
BLUE_LO <- "#DEEBF7"
BLUE_HI <- "#08306B"
plist_b <- lapply(seq_along(score_cols), function(i) {
  feat_plot(df_b, score_cols[i], display_titles[i], italic_title = FALSE,
            low = BLUE_LO, high = BLUE_HI)
})
pB <- wrap_plots(plist_b, ncol = 2) +
  plot_annotation(title = "Microglia gene-set signatures on Stroke UMAP") &
  theme(plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0.5, family = FAM))

# -----------------------------------------------------------------------------
# Panel c — DimPlot n vs p (n1+n2, p1+p2)
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
df_c <- read.csv(file.path(SRC, "FigS7C_DimPlot_np.csv"), stringsAsFactors = FALSE)
df_c$group <- factor(df_c$group, levels = c("n", "p"))
df_c$cluster <- factor(as.character(df_c$cluster), levels = as.character(0:9))
lab_c <- df_c %>%
  group_by(group, cluster) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
# display subsample for speed while keeping all clusters
df_c_plot <- df_c %>%
  group_by(group, cluster) %>%
  group_modify(~ .x[sample.int(nrow(.x), min(nrow(.x), 3500)), , drop = FALSE]) %>%
  ungroup()

pC <- ggplot(df_c_plot, aes(umap_1, umap_2, colour = cluster)) +
  geom_point(size = 0.18, stroke = 0, shape = 16, alpha = 0.9) +
  geom_text(data = lab_c, aes(umap_1, umap_2, label = cluster),
            colour = "black", size = 2.3, fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = cluster_cols, guide = "none") +
  facet_wrap(~ group, nrow = 1) +
  coord_equal(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Stroke immune cells (n vs p)") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# -----------------------------------------------------------------------------
# Panel d — cluster proportions n vs p
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
prop <- read.csv(file.path(SRC, "FigS7D_cluster_proportion_np.csv"), stringsAsFactors = FALSE)
prop$cluster <- factor(as.character(prop$cluster), levels = as.character(0:9))
prop$group <- factor(prop$group, levels = c("n", "p"))

pD <- ggplot(prop, aes(group, proportion, fill = cluster)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion", title = "Cluster proportions (n vs p)") +
  theme(
    legend.position = "right",
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(2.4, "mm")
  )

# -----------------------------------------------------------------------------
# Panel e — KEGG (was c)
# -----------------------------------------------------------------------------
cat("Panel e ...\n")
kegg20 <- read.csv(file.path(SRC, "FigS7E_KEGG_signaling.csv"), stringsAsFactors = FALSE)
kegg20$Description <- factor(kegg20$Description, levels = kegg20$Description[order(kegg20$NES)])
pE <- ggplot(kegg20, aes(NES, Description)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(size = setSize, colour = neglog10_padj), shape = 16) +
  scale_colour_gradient(low = "#3B4CC0", high = "#B40426",
                        name = expression(-log[10](italic(q)))) +
  scale_size_area(max_size = 4.0, name = "setSize") +
  labs(x = "NES (cluster 5 vs 0)", y = NULL, title = "KEGG") +
  theme(
    axis.text.y = element_text(size = BASE - 0.9),
    legend.position = "right",
    plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0, family = FAM)
  )

# -----------------------------------------------------------------------------
# Panel f — WikiPathways (was d)
# -----------------------------------------------------------------------------
cat("Panel f ...\n")
wiki20 <- read.csv(file.path(SRC, "FigS7F_WikiPathways_top20.csv"), stringsAsFactors = FALSE)
wiki20$Description <- factor(wiki20$Description, levels = wiki20$Description[order(wiki20$NES)])
pF <- ggplot(wiki20, aes(NES, Description)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(size = setSize, colour = neglog10_padj), shape = 16) +
  scale_colour_gradient(low = "#3B4CC0", high = "#B40426",
                        name = expression(-log[10](italic(q)))) +
  scale_size_area(max_size = 4.0, name = "setSize") +
  labs(x = "NES (cluster 5 vs 0)", y = NULL, title = "WikiPathways") +
  theme(
    axis.text.y = element_text(size = BASE - 0.9),
    legend.position = "right",
    plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0, family = FAM)
  )

# -----------------------------------------------------------------------------
# Assemble a–f
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
fig <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  (wrap_elements(full = pC) | wrap_elements(full = pD)) /
  (pE | pF) +
  plot_layout(heights = c(1.05, 1.15, 1.05, 1.55)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.98)
  )

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS7")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig)
dev.off()

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM)
print(fig)
dev.off()

cat("Writing PNG ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig)
dev.off()

writeLines(c(
  "Fig. S7 | Lineage markers, microglia signatures, n/p composition, and pathway enrichment.",
  "",
  "a, Feature plots of canonical lineage markers on Stroke immune cells.",
  "b, Feature plots of AddModuleScore signatures (macrophage/microglia activation, ischemic-associated microglia, MHC genes, lipid metabolism) on the Stroke UMAP (blue scale).",
  "c, UMAP DimPlot with samples merged as n (n1+n2) versus p (p1+p2).",
  "d, Cluster proportions within n versus p.",
  "e, Top 18 KEGG GSEA terms ending in 'signaling pathway' for cluster 5 versus cluster 0. Dot size, setSize; colour, -log10(FDR q).",
  "f, Top 20 WikiPathways terms (signaling-pathway filter; 10 positive and 10 negative NES). Dot size, setSize; colour, -log10(FDR q).",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS7_legend.txt"))

writeLines(c(
  "Figure S7 | Nature-style native ggplot package",
  paste0("Generated: ", Sys.time()),
  paste0("Final size: ", FIG_W_MM, " x ", FIG_H_MM, " mm"),
  "Exports: FigureS7.pdf / .svg / .png (600 dpi)",
  "Backend: R 4.4.1 exclusive; native ggplot2 geom_point (no scattermore)",
  "",
  "Panel b colour: light blue-grey (#DEEBF7) -> deep blue (#08306B)",
  "Panel c/d: n1+n2 -> n; p1+p2 -> p (Dimplot.pdf logic)",
  "Panel c display subsample: <=3500 cells / group / cluster (labels from full medians)",
  "Former panels c/d renumbered to e/f (KEGG / WikiPathways)",
  "New folder FigS7_NC_np_blue/ — does not overwrite prior packages."
), file.path(OUTDIR, "FigureS7_QA.txt"))

file.copy(file.path(OUTDIR, "scripts", "make_figureS7.R"),
          file.path(OUTDIR, "make_figureS7.R"), overwrite = TRUE)

writeLines(c(
  "FigS7_NC_np_blue",
  "Added DimPlot n/p (c) and proportions (d); pathways moved to e/f;",
  "panel b module scores use blue colour scale."
), file.path(OUTDIR, "README.txt"))

cat("DONE ->", OUTDIR, "\n")
