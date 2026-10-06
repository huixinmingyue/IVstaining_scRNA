#!/usr/bin/env Rscript
# =============================================================================
# Fig. S7 v3 — Nature SI multipanel (R / ggplot2 exclusive)
#
# Core conclusion
#   After marker-defined identity of CNS immune clusters, CellChat
#   counts redistribute between POS and NEG; T/B KEGG signalling shifts
#   (T exploratory; B IL-17/TNF FDR-significant); locked T/B genes mark
#   the unfiltered POS versus NEG contrast.
#
# Evidence chain
#   a  Top marker heatmap
#   b  Canonical lineage markers
#   c  CellChat interaction-number networks  (was v2 g)
#   d  T-cell KEGG signalling GSEA
#   e  B-cell KEGG signalling GSEA
#   f  T-cell locked genes, all cells only
#   g  B-cell locked genes, all cells only
#
# Archetype : asymmetric mixed-modality
# Size      : 183 x 473 mm
# New folder — does not overwrite FigS7_Nature_v1 or v2
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS7_Nature_v3")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183  # fig_width_mm = 183
FIG_H_MM <- 473
FAM      <- "Arial"
BASE     <- 6.5
PT5      <- 5.2
LABEL_PT <- PT5 / ggplot2::.pt

CT_ORD <- c("Bcell", "Tcell", "Mac", "Neu1", "Neu2", "Mg", "PC", "DC", "UK")
CT_COL <- c(
  Bcell = "#E64B35", Tcell = "#4DBBD5", Mac = "#00A087", Neu1 = "#3C5488",
  Neu2 = "#F39B7F", Mg = "#8491B4", PC = "#91D1C2", DC = "#DC0000", UK = "#7E6148"
)
dir_cols <- c("Up in POS" = "#C45C4A", "Up in NEG" = "#3C5488")

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

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title = element_text(size = BASE - 0.2, colour = "black"),
    axis.text = element_text(size = BASE - 0.6, colour = "black"),
    legend.title = element_text(size = BASE - 0.4),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5, colour = "grey30"),
    plot.tag = element_text(size = 8, face = "bold", family = FAM,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    plot.margin = margin(1.5, 2, 1, 1),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = BASE - 0.6, face = "bold", family = FAM)
  )
theme_set(theme_pub)
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = FAM))

kegg_display <- function(x) {
  x <- gsub(" - Mus musculus.*|house mouse", "", x)
  x <- gsub("[()]", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  x <- sub("\\s+signaling pathway\\b", " signaling", x, ignore.case = TRUE)
  x <- gsub("TGF-beta", "TGF-\u03b2", x)
  x <- gsub("NF-kappa B", "NF-\u03baB", x)
  x
}

# -----------------------------------------------------------------------------
# a  heatmap
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
hm <- readRDS(file.path(SRC, "FigS7a_heatmap.rds"))
hm_long <- hm$hm_long
gene_order <- hm$gene_order
zlim <- hm$zlim
cl_lab <- setNames(sprintf("%s\n%s", as.character(0:8), CT_ORD), as.character(0:8))
hm_long$gene <- factor(hm_long$gene, levels = rev(gene_order))
hm_long$cluster <- factor(hm_long$cluster, levels = as.character(0:8))
hm_long$cluster_lab <- factor(cl_lab[as.character(hm_long$cluster)], levels = cl_lab)
stopifnot(length(gene_order) == 88L)

pA <- ggplot(hm_long, aes(cell, gene, fill = z)) +
  geom_tile(width = 1, height = 1) +
  facet_grid(~ cluster_lab, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-zlim, zlim), oob = squish,
    breaks = c(-zlim, 0, zlim),
    labels = c(sprintf("%.1f", -zlim), "0", sprintf("%.1f", zlim)),
    name = "Scaled\nexpr.",
    guide = guide_colorbar(barwidth = unit(2.4, "mm"), barheight = unit(28, "mm"),
                           ticks.colour = "black", frame.colour = "black",
                           frame.linewidth = 0.2)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Top marker genes by cluster",
       subtitle = sprintf(
         "Display downsample %d cells / cluster (seed 42); full n = %s",
         hm$max_per_cl, format(hm$n_cells, big.mark = ",")),
       tag = "a") +
  theme_pub +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(size = PT5, face = "italic", colour = "black",
                               family = FAM),
    strip.text = element_text(size = 5.4, face = "bold", lineheight = 0.9),
    strip.background = element_rect(fill = "grey95", colour = "grey80",
                                    linewidth = 0.3),
    panel.spacing.x = unit(0.45, "mm"),
    panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.25),
    legend.position = "right",
    plot.subtitle = element_text(size = 5.2)
  )

# -----------------------------------------------------------------------------
# b  dot plot
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
dot <- read.csv(file.path(SRC, "FigS4b_dotplot.csv"), stringsAsFactors = FALSE)
feature8 <- c("Cd79a", "Ighd", "Trbc2", "Cd3g", "Plac8", "Ccr2", "Csf3r", "Cxcr2",
              "S100a8", "S100a9", "Hexb", "P2ry12", "Jchain", "Ighg1", "Runx2", "Irf8")
dot$gene <- factor(dot$gene, levels = feature8)
ylab_levels <- sprintf("%s (%s)", as.character(0:8), CT_ORD)
dot$ylab <- factor(sprintf("%s (%s)", dot$cluster, dot$celltype), levels = ylab_levels)
zlim_d <- 2.5
dot$avg_plot <- pmax(pmin(dot$avg_scaled, zlim_d), -zlim_d)

pB <- ggplot(dot, aes(gene, ylab)) +
  geom_point(aes(size = pct, colour = avg_plot), shape = 16) +
  scale_size_area(max_size = 4.0, breaks = c(10, 30, 50, 75),
                  name = "% expressing", limits = c(0, 100)) +
  scale_colour_gradientn(
    colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(-zlim_d, zlim_d), name = "Scaled mean",
    guide = guide_colorbar(barwidth = unit(2.4, "mm"), barheight = unit(16, "mm"),
                           ticks.colour = "black", frame.colour = "black",
                           frame.linewidth = 0.2)
  ) +
  labs(x = NULL, y = NULL, title = "Canonical lineage markers", tag = "b") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", size = PT5),
    axis.text.y = element_text(size = PT5),
    legend.position = "right",
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line = element_blank()
  )

# -----------------------------------------------------------------------------
# c  CellChat (was v2 g)
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
int_long <- read.csv(file.path(SRC, "FigS4g_interaction_counts.csv"),
                     stringsAsFactors = FALSE)
to_mat <- function(g) {
  sub <- int_long[int_long$group == g, ]
  m <- matrix(0, length(CT_ORD), length(CT_ORD), dimnames = list(CT_ORD, CT_ORD))
  for (i in seq_len(nrow(sub))) {
    if (sub$source[i] %in% CT_ORD && sub$target[i] %in% CT_ORD) {
      m[sub$source[i], sub$target[i]] <- sub$count[i]
    }
  }
  m
}
net_pos <- to_mat("EAE_POS")
net_neg <- to_mat("EAE_NEG")
emax <- max(c(as.numeric(net_pos), as.numeric(net_neg)))
n_int <- c(EAE_POS = sum(net_pos), EAE_NEG = sum(net_neg))

circle_ggplot <- function(net, title) {
  labsn <- rownames(net)
  n <- length(labsn)
  theta <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)] + pi / 2
  nodes <- data.frame(
    celltype = factor(labsn, levels = CT_ORD),
    x = cos(theta), y = sin(theta),
    deg = as.numeric(rowSums(net) + colSums(net)),
    stringsAsFactors = FALSE
  )
  rlab <- 1.36
  nodes$lx <- rlab * cos(theta)
  nodes$ly <- rlab * sin(theta)
  edges <- which(net != 0, arr.ind = TRUE)
  ed <- data.frame(
    from = labsn[edges[, 1]], to = labsn[edges[, 2]],
    w = as.numeric(net[edges]), stringsAsFactors = FALSE
  )
  ed <- ed %>%
    left_join(nodes %>% select(celltype, x, y), by = c("from" = "celltype")) %>%
    left_join(nodes %>% select(celltype, xend = x, yend = y),
              by = c("to" = "celltype"))
  ed$col <- unname(CT_COL[ed$from])
  ed_ns <- ed[ed$from != ed$to, , drop = FALSE]
  ggplot() +
    geom_curve(
      data = ed_ns,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = w),
      colour = ed_ns$col, alpha = 0.55, curvature = 0.25,
      arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
      lineend = "round"
    ) +
    scale_linewidth_continuous(range = c(0.2, 2.2), limits = c(0, emax),
                               guide = "none") +
    geom_point(data = nodes, aes(x, y, fill = celltype, size = deg),
               shape = 21, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = CT_COL, guide = "none") +
    scale_size_area(max_size = 10, guide = "none") +
    geom_text(data = nodes, aes(lx, ly, label = celltype),
              size = LABEL_PT, family = FAM, colour = "black") +
    coord_equal(xlim = c(-1.68, 1.68), ylim = c(-1.68, 1.68), expand = FALSE) +
    labs(title = title) +
    theme_void(base_family = FAM) +
    theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0.5,
                                    family = FAM),
          plot.margin = margin(1, 1, 1, 1, "mm"))
}

pC_pos <- circle_ggplot(net_pos, sprintf("EAE_POS (n = %d)", n_int[["EAE_POS"]]))
pC_neg <- circle_ggplot(net_neg, sprintf("EAE_NEG (n = %d)", n_int[["EAE_NEG"]]))
pC <- (pC_pos | pC_neg) +
  plot_annotation(
    title = "CellChat interaction number",
    subtitle = "Edge width = interaction count (shared scale); node size = total degree",
    theme = theme(
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
      plot.subtitle = element_text(size = 5.2, hjust = 0.5, colour = "grey30",
                                   family = FAM)
    )
  )
pC <- wrap_elements(full = pC) + labs(tag = "c") + tag_theme

# -----------------------------------------------------------------------------
# d / e  KEGG signalling
# -----------------------------------------------------------------------------
make_kegg <- function(path, title, tag) {
  k <- read.csv(path, stringsAsFactors = FALSE)
  k$lab <- kegg_display(as.character(k$Description))
  k$lab <- factor(k$lab, levels = k$lab[order(k$NES)])
  k$neglog10_q <- -log10(pmax(k$p.adjust, 1e-300))
  n_fdr <- sum(k$fdr_sig)
  sub <- if (n_fdr > 0) {
    sprintf("Black outline: FDR < 0.05 (n = %d / 10)", n_fdr)
  } else {
    "None reach FDR < 0.05 (exploratory)"
  }
  ggplot(k, aes(NES, lab)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(aes(size = setSize, fill = neglog10_q),
               shape = 21, colour = "grey55", stroke = 0.35) +
    geom_point(data = k[k$fdr_sig, , drop = FALSE],
               aes(NES, lab, size = setSize),
               shape = 21, fill = NA, colour = "black", stroke = 0.7,
               show.legend = FALSE) +
    scale_fill_gradient(
      low = "#C6DBEF", high = "#08519C", name = "-log10 FDR",
      guide = guide_colorbar(barwidth = unit(2.4, "mm"), barheight = unit(16, "mm"),
                             ticks.colour = "black", frame.colour = "black",
                             frame.linewidth = 0.2, order = 1)
    ) +
    scale_size_area(max_size = 4.4, name = "Set size",
                    guide = guide_legend(order = 2)) +
    labs(x = "NES", y = NULL, title = title, subtitle = sub, tag = tag) +
    theme_pub +
    theme(
      axis.text.y = element_text(size = PT5),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(1.5, "mm"),
      legend.key.size = unit(2.2, "mm"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      axis.line = element_blank(),
      plot.subtitle = element_text(size = 5.2),
      plot.margin = margin(2, 4, 2, 2)
    )
}

cat("Panels d and e ...\n")
pD <- make_kegg(file.path(SRC, "FigS4_KEGG_Tcell.csv"), "T cell KEGG signalling", "d")
pE <- make_kegg(file.path(SRC, "FigS4_KEGG_Bcell.csv"), "B cell KEGG signalling", "e")

# -----------------------------------------------------------------------------
# f / g  locked genes, all cells only
# -----------------------------------------------------------------------------
cat("Panels f and g ...\n")
lock <- read.csv(file.path(SRC, "locked_genes_cellfilter.csv"),
                 stringsAsFactors = FALSE)
dropn <- read.csv(file.path(SRC, "deg_drop_counts.csv"), stringsAsFactors = FALSE)
n_all <- function(ct, grp) {
  sum(dropn$n[dropn$celltype == ct & dropn$IV_label == grp])
}

T_ORD <- c("Pdcd1", "Cxcr6", "Ctla4", "Ifng", "Ccr7", "Il7r",
           "Lef1", "Sell", "S1pr1", "Klf2", "S100a9", "S100a8")
B_ORD <- c("Retnlg", "Il1b", "Cd72", "Gpr183", "Ciita", "Cr2",
           "Cxcl2", "S100a8", "S100a9")

make_locked <- function(ct, gene_ord, title, tag, show_legend = TRUE) {
  d <- lock[lock$celltype == ct, , drop = FALSE]
  stopifnot(nrow(d) == length(gene_ord), setequal(d$gene, gene_ord))
  d$gene <- factor(d$gene, levels = rev(gene_ord))
  d$direction <- ifelse(d$all_log2FC >= 0, "Up in POS", "Up in NEG")
  d$direction <- factor(d$direction, levels = c("Up in POS", "Up in NEG"))
  npos <- n_all(ct, "EAE_POS")
  nneg <- n_all(ct, "EAE_NEG")
  cat("  ", ct, "all cells POS/NEG", npos, nneg, "\n")
  ggplot(d, aes(all_log2FC, gene, colour = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 1.7, stroke = 0) +
    scale_colour_manual(values = dir_cols, name = NULL, drop = FALSE) +
    labs(x = "log2 FC (POS / NEG)", y = NULL, title = title, tag = tag) +
    theme_pub +
    theme(
      axis.text.y = element_text(size = PT5, face = "italic", family = FAM),
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      axis.line = element_blank(),
      plot.margin = margin(2, 3, 3, 2)
    )
}

pF <- make_locked("Tcell", T_ORD, "T cell locked genes, all cells", "f", TRUE)
pG <- make_locked("Bcell", B_ORD, "B cell locked genes, all cells", "g", FALSE)

# -----------------------------------------------------------------------------
# Assemble: a / b / c CellChat / d|e KEGG / f|g locked
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
seal <- function(x) wrap_elements(full = x)
fig <- seal(pA) / seal(pB) / pC /
  ((pD | pE) + plot_layout(widths = c(1, 1))) /
  ((pF | pG) + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = grid::unit.c(
    grid::unit(175, "mm"),
    grid::unit(52, "mm"),
    grid::unit(90, "mm"),
    grid::unit(68, "mm"),
    grid::unit(88, "mm")
  ))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "FigureS7")

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

write.csv(hm_long, file.path(PLOTTED, "SourceData_FigS7a_heatmap.csv"),
          row.names = FALSE)
write.csv(dot, file.path(PLOTTED, "SourceData_FigS7b_dotplot.csv"),
          row.names = FALSE)
write.csv(int_long, file.path(PLOTTED, "SourceData_FigS7c_interactions.csv"),
          row.names = FALSE)
write.csv(read.csv(file.path(SRC, "FigS4_KEGG_Tcell.csv"), stringsAsFactors = FALSE),
          file.path(PLOTTED, "SourceData_FigS7d_KEGG_Tcell.csv"), row.names = FALSE)
write.csv(read.csv(file.path(SRC, "FigS4_KEGG_Bcell.csv"), stringsAsFactors = FALSE),
          file.path(PLOTTED, "SourceData_FigS7e_KEGG_Bcell.csv"), row.names = FALSE)
lock_tb <- lock[lock$celltype %in% c("Tcell", "Bcell"),
                c("celltype", "gene", "all_log2FC", "all_padj", "all_sig")]
write.csv(lock_tb, file.path(PLOTTED, "SourceData_FigS7fg_locked_allcells.csv"),
          row.names = FALSE)

k_t <- read.csv(file.path(SRC, "FigS4_KEGG_Tcell.csv"), stringsAsFactors = FALSE)
k_b <- read.csv(file.path(SRC, "FigS4_KEGG_Bcell.csv"), stringsAsFactors = FALSE)
cat("T KEGG FDR<0.05:", sum(k_t$fdr_sig), " B:", sum(k_b$fdr_sig), "\n")
cat("c interactions POS/NEG:", n_int[["EAE_POS"]], n_int[["EAE_NEG"]], "\n")
cat("DONE ->", OUTDIR, "\n")
