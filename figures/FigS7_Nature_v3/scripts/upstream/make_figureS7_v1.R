#!/usr/bin/env Rscript
# =============================================================================
# Fig. S7 v1 — Nature SI multipanel (R / ggplot2 exclusive)
#
# Core conclusion
#   After marker-defined identity of CNS immune clusters, EAE_POS versus
#   EAE_NEG T and B cells shift KEGG signalling programmes, and DEGs that
#   belong to those displayed pathways mark the contrast; CellChat counts
#   redistribute across the same lineages.
#
# Evidence chain
#   a  Top marker heatmap (cluster identity)
#   b  Canonical lineage markers
#   c  T-cell volcano; labels = signalling genes from d
#   d  T-cell KEGG signalling GSEA
#   e  B-cell volcano; labels = signalling genes from f
#   f  B-cell KEGG signalling GSEA
#   g  CellChat interaction-number networks
#
# Archetype : asymmetric mixed-modality
# Size      : 183 x 510 mm
# New folder — does not overwrite any prior package
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
  library(ggrepel)
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
  file.path(IV_FIG_ROOT, "FigS7_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183  # fig_width_mm = 183
FIG_H_MM <- 510
FAM      <- "Arial"
BASE     <- 6.5
PT5      <- 5.2
LABEL_PT <- PT5 / ggplot2::.pt

CT_ORD <- c("Bcell", "Tcell", "Mac", "Neu1", "Neu2", "Mg", "PC", "DC", "UK")
CT_COL <- c(
  Bcell = "#E64B35", Tcell = "#4DBBD5", Mac = "#00A087", Neu1 = "#3C5488",
  Neu2 = "#F39B7F", Mg = "#8491B4", PC = "#91D1C2", DC = "#DC0000", UK = "#7E6148"
)
dir_cols <- c("Up in POS" = "#C45C4A", "Up in NEG" = "#3C5488", "NS" = "grey78")

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
# a  heatmap (NC-submission geometry: y-axis gene labels)
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
# c / e  volcanoes: all genes plotted; labels = paired KEGG core DEGs
# -----------------------------------------------------------------------------
make_volcano <- function(d, title, tag) {
  d$direction <- factor(as.character(d$direction),
                        levels = c("Up in POS", "Up in NEG", "NS"))
  lab <- d[d$label_me, , drop = FALSE]
  ymax <- max(d$neglog10_padj, na.rm = TRUE)
  xmax <- max(abs(d$log2FC), na.rm = TRUE)
  ggplot(d, aes(log2FC, neglog10_padj)) +
    geom_point(aes(colour = direction), data = d[!d$label_me, ],
               size = 0.42, stroke = 0, alpha = 0.55) +
    geom_point(aes(colour = direction), data = lab,
               size = 1.25, stroke = 0, alpha = 1) +
    geom_vline(xintercept = c(-0.5, 0.5), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.001), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_text_repel(
      data = lab, aes(label = gene),
      size = LABEL_PT, colour = "black", family = FAM, fontface = "italic",
      max.overlaps = Inf, box.padding = 0.22, point.padding = 0.14,
      min.segment.length = 0, segment.colour = "grey40",
      segment.size = 0.22, seed = 42, show.legend = FALSE,
      nudge_x = ifelse(lab$log2FC >= 0, 0.42, -0.42),
      direction = "y", max.iter = 20000, max.time = 2,
      ylim = c(6, NA)
    ) +
    scale_colour_manual(
      values = dir_cols, name = NULL, drop = FALSE,
      guide = guide_legend(nrow = 1, override.aes = list(size = 2.2, alpha = 1))
    ) +
    coord_cartesian(xlim = c(-xmax * 1.42, xmax * 1.32),
                    ylim = c(0, ymax * 1.08)) +
    labs(x = "log2 FC (POS / NEG)", y = "-log10 adjusted P",
         title = title,
         subtitle = sprintf(
           "Labels: KEGG signalling DEGs from the paired panel (n = %d)",
           nrow(lab)),
         tag = tag) +
    theme_pub +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-2, 0, 0, 0),
      legend.key.size = unit(2.2, "mm"),
      plot.subtitle = element_text(size = 5.2, lineheight = 1.05),
      plot.margin = margin(2, 3, 4, 2)
    )
}

cat("Panels c and e ...\n")
vol_t <- read.csv(file.path(SRC, "FigS7c_volcano_Tcell.csv"), stringsAsFactors = FALSE)
vol_b <- read.csv(file.path(SRC, "FigS7e_volcano_Bcell.csv"), stringsAsFactors = FALSE)
stopifnot("label_me" %in% names(vol_t), "label_me" %in% names(vol_b))
pC <- make_volcano(vol_t, "T cell POS versus NEG", "c")
pE <- make_volcano(vol_b, "B cell POS versus NEG", "e") +
  theme(legend.position = "none")
cat("  T labelled", sum(vol_t$label_me), " B labelled", sum(vol_b$label_me), "\n")

# -----------------------------------------------------------------------------
# d / f  KEGG signalling
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

cat("Panels d and f ...\n")
pD <- make_kegg(file.path(SRC, "FigS4_KEGG_Tcell.csv"), "T cell KEGG signalling", "d")
pF <- make_kegg(file.path(SRC, "FigS4_KEGG_Bcell.csv"), "B cell KEGG signalling", "f")

# -----------------------------------------------------------------------------
# g  CellChat interaction number
# -----------------------------------------------------------------------------
cat("Panel g ...\n")
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

pG_pos <- circle_ggplot(net_pos, sprintf("EAE_POS (n = %d)", n_int[["EAE_POS"]]))
pG_neg <- circle_ggplot(net_neg, sprintf("EAE_NEG (n = %d)", n_int[["EAE_NEG"]]))
pG <- (pG_pos | pG_neg) +
  plot_annotation(
    title = "CellChat interaction number",
    subtitle = "Edge width = interaction count (shared scale); node size = total degree",
    theme = theme(
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
      plot.subtitle = element_text(size = 5.2, hjust = 0.5, colour = "grey30",
                                   family = FAM)
    )
  )
pG <- wrap_elements(full = pG) + labs(tag = "g") + tag_theme

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
seal <- function(x) wrap_elements(full = x)
row1 <- seal(pA)
row2 <- seal(pB)
row3 <- (pC | pD) + plot_layout(widths = c(1.22, 1))
row4 <- (pE | pF) + plot_layout(widths = c(1.22, 1))
row5 <- pG

fig <- row1 / row2 / row3 / row4 / row5 +
  plot_layout(heights = grid::unit.c(
    grid::unit(175, "mm"),
    grid::unit(52, "mm"),
    grid::unit(100, "mm"),
    grid::unit(88, "mm"),
    grid::unit(95, "mm")
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
write.csv(vol_t, file.path(PLOTTED, "SourceData_FigS7c_volcano_Tcell.csv"),
          row.names = FALSE)
write.csv(vol_t[vol_t$label_me, ],
          file.path(PLOTTED, "SourceData_FigS7c_labelled.csv"), row.names = FALSE)
write.csv(read.csv(file.path(SRC, "FigS4_KEGG_Tcell.csv"), stringsAsFactors = FALSE),
          file.path(PLOTTED, "SourceData_FigS7d_KEGG_Tcell.csv"), row.names = FALSE)
write.csv(vol_b, file.path(PLOTTED, "SourceData_FigS7e_volcano_Bcell.csv"),
          row.names = FALSE)
write.csv(vol_b[vol_b$label_me, ],
          file.path(PLOTTED, "SourceData_FigS7e_labelled.csv"), row.names = FALSE)
write.csv(read.csv(file.path(SRC, "FigS4_KEGG_Bcell.csv"), stringsAsFactors = FALSE),
          file.path(PLOTTED, "SourceData_FigS7f_KEGG_Bcell.csv"), row.names = FALSE)
write.csv(int_long, file.path(PLOTTED, "SourceData_FigS7g_interactions.csv"),
          row.names = FALSE)

k_t <- read.csv(file.path(SRC, "FigS4_KEGG_Tcell.csv"), stringsAsFactors = FALSE)
k_b <- read.csv(file.path(SRC, "FigS4_KEGG_Bcell.csv"), stringsAsFactors = FALSE)
cat("T KEGG FDR<0.05:", sum(k_t$fdr_sig), " B:", sum(k_b$fdr_sig), "\n")
cat("g interactions POS/NEG:", n_int[["EAE_POS"]], n_int[["EAE_NEG"]], "\n")
cat("DONE ->", OUTDIR, "\n")
