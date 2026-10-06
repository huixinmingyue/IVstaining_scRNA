#!/usr/bin/env Rscript
# =============================================================================
# Fig. S1 — QC, marker identities, B / true-CD4 / true-CD8 DEG, plasma isotype
# New folder only. Does not overwrite FigS1_Nature_v2 or other S1 packages.
#
# Core conclusion
#   IV compartments have comparable library quality; canonical markers support
#   the RNA-rule identities in Fig. 1; B-cell Intra DEGs overlap MZ versus Fo;
#   residual Intra DEGs in hard-gated true CD4 and true CD8 are restricted;
#   plasma IV+ versus IV- follows IgM versus IgG, not a GC programme.
#
# Archetype : quantitative grid
# Backend   : R exclusive
# Size      : 183 x 270 mm
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
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
ROOT <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  normalizePath(getwd())
}

SRC <- file.path(ROOT, "source_data")
SDP <- file.path(ROOT, "source_data_plotted")
QA  <- file.path(ROOT, "qa")
dir.create(SDP, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 270
width_mm = 183
height_mm = 270
BASE <- 6.5
FAM  <- "Arial"
fc_cut <- 0.5
padj_cut <- 0.001
zlim <- 2.5

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) cat("font:", FAM, "->", fm$path, "\n")

IV_COL <- c("IV\u2212" = "#B8C0C8", "IV+" = "#3D4F5F")
DIR_COL <- c("IntraNeg up" = "#2F6F8F", "NS" = "#B0B0B0", "IntraPos up" = "#B2473E")

iv_short <- function(x) {
  ifelse(grepl("intravascular", x) & !grepl("extravascular", x),
         "IV+", "IV\u2212")
}

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = BASE),
      axis.text = element_text(size = BASE - 0.5, colour = "black"),
      legend.title = element_text(size = BASE - 0.3),
      legend.text = element_text(size = BASE - 0.7),
      legend.key.size = unit(2.6, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold"),
      plot.title = element_text(size = BASE + 0.4, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.margin = margin(3, 3, 2, 5),
      panel.grid = element_blank()
    )
)

save_pub <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot); invisible(dev.off())
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h,
                       family = FAM, bg = "white")
  print(plot); invisible(dev.off())
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in",
                res = dpi, background = "white")
  print(plot); invisible(dev.off())
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in",
                 res = dpi, compression = "lzw", background = "white")
  print(plot); invisible(dev.off())
}

# =============================================================================
# a  QC by compartment
# =============================================================================
cells <- read.csv(file.path(SRC, "cell_audit_table.csv"), stringsAsFactors = FALSE)
n_before_a <- nrow(cells)
qc <- cells %>%
  transmute(
    IV = factor(iv_short(IV_lab), levels = c("IV\u2212", "IV+")),
    `nCount RNA` = nCount_RNA,
    `nFeature RNA` = nFeature_RNA,
    `Mitochondrial %` = percent.mt,
    `Stress module` = stress
  )
n_after_a <- nrow(qc)
stopifnot(n_before_a == n_after_a, n_after_a == 15271L)

qc_long <- qc %>%
  pivot_longer(-IV, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("nCount RNA", "nFeature RNA",
                                            "Mitochondrial %", "Stress module")))
write.csv(qc_long, file.path(SDP, "SourceData_FigS1a_QC.csv"), row.names = FALSE)

pA <- ggplot(qc_long, aes(IV, value, fill = IV)) +
  geom_violin(linewidth = 0.25, colour = "black", scale = "width") +
  geom_boxplot(width = 0.18, outlier.size = 0.2, linewidth = 0.25, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = IV_COL, name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none")

# =============================================================================
# b  Marker dot plot
# =============================================================================
ord <- read.csv(file.path(SRC, "lineage_order.csv"), stringsAsFactors = FALSE)
ct_order <- as.character(ord$lineage)
gene_blocks <- list(
  "CD4 T"       = c("Il7r", "Cd4"),
  "CD8 T"       = c("Cd8b1", "Cd8a", "Cd27"),
  "NKT-like"    = c("Klrb1c"),
  "T other"     = c("Cd3e"),
  "B cell"      = c("Cd19", "Ms4a1", "Cd79a"),
  "cDC2"        = c("Klf4", "Clec9a", "Cd68", "Itgax"),
  "NK"          = c("Ncr1", "Cd160", "Ccl3", "Gzmb", "Nkg7"),
  "Plasma cell" = c("Jchain", "Mzb1"),
  "Neutrophil"  = c("S100a8", "S100a9")
)
block_names <- ct_order[ct_order %in% names(gene_blocks)]
gene_blocks <- gene_blocks[block_names]
gene_order  <- unlist(gene_blocks, use.names = FALSE)

dot <- read.csv(file.path(SRC, "FigS1b_markers_dotplot.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(gene_order, unique(dot$gene)))
dot$avg_plot <- pmax(pmin(dot$avg_scaled, zlim), -zlim)
dot$gene     <- factor(dot$gene, levels = gene_order)
dot$celltype <- factor(dot$celltype, levels = rev(ct_order))
dot <- droplevels(dot[!is.na(dot$gene) & !is.na(dot$celltype), ])
block_edges <- head(cumsum(vapply(gene_blocks, length, integer(1))), -1) + 0.5
write.csv(dot[, c("celltype", "gene", "avg", "avg_scaled", "avg_plot", "pct")],
          file.path(SDP, "SourceData_FigS1b_markers.csv"), row.names = FALSE)

pB <- ggplot(dot, aes(gene, celltype)) +
  geom_vline(xintercept = block_edges, colour = "grey88", linewidth = 0.25) +
  geom_point(aes(size = pct, colour = avg_plot), shape = 16) +
  scale_size_area(max_size = 3.6, breaks = c(10, 30, 50, 75),
                  limits = c(0, 100), name = "Expressing\ncells (%)") +
  scale_colour_gradientn(
    colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(-zlim, zlim), breaks = c(-2, 0, 2),
    name = "Scaled mean\nexpression",
    guide = guide_colorbar(barwidth = unit(2.2, "mm"), barheight = unit(16, "mm"),
                           ticks.colour = "black", frame.colour = "black",
                           frame.linewidth = 0.2, order = 1)
  ) +
  guides(size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", size = BASE - 0.8),
    axis.text.y = element_text(size = BASE - 0.4, colour = "black"),
    legend.position = "right",
    legend.box.spacing = unit(1.2, "mm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.25)
  )

# =============================================================================
# c–e  Volcanoes
# =============================================================================
prep_volc <- function(d) {
  if (!"neglog10_padj" %in% names(d)) {
    d$neglog10_padj <- pmin(-log10(pmax(d$p_val_adj, 1e-300)), 50)
  }
  d$direction <- factor(d$direction, levels = c("IntraNeg up", "NS", "IntraPos up"))
  d
}

volcano_sidestack <- function(d, title, lab_genes, show_legend = FALSE) {
  d <- prep_volc(d)
  lab <- d %>% filter(gene %in% lab_genes) %>% arrange(p_val_adj)
  n_left <- sum(as.character(lab$direction) == "IntraNeg up")
  n_right <- sum(as.character(lab$direction) == "IntraPos up")
  n_side_max <- max(n_left, n_right, 1)
  ncols_guess <- if (n_side_max > 70) 3L else if (n_side_max > 35) 2L else 1L
  xmax_data <- max(3, quantile(abs(d$avg_log2FC), 0.999, na.rm = TRUE) * 1.05, fc_cut * 1.2)
  xmax <- xmax_data * (1.18 + 0.62 * ncols_guess)
  ymax <- max(5, max(d$neglog10_padj, na.rm = TRUE) * 1.05)
  lab_size <- 5.3 / .pt

  p <- ggplot(d, aes(avg_log2FC, neglog10_padj, colour = direction)) +
    geom_point(size = 0.4, stroke = 0, alpha = 0.7) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(padj_cut), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    scale_colour_manual(
      values = DIR_COL, name = NULL, drop = FALSE,
      labels = c("IntraNeg up" = "IV\u2212 up", "NS" = "NS",
                 "IntraPos up" = "IV+ up")
    ) +
    coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), clip = "off") +
    labs(x = "log2FC (IV+ / IV-)",
         y = "-log10 adjusted P",
         title = title) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.margin = margin(2, 4, 2, 2)
    )
  if (show_legend) {
    p <- p + guides(colour = guide_legend(
      override.aes = list(size = 2, alpha = 1), nrow = 1
    ))
  }

  stack_side <- function(dat, side = c("left", "right")) {
    side <- match.arg(side)
    n <- nrow(dat)
    if (!n) return(dat[0, ])
    dat <- dat %>% arrange(desc(neglog10_padj), desc(abs(avg_log2FC)))
    ncols <- if (n > 70) 3L else if (n > 35) 2L else 1L
    n_row <- ceiling(n / ncols)
    dat$col <- rep(seq_len(ncols), each = n_row)[seq_len(n)]
    dat$row <- ave(seq_len(n), dat$col, FUN = seq_along)
    y_top <- ymax * 0.985
    y_bot <- ymax * 0.015
    gap <- (xmax - xmax_data * 1.05) / max(ncols, 1)
    if (side == "left") {
      dat$lab_x <- -xmax + 0.06 * xmax_data + (dat$col - 1) * gap
    } else {
      dat$lab_x <- xmax - 0.06 * xmax_data - (dat$col - 1) * gap
    }
    dat$lab_y <- if (n_row == 1) (y_top + y_bot) / 2 else
      y_top - (dat$row - 1) * (y_top - y_bot) / (n_row - 1)
    dat
  }

  lab_L <- stack_side(lab %>% filter(as.character(direction) == "IntraNeg up"), "left")
  lab_R <- stack_side(lab %>% filter(as.character(direction) == "IntraPos up"), "right")
  if (nrow(lab_L)) {
    p <- p + geom_text(
      data = lab_L, aes(lab_x, lab_y, label = gene, colour = direction),
      size = lab_size, fontface = "italic", family = FAM,
      hjust = 0, vjust = 0.5, show.legend = FALSE, inherit.aes = FALSE
    )
  }
  if (nrow(lab_R)) {
    p <- p + geom_text(
      data = lab_R, aes(lab_x, lab_y, label = gene, colour = direction),
      size = lab_size, fontface = "italic", family = FAM,
      hjust = 1, vjust = 0.5, show.legend = FALSE, inherit.aes = FALSE
    )
  }
  p
}

pick_lab <- function(d, n_lab) {
  d <- prep_volc(d)
  d %>%
    filter(direction != "NS", !tech_or_ribo) %>%
    group_by(direction) %>%
    slice_min(order_by = p_val_adj, n = n_lab, with_ties = FALSE) %>%
    ungroup()
}

deg_b <- read.csv(file.path(SRC, "DEG_Bcell_IntraPos_vs_IntraNeg.csv"),
                  stringsAsFactors = FALSE)
hit_c <- read.csv(file.path(SRC, "S1c_MZvsFo_x_BcellDEG.csv"),
                  stringsAsFactors = FALSE)
lab_c_genes <- unique(as.character(hit_c$gene))
cat("Panel c labelled genes:", length(lab_c_genes), "\n")

deg_cd4 <- read.csv(file.path(SRC, "DEG_true_CD4_IntraPos_vs_IntraNeg.csv"),
                    stringsAsFactors = FALSE)
deg_cd8 <- read.csv(file.path(SRC, "DEG_true_CD8_IntraPos_vs_IntraNeg.csv"),
                    stringsAsFactors = FALSE)
lab_d <- pick_lab(deg_cd4, 12)
lab_e <- pick_lab(deg_cd8, 10)

pC <- volcano_sidestack(deg_b, "B cell", lab_c_genes)
pD <- volcano_sidestack(deg_cd4, "True CD4", lab_d$gene)
pE <- volcano_sidestack(deg_cd8, "True CD8", lab_e$gene)
leg_volc <- cowplot::get_legend(
  volcano_sidestack(deg_cd8, "True CD8", lab_e$gene, show_legend = TRUE)
)

lab_c_out <- hit_c %>%
  transmute(panel = "c", celltype = "B cell", gene, direction,
            avg_log2FC, p_val_adj, MZ_Fo_direction)
write.csv(lab_c_out, file.path(SDP, "SourceData_FigS1c_labelled.csv"), row.names = FALSE)

lab_de <- bind_rows(
  lab_d %>% mutate(panel = "d", celltype = "True CD4"),
  lab_e %>% mutate(panel = "e", celltype = "True CD8")
) %>%
  select(panel, celltype, gene, direction, avg_log2FC, p_val_adj, pct.1, pct.2)
write.csv(lab_de, file.path(SDP, "SourceData_FigS1de_labelled.csv"), row.names = FALSE)

# =============================================================================
# f  Plasma genes; Aicda dropped
# =============================================================================
pc <- read.csv(file.path(SRC, "plasma_gene_expr.csv"), stringsAsFactors = FALSE)
plasma_genes_before <- c("Ighm", "Ighg1", "Aicda", "Mki67", "Sdc1", "Cxcr4")
plasma_genes_after  <- c("Ighm", "Ighg1", "Mki67", "Sdc1", "Cxcr4")
n_genes_before <- length(unique(pc$gene))
n_cells_f <- nrow(pc) / n_genes_before
pc_plot <- pc %>%
  filter(gene != "Aicda") %>%
  mutate(
    IV = factor(iv_short(IV_lab), levels = c("IV\u2212", "IV+")),
    gene = factor(gene, levels = plasma_genes_after)
  )
n_genes_after <- nlevels(droplevels(pc_plot$gene))
stopifnot(n_genes_before == 6L, n_genes_after == 5L)
cat(sprintf("Panel f Aicda exclusion: n_genes_before=%d n_genes_after=%d n_cells=%.0f\n",
            n_genes_before, n_genes_after, n_cells_f))
write.csv(pc_plot[, c("IV", "gene", "expr")],
          file.path(SDP, "SourceData_FigS1f_plasma.csv"), row.names = FALSE)

pF <- ggplot(pc_plot, aes(IV, expr, fill = IV)) +
  geom_violin(scale = "width", linewidth = 0.25, colour = "black") +
  geom_boxplot(width = 0.15, outlier.size = 0.25, linewidth = 0.25, fill = "white") +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = IV_COL, name = NULL) +
  labs(x = NULL, y = "Log-normalised expression") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = BASE, face = "italic")
  )

pA <- pA + labs(tag = "a")
pB <- pB + labs(tag = "b")
pC <- pC + labs(tag = "c")
pD <- pD + labs(tag = "d")
pE <- pE + labs(tag = "e")
pF <- pF + labs(tag = "f")

fig <- pA / pB / (pC | pD | pE) / wrap_elements(full = leg_volc) / pF +
  plot_layout(heights = c(0.80, 1.18, 1.38, 0.11, 0.80))

base <- file.path(ROOT, "FigureS1")
cat("Writing FigureS1 ...\n")
save_pub(fig, base, FIG_W_MM, FIG_H_MM)

cat("DONE ->", ROOT, "\n")
cat("Size:", FIG_W_MM, "x", FIG_H_MM, "mm\n")
