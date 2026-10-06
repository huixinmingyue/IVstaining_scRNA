#!/usr/bin/env Rscript
# =============================================================================
# Fig. S2 — B-cell subcluster identity, ImmGen MZ/Fo recovery, and
# nCount-adjusted within-cluster IV GSEA (all KEGG vs signaling subset).
# New folder only. Does not overwrite any previous FigS2 or Fig2 package.
#
# Core conclusion
#   Cluster 1 versus 0 identity genes recover ImmGen MZ versus Fo; cluster 2
#   is not T1, GC or plasmablast; nCount-adjusted within-cluster IV GSEA is
#   housekeeping/metabolic; the signaling subset of the same GSEA has few
#   FDR-significant sets (cluster 1 n = 1; cluster 0 n = 4).
#
# Archetype : quantitative grid
# Backend   : R exclusive
# Size      : 183 x 320 mm
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(stringr)
  library(ggrepel)
  library(svglite)
  library(ragg)
  library(systemfonts)
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
FIG_H_MM <- 320
width_mm = 183
height_mm = 320
BASE <- 6.5
FAM  <- "Arial"

# Helvetica stand-in if any leftover pdf() device asks for Arial.
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM
    do.call(reg, args)
  }
}


fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) cat("font:", FAM, "->", fm$path, "\n")

CL_COL <- c("0" = "#4C6E8A", "1" = "#D1495B", "2" = "#E9A03B")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"
DIV_LO <- "#3B6FA0"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#C0392B"

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
      strip.text = element_text(size = BASE - 0.3, face = "italic"),
      plot.title = element_text(size = BASE + 0.4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE - 0.7, hjust = 0.5, colour = "#333333"),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.margin = margin(2, 3, 2, 5),
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

p3_genes <- c("Cr2", "Fcer2a", "Cd93", "Cd24a", "Mki67", "Aicda", "Bcl6",
              "Sdc1", "Jchain", "Xbp1", "Cd1d1", "Ighd", "Ighm", "Sell", "Adgre1")
p6_genes <- c("Adgre5", "S1pr1", "Cxcr5", "Itgb1", "Itgb2", "Notch2",
              "Cr2", "Fcer2a", "Cd1d1", "Dtx1")
p6_only  <- setdiff(p6_genes, p3_genes)

# =============================================================================
# a  Query markers (P3a) + unique positioning genes (P6a)
# =============================================================================
vlong <- read.csv(file.path(SRC, "FigS2a_violin_expr.csv"),
                  stringsAsFactors = FALSE)
n_cell_a <- n_distinct(vlong$cell)
n_gene_a <- n_distinct(vlong$gene)
stopifnot(n_cell_a == 5975L, n_gene_a == 22L)
stopifnot(nrow(vlong) == n_cell_a * n_gene_a)
vlong$seurat_clusters <- factor(as.character(vlong$seurat_clusters),
                                levels = c("0", "1", "2"))
cl_n <- vlong %>% distinct(cell, seurat_clusters) %>% count(seurat_clusters)
cat("S2a cells by cluster:\n"); print(cl_n)

v_query <- vlong %>%
  filter(block == "Query markers") %>%
  mutate(gene = factor(gene, levels = p3_genes))
v_pos <- vlong %>%
  filter(block == "Positioning genes") %>%
  mutate(gene = factor(gene, levels = p6_only))
stopifnot(n_distinct(v_query$gene) == 15L, n_distinct(v_pos$gene) == 7L)
# Shared Cr2 / Fcer2a / Cd1d1 are kept in the query block only.
stopifnot(!any(c("Cr2", "Fcer2a", "Cd1d1") %in% as.character(v_pos$gene)))

violin_panel <- function(d, ttl, ncol = 5, show_legend = FALSE) {
  ggplot(d, aes(seurat_clusters, expr, fill = seurat_clusters)) +
    geom_violin(scale = "width", linewidth = 0.15, colour = "grey20") +
    stat_summary(fun = median, geom = "point", size = 0.55, colour = "black") +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = CL_COL, name = "Subcluster") +
    labs(x = "B-cell subcluster", y = "log-normalised expression", title = ttl) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      strip.text = element_text(size = BASE - 0.4, face = "italic"),
      axis.text.x = element_text(size = BASE - 0.4),
      axis.text.y = element_text(size = BASE - 0.8)
    )
}

pA_query <- violin_panel(v_query, "Query markers for cluster 2", ncol = 5) +
  labs(tag = "a") +
  theme(plot.margin = margin(8, 4, 2, 10))
pA_pos <- violin_panel(v_pos, "Positioning genes", ncol = 7, show_legend = TRUE)

write.csv(vlong, file.path(SDP, "SourceData_FigS2a_violin.csv"), row.names = FALSE)

# =============================================================================
# b  ImmGen identity-gene heatmap (P2 panel c)
# =============================================================================
zdf <- read.csv(file.path(SRC, "P2_ImmGen_zscore_identity_genes.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
id_genes <- read.csv(file.path(SRC, "P2_identity_gene_set.csv"),
                     stringsAsFactors = FALSE)
expect <- read.csv(file.path(SRC, "P2_expected_genes.csv"),
                   stringsAsFactors = FALSE)$gene

zmat <- as.matrix(zdf[, -1])
rownames(zmat) <- zdf$gene
n_gene_b <- nrow(zmat)
stopifnot(n_gene_b == 85L)
lab_map <- c("B.T1.Sp" = "T1", "B.T2.Sp" = "T2", "B.T3.Sp" = "T3",
             "B.Fo.Sp" = "Fo", "B.MZ.Sp" = "MZ", "B.mem.Sp" = "Mem",
             "B.GC.CB.Sp" = "GC-CB", "B.GC.CC.Sp" = "GC-CC")
stopifnot(all(colnames(zmat) %in% names(lab_map)))
colnames(zmat) <- unname(lab_map[colnames(zmat)])

memb <- ifelse(id_genes$avg_log2FC[match(rownames(zmat), id_genes$gene)] > 0,
               "cluster 1 up", "cluster 0 up")
names(memb) <- rownames(zmat)
memb[is.na(memb)] <- "unassigned"
n_unassigned <- sum(memb == "unassigned")
mz_vec <- as.numeric(zmat[, "MZ"])
z_ord <- zmat[order(memb[rownames(zmat)], -mz_vec), , drop = FALSE]
# ggplot tile heatmap (same z matrix, split and expected-gene marks as the
# source ComplexHeatmap). wrap_elements(grid.grabExpr) clipped column names
# and collapsed row height in the assembled patchwork.
gene_levels <- rev(rownames(z_ord))
expect_ord <- gene_levels[gene_levels %in% expect]
z_tile <- as.data.frame(z_ord)
z_tile$gene <- factor(rownames(z_ord), levels = gene_levels)
z_tile$set <- factor(unname(memb[rownames(z_ord)]),
                     levels = c("cluster 1 up", "cluster 0 up"))
z_long <- z_tile %>%
  pivot_longer(-c(gene, set), names_to = "subset", values_to = "z") %>%
  mutate(
    subset = factor(subset, levels = c("T1", "T2", "T3", "Fo", "MZ",
                                       "Mem", "GC-CB", "GC-CC")),
    z = pmin(pmax(as.numeric(z), -2.2), 2.2)
  )

lab_df <- z_long %>%
  filter(as.character(gene) %in% expect, subset == "GC-CC") %>%
  distinct(gene, set, .keep_all = TRUE)

pB <- ggplot(z_long, aes(subset, gene, fill = z)) +
  geom_tile() +
  geom_text_repel(
    data = lab_df, aes(x = 8.15, y = gene, label = gene),
    inherit.aes = FALSE, size = 1.8, family = FAM, fontface = "italic",
    direction = "y", nudge_x = 0.7, segment.size = 0.2,
    segment.colour = "#666666", min.segment.length = 0,
    box.padding = 0.12, point.padding = 0.05, max.overlaps = 20,
    show.legend = FALSE
  ) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = 0,
    limits = c(-2.2, 2.2), name = "z",
    guide = guide_colourbar(barheight = unit(22, "mm"),
                            barwidth = unit(2.4, "mm"),
                            title.position = "top")
  ) +
  scale_y_discrete(breaks = NULL, expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "Identity genes, ImmGen splenic B subsets",
       tag = "b") +
  theme(
    axis.text.x = element_text(size = BASE - 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    strip.text.y.left = element_text(size = BASE - 0.4, face = "bold",
                                     angle = 90),
    strip.placement = "outside",
    legend.position = "right",
    plot.margin = margin(4, 14, 2, 8)
  )

write.csv(z_long, file.path(SDP, "SourceData_FigS2b_ImmGen_z.csv"),
          row.names = FALSE)

# =============================================================================
# c / d  Within-cluster IV KEGG GSEA after nCount regression
# c = top 15 of all sets; d = top 15 of the signaling subset
# =============================================================================
kegg_all <- read.csv(file.path(SRC, "P4_gseKEGG_all.csv"), stringsAsFactors = FALSE)
kegg_all$cluster <- as.character(kegg_all$cluster)
kegg_all$is_signaling <- grepl("signaling", kegg_all$Description, ignore.case = TRUE)
n_kegg_c <- nrow(kegg_all)
n_sig_c1 <- sum(kegg_all$cluster == "1" & kegg_all$p.adjust < 0.05, na.rm = TRUE)
n_sig_c0 <- sum(kegg_all$cluster == "0" & kegg_all$p.adjust < 0.05, na.rm = TRUE)
n_sig_d1 <- sum(kegg_all$cluster == "1" & kegg_all$is_signaling &
                  kegg_all$p.adjust < 0.05, na.rm = TRUE)
n_sig_d0 <- sum(kegg_all$cluster == "0" & kegg_all$is_signaling &
                  kegg_all$p.adjust < 0.05, na.rm = TRUE)
n_sig_n <- c(cluster1 = sum(kegg_all$cluster == "1" & kegg_all$is_signaling),
             cluster0 = sum(kegg_all$cluster == "0" & kegg_all$is_signaling))

prep_dot <- function(df, cl, signaling_only = FALSE, strip_signaling = FALSE) {
  d <- df %>%
    filter(cluster == cl, is.finite(NES), is.finite(pvalue))
  if (signaling_only) d <- filter(d, is_signaling)
  d <- d %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(
      raw = gsub(" - Mus musculus.*", "", Description),
      label = raw,
      sig = p.adjust < 0.05,
      neglogp = -log10(pmax(pvalue, 1e-12))
    )
  if (strip_signaling) {
    d$label <- trimws(sub(" signaling pathway$", "", d$label, ignore.case = TRUE))
    d$label <- trimws(gsub("\\s+", " ", d$label))
    d$label <- gsub("TGF-beta", "TGF-\u03b2", d$label)
    d$label <- ifelse(d$label == "", d$raw, d$label)
  }
  d$label <- stringr::str_trunc(d$label, 55)
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  d
}

d_c1 <- prep_dot(kegg_all, "1")
d_c0 <- prep_dot(kegg_all, "0")
d_d1 <- prep_dot(kegg_all, "1", signaling_only = TRUE, strip_signaling = TRUE)
d_d0 <- prep_dot(kegg_all, "0", signaling_only = TRUE, strip_signaling = TRUE)
stopifnot(nrow(d_c1) == 15L, nrow(d_c0) == 15L,
          nrow(d_d1) == 15L, nrow(d_d0) == 15L)

gsea_iv <- function(d, ttl, n_sig_all, fill_lim, size_lim,
                    tag = NULL, show_legend = TRUE) {
  p <- ggplot(d, aes(NES, label)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "#999999") +
    geom_segment(aes(x = 0, xend = NES, y = label, yend = label),
                 linewidth = 0.2, colour = "#CFD6DA") +
    geom_point(aes(size = setSize, fill = neglogp),
               shape = 21, colour = "black", stroke = 0.2) +
    geom_point(data = subset(d, sig), aes(size = setSize),
               shape = 21, fill = NA, colour = "black", stroke = 0.7,
               show.legend = FALSE) +
    scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, name = "-log10 P",
                        limits = fill_lim) +
    scale_size_continuous(range = c(0.8, 2.6), name = "Set size",
                          limits = size_lim) +
    labs(x = "NES", y = NULL, title = ttl,
         subtitle = sprintf("minGSSize 15; ring = FDR q < 0.05 (n sig = %d)",
                            n_sig_all)) +
    theme(
      axis.text.y = element_text(size = 5.8, face = "plain"),
      panel.grid.major.y = element_line(linewidth = 0.15, colour = "#F2F4F5"),
      legend.position = if (show_legend) "right" else "none"
    )
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

fill_lim_c <- range(c(d_c1$neglogp, d_c0$neglogp))
size_lim_c <- range(c(d_c1$setSize, d_c0$setSize))
fill_lim_d <- range(c(d_d1$neglogp, d_d0$neglogp))
size_lim_d <- range(c(d_d1$setSize, d_d0$setSize))

pC1 <- gsea_iv(d_c1, "KEGG, cluster 1 IV+ vs IV\u2212", n_sig_c1,
               fill_lim_c, size_lim_c, tag = "c", show_legend = FALSE)
pC0 <- gsea_iv(d_c0, "KEGG, cluster 0 IV+ vs IV\u2212", n_sig_c0,
               fill_lim_c, size_lim_c, show_legend = TRUE)
pD1 <- gsea_iv(d_d1, "KEGG signaling, cluster 1 IV+ vs IV\u2212", n_sig_d1,
               fill_lim_d, size_lim_d, tag = "d", show_legend = FALSE)
pD0 <- gsea_iv(d_d0, "KEGG signaling, cluster 0 IV+ vs IV\u2212", n_sig_d0,
               fill_lim_d, size_lim_d, show_legend = TRUE)

write.csv(bind_rows(d_c1, d_c0, .id = "panel_col"),
          file.path(SDP, "SourceData_FigS2c_withinCluster_GSEA_top15.csv"),
          row.names = FALSE)
write.csv(bind_rows(d_d1, d_d0, .id = "panel_col"),
          file.path(SDP, "SourceData_FigS2d_signaling_GSEA_top15.csv"),
          row.names = FALSE)

# =============================================================================
# Assemble
# =============================================================================
row_c <- (pC1 | pC0) + plot_layout(widths = c(1, 1.12))
row_d <- (pD1 | pD0) + plot_layout(widths = c(1, 1.12))

fig <- (pA_query / pA_pos / pB / row_c / row_d) +
  plot_layout(heights = c(1.28, 0.62, 1.55, 1.12, 1.28)) &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FAM))

out <- file.path(ROOT, "FigureS2")
save_pub(fig, out, FIG_W_MM, FIG_H_MM)
cat("wrote", out, "pdf/svg/png/tiff\n")

# QA numbers
qa <- file.path(QA, "make_figureS2_stats.txt")
sink(qa)
cat("Fig. S2 stats\n")
cat("cells S2a:", n_cell_a, "\n")
print(cl_n)
cat("genes S2a unique:", n_gene_a, "(P3 15 + P6-only 7; Cr2/Fcer2a/Cd1d1 not duplicated)\n")
cat("ImmGen genes S2b:", n_gene_b, "unassigned:", n_unassigned, "\n")
cat("expect marks:", paste(expect_ord, collapse = ", "), "\n")
cat("P4 GSEA rows:", n_kegg_c, " FDR q<0.05 cluster1:", n_sig_c1,
    "cluster0:", n_sig_c0, "\n")
cat("P4 signaling n cluster1:", n_sig_n[["cluster1"]],
    "cluster0:", n_sig_n[["cluster0"]],
    " FDR q<0.05 cluster1:", n_sig_d1, "cluster0:", n_sig_d0, "\n")
sink()
cat("wrote", qa, "\n")
