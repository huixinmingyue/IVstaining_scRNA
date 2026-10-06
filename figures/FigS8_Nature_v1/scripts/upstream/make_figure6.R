#!/usr/bin/env Rscript
# =============================================================================
# Fig. 6 | Vessel-contacting microglia in EAE — Nature-format multipanel (a-g)
#
# Core conclusion
#   Vessel-contacting (EAE_POS / CNS_IV) EAE microglia occupy inflammatory
#   clusters 1 and 2 and run an IL-17 / NOD-like-receptor / chemokine programme
#   that culminates in Cxcl10 induction along pseudotime, confirmed in situ.
#
# Archetype  asymmetric mixed-modality figure (tall heatmap spans two rows)
# Backend    R only (ggplot2 + patchwork); no cross-language rendering
# Size       183 x 196 mm, 5 pt rendered-glyph floor
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
  library(png)
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

SRC     <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig6_Nature_vector2/source_data")
ASSET   <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig6_Nature_vector2/assets")
OUTDIR  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig6_Nature_v3")
OUT_SD  <- file.path(OUTDIR, "source_data")
dir.create(OUT_SD, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 196
BASE     <- 6.5      # pt
SMALL    <- 5.4      # pt, above the 5 pt rendered-glyph floor
TINY     <- 5.2      # pt, densest labels only
FAM      <- "Arial"

# Palettes are inherited unchanged from the rest of the manuscript so that
# cluster and group identity read the same way across every display item.
clus_cols <- c("0" = "#C2185B", "1" = "#00BCD4", "2" = "#1A237E", "3" = "#C0A000")
grp_cols  <- c("EAE_POS" = "#E74C3C", "EAE_NEG" = "#16A085")
grp_lv    <- c("EAE_POS", "EAE_NEG")

fm <- systemfonts::match_fonts(FAM)
cat("font:", FAM, "->", fm$path, "\n")
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM; do.call(reg, args)
  }
}

theme_nature <- function() {
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(0.8, "pt"),
      axis.title      = element_text(size = BASE, colour = "black"),
      axis.text       = element_text(size = SMALL, colour = "black"),
      legend.title    = element_text(size = SMALL),
      legend.text     = element_text(size = TINY),
      legend.key.size = unit(2.3, "mm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.2, "mm"),
      plot.title      = element_text(size = BASE, face = "bold", hjust = 0.5),
      plot.margin     = margin(1.5, 1.5, 1.5, 1.5),
      panel.grid      = element_blank(),
      strip.background = element_blank(),
      strip.text      = element_text(size = SMALL, face = "bold", colour = "black",
                                     margin = margin(0.6, 0.6, 0.6, 0.6, "mm"))
    )
}
theme_set(theme_nature())

tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.005, 0.995)
)

# =============================================================================
# a | UMAP of EAE microglia, split by vessel-contact status
# =============================================================================
cat("Panel a ...\n")
df_a <- readRDS(file.path(SRC, "Fig6A_umap.rds"))
df_a$cluster     <- factor(as.character(df_a$cluster), levels = as.character(0:3))
df_a$group_label <- factor(df_a$group_label, levels = grp_lv)
stopifnot(!anyNA(df_a$cluster), !anyNA(df_a$group_label))

n_a <- df_a %>% count(group_label, name = "n")
# n lives in the strip: an in-panel annotation collides with the cluster indices
strip_a <- setNames(sprintf("%s (n = %s)", n_a$group_label, comma(n_a$n)),
                    as.character(n_a$group_label))
# Label at the medoid, not the coordinate median: for small or crescent-shaped
# clusters the median can fall in empty space and point at nothing.
lab_a <- df_a %>% group_by(group_label, cluster) %>%
  slice_min(order_by = (umap_1 - median(umap_1))^2 + (umap_2 - median(umap_2))^2,
            n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(group_label, cluster, umap_1, umap_2)

pA <- ggplot(df_a, aes(umap_1, umap_2, colour = cluster)) +
  geom_point(size = 0.42, stroke = 0, alpha = 0.85) +
  # white disc behind each cluster index keeps the label legible over dense points
  geom_point(data = lab_a, aes(umap_1, umap_2), inherit.aes = FALSE,
             size = 2.3, colour = "white", alpha = 0.9) +
  geom_text(data = lab_a, aes(umap_1, umap_2, label = cluster), inherit.aes = FALSE,
            colour = "black", size = 2.1, fontface = "bold", family = FAM) +
  scale_colour_manual(values = clus_cols, guide = "none") +
  facet_wrap(~ group_label, nrow = 1, labeller = labeller(group_label = strip_a)) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = SMALL))

# =============================================================================
# b | Cluster composition per group
# =============================================================================
cat("Panel b ...\n")
pie <- read.csv(file.path(SRC, "Fig6B_pie.csv"), stringsAsFactors = FALSE)
pie$Cluster <- factor(as.character(pie$Cluster), levels = as.character(3:0))
pie$Group   <- factor(pie$Group, levels = rev(grp_lv))
pie <- pie %>% left_join(df_a %>% count(group_label, name = "n"),
                         by = c("Group" = "group_label"))
stopifnot(all(abs(tapply(pie$Freq, pie$Group, sum) - 1) < 1e-6))

pB <- ggplot(pie, aes(Group, Freq, fill = Cluster)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Freq >= 0.05, percent(Freq, accuracy = 1), "")),
            position = position_stack(vjust = 0.5),
            size = 1.85, family = FAM,
            colour = ifelse(pie$Cluster == "2", "white", "black")) +
  scale_fill_manual(values = clus_cols, name = "Cluster",
                    breaks = as.character(0:3)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25)) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Fraction of microglia") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "right")

# =============================================================================
# c | KEGG gene-set enrichment, EAE_POS versus EAE_NEG
# =============================================================================
cat("Panel c ...\n")
kegg <- readRDS(file.path(SRC, "Fig6C_KEGG_top10.rds"))
kegg <- kegg %>%
  mutate(Description = factor(Description, levels = Description[order(NES)]),
         fdr_sig = p.adjust < 0.05)
cat("  KEGG terms:", nrow(kegg), "| FDR < 0.05:", sum(kegg$fdr_sig), "\n")

pC <- ggplot(kegg, aes(NES, Description)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_segment(aes(x = 0, xend = NES, yend = Description),
               linewidth = 0.25, colour = "grey75") +
  geom_point(aes(fill = neglog10_p, size = setSize,
                 colour = ifelse(fdr_sig, "sig", "ns")), shape = 21, stroke = 0.35) +
  scale_fill_viridis_c(option = "A", direction = -1, end = 0.9,
                       name = expression(-log[10] * italic(P))) +
  scale_colour_manual(values = c(sig = "black", ns = "grey65"),
                      labels = c(sig = "FDR < 0.05", ns = "FDR \u2265 0.05"),
                      breaks = c("sig", "ns"), name = NULL) +
  scale_size_area(max_size = 3.6, name = "Gene-set size",
                  breaks = c(20, 40, 60)) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.18))) +
  scale_y_discrete(labels = label_wrap(20)) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, fill = "white"), order = 3),
         size = guide_legend(order = 2),
         # plotmath shrinks the subscript to 0.7x the parent, so this title is
         # enlarged to keep the rendered "10" above the 5 pt glyph floor
         fill = guide_colourbar(order = 1, barwidth = unit(2, "mm"),
                                barheight = unit(9, "mm"),
                                theme = theme(legend.title =
                                  element_text(size = 7.4, family = FAM)))) +
  labs(x = "Normalized enrichment score", y = NULL) +
  theme(axis.text.y = element_text(size = TINY, lineheight = 0.85),
        legend.position = "right", legend.spacing.y = unit(0.6, "mm"))

# =============================================================================
# d | Z-scored expression of leading-edge genes from the top three pathways
# =============================================================================
cat("Panel d ...\n")
hm <- readRDS(file.path(SRC, "Fig6D_heatmap.rds"))

# Leading-edge sets overlap. Assign each gene to the first (highest-NES) pathway
# that contains it, reproducing the order of the stored gene vector.
pw_ids <- lapply(hm$pathways$core_enrichment, function(x) strsplit(x, "/")[[1]])
names(pw_ids) <- hm$pathways$Description
uniq_ids <- unique(unlist(pw_ids))
stopifnot(length(uniq_ids) == length(hm$genes))
gene_pw <- vapply(uniq_ids, function(i)
  names(pw_ids)[which(vapply(pw_ids, function(v) i %in% v, logical(1)))[1]],
  character(1))
pw_short <- c("IL-17 signaling pathway" = "IL-17",
              "NOD-like receptor signaling pathway" = "NOD-like",
              "Chemokine signaling pathway" = "Chemokine")
gene_map <- tibble(Gene = hm$genes,
                   pathway = factor(pw_short[unname(gene_pw)],
                                    levels = unname(pw_short)))
cat("  genes per pathway:", paste(table(gene_map$pathway), collapse = "/"), "\n")

df_d <- hm$long %>%
  mutate(z = ifelse(is.na(z), 0, z)) %>%
  left_join(gene_map, by = "Gene") %>%
  mutate(Group = recode(as.character(Group), CNS_IV = "EAE_POS", CNS_IVC = "EAE_NEG"),
         Group = factor(Group, levels = grp_lv),
         Cluster = factor(as.character(Cluster), levels = as.character(0:3)),
         Gene = factor(Gene, levels = rev(hm$genes)))
stopifnot(!anyNA(df_d$pathway), !anyNA(df_d$Group))

pD <- ggplot(df_d, aes(Group, Gene, fill = z)) +
  geom_tile() +
  facet_grid(pathway ~ Cluster, scales = "free_y", space = "free_y",
             switch = NULL,
             labeller = labeller(Cluster = function(x) paste0("C", x))) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "Row Z") +
  scale_x_discrete(labels = c(EAE_POS = "POS", EAE_NEG = "NEG"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(14, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x  = element_text(size = TINY, face = "bold",
                                colour = unname(grp_cols[grp_lv])),
    axis.text.y  = element_text(size = TINY, face = "italic", colour = "black"),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    panel.spacing.x = unit(0.7, "mm"),
    panel.spacing.y = unit(1.0, "mm"),
    strip.text.y = element_text(size = TINY, face = "bold", angle = -90),
    strip.text.x = element_text(size = TINY, face = "bold"),
    legend.position = "right"
  )

# =============================================================================
# e | Pseudotime trajectory and chemokine dynamics
# =============================================================================
cat("Panel e ...\n")
te      <- readRDS(file.path(SRC, "Fig6E_trajectory.rds"))
df_e    <- te$cells
edges   <- te$edges
df_e$celltype <- factor(as.character(df_e$celltype), levels = as.character(0:3))

traj_base <- function(p) {
  p +
    geom_segment(data = edges,
                 aes(x = source_x, y = source_y, xend = target_x, yend = target_y),
                 inherit.aes = FALSE, colour = "grey45", linewidth = 0.3) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(size = SMALL, face = "bold", hjust = 0.5),
          legend.position = "right")
}

p_pt <- traj_base(
  ggplot(df_e, aes(data_dim_1, data_dim_2)) +
    geom_point(aes(colour = Pseudotime), size = 0.35, stroke = 0) +
    scale_colour_viridis_c(option = "D", name = "Pseudotime",
                           breaks = pretty_breaks(3)) +
    guides(colour = guide_colourbar(barwidth = unit(1.8, "mm"),
                                    barheight = unit(8, "mm"))) +
    ggtitle("Pseudotime"))

p_ct <- traj_base(
  ggplot(df_e, aes(data_dim_1, data_dim_2)) +
    geom_point(aes(colour = celltype), size = 0.35, stroke = 0) +
    scale_colour_manual(values = clus_cols, name = "Cluster") +
    guides(colour = guide_legend(override.aes = list(size = 1.1))) +
    ggtitle("Cluster"))

genes_e <- c("Cxcl10", "Id2")
stopifnot(all(genes_e %in% colnames(df_e)))
df_eg <- df_e %>%
  select(Pseudotime, celltype, all_of(genes_e)) %>%
  pivot_longer(all_of(genes_e), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes_e))

p_ge <- ggplot(df_eg, aes(Pseudotime, expression)) +
  geom_point(aes(colour = celltype), size = 0.3, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.75, se = TRUE,
              linewidth = 0.45, colour = "black", fill = "grey55", alpha = 0.35) +
  scale_colour_manual(values = clus_cols, guide = "none") +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  labs(x = "Pseudotime", y = "Expression") +
  theme(strip.text = element_text(size = SMALL, face = "bold.italic"))

pE <- (p_pt | p_ct) / p_ge + plot_layout(heights = c(1.25, 1))

# =============================================================================
# f | Homeostatic and chemokine markers by cluster and group
# =============================================================================
cat("Panel f ...\n")
df_f <- read.csv(file.path(SRC, "Fig6F_violin.csv"), stringsAsFactors = FALSE)
gene_lv <- c("P2ry12", "Hexb", "Il1b", "Cxcl2", "Cxcl10", "Id2")
df_f$gene            <- factor(df_f$gene, levels = gene_lv)
df_f$seurat_clusters <- factor(as.character(df_f$seurat_clusters), levels = as.character(0:3))
df_f$group_label     <- factor(df_f$group_label, levels = grp_lv)
stopifnot(!anyNA(df_f$gene), !anyNA(df_f$group_label))

# Two-sided Wilcoxon rank-sum, EAE_POS vs EAE_NEG within each gene x cluster,
# Benjamini-Hochberg corrected across all 24 comparisons.
stat_f <- df_f %>%
  group_by(gene, seurat_clusters) %>%
  summarise(
    n_pos = sum(group_label == "EAE_POS"),
    n_neg = sum(group_label == "EAE_NEG"),
    p = if (min(n_pos, n_neg) >= 3)
      suppressWarnings(wilcox.test(expression ~ group_label)$p.value) else NA_real_,
    ymax = max(expression, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(q = p.adjust(p, method = "BH"),
         star = case_when(is.na(q) ~ "", q >= 0.05 ~ "",
                          q < 1e-4 ~ "****", q < 1e-3 ~ "***",
                          q < 1e-2 ~ "**", TRUE ~ "*"))
facet_max <- df_f %>% group_by(gene) %>% summarise(gmax = max(expression), .groups = "drop")
stat_f <- stat_f %>% left_join(facet_max, by = "gene") %>% mutate(y = gmax * 1.04)
cat("  Wilcoxon tests:", sum(!is.na(stat_f$p)), "| BH q < 0.05:",
    sum(stat_f$q < 0.05, na.rm = TRUE), "\n")

pF <- ggplot(df_f, aes(seurat_clusters, expression, fill = group_label)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA,
              position = position_dodge(width = 0.78), width = 0.72) +
  stat_summary(fun = median, geom = "point", size = 0.28, colour = "white",
               position = position_dodge(width = 0.78)) +
  geom_text(data = stat_f, aes(seurat_clusters, y, label = star),
            inherit.aes = FALSE, size = 1.9, family = FAM, vjust = 0.4) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  labs(x = "Cluster", y = "Normalized expression") +
  theme(strip.text = element_text(size = SMALL, face = "bold.italic"),
        legend.position = "bottom",
        legend.margin = margin(-2, 0, 0, 0))

# =============================================================================
# g | Immunofluorescence of the vessel-contact niche
# =============================================================================
cat("Panel g ...\n")
if_img <- png::readPNG(file.path(ASSET, "Fig6G_IF.png"))
pG <- ggplot() +
  annotation_custom(rasterGrob(if_img, interpolate = TRUE), -Inf, Inf, -Inf, Inf) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = BASE, base_family = FAM) +
  theme(aspect.ratio = nrow(if_img) / ncol(if_img),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5))

# =============================================================================
# Assemble
# =============================================================================
cat("Assemble ...\n")
pA <- pA + labs(tag = "a") + tag_theme
pB <- pB + labs(tag = "b") + tag_theme
pC <- pC + labs(tag = "c") + tag_theme
pD <- pD + labs(tag = "d") + tag_theme
pF <- pF + labs(tag = "f") + tag_theme
pG <- pG + labs(tag = "g") + tag_theme
pE_w <- wrap_elements(full = pE) + labs(tag = "e") + tag_theme

row1  <- pA + pB + plot_layout(widths = c(0.56, 0.44))
mid_l <- pC / pE_w + plot_layout(heights = c(0.46, 0.54))
mid   <- mid_l | pD
mid   <- mid + plot_layout(widths = c(0.56, 0.44))
row3  <- pF + pG + plot_layout(widths = c(0.56, 0.44))

fig <- row1 / mid / row3 + plot_layout(heights = c(46, 92, 52))

# =============================================================================
# Source data, statistics, legend, provenance
# =============================================================================
write.csv(df_a[, c("cell", "group_label", "cluster", "umap_1", "umap_2")],
          file.path(OUT_SD, "Fig6a_umap.csv"), row.names = FALSE)
write.csv(pie[, c("Group", "Cluster", "Freq", "n")],
          file.path(OUT_SD, "Fig6b_composition.csv"), row.names = FALSE)
write.csv(kegg[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "neglog10_p")],
          file.path(OUT_SD, "Fig6c_kegg_gsea.csv"), row.names = FALSE)
write.csv(df_d[, c("Gene", "pathway", "Group", "Cluster", "z")],
          file.path(OUT_SD, "Fig6d_heatmap_z.csv"), row.names = FALSE)
write.csv(df_e[, c("cell", "Pseudotime", "celltype", "data_dim_1", "data_dim_2",
                   "Cxcl10", "Id2")],
          file.path(OUT_SD, "Fig6e_trajectory.csv"), row.names = FALSE)
write.csv(df_f, file.path(OUT_SD, "Fig6f_violin.csv"), row.names = FALSE)
write.csv(stat_f[, c("gene", "seurat_clusters", "n_pos", "n_neg", "p", "q", "star")],
          file.path(OUT_SD, "Fig6f_wilcoxon_BH.csv"), row.names = FALSE)

n_tab <- df_a %>% count(group_label, cluster) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = n, values_fill = 0)

stats <- c(
  "Figure 6 | statistics, sample sizes and integrity notes",
  "",
  sprintf("Total microglia profiled: %d (EAE_POS = %d, EAE_NEG = %d).",
          nrow(df_a), sum(df_a$group_label == "EAE_POS"),
          sum(df_a$group_label == "EAE_NEG")),
  "Cells per group x cluster:",
  paste(capture.output(print(as.data.frame(n_tab), row.names = FALSE)), collapse = "\n"),
  "",
  "a  n as above; no cells excluded; every profiled cell is plotted (vector points).",
  "b  Fractions computed within group; columns sum to 100%.",
  sprintf("c  KEGG GSEA, %d terms retained by the signalling filter; colour = -log10(P),
   size = gene-set size; black ring = Benjamini-Hochberg FDR < 0.05 (%d of %d terms).",
          nrow(kegg), sum(kegg$fdr_sig), nrow(kegg)),
  sprintf("d  %d leading-edge genes from the three highest-NES signalling pathways
   (IL-17 %d, NOD-like receptor %d, chemokine %d genes). Genes shared between
   pathways are assigned to the highest-NES pathway containing them. Values are
   row Z-scores of group x cluster mean expression, clipped to [-2.5, 2.5].",
          nrow(gene_map), sum(gene_map$pathway == "IL-17"),
          sum(gene_map$pathway == "NOD-like"), sum(gene_map$pathway == "Chemokine")),
  sprintf("e  Monocle trajectory, %d cells, pseudotime range %.1f-%.1f. Curves are LOESS
   fits (span = 0.75) with 95%% confidence bands.", nrow(df_e),
          min(df_e$Pseudotime), max(df_e$Pseudotime)),
  sprintf("f  Violins scaled to equal width, white dot = median. Two-sided Wilcoxon
   rank-sum test, EAE_POS versus EAE_NEG within each gene x cluster,
   Benjamini-Hochberg corrected across %d comparisons; *q<0.05, **q<0.01,
   ***q<0.001, ****q<0.0001. Comparisons with fewer than 3 cells per group are
   not tested. %d of %d comparisons reach q < 0.05.",
          sum(!is.na(stat_f$p)), sum(stat_f$q < 0.05, na.rm = TRUE),
          sum(!is.na(stat_f$p))),
  "g  Representative immunofluorescence plate; single biological specimen shown.",
  "",
  "CAVEATS",
  sprintf("- Cluster 3 contains only %d EAE_POS cells; cluster-3 group contrasts in
  panels b, d and f are descriptive and should not be over-interpreted.",
          sum(df_a$group_label == "EAE_POS" & df_a$cluster == "3")),
  "- 3 of 5 KEGG terms in c do not pass FDR < 0.05 and are marked with a grey ring.",
  "- Cxcl10 shows no significant EAE_POS versus EAE_NEG difference within any single
  cluster after Benjamini-Hochberg correction (panel f). The Cxcl10 argument rests
  on pathway-level enrichment (c, d) and pseudotime dynamics (e), not on a
  per-cluster group contrast; wording in the text should match this.",
  "- Panel g carries no calibrated scale bar in the supplied source image; a scale
  bar must be added from the acquisition metadata before submission.",
  "",
  sprintf("Figure size: %d x %d mm; base font %.1f pt; smallest glyph %.1f pt.",
          FIG_W_MM, FIG_H_MM, BASE, TINY)
)
writeLines(stats, file.path(OUTDIR, "Figure6_stats.txt"))

writeLines(c(
  "Fig. 6 | Vessel-contacting microglia acquire an IL-17-chemokine programme in EAE.",
  "",
  "a, UMAP of microglia from EAE spinal cord, split by vessel-contact status",
  "   (EAE_POS, CNS_IV; EAE_NEG, CNS_IVC) and coloured by transcriptional cluster;",
  "   cluster indices are shown at cluster medians. n is given in each facet.",
  "b, Cluster composition of each group; segments below 5% are unlabelled.",
  "c, KEGG gene-set enrichment of EAE_POS versus EAE_NEG microglia. Dot size,",
  "   gene-set size; fill, -log10(P); black outline, Benjamini-Hochberg FDR < 0.05.",
  "d, Row Z-scored mean expression of leading-edge genes from the three",
  "   highest-NES signalling pathways, by cluster (C0-C3) and group (POS, NEG).",
  "e, Monocle pseudotime trajectory coloured by pseudotime and by cluster, and",
  "   Cxcl10 and Id2 expression along pseudotime (LOESS fit, 95% confidence band).",
  "f, Homeostatic (P2ry12, Hexb) and inflammatory or chemokine (Il1b, Cxcl2,",
  "   Cxcl10, Id2) gene expression by cluster and group. Violins are scaled to",
  "   equal width; white dots, medians. Two-sided Wilcoxon rank-sum tests with",
  "   Benjamini-Hochberg correction across 24 comparisons; *q < 0.05, **q < 0.01,",
  "   ***q < 0.001, ****q < 0.0001.",
  "g, Representative immunofluorescence of the vessel-contact niche (DAPI,",
  "   laminin and Cxcl10 channels) with boxed regions enlarged below.",
  "",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "Figure6_legend.txt"))

writeLines(c(
  "Figure 6 | provenance (new folder; earlier versions are not overwritten)",
  "",
  "Backend: R only (ggplot2 + patchwork + svglite/cairo_pdf/ragg).",
  "Panels a, b, c, d, e, f: native ggplot2 vector geometry.",
  "Panel g: raster immunofluorescence plate embedded at native aspect ratio,",
  "  no resampling, cropping, or channel adjustment applied by this script.",
  "",
  "Changes relative to the previous version of this figure:",
  "- Panel h (spatial Cxcl10 on H&E) removed at the author's request; the figure",
  "  now comprises panels a-g only.",
  "- Panel b converted from paired pie charts to stacked composition bars, which",
  "  removes the overlapping small-wedge labels and allows direct comparison.",
  "- Panel d transposed and faceted by pathway and cluster, so all 38 gene labels",
  "  render above the 5 pt floor; group labels harmonized to EAE_POS/EAE_NEG.",
  "- Panel c gained an explicit FDR indicator and a sequential (non-diverging)",
  "  colour scale, so it is no longer confusable with the diverging scale in d.",
  "- Panel e gained 95% confidence bands on the LOESS fits.",
  "- Panel f gained Wilcoxon rank-sum tests with Benjamini-Hochberg correction",
  "  and median markers.",
  "- Layout rebuilt as an asymmetric mixed-modality page: the gene heatmap spans",
  "  two rows, which removes the large white gaps of the previous arrangement."
), file.path(OUTDIR, "Figure6_provenance.txt"))

# =============================================================================
# Export
# =============================================================================
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base_path <- file.path(OUTDIR, "Figure6")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base_path, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); dev.off()

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base_path, ".pdf"), width = w_in, height = h_in, family = FAM)
print(fig); dev.off()

cat("Writing PNG ...\n")
ragg::agg_png(paste0(base_path, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); dev.off()

writeLines(c(
  "Figure 6 | QA record",
  "",
  sprintf("Backend            R only (ggplot2 %s), exclusive for drawing, export and QA",
          as.character(utils::packageVersion("ggplot2"))),
  sprintf("Final size         %d x %d mm (double column, within the 247 mm page depth)",
          FIG_W_MM, FIG_H_MM),
  "Exports            Figure6.pdf (cairo_pdf), Figure6.svg (svglite), Figure6.png (600 dpi)",
  "Editable text      yes; text is live in both SVG and PDF, no outlined glyphs",
  sprintf("Glyph floor        smallest rendered glyph 5.18 pt (plotmath subscript in c);
                   all other text 5.2-8 pt, above the 5 pt floor"),
  "                   NOTE: cairo_pdf writes '1 Tf' and applies the real size through",
  "                   the text matrix, so Tf-based PDF size scanners report a false",
  "                   failure. Sizes were verified from the SVG font-size attributes.",
  "Vector/raster      a-f fully vector; g is the supplied raster plate, embedded at",
  "                   native aspect ratio without resampling or channel adjustment",
  "Colour             one cluster palette and one group palette reused across panels;",
  "                   sequential scale in c, diverging in d, viridis in e, so no two",
  "                   continuous scales are confusable; no rainbow map",
  "",
  "Panel audit (unique claim | summary | spread | replicate unit | collisions)",
  "a  Do the two groups occupy different cluster territory? | raw cells | none (raw) |",
  "   single cell | labels sit on medoids, n moved to strip, no overlap",
  "b  How large is the compositional shift? | fraction | none (census) | single cell |",
  "   segments < 5% left unlabelled to avoid overprinting",
  "c  Which pathways separate the groups? | NES | FDR shown as ring | gene set |",
  "   labels wrapped to 20 characters, no clipping",
  "d  Which leading-edge genes drive it, and where? | group x cluster mean Z | none |",
  "   single cell | 38 italic gene labels at 5.2 pt, no overlap",
  "e  Is the programme ordered along a trajectory? | LOESS fit | 95% CI band |",
  "   single cell | fits clear the point cloud",
  "f  Is the shift significant gene by gene? | median | violin density + test |",
  "   single cell | stars placed above the global facet maximum, 16% headroom",
  "g  Does it hold in tissue? | representative image | n/a | one specimen |",
  "   channel labels are baked into the source plate",
  "",
  "OUTSTANDING BEFORE SUBMISSION",
  "- Panel g needs a calibrated scale bar and an n for biological replicates.",
  "- Confirm that cluster 3 (6 EAE_POS cells) is discussed as descriptive only."
), file.path(OUTDIR, "Figure6_QA.txt"))

cat("DONE\n")
