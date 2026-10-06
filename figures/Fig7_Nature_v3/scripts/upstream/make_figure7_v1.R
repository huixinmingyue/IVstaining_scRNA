#!/usr/bin/env Rscript
# =============================================================================
# Fig. 7 v1 — Nature multipanel (Fig7_Nature_v1)
# New folder; does not overwrite any earlier package.
#
# Source mapping (requested)
#   a-c  previous microglia Fig. 6 panels a-c
#   d    IFN-γ / IL-17-axis module scores, All cells only
#   e-h  previous microglia Fig. 6 panels d-g
#
# Core conclusion
#   Vessel-contacting (POS) EAE microglia occupy inflammatory clusters 1-2,
#   carry an IFN-γ / IL-17-axis programme already in all cells, and run an
#   IL-17 / NOD-like / chemokine leading-edge ordered along pseudotime.
#
# Evidence chain
#   a  Where do POS and NEG microglia sit in cluster space?         UMAP
#   b  How large is the compositional shift?                        bars
#   c  Which signalling pathways separate the groups?               GSEA
#   d  Do IFN-γ and IL-17 axes already differ in all cells?         scores
#   e  Which leading-edge genes drive the programme, and where?     heatmap
#   f  Is the programme ordered along a trajectory?                 PT
#   g  Are marker shifts significant gene by gene?                  violins
#   h  Does the niche hold in tissue?                               IF
#
# Archetype : asymmetric mixed-modality
# Backend   : R exclusive
# Size      : 183 x 230 mm
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "Fig7_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
SRCFIG  <- file.path(OUTDIR, "source_figures")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 247
width_mm = 183
height_mm = 247
BASE  <- 6.5
SMALL <- 5.4
TINY  <- 5.2
FAM   <- "Arial"
PT5   <- 5.2

clus_cols <- c("0" = "#C2185B", "1" = "#00BCD4", "2" = "#1A237E", "3" = "#C0A000")
grp_cols  <- c(POS = "#E74C3C", NEG = "#16A085")
grp_lv    <- c("POS", "NEG")

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

to_posneg <- function(x) {
  x <- as.character(x)
  x[x %in% c("EAE_POS", "CNS_IV", "POS")] <- "POS"
  x[x %in% c("EAE_NEG", "CNS_IVC", "NEG")] <- "NEG"
  factor(x, levels = grp_lv)
}

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(0.8, "pt"),
      axis.title = element_text(size = BASE, colour = "black"),
      axis.text = element_text(size = SMALL, colour = "black"),
      legend.title = element_text(size = SMALL),
      legend.text = element_text(size = TINY),
      legend.key.size = unit(2.3, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.2, "mm"),
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = TINY, hjust = 0.5),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = SMALL, face = "bold", colour = "black",
                                margin = margin(0.6, 0.6, 0.6, 0.6, "mm"))
    )
)

tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.005, 0.995)
)

# -----------------------------------------------------------------------------
# a — UMAP
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
df_a <- read.csv(file.path(SRC, "Fig6A_umap.csv"), stringsAsFactors = FALSE)
n_a_before <- nrow(df_a)
df_a <- df_a[!is.na(df_a$umap_1) & !is.na(df_a$umap_2) &
               !is.na(df_a$cluster) & !is.na(df_a$group_label), ]
n_a_after <- nrow(df_a)
df_a$cluster <- factor(as.character(df_a$cluster), levels = as.character(0:3))
df_a$group_show <- to_posneg(df_a$group_label)
stopifnot(!anyNA(df_a$cluster), !anyNA(df_a$group_show))
n_a_pos <- sum(df_a$group_show == "POS")
n_a_neg <- sum(df_a$group_show == "NEG")
n_a_c3_pos <- sum(df_a$group_show == "POS" & df_a$cluster == "3")

n_a <- df_a %>% count(group_show, name = "n")
strip_a <- setNames(sprintf("%s (n = %s)", n_a$group_show, comma(n_a$n)),
                    as.character(n_a$group_show))
lab_a <- df_a %>% group_by(group_show, cluster) %>%
  slice_min(order_by = (umap_1 - median(umap_1))^2 + (umap_2 - median(umap_2))^2,
            n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(group_show, cluster, umap_1, umap_2)

write.csv(df_a[, c("cell", "group_label", "group_show", "cluster", "umap_1", "umap_2")],
          file.path(PLOTTED, "SourceData_Fig7a_umap.csv"), row.names = FALSE)
write.csv(lab_a, file.path(PLOTTED, "SourceData_Fig7a_labels.csv"), row.names = FALSE)

pA <- ggplot(df_a, aes(umap_1, umap_2, colour = cluster)) +
  geom_point(size = 0.42, stroke = 0, alpha = 0.85) +
  geom_point(data = lab_a, aes(umap_1, umap_2), inherit.aes = FALSE,
             size = 2.3, colour = "white", alpha = 0.9) +
  geom_text(data = lab_a, aes(umap_1, umap_2, label = cluster), inherit.aes = FALSE,
            colour = "black", size = 2.1, fontface = "bold", family = FAM) +
  scale_colour_manual(values = clus_cols, guide = "none") +
  facet_wrap(~ group_show, nrow = 1, labeller = labeller(group_show = strip_a)) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "a") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = SMALL)) +
  tag_theme

# -----------------------------------------------------------------------------
# b — composition
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
pie <- read.csv(file.path(SRC, "Fig6b_composition.csv"), stringsAsFactors = FALSE)
n_b_before <- nrow(pie)
pie <- pie[!is.na(pie$Freq), ]
n_b_after <- nrow(pie)
pie$Cluster <- factor(as.character(pie$Cluster), levels = as.character(3:0))
pie$group_show <- to_posneg(pie$Group)
pie$group_show <- factor(pie$group_show, levels = rev(grp_lv))
stopifnot(all(abs(tapply(pie$Freq, pie$group_show, sum) - 1) < 1e-6))
write.csv(pie, file.path(PLOTTED, "SourceData_Fig7b_composition.csv"),
          row.names = FALSE)

pB <- ggplot(pie, aes(group_show, Freq, fill = Cluster)) +
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
  labs(x = NULL, y = "Fraction of microglia", tag = "b") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "right") +
  tag_theme

# -----------------------------------------------------------------------------
# c — KEGG GSEA
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
kegg <- read.csv(file.path(SRC, "Fig6c_kegg_gsea.csv"), stringsAsFactors = FALSE)
n_c_before <- nrow(kegg)
kegg <- kegg %>%
  mutate(Description = factor(Description, levels = Description[order(NES)]),
         fdr_sig = p.adjust < 0.05)
n_c_fdr <- sum(kegg$fdr_sig)
cat("  KEGG terms:", nrow(kegg), "| FDR < 0.05:", n_c_fdr, "\n")
write.csv(kegg, file.path(PLOTTED, "SourceData_Fig7c_kegg_gsea.csv"),
          row.names = FALSE)

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
         fill = guide_colourbar(order = 1, barwidth = unit(2, "mm"),
                                barheight = unit(9, "mm"),
                                theme = theme(legend.title =
                                  element_text(size = 7.4, family = FAM)))) +
  labs(x = "Normalized enrichment score", y = NULL, tag = "c") +
  theme(axis.text.y = element_text(size = TINY, lineheight = 0.85),
        legend.position = "right", legend.spacing.y = unit(0.6, "mm")) +
  tag_theme

# -----------------------------------------------------------------------------
# d — IFN-γ / IL-17 axes, All cells only
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
mod <- read.csv(file.path(SRC, "module_score_summary.csv"), stringsAsFactors = FALSE)
n_mod_before <- nrow(mod)
id_cell <- read.csv(file.path(SRC, "identity_ifn_scores.csv"), stringsAsFactors = FALSE)
n_id_before <- nrow(id_cell)
id_cell <- id_cell[!is.na(id_cell$mod_ifn) & !is.na(id_cell$mod_il17r) &
                     !is.na(id_cell$mod_il17e), ]
n_id_after <- nrow(id_cell)

# All cells = every profiled microglia; Kept series omitted by request
n_mod_all <- sum(mod$set == "all")
n_mod_kept <- sum(mod$set == "kept")
mod_keep <- c("mod_ifn", "mod_il17r", "mod_il17e")
n_d_series_drop <- 1L

d_long <- id_cell %>%
  transmute(cluster = factor(as.character(cluster), levels = as.character(0:3)),
            group_show = to_posneg(group_label),
            ifn = mod_ifn,
            il17r = mod_il17r,
            il17e = mod_il17e) %>%
  pivot_longer(c(ifn, il17r, il17e), names_to = "module", values_to = "score")
mod_lab <- c(ifn = paste0("IFN-", intToUtf8(0x03B3)),
             il17r = "IL-17 receptor",
             il17e = "IL-17 effector")
d_long$module <- factor(mod_lab[d_long$module],
                        levels = unname(mod_lab))
d_sum <- d_long %>%
  group_by(cluster, group_show, module) %>%
  summarise(n = n(),
            mean = mean(score),
            se = sd(score) / sqrt(n()),
            .groups = "drop")
stopifnot(nrow(d_sum) == 3 * 4 * 2)

write.csv(mod, file.path(PLOTTED, "SourceData_Fig7d_module_score_summary_allsets.csv"),
          row.names = FALSE)
write.csv(d_sum, file.path(PLOTTED, "SourceData_Fig7d_IFN_IL17_allcells.csv"),
          row.names = FALSE)

pD <- ggplot(d_sum, aes(cluster, mean, fill = group_show)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_col(position = position_dodge(width = 0.72), width = 0.68,
           colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.72),
                width = 0.22, linewidth = 0.25) +
  facet_wrap(~ module, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, name = NULL, drop = FALSE) +
  labs(x = "Cluster", y = "Module score",
       title = paste0("IFN-", intToUtf8(0x03B3), " / IL-17 axes"),
       subtitle = "All cells",
       tag = "d") +
  theme(legend.position = "bottom",
        legend.margin = margin(-1, 0, 0, 0),
        strip.text = element_text(size = TINY, face = "bold"),
        plot.subtitle = element_text(size = TINY, hjust = 0.5)) +
  tag_theme

# -----------------------------------------------------------------------------
# e — leading-edge heatmap (previous d)
# -----------------------------------------------------------------------------
cat("Panel e ...\n")
df_e <- read.csv(file.path(SRC, "Fig6d_heatmap_z.csv"), stringsAsFactors = FALSE)
gene_ord <- read.csv(file.path(SRC, "Fig6D_gene_order.csv"), stringsAsFactors = FALSE)
pw_lv <- c("IL-17", "NOD-like", "Chemokine")
df_e$z <- ifelse(is.na(df_e$z), 0, df_e$z)
df_e$group_show <- to_posneg(df_e$Group)
df_e$Cluster <- factor(as.character(df_e$Cluster), levels = as.character(0:3))
df_e$pathway <- factor(as.character(df_e$pathway), levels = pw_lv)
df_e$Gene <- factor(df_e$Gene, levels = rev(gene_ord$gene))
stopifnot(!anyNA(df_e$pathway), !anyNA(df_e$group_show))
n_gene_pw <- table(distinct(df_e, Gene, pathway)$pathway)
write.csv(df_e, file.path(PLOTTED, "SourceData_Fig7e_heatmap_z.csv"),
          row.names = FALSE)

pE <- ggplot(df_e, aes(group_show, Gene, fill = z)) +
  geom_tile() +
  facet_grid(pathway ~ Cluster, scales = "free_y", space = "free_y",
             labeller = labeller(Cluster = function(x) paste0("C", x))) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "Row Z") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"),
                                barheight = unit(14, "mm"))) +
  labs(x = NULL, y = NULL, tag = "e") +
  theme(
    axis.text.x = element_text(size = TINY, face = "bold",
                               colour = unname(grp_cols[grp_lv])),
    axis.text.y = element_text(size = TINY, face = "italic", colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing.x = unit(0.7, "mm"),
    panel.spacing.y = unit(1.0, "mm"),
    strip.text.y = element_text(size = TINY, face = "bold", angle = -90),
    strip.text.x = element_text(size = TINY, face = "bold"),
    legend.position = "right"
  ) + tag_theme

# -----------------------------------------------------------------------------
# f — pseudotime (previous e)
# -----------------------------------------------------------------------------
cat("Panel f ...\n")
df_f <- read.csv(file.path(SRC, "Fig6E_cells.csv"), stringsAsFactors = FALSE)
edges <- read.csv(file.path(SRC, "Fig6E_edges.csv"), stringsAsFactors = FALSE)
n_f_before <- nrow(df_f)
df_f <- df_f[!is.na(df_f$Pseudotime) & !is.na(df_f$data_dim_1), ]
n_f_after <- nrow(df_f)
df_f$celltype <- factor(as.character(df_f$celltype), levels = as.character(0:3))
genes_f <- c("Cxcl10", "Id2")
stopifnot(all(genes_f %in% colnames(df_f)))
df_fg <- df_f %>%
  select(Pseudotime, celltype, all_of(genes_f)) %>%
  pivot_longer(all_of(genes_f), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes_f))
write.csv(df_f[, c("cell", "Pseudotime", "celltype", "data_dim_1", "data_dim_2",
                   "Cxcl10", "Id2")],
          file.path(PLOTTED, "SourceData_Fig7f_trajectory.csv"), row.names = FALSE)
write.csv(edges, file.path(PLOTTED, "SourceData_Fig7f_edges.csv"), row.names = FALSE)

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
  ggplot(df_f, aes(data_dim_1, data_dim_2)) +
    geom_point(aes(colour = Pseudotime), size = 0.35, stroke = 0) +
    scale_colour_viridis_c(option = "D", name = "Pseudotime",
                           breaks = pretty_breaks(3)) +
    guides(colour = guide_colourbar(barwidth = unit(1.8, "mm"),
                                    barheight = unit(8, "mm"))) +
    ggtitle("Pseudotime"))

p_ct <- traj_base(
  ggplot(df_f, aes(data_dim_1, data_dim_2)) +
    geom_point(aes(colour = celltype), size = 0.35, stroke = 0) +
    scale_colour_manual(values = clus_cols, name = "Cluster") +
    guides(colour = guide_legend(override.aes = list(size = 1.1))) +
    ggtitle("Cluster"))

p_ge <- ggplot(df_fg, aes(Pseudotime, expression)) +
  geom_point(aes(colour = celltype), size = 0.3, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.75, se = TRUE,
              linewidth = 0.45, colour = "black", fill = "grey55", alpha = 0.35) +
  scale_colour_manual(values = clus_cols, guide = "none") +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  labs(x = "Pseudotime", y = "Expression") +
  theme(strip.text = element_text(size = SMALL, face = "bold.italic"))

p_pt <- p_pt + labs(tag = "f") + tag_theme
pF <- (p_pt | p_ct) / p_ge + plot_layout(heights = c(1.2, 1))

# -----------------------------------------------------------------------------
# g — violins (previous f)
# -----------------------------------------------------------------------------
cat("Panel g ...\n")
df_g <- read.csv(file.path(SRC, "Fig6F_violin.csv"), stringsAsFactors = FALSE)
n_g_before <- nrow(df_g)
gene_lv <- c("P2ry12", "Hexb", "Il1b", "Cxcl2", "Cxcl10", "Id2")
df_g$gene <- factor(df_g$gene, levels = gene_lv)
df_g$seurat_clusters <- factor(as.character(df_g$seurat_clusters),
                               levels = as.character(0:3))
df_g$group_show <- to_posneg(df_g$group_label)
df_g <- df_g[!is.na(df_g$expression) & !is.na(df_g$gene), ]
n_g_after <- nrow(df_g)
stopifnot(!anyNA(df_g$gene), !anyNA(df_g$group_show))

stat_g <- df_g %>%
  group_by(gene, seurat_clusters) %>%
  summarise(
    n_pos = sum(group_show == "POS"),
    n_neg = sum(group_show == "NEG"),
    p = if (min(n_pos, n_neg) >= 3)
      suppressWarnings(wilcox.test(expression ~ group_show)$p.value) else NA_real_,
    ymax = max(expression, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(q = p.adjust(p, method = "BH"),
         star = case_when(is.na(q) ~ "", q >= 0.05 ~ "",
                          q < 1e-4 ~ "****", q < 1e-3 ~ "***",
                          q < 1e-2 ~ "**", TRUE ~ "*"))
facet_max <- df_g %>% group_by(gene) %>% summarise(gmax = max(expression), .groups = "drop")
stat_g <- stat_g %>% left_join(facet_max, by = "gene") %>% mutate(y = gmax * 1.04)
n_g_tests <- sum(!is.na(stat_g$p))
n_g_sig <- sum(stat_g$q < 0.05, na.rm = TRUE)
cat("  Wilcoxon tests:", n_g_tests, "| BH q < 0.05:", n_g_sig, "\n")
write.csv(df_g, file.path(PLOTTED, "SourceData_Fig7g_violin.csv"), row.names = FALSE)
write.csv(stat_g[, c("gene", "seurat_clusters", "n_pos", "n_neg", "p", "q", "star")],
          file.path(PLOTTED, "SourceData_Fig7g_wilcoxon_BH.csv"), row.names = FALSE)

pG <- ggplot(df_g, aes(seurat_clusters, expression, fill = group_show)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA,
              position = position_dodge(width = 0.78), width = 0.72) +
  stat_summary(fun = median, geom = "point", size = 0.28, colour = "white",
               position = position_dodge(width = 0.78)) +
  geom_text(data = stat_g, aes(seurat_clusters, y, label = star),
            inherit.aes = FALSE, size = 1.9, family = FAM, vjust = 0.4) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  labs(x = "Cluster", y = "Normalized expression", tag = "g") +
  theme(strip.text = element_text(size = SMALL, face = "bold.italic"),
        legend.position = "bottom",
        legend.margin = margin(-2, 0, 0, 0)) +
  tag_theme

# -----------------------------------------------------------------------------
# h — IF plate (previous g)
# -----------------------------------------------------------------------------
cat("Panel h ...\n")
if_png <- file.path(SRCFIG, "Fig6G_IF.png")
stopifnot(file.exists(if_png))
if_img <- png::readPNG(if_png)
n_h_px <- c(ncol(if_img), nrow(if_img))
pH <- ggplot() +
  annotation_custom(rasterGrob(if_img, interpolate = TRUE), -Inf, Inf, -Inf, Inf) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = BASE, base_family = FAM) +
  labs(tag = "h") +
  theme(aspect.ratio = nrow(if_img) / ncol(if_img),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5)) +
  tag_theme

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- pA + pB + plot_layout(widths = c(0.56, 0.44))
mid_l <- pC / pD + plot_layout(heights = c(0.92, 1.18))
mid <- mid_l | pE
mid <- mid + plot_layout(widths = c(0.52, 0.48))
row3 <- pF
row4 <- pG + pH + plot_layout(widths = c(0.56, 0.44))
fig <- row1 / mid / row3 / row4 +
  plot_layout(heights = c(40, 88, 72, 50))

# -----------------------------------------------------------------------------
# Stats / legend / provenance
# -----------------------------------------------------------------------------
n_tab <- df_a %>% count(group_show, cluster) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = n, values_fill = 0)

stats <- c(
  "Figure 7 statistics / integrity (Fig7_Nature_v1)",
  sprintf("Total microglia: %d (POS = %d, NEG = %d).",
          n_a_after, n_a_pos, n_a_neg),
  "Cells per group x cluster:",
  paste(capture.output(print(as.data.frame(n_tab), row.names = FALSE)), collapse = "\n"),
  sprintf("Panel a: before %d after NA drop %d; all retained", n_a_before, n_a_after),
  sprintf("Panel b: %d rows; columns sum to 100%% within group", n_b_after),
  sprintf("Panel c: %d KEGG terms; FDR < 0.05: %d", nrow(kegg), n_c_fdr),
  sprintf("Panel d module-score table: %d rows (%d all, %d kept); Kept series omitted by request (%d series dropped)",
          n_mod_before, n_mod_all, n_mod_kept, n_d_series_drop),
  sprintf("Panel d cell-level scores: before %d after NA drop %d; All cells plotted; bars = mean \u00b1 s.e.m.",
          n_id_before, n_id_after),
  sprintf("Panel e: %d genes (IL-17 %d, NOD-like %d, Chemokine %d); row Z clipped to [-2.5, 2.5]",
          nlevels(df_e$Gene), n_gene_pw["IL-17"], n_gene_pw["NOD-like"], n_gene_pw["Chemokine"]),
  sprintf("Panel f: %d cells, pseudotime %.1f-%.1f; LOESS span 0.75 with 95%% CI",
          n_f_after, min(df_f$Pseudotime), max(df_f$Pseudotime)),
  sprintf("Panel g: Wilcoxon POS vs NEG within gene x cluster; BH across %d tests; q < 0.05: %d",
          n_g_tests, n_g_sig),
  sprintf("Panel h IF raster: %d x %d px; no contrast change", n_h_px[1], n_h_px[2]),
  sprintf("Cluster 3 POS n = %d (descriptive)", n_a_c3_pos),
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
)
writeLines(stats, file.path(OUTDIR, "Figure7_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. 7 | Vessel-contacting microglia acquire an IFN-\u03B3 / IL-17-chemokine programme in EAE.",
  "",
  sprintf("a, UMAP of microglia from EAE spinal cord, split by vessel-contact status (POS, n = %s; NEG, n = %s) and coloured by transcriptional cluster. Cluster indices mark cluster medoids. All cells are plotted.",
          format(n_a_pos, big.mark = ","), format(n_a_neg, big.mark = ",")),
  "b, Cluster composition of each group; segments below 5% are unlabelled.",
  sprintf("c, KEGG gene-set enrichment of POS versus NEG microglia. Dot size, gene-set size; fill, -log10(P); black outline, Benjamini-Hochberg FDR < 0.05 (%d of %d terms).",
          n_c_fdr, nrow(kegg)),
  "d, Mean AddModuleScore (\u00b1 s.e.m.) for IFN-\u03B3 (Stat1, Ifngr1, Irf1, Cxcl9, Cxcl10, Gbp2, Isg15), IL-17 receptor (Il17ra, Il17rc, Traf3ip2) and IL-17 effector (Cxcl1, Cxcl2, Csf2, Il6, Lcn2) programmes, computed on all cells. The cell-filtered (Kept) series is omitted.",
  "e, Row Z-scored mean expression of leading-edge genes from the three highest-NES signalling pathways, by cluster (C0-C3) and group.",
  "f, Monocle pseudotime trajectory coloured by pseudotime and by cluster, and Cxcl10 and Id2 expression along pseudotime (LOESS fit, 95% confidence band).",
  sprintf("g, Homeostatic (P2ry12, Hexb) and inflammatory or chemokine (Il1b, Cxcl2, Cxcl10, Id2) gene expression by cluster and group. Violins are scaled to equal width; white dots, medians. Two-sided Wilcoxon rank-sum tests with Benjamini-Hochberg correction across %d comparisons; *q < 0.05, **q < 0.01, ***q < 0.001, ****q < 0.0001.",
          n_g_tests),
  "h, Representative immunofluorescence of the vessel-contact niche (DAPI, laminin and Cxcl10) with boxed regions enlarged below.",
  "",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "Figure7_legend.txt"))

writeLines(c(
  "Figure 7 | new folder Fig7_Nature_v1 (does not overwrite earlier packages).",
  "Reuse level: structural adaptation of the microglia Fig. 6 a-g plus All-cells IFN/IL-17 module scores.",
  "Panels a-c: same observations as the microglia Fig. 6 a-c.",
  "Panel d: All-cells module scores; Kept series omitted by request.",
  "Panels e-h: same observations as the microglia Fig. 6 d-g (heatmap, pseudotime, violins, IF).",
  "Panel h: IF raster embedded at native aspect; no resampling or channel adjustment.",
  "Backend: R exclusive."
), file.path(OUTDIR, "Figure7_provenance.txt"))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
base <- file.path(OUTDIR, "Figure7")

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
  "Build notes (full visual QA is in Figure7_QA.txt)",
  "- Backend: R exclusive",
  "- a-g: native ggplot vector (UMAP/trajectory marks are vector points)",
  "- d: All-cells IFN-\u03B3 / IL-17 scores; Kept omitted",
  "- h: supplied IF raster, native aspect, no contrast change",
  sprintf("- Exports: PDF, SVG, PNG, TIFF at %d x %d mm, raster 600 dpi",
          FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "Figure7_build_notes.txt"))

cat("DONE\n")
