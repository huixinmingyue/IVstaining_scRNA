#!/usr/bin/env Rscript
# =============================================================================
# 04  Figures for hard-gated CD4_T reannotation (R exclusive)
# Writes into FigS1b_reannotate_ImmGen_v1 only. Does not overwrite old S1/Fig1.
# =============================================================================
options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(scales)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS1b_reannotate_ImmGen_v1")
SRC <- file.path(OUT, "source_data")
DEGDIR <- file.path(OUT, "DEG_tables")
SDP <- file.path(OUT, "source_data_plotted")
QA <- file.path(OUT, "qa")
dir.create(SDP, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FAM <- "Arial"
BASE <- 6.5
FIG_W <- 183
COL <- list(
  true_CD4 = "#2F6F8F", true_CD8 = "#5B8F4E", NKT = "#B2473E",
  DP = "#C4A35A", unassigned = "#9A9A9A", doublet_c2 = "#5A5A5A",
  IntraNeg = "#2F6F8F", IntraPos = "#B2473E", NS = "#B0B0B0"
)
CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
CLASS_LAB <- c(true_CD4 = "True CD4", true_CD8 = "True CD8", NKT = "NKT",
               DP = "DP", unassigned = "Unassigned", doublet_c2 = "Doublet (c2)")

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) cat("font:", FAM, "->", fm$path, "\n")

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = BASE),
      axis.text = element_text(size = BASE - 0.5, colour = "black"),
      plot.title = element_text(size = BASE + 0.4, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid = element_blank(),
      legend.key.size = unit(2.4, "mm"),
      legend.title = element_text(size = BASE - 0.3),
      legend.text = element_text(size = BASE - 0.7)
    )
)

save_pub <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = FAM)
  print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in",
                 res = dpi, compression = "lzw")
  print(plot); dev.off()
}

md <- read.csv(file.path(SRC, "cell_metadata_reannotated.csv"), stringsAsFactors = FALSE)
md$t_class <- factor(md$t_class, levels = CLASS_LV)
md$IV_group <- factor(md$IV_group, levels = c("IntraNeg", "IntraPos"))
comp <- read.csv(file.path(SRC, "n_by_class_IV.csv"), stringsAsFactors = FALSE)
comp$t_class <- factor(comp$t_class, levels = CLASS_LV)
comp$IV_lab <- factor(ifelse(comp$IV_group == "IntraPos", "IV+", "IV-"),
                      levels = c("IV-", "IV+"))
fish <- read.csv(file.path(SRC, "fisher_class_IV.csv"), stringsAsFactors = FALSE)
summ <- read.csv(file.path(SRC, "DEG_run_summary.csv"), stringsAsFactors = FALSE)
hom <- read.csv(file.path(SRC, "homing_genes_raw.csv"), stringsAsFactors = FALSE)

# ----- a UMAP (all 6,098 cells) -----
p_a <- ggplot(md, aes(UMAP_1, UMAP_2, colour = t_class)) +
  geom_point(size = 0.18, stroke = 0, alpha = 0.75, shape = 16) +
  scale_colour_manual(values = unlist(COL[CLASS_LV]), labels = CLASS_LAB, name = NULL) +
  coord_equal() +
  labs(title = "Hard-gated identity", tag = "a", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = BASE - 0.8),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 2))

# ----- b composition -----
p_b <- ggplot(comp, aes(IV_lab, pct, fill = t_class)) +
  geom_col(width = 0.72, colour = NA) +
  scale_fill_manual(values = unlist(COL[CLASS_LV]), labels = CLASS_LAB, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Compartment composition", tag = "b",
       x = NULL, y = "% of CD4_T object") +
  theme(legend.position = "none")

# ----- c DEG collapse -----
deg_n <- summ %>%
  filter(label %in% c("original_CD4_T", "true_CD4", "true_CD8", "NKT",
                      "trueCD4_1a", "trueCD4_1b")) %>%
  transmute(
    label = recode(label,
                   original_CD4_T = "Original CD4_T",
                   true_CD4 = "True CD4",
                   true_CD8 = "True CD8",
                   NKT = "NKT",
                   trueCD4_1a = "True CD4 1a",
                   trueCD4_1b = "True CD4 1b"),
    IntraPos = n_IntraPos_up_filt,
    IntraNeg = n_IntraNeg_up_filt
  ) %>%
  pivot_longer(c(IntraPos, IntraNeg), names_to = "direction", values_to = "n")
deg_n$label <- factor(deg_n$label, levels = c("Original CD4_T", "True CD4", "True CD8",
                                              "NKT", "True CD4 1a", "True CD4 1b"))
p_c <- ggplot(deg_n, aes(label, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68, colour = NA) +
  scale_fill_manual(values = c(IntraPos = COL$IntraPos, IntraNeg = COL$IntraNeg),
                    labels = c(IntraPos = "IV+ up", IntraNeg = "IV- up"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Intra DEG count (tech-filtered)", tag = "c",
       x = NULL, y = "Genes |log2FC|\u22650.5, padj<0.001") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

# ----- volcano helper -----
prep_deg <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$neglog10_padj <- pmin(-log10(pmax(d$p_val_adj, 1e-300)), 50)
  d$direction <- factor(d$direction, levels = c("IntraNeg up", "NS", "IntraPos up"))
  d
}
volcano <- function(d, title, tag, n_lab = 12) {
  lab <- d %>%
    filter(direction != "NS", !tech_or_ribo) %>%
    group_by(direction) %>%
    slice_min(order_by = p_val_adj, n = n_lab, with_ties = FALSE) %>%
    ungroup()
  xmax <- max(2.8, as.numeric(quantile(abs(d$avg_log2FC), 0.995, na.rm = TRUE)) * 1.15)
  ymax <- max(8, max(d$neglog10_padj, na.rm = TRUE) * 1.05)
  ggplot(d, aes(avg_log2FC, neglog10_padj, colour = direction)) +
    geom_point(size = 0.45, stroke = 0, alpha = 0.7, shape = 16) +
    geom_vline(xintercept = c(-0.5, 0.5), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.001), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(
      data = lab, aes(label = gene), size = 1.85, fontface = "italic",
      family = FAM, max.overlaps = 40, min.segment.length = 0.05,
      segment.size = 0.15, show.legend = FALSE, seed = 1
    ) +
    scale_colour_manual(values = c("IntraNeg up" = COL$IntraNeg, NS = COL$NS,
                                   "IntraPos up" = COL$IntraPos), drop = FALSE, name = NULL) +
    coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), clip = "off") +
    labs(title = title, tag = tag,
         x = expression(log[2]*"FC (IV+ / IV-)"),
         y = expression(-log[10]*" adjusted "*italic(P))) +
    theme(legend.position = "none")
}

d_cd4 <- prep_deg(file.path(DEGDIR, "DEG_true_CD4_IntraPos_vs_IntraNeg.csv"))
d_cd8 <- prep_deg(file.path(DEGDIR, "DEG_true_CD8_IntraPos_vs_IntraNeg.csv"))
p_d <- volcano(d_cd4, "True CD4", "d", 12)
p_e <- volcano(d_cd8, "True CD8", "e", 10)

# ----- f homing: always report the six genes -----
H6 <- c("Ccr7", "Sell", "S1pr1", "Cd69", "Cxcr5", "Foxp3")
h <- hom %>% filter(label == "true_CD4", gene %in% H6)
h$gene <- factor(h$gene, levels = rev(H6))
h$sig <- ifelse(!is.na(h$p_val) & h$p_val < 0.05, "Uncorrected P < 0.05", "NS")
h$pct_lab <- sprintf("%.0f / %.0f", 100 * h$pct.1, 100 * h$pct.2)
write.csv(h, file.path(SDP, "SourceData_homing_trueCD4.csv"), row.names = FALSE)
p_f <- ggplot(h, aes(avg_log2FC, gene, colour = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(size = 2.2, stroke = 0) +
  geom_text(aes(label = pct_lab), hjust = ifelse(h$avg_log2FC >= 0, -0.25, 1.25),
            size = 1.85, family = FAM, colour = "grey20") +
  scale_colour_manual(values = c("Uncorrected P < 0.05" = COL$IntraPos, NS = "#555555"),
                      name = NULL) +
  scale_x_continuous(limits = c(-2.4, 1.4), expand = expansion(mult = 0.02)) +
  labs(title = "Homing genes, true CD4", tag = "f",
       x = expression(log[2]*"FC (IV+ / IV-)  (labels: %IV+ / %IV-)"), y = NULL) +
  theme(legend.position = "bottom") +
  coord_cartesian(clip = "off")

fig1 <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_layout(heights = c(1, 1.15))
cat("Saving Figure_reannotate ...\n")
save_pub(fig1, file.path(OUT, "Figure_reannotate"), FIG_W, 210)

# ----- gate diagnostic (Fig. 1f diagonal should vanish) -----
det_long <- md %>%
  select(t_class, det_Cd4, det_Cd8a, det_Cd8b1, det_Klrb1c, det_Cd3e, det_Ms4a1, det_Ncr1) %>%
  pivot_longer(-t_class, names_to = "gene", values_to = "det") %>%
  mutate(gene = sub("^det_", "", gene)) %>%
  group_by(t_class, gene) %>%
  summarise(pct = 100 * mean(det > 0), n = n(), .groups = "drop")
det_long$t_class <- factor(det_long$t_class, levels = CLASS_LV)
det_long$gene <- factor(det_long$gene, levels = c("Cd4", "Cd8a", "Cd8b1", "Klrb1c",
                                                   "Cd3e", "Ncr1", "Ms4a1"))
write.csv(det_long, file.path(SDP, "SourceData_lineage_detection.csv"), row.names = FALSE)
p_gate <- ggplot(det_long, aes(gene, t_class, fill = pct)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f", pct)), size = 1.85, family = FAM, colour = "grey15") +
  scale_fill_gradient(low = "white", high = "#2F6F8F", name = "% detected",
                      limits = c(0, 100)) +
  scale_y_discrete(labels = CLASS_LAB) +
  labs(title = "RNA detection after the hard gate", tag = "a",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank())
save_pub(p_gate, file.path(OUT, "Figure_gate_detection"), 120, 70)

# ----- ImmGen ggplot heatmap: true CD4 genes only -----
z <- as.matrix(read.csv(file.path(SRC, "ImmGen_zscore_matrix.csv"),
                        row.names = 1, check.names = FALSE))
lab4 <- read.csv(file.path(SRC, "S1b_true_CD4_sig_labelled.csv"), stringsAsFactors = FALSE)
genes4 <- intersect(lab4$gene, rownames(z))
z4 <- z[genes4, , drop = FALSE]
dir4 <- lab4$direction[match(rownames(z4), lab4$gene)]
fc4 <- lab4$avg_log2FC[match(rownames(z4), lab4$gene)]
ord <- order(dir4, ifelse(dir4 == "IntraNeg up", fc4, -fc4))
z4 <- z4[ord, , drop = FALSE]
dir4 <- dir4[ord]
col_lab <- c(CD4_naive = "Naive CD4", CD4_total = "Total CD4", Treg = "Treg",
             CD8_naive = "Naive CD8", NKT = "NKT", Tgd = "gdT",
             NK_immature = "NK immature", NK_int = "NK int", NK_mature = "NK mature")
z_long <- as.data.frame(z4, check.names = FALSE)
z_long$gene <- factor(rownames(z4), levels = rev(rownames(z4)))
z_long$direction <- dir4
z_long <- z_long %>%
  pivot_longer(all_of(colnames(z4)), names_to = "subset", values_to = "z")
z_long$subset <- factor(z_long$subset, levels = colnames(z4),
                        labels = unname(col_lab[colnames(z4)]))
z_long$direction <- factor(z_long$direction, levels = c("IntraNeg up", "IntraPos up"),
                           labels = c("IV- up", "IV+ up"))
write.csv(z_long, file.path(SDP, "SourceData_ImmGen_trueCD4.csv"), row.names = FALSE)

p_ht <- ggplot(z_long, aes(subset, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2F6F8F", mid = "white", high = "#B2473E",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                       name = "z") +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  labs(title = "ImmGen spleen T/NK (true-CD4 Intra DEGs)", tag = "a",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 0.6),
        axis.text.y = element_text(face = "italic", size = 5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = BASE - 0.4, face = "bold"),
        legend.position = "right")

res <- deg_n %>% filter(label %in% c("True CD4", "True CD4 1a", "True CD4 1b"))
p_res <- ggplot(res, aes(label, n, fill = direction)) +
  geom_col(position = position_dodge(0.72), width = 0.68, colour = NA) +
  scale_fill_manual(values = c(IntraPos = COL$IntraPos, IntraNeg = COL$IntraNeg),
                    labels = c(IntraPos = "IV+ up", IntraNeg = "IV- up"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Residual DE inside true CD4", tag = "b",
       x = NULL, y = "Tech-filtered Intra DEGs") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

# spatial maps: all bins, y flipped to match Fig. 4
spat_f <- file.path(SRC, "spatial_module_scores.csv")
has_spat <- file.exists(spat_f)
if (has_spat) {
  sp <- read.csv(spat_f, stringsAsFactors = FALSE)
  sp$xp <- as.numeric(sp$x)
  sp$yp <- -as.numeric(sp$y)
  clip01 <- function(v) {
    q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
    pmin(pmax(v, q[1]), q[2])
  }
  sp$TfhTreg_c <- clip01(sp$TfhTreg)
  sp$Cytotoxic_c <- clip01(sp$Cytotoxic)
  xlm <- range(sp$xp, na.rm = TRUE)
  ylm <- range(sp$yp, na.rm = TRUE)
  sb_len <- 1000
  sb_x1 <- xlm[2] - 0.06 * diff(xlm)
  sb_x0 <- sb_x1 - sb_len
  sb_y <- ylm[1] + 0.045 * diff(ylm)
  spat_geom <- if (requireNamespace("ggrastr", quietly = TRUE)) {
    function(...) ggrastr::geom_point_rast(..., raster.dpi = 300)
  } else {
    function(...) geom_point(...)
  }
  spat_map <- function(df, zcol, title, tag, pal) {
    ggplot(df, aes(xp, yp, colour = .data[[zcol]])) +
      spat_geom(size = 0.12, stroke = 0, shape = 16) +
      annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
               linewidth = 0.55, colour = "grey20") +
      annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y - 0.03 * diff(ylm),
               label = "1 mm", size = 5 / .pt, family = FAM, colour = "grey20") +
      scale_colour_gradientn(colours = pal, name = "Score") +
      coord_fixed(xlim = xlm, ylim = ylm, expand = FALSE) +
      labs(title = title, tag = tag, x = NULL, y = NULL) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "right",
            legend.key.height = unit(8, "mm"),
            panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.25))
  }
  pal_tfh <- c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
  pal_cyt <- c("#F7F7F7", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704")
  p_c_sp <- spat_map(sp, "TfhTreg_c", "Tfh/Treg module", "c", pal_tfh)
  p_d_sp <- spat_map(sp, "Cytotoxic_c", "Cytotoxic module", "d", pal_cyt)

  if ("shown_in_g" %in% colnames(sp)) {
    sp$niche <- ifelse(!is.na(sp$shown_in_g) & as.logical(sp$shown_in_g),
                       "B/T follicular", "Rest of section")
    sp$cd4lab <- ifelse("in_CD4_niche" %in% colnames(sp) & sp$in_CD4_niche %in% TRUE,
                        "CD4 bins", "Other bins")
    sp_box <- bind_rows(
      sp %>% mutate(split = "Section niche") %>%
        transmute(split, group = niche, TfhTreg, Cytotoxic),
      sp %>% mutate(split = "CD4 prediction") %>%
        transmute(split, group = cd4lab, TfhTreg, Cytotoxic)
    ) %>%
      pivot_longer(c(TfhTreg, Cytotoxic), names_to = "module", values_to = "score")
    sp_box$module <- recode(sp_box$module, TfhTreg = "Tfh/Treg", Cytotoxic = "Cytotoxic")
    sp_box$split <- factor(sp_box$split, levels = c("Section niche", "CD4 prediction"))
    write.csv(sp_box, file.path(SDP, "SourceData_spatial_module_by_niche.csv"),
              row.names = FALSE)
    p_e_sp <- ggplot(sp_box, aes(group, score, fill = group)) +
      geom_boxplot(outlier.size = 0.12, outlier.alpha = 0.12, linewidth = 0.3,
                   width = 0.62) +
      facet_grid(module ~ split, scales = "free") +
      scale_fill_manual(values = c("B/T follicular" = "#6BAED6",
                                   "Rest of section" = "#FD8D3C",
                                   "CD4 bins" = "#2F6F8F",
                                   "Other bins" = "#B0B0B0"), guide = "none") +
      labs(title = "Module scores by niche", tag = "e",
           x = NULL, y = "AddModuleScore") +
      theme(axis.text.x = element_text(angle = 18, hjust = 1),
            strip.text = element_text(size = BASE - 0.6, face = "bold"))
  } else {
    p_e_sp <- ggplot() + labs(title = "Niche scores unavailable", tag = "e")
  }
} else {
  p_c_sp <- ggplot() + labs(title = "Spatial Tfh/Treg (not run)", tag = "c")
  p_d_sp <- ggplot() + labs(title = "Spatial cytotoxic (not run)", tag = "d")
  p_e_sp <- ggplot() + labs(title = "Spatial niche (not run)", tag = "e")
}

fig2 <- (p_ht | (p_res / p_e_sp)) / (p_c_sp | p_d_sp) +
  plot_layout(heights = c(1.35, 1))
cat("Saving Figure_ImmGen_spatial ...\n")
save_pub(fig2, file.path(OUT, "Figure_ImmGen_spatial"), FIG_W, 240)

# plotted source extracts
write.csv(comp, file.path(SDP, "SourceData_class_composition.csv"), row.names = FALSE)
write.csv(deg_n, file.path(SDP, "SourceData_DEG_counts.csv"), row.names = FALSE)
write.csv(fish, file.path(SDP, "SourceData_fisher_OR.csv"), row.names = FALSE)

lg1 <- paste0(
  "Fig. Sx | Hard-gated reannotation of the original CD4 T object.\n",
  "a, UMAP of 6,098 cells from the original CD4_T annotation, coloured by an RNA hard gate ",
  "(counts > 0): true CD4 (Cd4+ Cd8a- Cd8b1- Klrb1c-; n = 1,185), true CD8 (Cd8a+ or Cd8b1+, Cd4- Klrb1c-; n = 2,887), ",
  "NKT (Klrb1c+; n = 473), DP (Cd4+ and Cd8a/Cd8b1+), unassigned (lineage genes undetected), ",
  "and cluster-2 doublets (removed from DEG). b, Class composition within each IV compartment ",
  "(one pooled sample per compartment). True CD4 is enriched in IV- (Fisher OR = 0.31, Holm P = 2.0e-61); ",
  "NKT in IV+ (OR = 3.39, Holm P = 4.9e-34); true CD8 in IV+ (OR = 1.60, Holm P = 5.4e-19). ",
  "c, Number of IntraPos- versus IntraNeg-up genes after removing ",
  "ribosomal, mitochondrial, Gm* and shared technical genes (Lars2, Cmss1, Gphn, Cdk8, Camk1d). ",
  "The original mixed CD4_T set collapses from 335 to 54 genes in bona fide CD4. ",
  "d,e, Wilcoxon volcanoes (IV+ versus IV-) in true CD4 (n = 925 IV- / 260 IV+) and true CD8 ",
  "(n = 1,477 / 1,410); labelled genes are the top 12 / 10 per direction after the tech filter. ",
  "f, Pre-specified homing-axis genes in true CD4, reported regardless of genome-wide significance. ",
  "Labels are detection rates in IV+ / IV-. Uncorrected P < 0.05 is marked because Bonferroni-adjusted P ",
  "is conservative for a six-gene set. Each IV compartment is one pooled sample; Milo/scCODA were not run. ",
  "Source data are provided as a Source Data file."
)
lg2 <- paste0(
  "Fig. Sy | ImmGen identity overlay and Stereo-seq module maps.\n",
  "a, ImmGen ULI (GSE109125) spleen T/NK z-scores for true-CD4 Intra DEGs. Columns include ",
  "naive CD4, total CD4, Treg, naive CD8, NKT, gdT and NK states. This table has no Tfh or memory CD8 reference. ",
  "Rows are split by IV direction and ordered by log2FC. b, Residual Intra DEG counts inside true-CD4 clusters 1a and 1b ",
  "after the hard gate. c,d, Stereo-seq AddModuleScore maps for a Tfh/Treg gene set ",
  "(Izumo1r, Tox, Cd27, Ctla4, Pdcd1, Icos, Maf, Btla, Foxp3, Cxcr5, Slamf6, Il6st) versus ",
  "a cytotoxic set (Gzma, Gzmb, Ccl5, Nkg7, Klrk1, Klrb1c, Klrd1, Id2, Zeb2, Ly6c2, Cxcr6, Ifng, Prf1). ",
  "Colour scales are clipped to the 1st-99th percentile of each module. Scale bar, 1 mm. ",
  "e, Module scores in Fig. 4 B/T follicular bins (Louvain 1, 4, 5, 8) versus the rest of the section, ",
  "and in CD4-predicted bins versus other bins. Boxes show median and IQR. ",
  "Both modules are higher in CD4 bins than in the remainder of the section; ",
  "relative to the B/T follicular mask, both scores are higher outside that mask (extrafollicular). ",
  "This does not recover a clean Tfh-in-PALS versus cytotoxic-in-red-pulp split at bin resolution. ",
  "Source data are provided as a Source Data file."
)
writeLines(lg1, file.path(OUT, "Figure_reannotate_legend.txt"))
writeLines(lg2, file.path(OUT, "Figure_ImmGen_spatial_legend.txt"))
cat("DONE 04 figures\n")
