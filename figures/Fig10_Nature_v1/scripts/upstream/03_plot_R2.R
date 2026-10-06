#!/usr/bin/env Rscript
# Nature-format reanalysis figures. R exclusive. New folder only.
options(warn = 1)
set.seed(42)

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig9_S8S9_depth_decontam_reanalysis")
SRC <- file.path(OUT, "source_data")
FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

FAM <- "Arial"
BASE <- 6.5
if (!FAM %in% systemfonts::system_fonts()$family) {
  stop("Arial not registered with systemfonts")
}
cowplot::set_null_device("agg")

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
      plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = BASE - 0.6, colour = "grey30"),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = BASE - 0.2, face = "bold")
    )
)

ARM <- c(`NEG-dominant` = "#E67E22", `POS-dominant` = "#9B59B6",
         NEG = "#E67E22", POS = "#9B59B6")
CL <- c(`0` = "#0072B2", `3` = "#009E73", `5` = "#A6301F")
GENE_ORD <- c("Lgals3", "Mki67", "Cxcl10", "P2ry12", "Prdx1")

save_pub <- function(plot, filename, width_mm = 183, height_mm = 160) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = FAM)
  print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 600)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(plot); dev.off()
}

fmt_q <- function(q) {
  ifelse(q < 1e-300, "q<1e-300",
         ifelse(q < 0.001, sprintf("q=%.0e", q), sprintf("q=%.3f", q)))
}

# =============================================================================
# Figure R1 — spatial depth + S9 all-spot + Cxcl10
# =============================================================================
d3 <- read.csv(file.path(SRC, "fig9_d3_depth_spots.csv"), stringsAsFactors = FALSE)
d3$dominance <- factor(d3$dominance, c("NEG-dominant", "POS-dominant"))
t9 <- read.csv(file.path(SRC, "fig9_d3_depth_tests.csv"), stringsAsFactors = FALSE)
sp_all <- read.csv(file.path(SRC, "spatial_depth_tests.csv"), stringsAsFactors = FALSE)
cx9 <- read.csv(file.path(SRC, "fig9_d3_cxcl10_partners.csv"), stringsAsFactors = FALSE)
ncs <- read.csv(file.path(SRC, "spatial_ncount_by_dominance.csv"), stringsAsFactors = FALSE)
win9 <- read.csv(file.path(SRC, "fig9_d3_cxcl10_windows.csv"), stringsAsFactors = FALSE)

pA <- ggplot(d3, aes(dominance, nCount, fill = dominance)) +
  geom_violin(colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.size = 0.15, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = ARM, guide = "none") +
  scale_y_continuous(trans = "log10", labels = label_comma()) +
  labs(x = NULL, y = "nCount (Spatial)", tag = "a",
       title = "Fig. 9 D3 library size",
       subtitle = sprintf("POS median %s vs NEG %s",
                          format(median(d3$nCount[d3$dominance == "POS-dominant"]), big.mark = ","),
                          format(median(d3$nCount[d3$dominance == "NEG-dominant"]), big.mark = ",")))

det <- t9 %>% filter(norm == "count") %>%
  mutate(gene = factor(gene, GENE_ORD)) %>%
  pivot_longer(c(det_pct_neg, det_pct_pos), names_to = "arm", values_to = "pct") %>%
  mutate(arm = recode(arm, det_pct_neg = "NEG-dominant", det_pct_pos = "POS-dominant"))
pB <- ggplot(det, aes(gene, pct, fill = arm)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  scale_fill_manual(values = ARM, name = NULL) +
  labs(x = NULL, y = "Spots with ≥1 count (%)", tag = "b",
       title = "Detection (raw counts)",
       subtitle = "Same split as Fig. 9h; not depth-corrected") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

sp9 <- t9 %>%
  filter(norm %in% c("count", "logcpm"), gene %in% GENE_ORD) %>%
  mutate(gene = factor(gene, GENE_ORD),
         norm = recode(norm, count = "Raw counts", logcpm = "logCPM"))
pC <- ggplot(sp9, aes(spearman_rho, gene, colour = norm)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbarh(aes(xmin = spearman_lo, xmax = spearman_hi),
                 height = 0.25, linewidth = 0.4,
                 position = position_dodge(0.55)) +
  geom_point(size = 1.6, position = position_dodge(0.55)) +
  scale_colour_manual(values = c("Raw counts" = "#767676", "logCPM" = "#A6301F"), name = NULL) +
  labs(x = "Spearman ρ vs POS fraction (bootstrap 95% CI)", y = NULL, tag = "c",
       title = "Fig. 9 D3, all 2,573 spots",
       subtitle = "Prdx1 reverses after logCPM; P2ry12 shrinks") +
  theme(legend.position = "top")

s9sp <- sp_all %>%
  filter(subset == "all_spots", norm == "sct",
         section %in% c("S9_section1_immune2", "S9_section2_immune4")) %>%
  mutate(
    gene = factor(gene, GENE_ORD),
    section = recode(section,
                     S9_section1_immune2 = "S9 section 1",
                     S9_section2_immune4 = "S9 section 2")
  )
pDE <- ggplot(s9sp, aes(spearman_rho, gene, colour = section)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbarh(aes(xmin = spearman_lo, xmax = spearman_hi),
                 height = 0.25, linewidth = 0.4,
                 position = position_dodge(0.55)) +
  geom_point(size = 1.6, position = position_dodge(0.55)) +
  scale_colour_manual(values = c("S9 section 1" = "#0072B2", "S9 section 2" = "#E67E22"),
                      name = NULL) +
  labs(x = "Spearman ρ vs Seurat POS fraction (SCT)", y = NULL, tag = "d",
       title = "S9 all spots, no top-100",
       subtitle = "Section 2 does not repeat the POS association") +
  theme(legend.position = "top")

cx <- cx9 %>% filter(norm == "logcpm") %>%
  mutate(partner = factor(partner, rev(c("Ccr2", "Plac8", "Ms4a7", "Gfap", "Aqp4",
                                         "Lgals3", "Mki67", "P2ry12", "Hexb", "Tmem119"))))
pF <- ggplot(cx, aes(rho, partner)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.35,
                 colour = "#3182BD") +
  geom_point(size = 1.5, colour = "#3182BD") +
  labs(x = "Spearman ρ vs Cxcl10 (logCPM)", y = NULL, tag = "e",
       title = "Fig. 9 D3 Cxcl10 partners",
       subtitle = "All |ρ| < 0.12; monocyte ≈ astrocyte ≈ microglia")

win9$class <- recode(win9$class,
                     astro_top10 = "Astrocyte",
                     mono_top10 = "Monocyte-like",
                     mg_homeo_top10 = "Homeostatic Mg")
pG <- ggplot(win9, aes(class, cx_pct)) +
  geom_col(fill = "#3182BD", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\nn=%d", cx_pct, n_hi)),
            vjust = -0.15, size = 1.8, family = FAM) +
  coord_cartesian(ylim = c(0, max(win9$cx_pct, na.rm = TRUE) * 1.25)) +
  labs(x = NULL, y = "Cxcl10+ in window (%)", tag = "f",
       title = "Cxcl10+ rate, lineage windows",
       subtitle = "Top 10% among spots with that lineage score > 0")

R1 <- (pA | pB) / (pC | pDE) / (pF | pG) +
  plot_layout(heights = c(1, 1.15, 1.1))
save_pub(R1, file.path(FIG, "Figure_R1_spatial_depth_S9_Cxcl10"), 183, 210)
cat("wrote Figure_R1\n")

# =============================================================================
# Figure R2 — identity, if source files exist
# =============================================================================
need <- c("composition_arm_before_after.csv", "contamination_by_cluster_contact.csv",
          "DEG_locked_genes.csv", "module_cluster_by_contact.csv",
          "module_POS_vs_NEG_tests.csv")
if (!all(file.exists(file.path(SRC, need)))) {
  cat("R2 skipped; scRNA tables not ready\n")
  quit(save = "no", status = 0)
}

arm <- read.csv(file.path(SRC, "composition_arm_before_after.csv"), stringsAsFactors = FALSE)
arm$cluster <- factor(arm$cluster, c("0", "3", "5"))
arm$contact <- factor(arm$contact, c("NEG", "POS"))
arm$cleaned <- factor(ifelse(arm$cleaned, "After cleaning", "Before cleaning"),
                      c("Before cleaning", "After cleaning"))
p2a <- ggplot(arm, aes(contact, pct, fill = cluster)) +
  geom_col(width = 0.7) +
  facet_wrap(~ cleaned) +
  scale_fill_manual(values = CL, name = "Cluster") +
  labs(x = NULL, y = "Microglia (%)", tag = "a",
       title = "State composition", subtitle = "Drop doublet or neu_high or bam_high")

enr <- read.csv(file.path(SRC, "contamination_by_cluster_contact.csv"), stringsAsFactors = FALSE)
enr$cluster <- factor(enr$cluster, c("0", "3", "5"))
enr$contact <- factor(enr$contact, c("NEG", "POS"))
enr_l <- enr %>%
  pivot_longer(c(pct_dbl, pct_neu, pct_bam), names_to = "flag", values_to = "pct") %>%
  mutate(flag = recode(flag, pct_dbl = "Doublet", pct_neu = "Neutrophil-high",
                       pct_bam = "BAM/monocyte-high"))
p2b <- ggplot(enr_l, aes(cluster, pct, fill = contact)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  facet_wrap(~ flag) +
  scale_fill_manual(values = ARM, name = NULL) +
  labs(x = "Cluster", y = "Flagged cells (%)", tag = "b",
       title = "Contamination flags in microglia")

lock <- read.csv(file.path(SRC, "DEG_locked_genes.csv"), stringsAsFactors = FALSE)
lock <- lock %>%
  filter(cleaned, grepl("^within_c", contrast),
         gene %in% c("Lgals3", "Ifit1", "Ifit3", "Isg15", "Cxcl10", "Spp1",
                     "P2ry12", "Hexb", "S100a8", "Ms4a7", "Mki67", "Trem2")) %>%
  mutate(
    cluster = sub("within_c([035])_POS_vs_NEG", "\\1", contrast),
    gene = factor(gene, c("Lgals3", "Spp1", "Ifit1", "Ifit3", "Isg15", "Cxcl10",
                          "Mki67", "Trem2", "P2ry12", "Hexb", "S100a8", "Ms4a7")),
    cluster = factor(cluster, c("0", "3", "5"))
  )
p2c <- ggplot(lock, aes(avg_log2FC, gene, colour = cluster)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(shape = sig), size = 1.7,
             position = position_dodge(0.55)) +
  scale_colour_manual(values = CL, name = "Cluster") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "padj<0.05\n|lfc|≥0.5") +
  labs(x = "log2FC (POS − NEG), cleaned", y = NULL, tag = "c",
       title = "Within-cluster POS vs NEG",
       subtitle = "Selection if empty; induction if cluster 0 already shifts")

mod <- read.csv(file.path(SRC, "module_cluster_by_contact.csv"), stringsAsFactors = FALSE)
mod <- mod %>%
  filter(set %in% c("DAM", "DAM_core", "ICAM", "SAMC", "IFN", "G2M")) %>%
  mutate(cluster = factor(as.character(cluster), c("0", "3", "5")),
         contact = factor(contact, c("NEG", "POS")),
         set = recode(set, DAM = "DAM-like", DAM_core = "DAM core"),
         set = factor(set, c("DAM-like", "DAM core", "ICAM", "SAMC", "IFN", "G2M")))
p2d <- ggplot(mod, aes(contact, set, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", mean)), size = 1.7, colour = "black") +
  facet_wrap(~ cluster, nrow = 1) +
  scale_fill_gradient2(low = "#0072B2", mid = "#F5F5F5", high = "#A6301F",
                       name = "Mean z") +
  labs(x = NULL, y = NULL, tag = "d",
       title = "Modules by cluster × contact",
       subtitle = "DAM core = Clec7a, Itgax, Cst7, Trem2")

ht <- read.csv(file.path(SRC, "module_POS_vs_NEG_tests.csv"), stringsAsFactors = FALSE)
ht <- ht %>%
  filter(set %in% c("H_glyc", "H_oxph", "H_ifna", "H_g2m", "H_glyc_minus_oxph",
                    "DAM", "DAM_core", "ICAM", "SAMC", "IFN", "G2M")) %>%
  mutate(set = recode(set,
                      H_glyc = "Hallmark glycolysis",
                      H_oxph = "Hallmark OXPHOS",
                      H_ifna = "Hallmark IFN-α",
                      H_g2m = "Hallmark G2M",
                      H_glyc_minus_oxph = "Glycolysis − OXPHOS",
                      DAM = "DAM-like", DAM_core = "DAM core",
                      ICAM = "ICAM", SAMC = "SAMC",
                      IFN = "IFN (custom)", G2M = "G2M (custom)"),
         set = factor(set, rev(c("Hallmark glycolysis", "Hallmark OXPHOS",
                                 "Glycolysis − OXPHOS", "Hallmark IFN-α",
                                 "Hallmark G2M", "DAM-like", "DAM core",
                                 "ICAM", "SAMC", "IFN (custom)", "G2M (custom)"))))
p2e <- ggplot(ht, aes(delta, set)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = delta > 0), size = 2) +
  geom_point(aes(x = delta_clean), shape = 1, size = 2.4, colour = "black") +
  scale_colour_manual(values = c(`TRUE` = "#9B59B6", `FALSE` = "#E67E22"), guide = "none") +
  labs(x = "Mean difference (POS − NEG)", y = NULL, tag = "e",
       title = "Gene-set scores, all Mg",
       subtitle = "Filled = all cells; open = after cleaning")

R2 <- (p2a | p2b) / (p2c | p2d) / p2e +
  plot_layout(heights = c(1, 1.15, 0.9))
save_pub(R2, file.path(FIG, "Figure_R2_decontam_modules_Hallmark"), 183, 220)
cat("wrote Figure_R2\n")
cat("DONE plot\n")
