#!/usr/bin/env Rscript
# Publication figures for the Fig.5/S4 decontamination reanalysis.
# R exclusive. New folder only.
options(warn = 1)
set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
  library(scattermore)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig5_S4_decontam_CellChat_Cxcl10")
SRC <- file.path(OUT, "source_data")
FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

PAL <- c(
  EAE_POS = "#C44E52",
  EAE_NEG = "#4C72B0",
  Mac = "#E28E2C",
  Mg = "#33B5A5",
  orig = "#767676",
  match = "#3182BD"
)
COL_CT <- c(
  Bcell = "#8DA0CB", Tcell = "#66C2A5", Mac = "#E28E2C",
  Neu1 = "#FC8D62", Neu2 = "#E78AC3", Mg = "#33B5A5",
  PC = "#A6D854", DC = "#FFD92F", UK = "#B3B3B3"
)

theme_set(
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      strip.text = element_text(size = 6.2, face = "bold"),
      plot.title = element_text(size = 7, face = "bold"),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
)

save_pub <- function(plot, filename, width_mm = 183, height_mm = 200, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h,
                       family = "Arial")
  print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h,
                units = "in", res = dpi, background = "white")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h,
                 units = "in", res = dpi, compression = "lzw")
  print(plot); dev.off()
}

enr <- read.csv(file.path(SRC, "doublet_hybrid_by_celltype_group.csv"),
                stringsAsFactors = FALSE)
enr$celltype <- factor(enr$celltype,
                       levels = c("Bcell", "Tcell", "Mac", "Neu1", "Neu2",
                                  "Mg", "PC", "DC", "UK"))
enr$IV_label <- factor(enr$IV_label, levels = c("EAE_POS", "EAE_NEG"))

p_a <- ggplot(enr, aes(celltype, pct_dbl, fill = IV_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = PAL[c("EAE_POS", "EAE_NEG")],
                    name = NULL) +
  labs(x = NULL, y = "scDblFinder doublets (%)",
       title = "Doublet call rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

p_b <- ggplot(enr, aes(celltype, pct_neu_high, fill = IV_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = PAL[c("EAE_POS", "EAE_NEG")],
                    name = NULL) +
  labs(x = NULL, y = "Neutrophil-module high (%)",
       title = "Hybrid score (non-Neu cells > P90)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

volcano_panel <- function(path, title, lock, lfc_cut) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[is.finite(d$avg_log2FC) & is.finite(d$p_val_adj) & d$pass_filter, ]
  d$neglog <- pmin(-log10(pmax(d$p_val_adj, 1e-300)), 300)
  d$cls <- "ns"
  d$cls[d$sig_mask & d$avg_log2FC > 0] <- "POS"
  d$cls[d$sig_mask & d$avg_log2FC < 0] <- "NEG"
  d$lab <- d$gene %in% lock & d$pass_filter
  ggplot(d, aes(avg_log2FC, neglog)) +
    geom_scattermore(aes(colour = cls), pointsize = 1.6, pixels = c(900, 700)) +
    geom_hline(yintercept = -log10(0.001), linetype = 2, linewidth = 0.25,
               colour = "grey40") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, linewidth = 0.25,
               colour = "grey40") +
    geom_text_repel(
      data = d[d$lab, ], aes(label = gene),
      size = 1.7, family = "Arial", max.overlaps = 40,
      min.segment.length = 0, segment.size = 0.2, seed = 42
    ) +
    scale_colour_manual(
      values = c(ns = "grey75", POS = PAL[["EAE_POS"]], NEG = PAL[["EAE_NEG"]]),
      guide = "none"
    ) +
    labs(x = "log2 FC (POS / NEG)", y = expression(-log[10](italic(P)[adj])),
         title = title)
}

mac_n <- read.csv(file.path(SRC, "DEG_Mac_clean.csv"), nrows = 1, stringsAsFactors = FALSE)
mg_n  <- read.csv(file.path(SRC, "DEG_Mg_clean.csv"), nrows = 1, stringsAsFactors = FALSE)
t_n   <- read.csv(file.path(SRC, "DEG_Tcell_clean.csv"), nrows = 1, stringsAsFactors = FALSE)
b_n   <- read.csv(file.path(SRC, "DEG_Bcell_clean.csv"), nrows = 1, stringsAsFactors = FALSE)
p_c <- volcano_panel(file.path(SRC, "DEG_Mac_clean.csv"),
                     sprintf("Mac, cell-filtered (n=%d/%d)",
                             mac_n$n_pos[1], mac_n$n_neg[1]),
                     c("Anxa1", "F13a1", "Cxcl9", "Cxcl10", "Cxcl16", "Ccl5",
                       "Cd74", "H2-Ab1", "Arg1", "Nos2"), 1)
p_d <- volcano_panel(file.path(SRC, "DEG_Mg_clean.csv"),
                     sprintf("Mg, cell-filtered (n=%d/%d)",
                             mg_n$n_pos[1], mg_n$n_neg[1]),
                     c("Cxcl13", "Cxcl2", "Cxcr2", "Ccl5", "Ccl2", "Ccl3",
                       "Ccl4", "Anxa1", "Csf1", "Cxcl10"), 1)
p_e <- volcano_panel(file.path(SRC, "DEG_Tcell_clean.csv"),
                     sprintf("T cell, cell-filtered (n=%d/%d)",
                             t_n$n_pos[1], t_n$n_neg[1]),
                     c("Klf2", "S1pr1", "Sell", "Lef1", "Il7r", "Ccr7",
                       "Ifng", "Ctla4", "Cxcr6", "Pdcd1"), 0.5)
p_f <- volcano_panel(file.path(SRC, "DEG_Bcell_clean.csv"),
                     sprintf("B cell, cell-filtered (n=%d/%d)",
                             b_n$n_pos[1], b_n$n_neg[1]),
                     c("Cr2", "Ciita", "Gpr183", "Cd72", "Cxcl2", "Il1b"), 0.5)

lock <- read.csv(file.path(SRC, "locked_genes_cellfilter.csv"),
                 stringsAsFactors = FALSE)
lock$gene <- factor(lock$gene, levels = unique(lock$gene))
lock$celltype <- factor(lock$celltype, levels = c("Mac", "Mg", "Tcell", "Bcell"))
lock_long <- lock %>%
  select(celltype, gene, all_log2FC, clean_log2FC) %>%
  pivot_longer(c(all_log2FC, clean_log2FC), names_to = "set",
               values_to = "log2FC") %>%
  mutate(set = recode(set, all_log2FC = "All cells",
                      clean_log2FC = "Cell-filtered"))
p_g <- ggplot(lock_long, aes(log2FC, gene, colour = set)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(size = 1.4, alpha = 0.9) +
  facet_wrap(~celltype, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c("All cells" = PAL[["orig"]],
                                 "Cell-filtered" = PAL[["match"]]),
                      name = NULL) +
  labs(x = "log2 FC (POS / NEG)", y = NULL,
       title = "Locked genes before vs after cell filter") +
  theme(legend.position = "top")

fig1 <- (p_a | p_b) /
  (p_c | p_d) /
  (p_e | p_f) /
  p_g +
  plot_layout(heights = c(0.9, 1.15, 1.15, 0.95)) +
  plot_annotation(tag_levels = "a")
save_pub(fig1, file.path(FIG, "Figure_decontam"), 183, 240)
cat("wrote Figure_decontam\n")

# ---- figure 2: CellChat + Cxcl10 (CellChat file may be missing) ----
cx <- read.csv(file.path(SRC, "Cxcl10_Mac_vs_Mg_scRNA.csv"),
               stringsAsFactors = FALSE)
cx$set <- factor(cx$set, levels = c("all", "clean"),
                 labels = c("All cells", "Cell-filtered"))
cx$IV_label <- factor(cx$IV_label, levels = c("EAE_POS", "EAE_NEG"))
p_cx <- ggplot(cx, aes(celltype, pct, fill = IV_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           colour = "black", linewidth = 0.2) +
  facet_wrap(~set) +
  scale_fill_manual(values = PAL[c("EAE_POS", "EAE_NEG")], name = NULL) +
  labs(x = NULL, y = "Cells with Cxcl10 > 0 (%)",
       title = "scRNA Cxcl10 detection") +
  theme(legend.position = "top")

p_cxm <- ggplot(cx, aes(celltype, mean, fill = IV_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           colour = "black", linewidth = 0.2) +
  facet_wrap(~set) +
  scale_fill_manual(values = PAL[c("EAE_POS", "EAE_NEG")], name = NULL) +
  labs(x = NULL, y = "Mean Cxcl10 (data layer)",
       title = "scRNA Cxcl10 mean") +
  theme(legend.position = "top")

vsum <- read.csv(file.path(SRC, "visium_cxcl10_by_MacMg.csv"),
                 stringsAsFactors = FALSE)
vsum$contact <- factor(vsum$contact, levels = c("POS", "NEG"),
                       labels = c("EAE_POS", "EAE_NEG"))
p_vs <- ggplot(vsum, aes(lineage, cxcl10_pct, fill = contact)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = PAL[c("EAE_POS", "EAE_NEG")], name = NULL) +
  labs(x = NULL, y = "Cxcl10+ spots (%)",
       title = "Visium, top-predicted Mac/Mg") +
  theme(legend.position = "top")

cors <- read.csv(file.path(SRC, "visium_cxcl10_marker_spearman.csv"),
                 stringsAsFactors = FALSE)
cors$cls <- ifelse(cors$gene %in% c("Ccr2", "Plac8"), "monocyte-like",
                   ifelse(cors$gene %in% c("P2ry12", "Tmem119", "Hexb"),
                          "microglia", "chemokine"))
cors$gene <- factor(cors$gene, levels = cors$gene[order(cors$spearman_rho)])
p_cor <- ggplot(cors, aes(spearman_rho, gene, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("monocyte-like" = PAL[["Mac"]],
                               microglia = PAL[["Mg"]],
                               chemokine = "grey50"),
                    name = NULL) +
  labs(x = "Spearman rho vs Cxcl10 (all spots)", y = NULL,
       title = "Visium Cxcl10 vs lineage markers") +
  theme(legend.position = "top")

spat <- readRDS(file.path(SRC, "visium_cxcl10_lineage.rds"))
set.seed(42)
# Display subsample of zeros plus every Cxcl10+ spot (seed recorded).
# geom_point: scattermore did not paint on cairo/ragg devices here.
cx_hit <- which(spat$Cxcl10 > 0)
bg <- sample.int(nrow(spat), min(6000L, nrow(spat)))
spat_sc <- spat[unique(c(bg, cx_hit)), ]
p_sc1 <- ggplot(spat_sc, aes(Ccr2, Cxcl10)) +
  geom_point(size = 0.35, alpha = 0.25, colour = PAL[["Mac"]], stroke = 0) +
  labs(x = "Ccr2", y = "Cxcl10", title = "Visium Cxcl10 vs Ccr2")
p_sc2 <- ggplot(spat_sc, aes(P2ry12, Cxcl10)) +
  geom_point(size = 0.35, alpha = 0.25, colour = PAL[["Mg"]], stroke = 0) +
  labs(x = "P2ry12", y = "Cxcl10", title = "Visium Cxcl10 vs P2ry12")

cc_file <- file.path(SRC, "cellchat_original_vs_matched.csv")
if (file.exists(cc_file)) {
  cmp <- read.csv(cc_file, stringsAsFactors = FALSE)
  keep_pw <- c("GALECTIN", "CCL", "MIF", "COMPLEMENT", "CXCL", "ANNEXIN",
               "VISFATIN", "TGFb", "GRN", "CSF")
  cmp <- cmp[cmp$pathway %in% keep_pw, ]
  cmp$pathway <- factor(cmp$pathway, levels = rev(keep_pw))
  long <- bind_rows(
    data.frame(pathway = cmp$pathway, source = "original",
               ratio = cmp$ratio_orig),
    data.frame(pathway = cmp$pathway, source = "composition-matched",
               ratio = cmp$ratio_match)
  )
  long$ratio[!is.finite(long$ratio) | long$ratio <= 0] <- NA
  long$source <- factor(long$source, levels = c("original", "composition-matched"),
                        labels = c("Original", "Composition-matched"))
  p_cc <- ggplot(long, aes(ratio, pathway, fill = source)) +
    geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey40") +
    geom_col(position = position_dodge(0.8), width = 0.7,
             colour = "black", linewidth = 0.2) +
    scale_x_continuous(trans = "log10") +
    scale_fill_manual(values = c(Original = PAL[["orig"]],
                                 "Composition-matched" = PAL[["match"]]),
                      name = NULL) +
    labs(x = "Outgoing POS / NEG (log10)", y = NULL,
         title = "CellChat outgoing after composition match") +
    theme(legend.position = "top")
  fig2 <- (p_cc | p_cx) /
    (p_cxm | p_vs) /
    (p_cor | p_sc1 | p_sc2) +
    plot_annotation(tag_levels = "a")
  save_pub(fig2, file.path(FIG, "Figure_cellchat_cxcl10"), 183, 210)
} else {
  fig2 <- (p_cx | p_cxm) /
    (p_vs | p_cor) /
    (p_sc1 | p_sc2) +
    plot_annotation(tag_levels = "a")
  save_pub(fig2, file.path(FIG, "Figure_cxcl10"), 183, 180)
  cat("CellChat comparison missing; wrote Figure_cxcl10 without panel a\n")
}
cat("DONE plot\n")
