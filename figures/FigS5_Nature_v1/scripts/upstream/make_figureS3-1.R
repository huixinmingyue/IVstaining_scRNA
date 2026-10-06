#!/usr/bin/env Rscript
# =============================================================================
# Fig. S3-1 revised — Nature SI (native ggplot2 vector)
# a labeled gel | b UMAP | c DotPlot (B markers + myeloid)
# d logistic | e PF>=6 vs <=1 split UMAP | f proportion | g C2 volcano
# New folder FigS3-1_NC_PFsplit — does not overwrite prior packages
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
  library(png)
  library(grid)
  library(cowplot)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-1_NC_PFsplit")
SRC <- file.path(OUTDIR, "source_data")
FIG_W_MM <- 183
FIG_H_MM <- 240
FAM <- "Arial"
BASE <- 6.5
FC_CUT <- 0.25

btype_cols <- c(
  "Follicular B cells (Gpr183+)" = "#F8766D",
  "Follicular B cells" = "#B79F00",
  "Marginal zone B cells" = "#00BA38",
  "Follicular B cells [3]" = "#00BFC4",
  "Myeloid-contaminated" = "#619CFF",
  "Immature / pre-B cells" = "#F564E3"
)
# noC4 lacks Myeloid-contaminated
btype_cols_noc4 <- btype_cols[names(btype_cols) != "Myeloid-contaminated"]

short_bt <- c(
  "Follicular B cells (Gpr183+)" = "Follicular B (Gpr183+)",
  "Follicular B cells" = "Follicular B",
  "Marginal zone B cells" = "Marginal zone B",
  "Follicular B cells [3]" = "Follicular B [3]",
  "Myeloid-contaminated" = "Myeloid-contaminated",
  "Immature / pre-B cells" = "Immature / pre-B"
)

label_genes <- c("Cd274", "Ezr", "Marcks", "Cr2", "Dtx1", "Klf2")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = BASE),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.5),
    legend.text = element_text(size = BASE - 1.1),
    legend.key.size = unit(2.4, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 1.3, hjust = 0.5, colour = "grey30"),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 3, 2, 2),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = BASE - 0.4, face = "bold", family = FAM)
  )

# ---------- a ----------
cat("Panel a ...\n")
img_a <- readPNG(file.path(SRC, "panel_a_schematic.png"))
pA <- ggdraw() +
  draw_grob(rasterGrob(img_a, interpolate = TRUE)) +
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.02, 0.99),
    plot.margin = margin(1, 2, 1, 1)
  ) +
  labs(tag = "a")

# ---------- b ----------
cat("Panel b ...\n")
df_b <- readRDS(file.path(SRC, "panel_b_umap.rds"))
df_b[["Bcelltype"]] <- factor(df_b[["Bcelltype"]], levels = names(btype_cols))
lab_b <- df_b %>%
  group_by(Bcelltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop") %>%
  mutate(lab = unname(short_bt[as.character(Bcelltype)]))

pB <- ggplot(df_b, aes(umap_1, umap_2, colour = Bcelltype)) +
  geom_point(size = 0.18, stroke = 0, shape = 16, alpha = 0.85) +
  geom_text_repel(
    data = lab_b, aes(umap_1, umap_2, label = lab),
    colour = "black", size = 1.7, family = FAM,
    max.overlaps = 30, box.padding = 0.18, segment.size = 0.18,
    min.segment.length = 0, seed = 42, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = btype_cols, name = NULL,
                      labels = short_bt[names(btype_cols)]) +
  coord_equal(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "B cells (Harmony, resolution 0.2)", tag = "b") +
  theme_pub +
  theme(legend.position = "right", legend.text = element_text(size = BASE - 1.4)) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))

# ---------- c ----------
cat("Panel c ...\n")
dot <- readRDS(file.path(SRC, "panel_c_dotplot.rds"))
ann_csv <- read.csv(
  file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_res0.2_PF_corr/data/Bcell_cluster_annotation.csv"),
  stringsAsFactors = FALSE
)
ylab <- paste0(as.character(ann_csv[["cluster"]]), ": ",
               unname(short_bt[as.character(ann_csv[["Bcelltype"]])]))
names(ylab) <- as.character(ann_csv[["cluster"]])
dot[["cluster_lab"]] <- factor(unname(ylab[as.character(dot[["cluster"]])]), levels = ylab)

pC <- ggplot(dot, aes(gene, cluster_lab)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_colour_gradient(low = "lightgrey", high = "#B2182B", name = "Avg. exp.\n(scaled)") +
  scale_size_continuous(name = "% exp.", range = c(0.35, 3.0), breaks = c(0, 25, 50, 75)) +
  labs(x = NULL, y = NULL,
       title = "B-lineage + myeloid markers", tag = "c") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = BASE - 1.7, face = "italic"),
    axis.text.y = element_text(size = BASE - 1.2),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.2)
  )

# ---------- d ----------
cat("Panel d ...\n")
pred <- readRDS(file.path(SRC, "panel_d_logistic.rds"))
pred[["Bcelltype"]] <- factor(pred[["Bcelltype"]], levels = names(btype_cols))

pD <- ggplot(pred, aes(Proteinfeatures, 100 * prob, colour = Bcelltype)) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = btype_cols, name = NULL,
                      labels = short_bt[names(btype_cols)]) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100)) +
  labs(x = "Proteinfeatures", y = "Predicted probability (%)",
       title = "Logistic P(Bcelltype) ~ log1p(PF)", tag = "d") +
  theme_pub +
  theme(legend.position = "right", legend.text = element_text(size = BASE - 1.4))

# ---------- e: split UMAP High vs Low ----------
cat("Panel e ...\n")
df_e <- readRDS(file.path(SRC, "panel_e_umap_split.rds"))
df_e[["Bcelltype"]] <- factor(df_e[["Bcelltype"]], levels = names(btype_cols_noc4))

pE <- ggplot(df_e, aes(umap_1, umap_2, colour = Bcelltype)) +
  geom_point(size = 0.16, stroke = 0, shape = 16, alpha = 0.85) +
  facet_wrap(~ PF_class, nrow = 1) +
  scale_colour_manual(values = btype_cols_noc4, name = NULL,
                      labels = short_bt[names(btype_cols_noc4)]) +
  coord_equal(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "no-C4 UMAP split by PF",
       subtitle = "High PF>=6 vs Low PF<=1", tag = "e") +
  theme_pub +
  theme(legend.position = "right", legend.text = element_text(size = BASE - 1.5)) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 6))

# ---------- f: proportion High vs Low ----------
cat("Panel f ...\n")
prop <- readRDS(file.path(SRC, "panel_f_proportion.rds"))
prop[["Bcelltype"]] <- factor(prop[["Bcelltype"]], levels = names(btype_cols_noc4))

pF <- ggplot(prop, aes(PF_class, pct, fill = Bcelltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = btype_cols_noc4, name = NULL,
                    labels = short_bt[names(btype_cols_noc4)]) +
  labs(x = NULL, y = "Proportion (%)",
       title = "no-C4 subtype composition",
       subtitle = "PF>=6 vs PF<=1", tag = "f") +
  theme_pub +
  theme(
    axis.text.x = element_text(size = BASE - 1.0),
    legend.position = "right",
    legend.text = element_text(size = BASE - 1.5)
  )

# ---------- g: volcano C2 ge6 le1, selected labels ----------
cat("Panel g ...\n")
deg <- read.csv(file.path(SRC, "panel_g_DEG_C2_ge6_le1.csv"), stringsAsFactors = FALSE)
if (file.exists(file.path(SRC, "panel_g_Cd274_added.csv"))) {
  cd <- read.csv(file.path(SRC, "panel_g_Cd274_added.csv"), stringsAsFactors = FALSE)
  if (!("Cd274" %in% deg[["gene"]])) deg <- bind_rows(deg, cd)
}
deg[["neglog10p"]] <- -log10(pmax(deg[["p_val_adj"]], 1e-300))
deg[["sig"]] <- ifelse(
  !is.na(deg[["p_val_adj"]]) & deg[["p_val_adj"]] < 0.05 & abs(deg[["avg_log2FC"]]) >= FC_CUT,
  "sig", "ns"
)
n_sig <- sum(deg[["sig"]] == "sig", na.rm = TRUE)
lab_g <- deg %>% filter(gene %in% label_genes)
cat("  labelling:", paste(lab_g[["gene"]], collapse = ", "), "\n")

pG <- ggplot(deg, aes(avg_log2FC, neglog10p)) +
  geom_point(aes(colour = sig), size = 0.75, alpha = 0.85) +
  geom_point(data = lab_g, shape = 21, size = 1.7, stroke = 0.5,
             colour = "#08306B", fill = NA) +
  scale_colour_manual(values = c(ns = "grey70", sig = "#B2182B"), guide = "none") +
  geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = 2,
             colour = "grey45", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = 2,
             colour = "grey45", linewidth = 0.35) +
  geom_text_repel(
    data = lab_g, aes(label = gene),
    size = 2.1, family = FAM, colour = "black", fontface = "italic",
    max.overlaps = Inf, box.padding = 0.3, point.padding = 0.15,
    segment.size = 0.22, segment.colour = "grey40",
    min.segment.length = 0, seed = 42
  ) +
  labs(
    x = expression(avg~log[2]*FC~(HighPF / LowPF)),
    y = expression(-log[10](adj~P)),
    title = "C2 MZ B: PF>=6 vs <=1",
    subtitle = sprintf("red: padj<0.05 & |FC|>=%.2f (n=%d)", FC_CUT, n_sig),
    tag = "g"
  ) +
  theme_pub

# ---------- Assemble ----------
cat("Assemble ...\n")
design <- "
AABB
AACC
DDEE
FFGG
"
fig <- pA + pB + pC + pD + pE + pF + pG +
  plot_layout(design = design, heights = c(1.05, 1.05, 1.0, 1.05)) &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FAM))

save_pub <- function(plot, stem, width_mm, height_mm) {
  # Nature raster export: dpi = 600
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  cat("  cairo_pdf\n")
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = FAM)
  print(plot); grDevices::dev.off()
  cat("  svglite\n")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); grDevices::dev.off()
  cat("  ragg png dpi = 600\n")
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600)
  print(plot); grDevices::dev.off()
}

stem <- file.path(OUTDIR, "FigureS3-1")
save_pub(fig, stem, FIG_W_MM, FIG_H_MM)

writeLines(c(
  "Fig. S3-1 | CITE-seq oligo QC, B-cell annotation, PF thresholds and MZ DEGs",
  "a, Oligo-CD45 conjugate gels with lane labels matching FigureS3 (M1/42, 30-F11 +/- Oligo; groups 1-2).",
  "b, UMAP of Harmony-reclustered B cells (resolution 0.2) by Bcelltype.",
  "c, Dot plot: B-lineage markers (Cr2..Rag1) plus myeloid markers (Slpi, Axl, Mafb, Hpgd, Fcna).",
  "d, Logistic P(Bcelltype) versus continuous Proteinfeatures.",
  "e, no-C4 UMAP split by High (PF>=6) vs Low (PF<=1).",
  "f, no-C4 Bcelltype composition for the same PF thresholds.",
  "g, C2 marginal-zone B volcano (PF>=6 vs <=1); labels: Cd274, Ezr, Marcks, Cr2, Dtx1, Klf2.",
  "Note: Cd274 was not in |FC|>=0.25 DEG export for ge6_le1; point added from Wilcoxon (FC~0.10, ns).",
  "Backend: R 4.4.1 exclusive; native ggplot vector for b-g."
), file.path(OUTDIR, "FigureS3-1_legend.txt"))

writeLines(c(
  "Figure S3-1 | Nature Communications SI revised package",
  paste0("Generated: ", Sys.time()),
  "Core conclusion: Oligo QC + B annotation; PF>=6 vs <=1 stratifies niches; MZ DEGs highlight selected genes.",
  "Archetype: asymmetric mixed-modality",
  "Backend: R exclusive (ggplot); user preferred Python but native ggplot requires R",
  paste0("Final size: ", FIG_W_MM, " x ", FIG_H_MM, " mm"),
  "E/F thresholds: High PF>=6, Low PF<=1 (no C4)",
  "G: former volcano-style panel; labels limited to six genes",
  "New folder FigS3-1_NC_PFsplit/; prior FigS3-1_NC_submission untouched"
), file.path(OUTDIR, "FigureS3-1_QA.txt"))

writeLines(c(
  "a: 图片2.tif + lane labels from Figures_Nature_R_20260824/FigureS3",
  "b/c/d: res0.2 annotated object",
  "e/f: noC4 object, PF>=6 vs <=1",
  "g: DEG_C2_ge6_le1 (+ Cd274 added)"
), file.path(OUTDIR, "FigureS3-1_provenance.txt"))

file.copy(file.path(OUTDIR, "scripts", "make_figureS3-1.R"),
          file.path(OUTDIR, "make_figureS3-1.R"), overwrite = TRUE)

cat("DONE ->", OUTDIR, "\n")
