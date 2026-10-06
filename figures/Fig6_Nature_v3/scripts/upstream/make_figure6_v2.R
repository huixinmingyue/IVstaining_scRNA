#!/usr/bin/env Rscript
# =============================================================================
# Fig. 6 v2 — Nature multipanel (Fig6_Nature_v2)
# New folder; does not overwrite Fig6_Nature_v1 or any earlier package.
#
# Core conclusion
#   IV CD45 gates separate vascular-contact (POS) from parenchyma-associated
#   (NEG) CNS immune cells. Original CellChat outgoing POS/NEG ratios for
#   GALECTIN to CSF mark a CXCL/ANNEXIN-high POS programme. Macrophage and
#   microglia DEGs enrich distinct KEGG signalling terms, and locked myeloid
#   genes keep direction after cell filtering.
#
# Evidence chain
#   a  How were IV-labelled CNS cells gated as NEG vs POS?          FACS
#   b  Do POS and NEG occupy different cell-type space?             UMAP
#   c  Does composition shift with IV contact?                      bars
#   d  What is the original outgoing POS/NEG ratio (GALECTIN-CSF)?  CellChat
#   e  What KEGG signalling programmes do Mac DEGs over-represent?  ORA
#   f  What KEGG signalling programmes do Mg DEGs over-represent?   ORA
#   g  Do locked Mac genes keep POS/NEG FC after cell filtering?    dots
#   h  Do locked Mg genes keep POS/NEG FC after cell filtering?     dots
#
# Changed relative to Fig6_Nature_v1
#   d  heatmap replaced by original-only POS/NEG outgoing ratios
#   e/f  previous Mac/Mg KEGG panels (were f/h)
#   g/h  locked-gene Mac/Mg from the decontamination figure (panel g)
#   previous volcanoes (e/g) omitted
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
  library(tibble)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(grid)
  library(png)
  library(ggrepel)
  library(scattermore)
  library(magick)
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
  file.path(IV_FIG_ROOT, "Fig6_Nature_v2")
}
SRC     <- file.path(OUTDIR, "source_data")
SRCFIG  <- file.path(OUTDIR, "source_figures")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 230
width_mm = 183
height_mm = 230
BASE <- 6.5
FAM  <- "Arial"
PT5  <- 5.2
LABEL_PT <- 5.2 / ggplot2::.pt
ARROW_COL <- "#2E75B6"
COL_ORIG <- "#767676"
COL_ALL  <- "#767676"
COL_CLEAN <- "#3182BD"

ct_levels <- c("Bcell", "Tcell", "Mac", "Neu1", "Neu2", "Mg", "PC", "DC", "UK")
ct_lab <- c(
  Bcell = "B cell", Tcell = "T cell", Mac = "Mac", Neu1 = "Neu1",
  Neu2 = "Neu2", Mg = "Mg", PC = "PC", DC = "DC", UK = "UK"
)
ct_cols <- c(
  Bcell = "#E64B35", Tcell = "#4DBBD5", Mac = "#00A087", Neu1 = "#3C5488",
  Neu2 = "#F39B7F", Mg = "#8491B4", PC = "#91D1C2", DC = "#DC0000", UK = "#7E6148"
)
names(ct_cols) <- ct_lab[names(ct_cols)]
grp_cols <- c(POS = "#E64B35", NEG = "#4DBBD5")
TOP10 <- c("GALECTIN", "CCL", "MIF", "COMPLEMENT", "CXCL",
           "ANNEXIN", "VISFATIN", "TGFb", "GRN", "CSF")

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

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = BASE),
      axis.text = element_text(size = PT5, colour = "black"),
      legend.title = element_text(size = PT5),
      legend.text = element_text(size = PT5),
      legend.key.size = unit(2.4, "mm"),
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PT5, hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold")
    )
)

tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.01, 0.98)
)

wrap_path <- function(x, width = 28) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

kegg_display <- function(x) {
  x <- sub("\\s+signaling pathway\\b", " signaling", x)
  wrap_path(x, width = 32)
}

# -----------------------------------------------------------------------------
# a — FACS plate
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
tif_a <- file.path(SRCFIG, "Fig6a_FACS_source.tif")
stopifnot(file.exists(tif_a))
im_raw <- magick::image_background(magick::image_read(tif_a), "white",
                                   flatten = TRUE)
im_crop <- magick::image_crop(im_raw, "429x405+75+69")
crop_png <- file.path(SRCFIG, "Fig6a_FACS_cropped.png")
magick::image_write(im_crop, crop_png, format = "png")
img_a <- png::readPNG(crop_png)
n_a_px <- c(ncol(img_a), nrow(img_a))
img_h <- 1 * nrow(img_a) / ncol(img_a)

pA <- ggplot() +
  annotation_raster(img_a, xmin = 0, xmax = 1, ymin = 0, ymax = img_h,
                    interpolate = TRUE) +
  annotate("text", x = 0.30, y = img_h + 0.085, label = "NEG",
           fontface = "bold", size = 7 / ggplot2::.pt, family = FAM,
           colour = "black") +
  annotate("text", x = 0.72, y = img_h + 0.085, label = "POS",
           fontface = "bold", size = 7 / ggplot2::.pt, family = FAM,
           colour = "black") +
  annotate(
    "segment", x = -0.09, xend = -0.09, y = 0.10, yend = img_h - 0.10,
    colour = ARROW_COL, linewidth = 0.45,
    arrow = arrow(length = unit(1.7, "mm"), type = "closed", angle = 20)
  ) +
  annotate("text", x = -0.19, y = img_h / 2, label = "Extra CD45", angle = 90,
           size = BASE / ggplot2::.pt, family = FAM, colour = "black") +
  annotate(
    "segment", x = 0.12, xend = 0.88, y = -0.08, yend = -0.08,
    colour = ARROW_COL, linewidth = 0.45,
    arrow = arrow(length = unit(1.7, "mm"), type = "closed", angle = 20)
  ) +
  annotate("text", x = 0.50, y = -0.18, label = "Intra CD45",
           size = BASE / ggplot2::.pt, family = FAM, colour = "black") +
  coord_fixed(xlim = c(-0.28, 1.04), ylim = c(-0.26, img_h + 0.14),
              expand = FALSE, clip = "off", ratio = 1) +
  theme_void(base_family = FAM) +
  labs(tag = "a") +
  theme(
    plot.margin = margin(2, 2, 4, 6),
    plot.background = element_rect(fill = "white", colour = NA)
  ) + tag_theme

# -----------------------------------------------------------------------------
# b — UMAP
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
df_b <- readRDS(file.path(SRC, "Fig4B_umap.rds"))
n_b_before <- nrow(df_b)
df_b <- df_b[!is.na(df_b$umap_1) & !is.na(df_b$umap_2) &
               !is.na(df_b$celltype) & !is.na(df_b$group_label), ]
n_b_after <- nrow(df_b)
df_b$celltype <- factor(as.character(df_b$celltype), levels = ct_levels)
df_b$cell_lab <- factor(ct_lab[as.character(df_b$celltype)],
                        levels = unname(ct_lab[ct_levels]))
df_b$group_show <- factor(
  ifelse(df_b$group_label == "EAE_POS", "POS", "NEG"),
  levels = c("POS", "NEG")
)
n_b_pos <- sum(df_b$group_show == "POS")
n_b_neg <- sum(df_b$group_show == "NEG")
lab_b <- df_b %>%
  group_by(group_show, cell_lab) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2),
            n = dplyr::n(), .groups = "drop")
write.csv(df_b, file.path(PLOTTED, "SourceData_Fig6b_umap.csv"), row.names = FALSE)
write.csv(lab_b, file.path(PLOTTED, "SourceData_Fig6b_labels.csv"), row.names = FALSE)

pB <- ggplot(df_b, aes(umap_1, umap_2, colour = cell_lab)) +
  geom_scattermore(pointsize = 1.6, pixels = c(1100, 900), interpolate = FALSE) +
  geom_text_repel(
    data = lab_b, aes(x = umap_1, y = umap_2, label = cell_lab),
    colour = "black", size = LABEL_PT, fontface = "bold", family = FAM,
    inherit.aes = FALSE, max.overlaps = Inf, seed = 42,
    box.padding = 0.28, point.padding = 0.15,
    min.segment.length = 0.12, segment.size = 0.2,
    segment.colour = "grey40", force = 2.2, force_pull = 0.4
  ) +
  scale_colour_manual(values = ct_cols, name = "Cell type") +
  facet_wrap(~ group_show, nrow = 1) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "CNS immune cells",
       subtitle = sprintf("n = %s POS; %s NEG",
                          format(n_b_pos, big.mark = ","),
                          format(n_b_neg, big.mark = ",")),
       tag = "b") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.size = unit(2.2, "mm"),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2))) +
  tag_theme

# -----------------------------------------------------------------------------
# c — composition
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
prop <- read.csv(file.path(SRC, "Fig4C_proportions.csv"), stringsAsFactors = FALSE)
n_c_before <- nrow(prop)
prop <- prop[!is.na(prop$proportion), ]
n_c_after <- nrow(prop)
prop$celltype <- factor(prop$celltype, levels = ct_levels)
prop$cell_lab <- factor(ct_lab[as.character(prop$celltype)],
                        levels = unname(ct_lab[ct_levels]))
prop$group_show <- factor(
  ifelse(prop$group == "EAE_POS", "POS", "NEG"),
  levels = c("POS", "NEG")
)
write.csv(prop, file.path(PLOTTED, "SourceData_Fig6c_proportions.csv"),
          row.names = FALSE)

pC <- ggplot(prop, aes(cell_lab, proportion, fill = group_show)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion", title = "Cell-type composition", tag = "c") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = PT5),
        legend.position = "top") +
  tag_theme

# -----------------------------------------------------------------------------
# d — CellChat original POS/NEG outgoing ratio (GALECTIN to CSF)
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
cc_cmp <- read.csv(file.path(SRC, "cellchat_original_vs_matched.csv"),
                   stringsAsFactors = FALSE)
n_cc_before <- nrow(cc_cmp)
cc_keep <- cc_cmp[cc_cmp$pathway %in% TOP10, ]
n_cc_after <- nrow(cc_keep)
stopifnot(nrow(cc_keep) == 10)
d_df <- data.frame(
  pathway = factor(cc_keep$pathway, levels = rev(TOP10)),
  POS = cc_keep$POS_orig,
  NEG = cc_keep$NEG_orig,
  ratio = cc_keep$ratio_orig,
  stringsAsFactors = FALSE
)
n_ratio_raw <- nrow(d_df)
d_plot <- d_df
d_plot$ratio[!is.finite(d_plot$ratio) | d_plot$ratio <= 0] <- NA
n_ratio_plotted <- sum(is.finite(d_plot$ratio))
write.csv(cc_cmp, file.path(PLOTTED, "SourceData_Fig6d_cellchat_orig_vs_matched.csv"),
          row.names = FALSE)
write.csv(d_df, file.path(PLOTTED, "SourceData_Fig6d_cellchat_original_top10.csv"),
          row.names = FALSE)

pD <- ggplot(d_plot, aes(ratio, pathway)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey40") +
  geom_col(width = 0.72, fill = COL_ORIG, colour = "black", linewidth = 0.2,
           na.rm = TRUE) +
  geom_point(size = 1.1, colour = "black", na.rm = TRUE) +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.2, 0.5, 1, 2, 5, 10, 20),
                     labels = c("0.2", "0.5", "1", "2", "5", "10", "20")) +
  labs(x = "Outgoing POS / NEG (log10)", y = NULL,
       title = "Outgoing signalling (original)",
       subtitle = "GALECTIN to CSF",
       tag = "d") +
  tag_theme

# -----------------------------------------------------------------------------
# e / f — KEGG signalling ORA (previous f / h)
# -----------------------------------------------------------------------------
signaling_only <- function(k, ct) {
  keep <- grepl("signaling", k$Description, ignore.case = TRUE)
  cat(ct, "signaling filter: KEGG terms", nrow(k), "-> displayed pool",
      sum(keep), "; FDR<0.05 among them", sum(k$FDR05[keep]), "\n")
  k[keep, , drop = FALSE]
}

make_kegg <- function(k, title) {
  top <- k %>% arrange(p.adjust) %>% slice_head(n = 10)
  n_fdr <- sum(top$FDR05)
  sub <- if (n_fdr > 0) {
    sprintf("Signaling terms only; black outline: FDR < 0.05 (n = %d)", n_fdr)
  } else {
    "Signaling terms only; none reach FDR < 0.05 (exploratory)"
  }
  top$Description <- kegg_display(as.character(top$Description))
  top$Description <- factor(top$Description, levels = rev(unique(top$Description)))
  ggplot(top, aes(GeneRatio_num, Description)) +
    geom_point(aes(size = Count, fill = neglog10_padj, colour = FDR05),
               shape = 21, stroke = 0.45) +
    scale_fill_gradient(low = "#C6DBEF", high = "#08519C",
                        name = "-log10(FDR)") +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey70"),
                        guide = "none") +
    scale_size_area(max_size = 4.4, name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.14))) +
    labs(x = "Gene ratio", y = "KEGG signalling pathway",
         title = title, subtitle = sub) +
    theme(axis.text.y = element_text(size = PT5),
          legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(1.5, "mm"),
          legend.key.size = unit(2.0, "mm"),
          plot.subtitle = element_text(size = PT5, hjust = 0),
          plot.margin = margin(4, 4, 2, 2))
}

cat("KEGG ...\n")
kegg_mac_all <- read.csv(file.path(SRC, "Fig5_KEGG_ORA_Mac.csv"),
                         stringsAsFactors = FALSE)
kegg_mg_all  <- read.csv(file.path(SRC, "Fig5_KEGG_ORA_Mg.csv"),
                         stringsAsFactors = FALSE)
kegg_mac <- signaling_only(kegg_mac_all, "Mac")
kegg_mg  <- signaling_only(kegg_mg_all, "Mg")
stopifnot(nrow(kegg_mac) >= 10, nrow(kegg_mg) >= 10)
n_fdr_mac <- sum(kegg_mac$FDR05)
n_fdr_mg  <- sum(kegg_mg$FDR05)
n_fdr_mac_top <- sum((kegg_mac %>% arrange(p.adjust) %>% slice_head(n = 10))$FDR05)
n_fdr_mg_top  <- sum((kegg_mg  %>% arrange(p.adjust) %>% slice_head(n = 10))$FDR05)

write.csv(kegg_mac_all, file.path(PLOTTED, "SourceData_Fig6e_KEGG_Mac_allterms.csv"),
          row.names = FALSE)
write.csv(kegg_mg_all,  file.path(PLOTTED, "SourceData_Fig6f_KEGG_Mg_allterms.csv"),
          row.names = FALSE)
write.csv(kegg_mac, file.path(PLOTTED, "SourceData_Fig6e_KEGG_Mac_signaling.csv"),
          row.names = FALSE)
write.csv(kegg_mg,  file.path(PLOTTED, "SourceData_Fig6f_KEGG_Mg_signaling.csv"),
          row.names = FALSE)

pE <- make_kegg(kegg_mac, "Mac KEGG") + labs(tag = "e") + tag_theme
pF <- make_kegg(kegg_mg,  "Mg KEGG")  + labs(tag = "f") + tag_theme

# -----------------------------------------------------------------------------
# g / h — locked genes Mac / Mg, all cells vs cell-filtered
# -----------------------------------------------------------------------------
cat("Panels g/h ...\n")
lock <- read.csv(file.path(SRC, "locked_genes_cellfilter.csv"),
                 stringsAsFactors = FALSE)
n_lock_before <- nrow(lock)
gene_lv <- unique(as.character(lock$gene))
lock_mm <- lock[lock$celltype %in% c("Mac", "Mg"), ]
n_lock_after <- nrow(lock_mm)
n_lock_dropped <- n_lock_before - n_lock_after
lock_long <- lock_mm %>%
  select(celltype, gene, all_log2FC, clean_log2FC) %>%
  pivot_longer(c(all_log2FC, clean_log2FC), names_to = "set",
               values_to = "log2FC") %>%
  mutate(set = recode(set, all_log2FC = "All cells",
                      clean_log2FC = "Cell-filtered"))
lock_long$set <- factor(lock_long$set, levels = c("All cells", "Cell-filtered"))
lock_long$gene <- factor(as.character(lock_long$gene), levels = gene_lv)
x_lock <- range(lock_long$log2FC, na.rm = TRUE)
x_lock <- c(min(-0.5, x_lock[1]), max(0.5, x_lock[2])) * c(1.08, 1.08)

write.csv(lock, file.path(PLOTTED, "SourceData_Fig6gh_locked_genes_alltypes.csv"),
          row.names = FALSE)
write.csv(lock_mm, file.path(PLOTTED, "SourceData_Fig6gh_locked_genes_MacMg.csv"),
          row.names = FALSE)
write.csv(lock_long, file.path(PLOTTED, "SourceData_Fig6gh_locked_genes_long.csv"),
          row.names = FALSE)

make_lock <- function(ct, title, tag) {
  d <- lock_long[lock_long$celltype == ct, ]
  d$gene <- droplevels(d$gene)
  ggplot(d, aes(log2FC, gene, colour = set)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 1.55, alpha = 0.95) +
    scale_colour_manual(values = c("All cells" = COL_ALL,
                                   "Cell-filtered" = COL_CLEAN),
                        name = NULL) +
    coord_cartesian(xlim = x_lock) +
    labs(x = "log2 FC (POS / NEG)", y = NULL, title = title, tag = tag) +
    theme(legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0)) +
    tag_theme
}

pG <- make_lock("Mac", "Mac locked genes", "g")
pH <- make_lock("Mg",  "Mg locked genes",  "h") +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- (pA | pB | pC) + plot_layout(widths = c(1.28, 1.42, 1.05))
row2 <- pD
row3 <- (pE | pF) + plot_layout(widths = c(1, 1))
row4 <- (pG | pH) + plot_layout(widths = c(1, 1))
fig <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(0.95, 0.92, 1.18, 1.08))

# -----------------------------------------------------------------------------
# Stats / legend / provenance
# -----------------------------------------------------------------------------
stats <- c(
  "Figure 6 statistics / integrity (Fig6_Nature_v2)",
  sprintf("Panel a raster after crop: %d x %d px; crop 429x405+75+69; no contrast change",
          n_a_px[1], n_a_px[2]),
  "Panel a vector overlays: NEG, POS, Extra CD45, Intra CD45. Gate percentages 58.3 and 15.4 remain native cytometer glyphs.",
  sprintf("Panel b cells: before %d after NA drop %d (POS %d, NEG %d); all retained",
          n_b_before, n_b_after, n_b_pos, n_b_neg),
  sprintf("Panel c rows: before %d after NA drop %d", n_c_before, n_c_after),
  sprintf("Panel d CellChat table: %d pathways; kept TOP10 = %s",
          n_cc_before, paste(TOP10, collapse = ", ")),
  sprintf("Panel d original ratios: %d rows; finite plotted %d; composition-matched series omitted by request",
          n_ratio_raw, n_ratio_plotted),
  sprintf("KEGG ORA Mac: %d terms tested, %d with FDR < 0.05 overall",
          nrow(kegg_mac_all), sum(kegg_mac_all$FDR05)),
  sprintf("KEGG ORA Mg: %d terms tested, %d with FDR < 0.05 overall",
          nrow(kegg_mg_all), sum(kegg_mg_all$FDR05)),
  sprintf("Panel e display pool: %d signaling terms, %d with FDR < 0.05; top 10 by FDR plotted",
          nrow(kegg_mac), n_fdr_mac),
  sprintf("Panel f display pool: %d signaling terms, %d with FDR < 0.05; top 10 by FDR plotted",
          nrow(kegg_mg), n_fdr_mg),
  sprintf("Panel e best signaling FDR: %.3g (%s)",
          min(kegg_mac$p.adjust),
          kegg_mac$Description[which.min(kegg_mac$p.adjust)]),
  sprintf("Panel f best signaling FDR: %.3g (%s)",
          min(kegg_mg$p.adjust),
          kegg_mg$Description[which.min(kegg_mg$p.adjust)]),
  sprintf("Locked-gene table: before %d rows; kept Mac/Mg %d; dropped T/B %d by request",
          n_lock_before, n_lock_after, n_lock_dropped),
  sprintf("Locked Mac genes: %s",
          paste(unique(lock_mm$gene[lock_mm$celltype == "Mac"]), collapse = ", ")),
  sprintf("Locked Mg genes: %s",
          paste(unique(lock_mm$gene[lock_mm$celltype == "Mg"]), collapse = ", ")),
  "KEGG: cached enrichKEGG ORA; BH FDR from the full KEGG test",
  "Signaling restriction in e/f is a display filter; FDR values were not recomputed on the subset",
  "Panel d log10 axis: ratios <= 0 or non-finite are not drawn (none in the TOP10 original set)",
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
)
writeLines(stats, file.path(OUTDIR, "Figure6_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. 6 | Vascular-contact CNS immune cells rewire original outgoing signalling and myeloid programmes.",
  "",
  "a, Flow-cytometric gates used to separate parenchyma-associated (NEG, 58.3%) from vascular-contact (POS, 15.4%) CNS immune cells. Extra CD45 is the intravenous surface label; Intra CD45 is the intracellular stain. Gate percentages are the native cytometer statistics on the plot.",
  sprintf("b, UMAP of CNS immune cells split by POS and NEG, coloured by cell type (n = %s POS; %s NEG). All cells are plotted; point marks are rasterised at 600 dpi with vector axes.",
          format(n_b_pos, big.mark = ","), format(n_b_neg, big.mark = ",")),
  "c, Cell-type proportions within each group.",
  "d, CellChat outgoing POS/NEG ratio for the ten pathways with the highest total outgoing strength (GALECTIN to CSF), computed from the original (not composition-matched) matrices. The x axis is log10-scaled; the vertical line marks a ratio of 1.",
  sprintf("e, KEGG over-representation of macrophage POS-versus-NEG DEGs (padj < 0.001 and |log2 FC| >= 1 after removal of Rpl/Rps, mt and Gm genes), restricted to signalling-pathway terms. Over-representation was tested against all %d KEGG terms; the ten most significant signalling terms are shown. BH-adjusted P values are from the full test. No macrophage signalling term reaches FDR < 0.05 (lowest FDR = %.2f), so panel e is exploratory.",
          nrow(kegg_mac_all), min(kegg_mac$p.adjust)),
  sprintf("f, KEGG over-representation of microglia DEGs, plotted and filtered as in e. Black outlines mark the %d displayed signalling terms that reach FDR < 0.05.",
          n_fdr_mg_top),
  "g, Locked macrophage genes: log2 fold change (POS / NEG) in all cells (grey) versus after doublet/hybrid cell filtering (blue).",
  "h, Locked microglia genes, plotted as in g.",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "Figure6_legend.txt"))

writeLines(c(
  "Figure 6 | new folder Fig6_Nature_v2 (does not overwrite Fig6_Nature_v1).",
  "Reuse level: structural adaptation of Fig6_Nature_v1 a-c and KEGG, plus original-only CellChat ratios and Mac/Mg locked-gene dots.",
  "Panel a: FACS TIFF crop geometry unchanged from the prior Fig. 6 build.",
  "Panel b: all cells from the UMAP table; scattermore rasterisation of marks only.",
  "Panel c: same observations as the source proportion table.",
  "Panel d: original POS/NEG outgoing ratios for GALECTIN through CSF; composition-matched series omitted by request.",
  "Panels e/f: cached KEGG ORA; signalling-name display filter after testing (previous Fig. 6 f/h).",
  "Panels g/h: locked-gene table restricted to Mac and Mg (T and B omitted by request).",
  "Backend: R exclusive."
), file.path(OUTDIR, "Figure6_provenance.txt"))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
base <- file.path(OUTDIR, "Figure6")

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
  "Build notes (full visual QA is in Figure6_QA.txt)",
  "- Backend: R exclusive",
  "- a: FACS TIFF cropped to the plot rectangle; vector NEG/POS/CD45 labels",
  "- b: scattermore raster marks; all cells retained",
  "- c/d/e/f/g/h: native ggplot vector",
  "- d: original CellChat POS/NEG ratios only; log10 x; TOP10 pathways",
  "- e/f: cached enrichKEGG, signalling display filter",
  "- g/h: locked Mac/Mg genes; all cells vs cell-filtered",
  sprintf("- Exports: PDF, SVG, PNG, TIFF at %d x %d mm, raster 600 dpi",
          FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "Figure6_build_notes.txt"))

cat("DONE\n")
