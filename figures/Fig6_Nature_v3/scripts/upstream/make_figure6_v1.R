#!/usr/bin/env Rscript
# =============================================================================
# Fig. 6 v1 — Nature multipanel (Fig6_Nature_v1)
# New folder; does not overwrite any earlier figure package.
#
# Core conclusion
#   Vascular-contact (POS) versus parenchyma-associated (NEG) CNS immune
#   cells, gated by Extra CD45 versus Intra CD45, differ in composition
#   and outgoing signalling. Macrophage and microglia DEGs that belong
#   to the top 10 outgoing pathways (GALECTIN to CSF) mark that shift,
#   and the same DEG lists enrich distinct KEGG signalling programmes.
#
# Evidence chain (one unique question per panel)
#   a  How were IV-labelled CNS cells gated as NEG vs POS?          FACS
#   b  Do POS and NEG occupy different cell-type space?             UMAP
#   c  Does composition shift with IV contact?                      bars
#   d  Which top-10 outgoing pathways (GALECTIN to CSF) dominate?   CellChat
#   e  Which Mac DEGs belong to those pathways?                     volcano
#   f  What KEGG signalling programmes do Mac DEGs over-represent?  ORA
#   g  Which Mg DEGs belong to those pathways?                      volcano
#   h  What KEGG signalling programmes do Mg DEGs over-represent?   ORA
#
# Changed relative to the prior CNS signalling figure
#   a  supplied FACS TIFF, cropped to the plot rectangle and annotated
#      NEG / POS / Extra CD45 / Intra CD45 as vector text
#   d  heatmap and side bars restricted to the dotted (top 10) rows;
#      side bars share the heatmap y-scale so each bar sits on its row
#
# Archetype : asymmetric mixed-modality
# Backend   : R exclusive (ggplot2 / patchwork / ggrepel / scattermore)
# Size      : 183 x 235 mm (Nature double column)
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
  file.path(IV_FIG_ROOT, "Fig6_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
SRCFIG  <- file.path(OUTDIR, "source_figures")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)
dir.create(SRCFIG, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 235
width_mm = 183
height_mm = 235
BASE <- 6.5
FAM  <- "Arial"
PT5  <- 5.2
LABEL_PT <- 5.2 / ggplot2::.pt
ARROW_COL <- "#2E75B6"

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
dir_cols <- c("Up in POS" = "#E64B35", "Up in NEG" = "#3C5488", "NS" = "grey75")
HM_COLS <- c("#FFFFCC", "#FD8D3C", "#800026")

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

drop_ribo_mt_gm <- function(genes) {
  grepl("^(Rpl|Rps)", genes) |
    grepl("^mt-", genes, ignore.case = TRUE) |
    grepl("^MT-", genes) |
    grepl("^Gm[0-9]", genes)
}

wrap_path <- function(x, width = 28) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

kegg_display <- function(x) {
  x <- sub("\\s+signaling pathway\\b", " signaling", x)
  wrap_path(x, width = 32)
}

# -----------------------------------------------------------------------------
# Top 10 outgoing pathways (GALECTIN → CSF) and CellChat gene sets
# -----------------------------------------------------------------------------
cc_long <- read.csv(file.path(SRC, "Fig4D_cellchat_outgoing_long.csv"),
                    stringsAsFactors = FALSE)
n_cc_before <- nrow(cc_long)
cc_long <- cc_long[!is.na(cc_long$strength), ]
n_cc_after <- nrow(cc_long)

path_tot <- cc_long %>%
  group_by(pathway) %>%
  summarise(total_outgoing = sum(strength),
            POS = sum(strength[group == "EAE_POS"]),
            NEG = sum(strength[group == "EAE_NEG"]),
            .groups = "drop") %>%
  arrange(desc(total_outgoing))
TOP10 <- head(path_tot$pathway, 10)
stopifnot(identical(TOP10[1], "GALECTIN"), identical(TOP10[10], "CSF"))

gs <- read.csv(file.path(SRC, "CellChat_pathway_gene_sets.csv"),
               stringsAsFactors = FALSE)
gs10 <- gs[gs$pathway %in% TOP10, ]
top10_genes <- unique(gs10$gene)
gene_to_pw <- gs10 %>%
  group_by(gene) %>%
  summarise(
    pathways = paste(unique(pathway), collapse = ";"),
    roles = paste(unique(role), collapse = ";"),
    .groups = "drop"
  )

write.csv(path_tot, file.path(PLOTTED, "SourceData_Fig6d_pathway_totals.csv"),
          row.names = FALSE)
write.csv(path_tot[path_tot$pathway %in% TOP10, ],
          file.path(PLOTTED, "SourceData_Fig6d_top10_pathway_totals.csv"),
          row.names = FALSE)
write.csv(gs10, file.path(PLOTTED, "SourceData_Fig6_top10_pathway_genes.csv"),
          row.names = FALSE)
cat("TOP10 outgoing:", paste(TOP10, collapse = " | "), "\n")
cat("TOP10 L/R genes:", length(top10_genes), "\n")

# -----------------------------------------------------------------------------
# a — FACS plate: crop supplied TIFF, vector NEG/POS and CD45 axis labels
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
tif_a <- file.path(SRCFIG, "Fig6a_FACS_source.tif")
stopifnot(file.exists(tif_a))
im_raw <- magick::image_background(magick::image_read(tif_a), "white",
                                   flatten = TRUE)
# Plot rectangle measured on the supplied TIFF (black margin + original
# arrows excluded). No contrast, gamma or recolouring.
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
# b — UMAP, all cells retained; rasterised marks only
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
# c — grouped bar proportions
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
# d — CellChat outgoing, TOP10 rows only; side bars aligned to heatmap rows
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
cc <- readRDS(file.path(SRC, "Fig4D_cellchat_outgoing.rds"))
path_ord_all <- names(sort(cc$side_pos + cc$side_neg, decreasing = TRUE))
stopifnot(identical(unname(path_ord_all[1:10]), unname(TOP10)))
path_ord <- TOP10
ct_ord <- cc$celltypes
ct_ord_lab <- unname(ct_lab[ct_ord])

subset_mat <- function(mat) {
  m <- mat[, path_ord, drop = FALSE]
  rownames(m) <- ct_lab[rownames(m)]
  colnames(m) <- colnames(m)
  m[ct_ord_lab, path_ord, drop = FALSE]
}
mat_pos <- subset_mat(cc$mat_pos)
mat_neg <- subset_mat(cc$mat_neg)
side_pos <- cc$side_pos[path_ord]
side_neg <- cc$side_neg[path_ord]
top_pos <- rowSums(mat_pos)[ct_ord_lab]
top_neg <- rowSums(mat_neg)[ct_ord_lab]
hm_max   <- max(c(mat_pos, mat_neg), na.rm = TRUE)
top_max  <- max(c(top_pos, top_neg), na.rm = TRUE)
side_max <- max(c(side_pos, side_neg), na.rm = TRUE)
n_path_dropped <- length(path_ord_all) - length(path_ord)
cat(sprintf("Panel d: kept %d / %d pathways (dropped %d below CSF)\n",
            length(path_ord), length(path_ord_all), n_path_dropped))
cat(sprintf("shared scales hm_max=%.3g top_max=%.3g side_max=%.3g\n",
            hm_max, top_max, side_max))

make_cc_panel <- function(mat, topv, sidev, title) {
  long <- as.data.frame(mat) %>%
    rownames_to_column("celltype") %>%
    pivot_longer(-celltype, names_to = "pathway", values_to = "strength")
  path_lab <- paste0("\u25CF ", path_ord)
  names(path_lab) <- path_ord
  long$celltype <- factor(long$celltype, levels = ct_ord_lab)
  long$pathway <- factor(long$pathway, levels = rev(path_ord))
  top_df <- data.frame(
    celltype = factor(names(topv), levels = ct_ord_lab),
    strength = as.numeric(topv)
  )
  side_df <- data.frame(
    pathway = factor(names(sidev), levels = rev(path_ord)),
    strength = as.numeric(sidev)
  )
  top_breaks <- pretty(c(0, top_max), n = 3)
  side_breaks <- pretty(c(0, side_max), n = 3)

  p_hm <- ggplot(long, aes(celltype, pathway, fill = strength)) +
    geom_tile(width = 1, height = 1, colour = NA) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(labels = path_lab, expand = c(0, 0)) +
    scale_fill_gradientn(
      colours = HM_COLS, limits = c(0, hm_max),
      breaks = pretty(c(0, hm_max), n = 4),
      name = "Strength", oob = scales::squish, guide = "none"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = BASE, base_family = FAM) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = PT5,
                                 family = FAM, colour = "black"),
      axis.text.y = element_text(size = PT5, family = FAM, colour = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5,
                                family = FAM),
      plot.margin = margin(1, 1, 1, 1)
    )

  p_top <- ggplot(top_df, aes(celltype, strength, fill = celltype)) +
    geom_col(width = 0.92, colour = NA) +
    scale_fill_manual(values = ct_cols, guide = "none") +
    scale_x_discrete(expand = c(0, 0), drop = FALSE) +
    scale_y_continuous(limits = c(0, top_max), breaks = top_breaks,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Outgoing") +
    theme_classic(base_size = PT5, base_family = FAM) +
    theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.y = element_text(size = PT5, family = FAM),
          axis.text.y = element_text(size = PT5, family = FAM),
          plot.margin = margin(1, 1, 0, 1)
    )

  p_side <- ggplot(side_df, aes(strength, pathway)) +
    geom_col(width = 0.92, fill = "#4D4D4D", colour = NA) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    scale_x_continuous(limits = c(0, side_max), breaks = side_breaks,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Outgoing", y = NULL) +
    theme_classic(base_size = PT5, base_family = FAM) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = PT5, family = FAM),
      axis.text.x = element_text(size = PT5, family = FAM),
      plot.margin = margin(1, 1, 1, 0)
    )

  # Top-right spacer matches p_top height so p_side shares the heatmap rows.
  p_top + plot_spacer() + p_hm + p_side +
    plot_layout(ncol = 2, nrow = 2, byrow = TRUE,
                heights = c(0.28, 1), widths = c(1, 0.28))
}

pD_pos <- make_cc_panel(mat_pos, top_pos, side_pos, "POS")
pD_neg <- make_cc_panel(mat_neg, top_neg, side_neg, "NEG")

p_cbar <- ggplot(data.frame(x = 0, y = 0, z = c(0, hm_max)), aes(x, y, fill = z)) +
  geom_point(shape = 21, alpha = 0, size = 0) +
  scale_fill_gradientn(
    colours = HM_COLS, limits = c(0, hm_max),
    breaks = pretty(c(0, hm_max), n = 4),
    name = "Strength",
    guide = guide_colourbar(
      barwidth = unit(2.2, "mm"), barheight = unit(22, "mm"),
      ticks.colour = "black", frame.colour = "black"
    )
  ) +
  theme_void(base_family = FAM) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = PT5, family = FAM),
    legend.text = element_text(size = PT5, family = FAM),
    plot.margin = margin(20, 2, 20, 0)
  )

pD <- wrap_elements(
  panel = (pD_pos | pD_neg | p_cbar) + plot_layout(widths = c(1, 1, 0.16))
) +
  labs(title = "Outgoing signalling (GALECTIN to CSF)", tag = "d") +
  theme(
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5,
                              family = FAM),
    plot.margin = margin(2, 2, 2, 2)
  ) + tag_theme

cc_long10 <- cc_long[cc_long$pathway %in% TOP10, ]
write.csv(cc_long10, file.path(PLOTTED, "SourceData_Fig6d_cellchat_top10_long.csv"),
          row.names = FALSE)
write.csv(cc_long, file.path(PLOTTED, "SourceData_Fig6d_cellchat_long.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Volcano helper — all genes plotted; labels = TOP10 pathway DEGs only
# -----------------------------------------------------------------------------
prep_volcano <- function(path, ct) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  n_before <- nrow(d)
  d <- d[!is.na(d$log2FC) & !is.na(d$padj), ]
  n_finite <- nrow(d)
  d <- d[!drop_ribo_mt_gm(d$gene), ]
  n_after <- nrow(d)
  d$neglog10_padj <- -log10(pmax(d$padj, 1e-300))
  d$sig <- d$padj < 0.001 & abs(d$log2FC) >= 1
  d$direction <- "NS"
  d$direction[d$sig & d$log2FC >= 1] <- "Up in POS"
  d$direction[d$sig & d$log2FC <= -1] <- "Up in NEG"
  d$direction <- factor(d$direction, levels = c("Up in POS", "Up in NEG", "NS"))
  d$label_me <- d$sig & d$gene %in% top10_genes
  d <- dplyr::left_join(d, gene_to_pw, by = "gene")
  d$celltype <- ct
  cat(ct, "volcano: before", n_before, " finite", n_finite,
      " after ribo/mt/Gm drop", n_after,
      " sig", sum(d$sig), " labelled", sum(d$label_me), "\n")
  attr(d, "n_before") <- n_before
  attr(d, "n_finite") <- n_finite
  attr(d, "n_after") <- n_after
  d
}

make_volcano <- function(d, title) {
  lab <- d[d$label_me, ]
  ymax <- max(d$neglog10_padj, na.rm = TRUE)
  xmax <- max(abs(d$log2FC), na.rm = TRUE)
  ggplot(d, aes(log2FC, neglog10_padj)) +
    geom_point(aes(colour = direction),
               data = d[!d$label_me, ],
               size = 0.45, stroke = 0, alpha = 0.7) +
    geom_point(aes(colour = direction),
               data = lab,
               size = 1.15, stroke = 0.25, alpha = 1) +
    geom_vline(xintercept = c(-1, 1), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.001), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_text_repel(
      data = lab, aes(label = gene),
      size = LABEL_PT, colour = "black", family = FAM,
      max.overlaps = Inf, box.padding = 0.18, point.padding = 0.12,
      min.segment.length = 0, segment.colour = "grey35",
      segment.size = 0.25, seed = 42, show.legend = FALSE
    ) +
    scale_colour_manual(values = dir_cols, name = NULL, drop = FALSE) +
    coord_cartesian(xlim = c(-xmax * 1.12, xmax * 1.12),
                    ylim = c(0, ymax * 1.10)) +
    labs(x = "log2 FC (POS / NEG)", y = "-log10(Padj)",
         title = title,
         subtitle = sprintf("Labels: GALECTIN-CSF DEGs (n = %d)",
                            nrow(lab))) +
    theme(legend.position = "bottom", legend.key.size = unit(2.0, "mm"))
}

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
          legend.key.size = unit(2.2, "mm"))
}

cat("Volcano + KEGG ...\n")
vol_mac <- prep_volcano(file.path(SRC, "Fig4_DEG_Mac_POS_vs_NEG.csv"), "Mac")
vol_mg  <- prep_volcano(file.path(SRC, "Fig4_DEG_Mg_POS_vs_NEG.csv"), "Mg")

kegg_mac_all <- read.csv(file.path(SRC, "Fig5_KEGG_ORA_Mac.csv"),
                         stringsAsFactors = FALSE)
kegg_mg_all  <- read.csv(file.path(SRC, "Fig5_KEGG_ORA_Mg.csv"),
                         stringsAsFactors = FALSE)
kegg_mac <- signaling_only(kegg_mac_all, "Mac")
kegg_mg  <- signaling_only(kegg_mg_all, "Mg")
stopifnot(nrow(kegg_mac) >= 10, nrow(kegg_mg) >= 10)

write.csv(kegg_mac_all, file.path(PLOTTED, "SourceData_Fig6f_KEGG_Mac_allterms.csv"),
          row.names = FALSE)
write.csv(kegg_mg_all,  file.path(PLOTTED, "SourceData_Fig6h_KEGG_Mg_allterms.csv"),
          row.names = FALSE)
write.csv(kegg_mac, file.path(PLOTTED, "SourceData_Fig6f_KEGG_Mac_signaling.csv"),
          row.names = FALSE)
write.csv(kegg_mg,  file.path(PLOTTED, "SourceData_Fig6h_KEGG_Mg_signaling.csv"),
          row.names = FALSE)
write.csv(vol_mac, file.path(PLOTTED, "SourceData_Fig6e_Mac_volcano.csv"),
          row.names = FALSE)
write.csv(vol_mg,  file.path(PLOTTED, "SourceData_Fig6g_Mg_volcano.csv"),
          row.names = FALSE)
write.csv(vol_mac[vol_mac$label_me, ],
          file.path(PLOTTED, "SourceData_Fig6e_Mac_labelled.csv"),
          row.names = FALSE)
write.csv(vol_mg[vol_mg$label_me, ],
          file.path(PLOTTED, "SourceData_Fig6g_Mg_labelled.csv"),
          row.names = FALSE)

pE <- make_volcano(vol_mac, "Mac DEG") + labs(tag = "e") + tag_theme
pF <- make_kegg(kegg_mac, "Mac KEGG") + labs(tag = "f") + tag_theme
pG <- make_volcano(vol_mg, "Mg DEG") +
  theme(legend.position = "none") + labs(tag = "g") + tag_theme
pH <- make_kegg(kegg_mg, "Mg KEGG") + labs(tag = "h") + tag_theme

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- (pA | pB | pC) + plot_layout(widths = c(1.28, 1.42, 1.05))
row2 <- pD
row3 <- wrap_elements((pE | pF) + plot_layout(widths = c(1.05, 1.10)))
row4 <- wrap_elements((pG | pH) + plot_layout(widths = c(1.05, 1.10)))
fig <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(1.00, 1.00, 1.14, 1.20))

# -----------------------------------------------------------------------------
# Stats / legend / provenance
# -----------------------------------------------------------------------------
n_lab_mac <- sum(vol_mac$label_me)
n_lab_mg  <- sum(vol_mg$label_me)
n_fdr_mac <- sum(kegg_mac$FDR05)
n_fdr_mg  <- sum(kegg_mg$FDR05)
n_fdr_mac_top <- sum((kegg_mac %>% arrange(p.adjust) %>% slice_head(n = 10))$FDR05)
n_fdr_mg_top  <- sum((kegg_mg  %>% arrange(p.adjust) %>% slice_head(n = 10))$FDR05)

stats <- c(
  "Figure 6 statistics / integrity (Fig6_Nature_v1)",
  sprintf("Panel a raster after crop: %d x %d px; crop 429x405+75+69 of the supplied TIFF; no contrast change",
          n_a_px[1], n_a_px[2]),
  "Panel a vector overlays: NEG, POS, Extra CD45, Intra CD45. Gate percentages 58.3 and 15.4 remain native cytometer glyphs in the TIFF.",
  sprintf("Panel b cells: before %d after NA drop %d (POS %d, NEG %d); all retained for plotting",
          n_b_before, n_b_after, n_b_pos, n_b_neg),
  sprintf("Panel c rows: before %d after NA drop %d", n_c_before, n_c_after),
  sprintf("Panel d pathways kept: %d / %d (TOP10 = %s); dropped %d unmarked rows below CSF",
          length(path_ord), length(path_ord_all),
          paste(TOP10, collapse = ", "), n_path_dropped),
  sprintf("Panel d shared scales (TOP10 subset): hm_max=%.3g top_max=%.3g side_max=%.3g",
          hm_max, top_max, side_max),
  "Panel d top bars are rowSums of the displayed (TOP10) matrix so they match the heatmap columns.",
  "Panel d side bars use the same pathway factor and expand=0 as the heatmap tiles.",
  sprintf("Panel d CellChat long table: before %d after NA drop %d; plotted rows %d",
          n_cc_before, n_cc_after, nrow(cc_long10)),
  sprintf("Volcano Mac: before %d finite %d after ribo/mt/Gm %d; sig %d; labelled %d",
          attr(vol_mac, "n_before"), attr(vol_mac, "n_finite"),
          attr(vol_mac, "n_after"), sum(vol_mac$sig), n_lab_mac),
  sprintf("Volcano Mg: before %d finite %d after ribo/mt/Gm %d; sig %d; labelled %d",
          attr(vol_mg, "n_before"), attr(vol_mg, "n_finite"),
          attr(vol_mg, "n_after"), sum(vol_mg$sig), n_lab_mg),
  sprintf("Mac labelled genes: %s",
          paste(sort(vol_mac$gene[vol_mac$label_me]), collapse = ", ")),
  sprintf("Mg labelled genes: %s",
          paste(sort(vol_mg$gene[vol_mg$label_me]), collapse = ", ")),
  sprintf("KEGG ORA Mac: %d terms tested, %d with FDR < 0.05 overall",
          nrow(kegg_mac_all), sum(kegg_mac_all$FDR05)),
  sprintf("KEGG ORA Mg: %d terms tested, %d with FDR < 0.05 overall",
          nrow(kegg_mg_all), sum(kegg_mg_all$FDR05)),
  sprintf("Panel f display pool: %d signaling terms, %d with FDR < 0.05; top 10 by FDR plotted",
          nrow(kegg_mac), n_fdr_mac),
  sprintf("Panel h display pool: %d signaling terms, %d with FDR < 0.05; top 10 by FDR plotted",
          nrow(kegg_mg), n_fdr_mg),
  sprintf("Panel f best signaling FDR: %.3g (%s)",
          min(kegg_mac$p.adjust),
          kegg_mac$Description[which.min(kegg_mac$p.adjust)]),
  sprintf("Panel h best signaling FDR: %.3g (%s)",
          min(kegg_mg$p.adjust),
          kegg_mg$Description[which.min(kegg_mg$p.adjust)]),
  "Volcano thresholds: padj < 0.001 and |log2FC| >= 1 after dropping Rpl/Rps/mt/Gm",
  "KEGG: cached enrichKEGG ORA of those DEGs; BH FDR from the full KEGG test",
  "Signaling restriction in f/h is a display filter; FDR values were not recomputed on the subset",
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
)
writeLines(stats, file.path(OUTDIR, "Figure6_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. 6 | Vascular-contact CNS immune cells rewire outgoing signalling and myeloid programmes.",
  "",
  "a, Flow-cytometric gates used to separate parenchyma-associated (NEG, 58.3%) from vascular-contact (POS, 15.4%) CNS immune cells. Extra CD45 is the intravenous surface label; Intra CD45 is the intracellular stain. Gate percentages are the native cytometer statistics on the plot.",
  sprintf("b, UMAP of CNS immune cells split by POS and NEG, coloured by cell type (n = %s POS; %s NEG). All cells are plotted; point marks are rasterised at 600 dpi with vector axes.",
          format(n_b_pos, big.mark = ","), format(n_b_neg, big.mark = ",")),
  "c, Cell-type proportions within each group.",
  "d, CellChat outgoing signalling strength for POS and NEG on shared colour and bar scales, restricted to the ten pathways with the highest total outgoing strength (GALECTIN to CSF). Rows are ordered by that total. Side bars use the same row order and discrete scale as the heatmap. Column bars are the per-cell-type sum of the displayed pathways. Ligand and receptor genes from these pathways are the only volcano labels in e and g.",
  "e, Macrophage differential expression (POS versus NEG). Points show every tested gene after removal of Rpl/Rps, mt and Gm genes. Dashed lines mark padj = 0.001 and |log2 FC| = 1. Only significant DEGs that belong to the panel-d pathways (GALECTIN to CSF) are labelled.",
  sprintf("f, KEGG over-representation of the macrophage DEGs defined in e, restricted to signalling-pathway terms (KEGG names containing \"signaling\"). Over-representation was tested against all %d KEGG terms returned for this gene set and the ten most significant signalling terms are displayed; the BH-adjusted P values shown are from the full test and were not recalculated on the displayed subset. Colour is -log10(FDR), point size is the number of DEGs in the term, and a black outline would mark FDR < 0.05. No macrophage signalling term reaches FDR < 0.05 (lowest FDR = %.2f), so panel f is exploratory.",
          nrow(kegg_mac_all), min(kegg_mac$p.adjust)),
  "g, Microglia differential expression plotted as in e, with the same GALECTIN-to-CSF outgoing-pathway label rule.",
  sprintf("h, KEGG over-representation of the microglia DEGs defined in g, plotted and filtered as in f. Black outlines mark the %d displayed signalling terms that reach FDR < 0.05.",
          n_fdr_mg_top),
  "Source data are provided as a Source Data file. The full unfiltered KEGG results behind f and h accompany the filtered tables."
), file.path(OUTDIR, "Figure6_legend.txt"))

writeLines(c(
  "Figure 6 | new folder Fig6_Nature_v1 (does not overwrite any earlier package).",
  "Reuse level: structural adaptation of the prior CNS CellChat figure, with a new FACS panel and a TOP10-only outgoing heatmap.",
  "Panel a: supplied flow-cytometry TIFF cropped to the plot rectangle (geometry 429x405+75+69); no recolouring, contrast or gamma change. NEG, POS, Extra CD45 and Intra CD45 are vector overlays. Gate percentages 58.3 and 15.4 remain in the raster.",
  "Panel b: all cells from the UMAP table; scattermore rasterisation of marks only.",
  "Panel c: same observations as the source proportion table.",
  "Panel d: same CellChat matrices; display restricted to GALECTIN through CSF; side bars aligned to heatmap rows; top bars recomputed on the displayed pathways.",
  "Panels e/g: same DEG tables; labels restricted to CellChat ligand/receptor genes of those ten pathways.",
  "Panels f/h: cached KEGG ORA of those DEGs; signalling-name display filter applied after testing.",
  "Backend: R exclusive."
), file.path(OUTDIR, "Figure6_provenance.txt"))

# -----------------------------------------------------------------------------
# Export PDF / SVG / PNG / TIFF
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
  "- Backend: R exclusive (ggplot2 / patchwork / ggrepel / scattermore / magick)",
  "- a: supplied FACS TIFF cropped to the plot rectangle; vector NEG/POS/CD45 labels",
  "- b: native axes; scattermore raster marks; all cells retained",
  "- c/d/e/f/g/h: native ggplot vector",
  "- d: TOP10 pathways only; shared hm_max / top_max / side_max; side bars aligned to heatmap rows",
  "- e/g: all genes plotted; labels restricted to TOP10 outgoing L/R DEGs; ggrepel, no opaque boxes",
  "- f/h: cached enrichKEGG ORA, signalling-name display filter, top 10 by FDR",
  "- Axis titles are plain text (no plotmath) so subscripts cannot fall below 5 pt",
  sprintf("- Exports: PDF, SVG, PNG, TIFF at %d x %d mm, raster 600 dpi",
          FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "Figure6_build_notes.txt"))

cat("DONE\n")
