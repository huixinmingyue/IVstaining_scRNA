#!/usr/bin/env Rscript
# =============================================================================
# Fig. S5 v1 - Oligo-CD45 conjugate QC and CITE-seq B-cell annotation
#
# Core conclusion
#   Oligo-CD45 conjugates are intact on gel, Harmony B cells split into
#   follicular, marginal-zone, immature and myeloid-contaminated islands,
#   and Proteinfeatures continuously shifts subset probability toward the
#   myeloid-contaminated island.
#
# Evidence chain
#   a  Oligo-CD45 conjugate gels (Coomassie and fluorescence)
#   b  Harmony UMAP of B cells (resolution 0.2)
#   c  B-lineage and myeloid marker dot plot
#   d  Logistic P(Bcelltype) versus Proteinfeatures
#
# Archetype : asymmetric mixed-modality
# Backend   : R exclusive
# Size      : 183 x 175 mm
# New folder — does not overwrite any prior package
# =============================================================================

options(warn = 1)
set.seed(42)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(png)
  library(grid)
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS5_Nature_v1")
}
SRC <- file.path(OUTDIR, "source_data")
SDP <- file.path(OUTDIR, "source_data_plotted")
QA  <- file.path(OUTDIR, "qa")
dir.create(SDP, showWarnings = FALSE)
dir.create(QA, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 175
width_mm <- 183
height_mm <- 175
BASE <- 6.5
FAM  <- "Arial"
PT <- function(x) x / .pt

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args <- list(reg()$Helvetica); names(args) <- FAM
      try(do.call(reg, args), silent = TRUE)
    }
  }
}

BT_LV <- c(
  "Follicular B cells (Gpr183+)",
  "Follicular B cells",
  "Marginal zone B cells",
  "Follicular B cells [3]",
  "Myeloid-contaminated",
  "Immature / pre-B cells"
)
BT_LAB <- c(
  "Follicular B cells (Gpr183+)" = "Follicular B (Gpr183+)",
  "Follicular B cells"           = "Follicular B",
  "Marginal zone B cells"        = "Marginal zone B",
  "Follicular B cells [3]"       = "Follicular B [3]",
  "Myeloid-contaminated"         = "Myeloid-contaminated",
  "Immature / pre-B cells"       = "Immature / pre-B"
)
BT_COL <- c(
  "Follicular B cells (Gpr183+)" = "#E07A5F",
  "Follicular B cells"           = "#C4A35A",
  "Marginal zone B cells"        = "#2A9D8F",
  "Follicular B cells [3]"       = "#5B8F9E",
  "Myeloid-contaminated"         = "#3D5A80",
  "Immature / pre-B cells"       = "#B2473E"
)

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title = element_text(size = BASE, colour = "black"),
    axis.text = element_text(size = BASE - 0.5, colour = "black"),
    legend.title = element_text(size = 5.6),
    legend.text = element_text(size = 5.4),
    legend.key.size = unit(2.4, "mm"),
    legend.background = element_blank(),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.margin = margin(1.5, 2, 1.5, 2, "mm"),
    panel.grid = element_blank()
  )
tag_it <- function(p, t) {
  p + labs(tag = t) +
    theme(plot.tag = element_text(size = 8, face = "bold", family = FAM,
                                  hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}

# -----------------------------------------------------------------------------
# a  oligo-CD45 conjugate gels
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
gel <- png::readPNG(file.path(SRC, "panel_a_raw.png"))
lane_lab <- c("M1/42-Oligo", "30-F11-Oligo", "M1/42", "30-F11",
              "M1/42-Oligo", "30-F11-Oligo")
# Leave a right gutter so 40-degree labels are not clipped by panel b.
gel_x0 <- 0.00
gel_x1 <- 0.84
xs <- gel_x0 + c(0.30, 0.415, 0.525, 0.630, 0.750, 0.865) * (gel_x1 - gel_x0)
pA <- ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x, y)) +
  geom_blank() +
  annotation_raster(gel, xmin = gel_x0, xmax = gel_x1, ymin = 0.00, ymax = 0.78,
                    interpolate = TRUE) +
  annotate("segment", x = xs[1] - 0.02, xend = xs[4] + 0.02,
           y = 0.955, yend = 0.955, linewidth = 0.4, colour = "black") +
  annotate("segment", x = xs[5] - 0.02, xend = xs[6] + 0.03,
           y = 0.955, yend = 0.955, linewidth = 0.4, colour = "black") +
  annotate("text", x = mean(xs[1:4]), y = 0.985, label = "1", size = PT(8),
           fontface = "bold", family = FAM, vjust = 1) +
  annotate("text", x = mean(xs[5:6]), y = 0.985, label = "2", size = PT(8),
           fontface = "bold", family = FAM, vjust = 1) +
  annotate("text", x = xs, y = 0.800, label = lane_lab, size = PT(5.2),
           family = FAM, angle = 40, hjust = 0, vjust = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void(base_family = FAM) +
  theme(plot.margin = margin(2, 1, 1, 1, "mm"))

# -----------------------------------------------------------------------------
# b  Harmony UMAP
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
df_b <- readRDS(file.path(SRC, "panel_b_umap.rds"))
stopifnot(nrow(df_b) > 0)
df_b$Bcelltype <- factor(df_b$Bcelltype, levels = BT_LV)
n_b <- nrow(df_b)
lab_b <- df_b %>%
  group_by(Bcelltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop") %>%
  mutate(lab = unname(BT_LAB[as.character(Bcelltype)]))

pB <- ggplot(df_b, aes(umap_1, umap_2, colour = Bcelltype)) +
  geom_point(size = 0.18, stroke = 0, shape = 16, alpha = 0.82) +
  geom_text_repel(
    data = lab_b, aes(umap_1, umap_2, label = lab),
    colour = "black", size = PT(5.2), family = FAM,
    max.overlaps = Inf, box.padding = 0.22, point.padding = 0.12,
    segment.size = 0.2, segment.colour = "grey40",
    min.segment.length = 0, seed = 42, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = BT_COL, guide = "none") +
  coord_equal(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_pub +
  theme(plot.margin = margin(2, 2, 1, 2, "mm"))
cat("  n cells =", n_b, "\n")

# -----------------------------------------------------------------------------
# c  marker dot plot
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
dot <- read.csv(file.path(SRC, "panel_c_dotplot.csv"), stringsAsFactors = FALSE)
ann <- unique(dot[, c("cluster", "Bcelltype")])
ann <- ann[order(as.integer(ann$cluster)), , drop = FALSE]
ylab <- paste0(ann$cluster, ": ", unname(BT_LAB[ann$Bcelltype]))
names(ylab) <- as.character(ann$cluster)
dot$cluster_lab <- factor(unname(ylab[as.character(dot$cluster)]), levels = ylab)
gene_lv <- unique(as.character(dot$gene))
# Preserve extraction order: B-lineage then myeloid.
gene_order <- c("Cr2", "Fcer2a", "Cd19", "Ms4a1", "Cd79a", "Vpreb3", "Igll1",
                "Cd93", "Sox4", "Cd24a", "Vpreb1", "Rag1",
                "Slpi", "Axl", "Mafb", "Hpgd", "Fcna")
gene_order <- gene_order[gene_order %in% gene_lv]
dot$gene <- factor(dot$gene, levels = gene_order)
stopifnot(!anyNA(dot$gene), !anyNA(dot$cluster_lab))

pC <- ggplot(dot, aes(gene, cluster_lab)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_colour_gradient(low = "#D9D9D9", high = "#B2182B",
                        name = "Scaled mean\nexpression",
                        limits = c(-2.5, 2.5), oob = squish) +
  scale_size_area(name = "Expressing\ncells (%)", max_size = 3.2,
                  breaks = c(10, 25, 50, 75), limits = c(0, 100)) +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", size = 5.4),
    axis.text.y = element_text(size = 5.5),
    legend.position = "right",
    legend.box.spacing = unit(1.0, "mm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line = element_blank()
  )

# -----------------------------------------------------------------------------
# d  logistic curves
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
pred <- read.csv(file.path(SRC, "panel_d_logistic.csv"), stringsAsFactors = FALSE)
pred$Bcelltype <- factor(pred$Bcelltype, levels = BT_LV)
stopifnot(!anyNA(pred$Bcelltype), all(is.finite(pred$prob)))

pD <- ggplot(pred, aes(Proteinfeatures, 100 * prob, colour = Bcelltype)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = BT_COL, labels = BT_LAB, name = NULL) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 2, 5, 10, 20, 50)) +
  labs(x = "Protein features", y = "Predicted probability (%)") +
  theme_pub +
  theme(legend.position = "right",
        legend.text = element_text(size = 5.4)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.1)))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
design <- "
AB
AC
DD
"
fig <- tag_it(pA, "a") + tag_it(pB, "b") + tag_it(pC, "c") + tag_it(pD, "d") +
  plot_layout(design = design, widths = c(0.92, 1.22),
              heights = c(1.05, 1.12, 0.95)) &
  theme(plot.margin = margin(1.4, 1.4, 1.2, 1.8, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS5")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())
cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in,
                     family = FAM, bg = "white")
print(fig); invisible(dev.off())
cat("Writing PNG ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())
cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(base, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); invisible(dev.off())

um_cols <- intersect(c("cell", "umap_1", "umap_2", "cluster", "Bcelltype",
                       "Proteinfeatures"), names(df_b))
write.csv(df_b[, um_cols], file.path(SDP, "SourceData_FigS5b_umap.csv"),
          row.names = FALSE)
write.csv(dot, file.path(SDP, "SourceData_FigS5c_dotplot.csv"), row.names = FALSE)
write.csv(pred, file.path(SDP, "SourceData_FigS5d_logistic.csv"), row.names = FALSE)
cat("DONE ->", OUTDIR, "n_b =", n_b, "\n")
