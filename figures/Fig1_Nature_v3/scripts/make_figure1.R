#!/usr/bin/env Rscript
# =============================================================================
# Fig. 1 v3 - Nature-style multipanel figure (R / ggplot2 exclusive)
#
# Core conclusion
#   Intravascular antibody labelling separates intravascular from
#   extravascular splenic immune cells, and scRNA-seq of the two compartments
#   shows a composition shift on RNA-rule CD4 T / CD8 T identities.
#
# Evidence chain
#   a  Are lymphoid lineages spatially partitioned in the spleen?     (IF plate)
#   b  Does the myeloid/vascular compartment show the same border?    (IF plate)
#   c  Can the two compartments be gated and sorted for scRNA-seq?    (flow)
#   d  Do both compartments contain the same transcriptional states?  (UMAP)
#   e  Do the compartments differ in composition?                     (proportions)
# Marker support for the identities is Fig. S1a (moved out of Fig. 1).
#
# Archetype : asymmetric mixed-modality figure (image plate + quantitative)
# Backend   : R exclusive (ggplot2 / patchwork / svglite / ragg)
# Exports   : PDF (cairo, editable text), SVG, PNG (600 dpi), TIFF (600 dpi)
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(png)
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

set.seed(42)

OUTDIR <- file.path(IV_FIG_ROOT, "Fig1_Nature_v3")
SRC    <- file.path(OUTDIR, "source_data")
ASSETS <- file.path(OUTDIR, "assets")
QA     <- file.path(OUTDIR, "qa")
dir.create(QA, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 230
BASE     <- 7
FAM      <- "Arial"
ROW1_MM  <- 88
ROW2_MM  <- 50

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

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(0.9, "pt"),
      axis.title      = element_text(size = BASE - 0.5, colour = "black"),
      axis.text       = element_text(size = BASE - 1.0, colour = "black"),
      legend.title    = element_text(size = BASE - 0.5, colour = "black"),
      legend.text     = element_text(size = BASE - 1.0, colour = "black"),
      legend.key.size = unit(2.6, "mm"),
      legend.margin   = margin(0, 0, 0, 1),
      legend.background = element_blank(),
      plot.title      = element_text(size = BASE - 0.5, hjust = 0.5, colour = "black"),
      plot.margin     = margin(1, 1, 1, 1),
      panel.grid      = element_blank(),
      strip.background = element_blank(),
      strip.text      = element_text(size = BASE - 0.5, colour = "black")
    )
)

CT_COL <- c(
  "CD4 T"        = "#377EB8",
  "CD8 T"        = "#4DAF4A",
  "NKT-like"     = "#1B9E77",
  "CD4/CD8 DP"   = "#E6AB02",
  "T other"      = "#A6CEE3",
  "B-T doublet"  = "#7A7A7A",
  "B cell"       = "#E41A1C",
  "cDC2"         = "#FF7F00",
  "NK"           = "#984EA3",
  "Plasma cell"  = "#F781BF",
  "Neutrophil"   = "#A65628"
)
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")

# -----------------------------------------------------------------------------
# Panels a-c : image plates (unchanged from Fig1_Nature_v2)
# -----------------------------------------------------------------------------
cat("Panels a-c (image plates) ...\n")

plate <- function(items, gap = 0.10, pad_l = 0.02, pad_t = 0.075,
                  fill_aspect = NULL) {
  imgs <- lapply(items, function(it) png::readPNG(it$path))
  asp  <- vapply(imgs, function(m) ncol(m) / nrow(m), numeric(1))
  hh   <- vapply(items, function(it) it$h, numeric(1))
  ww   <- hh * asp
  x0   <- cumsum(c(0, head(ww, -1) + gap))
  Htot <- max(hh)
  Wtot <- sum(ww) + gap * (length(ww) - 1)

  p <- ggplot()
  for (i in seq_along(imgs)) {
    y0 <- (Htot - hh[i]) / 2
    p <- p + annotation_raster(grDevices::as.raster(imgs[[i]]),
                               xmin = x0[i], xmax = x0[i] + ww[i],
                               ymin = y0,    ymax = y0 + hh[i],
                               interpolate = TRUE)
  }
  xlo <- -pad_l * Wtot
  xhi <- Wtot
  yhi <- Htot * (1 + pad_t)
  if (!is.null(fill_aspect)) {
    extra <- fill_aspect * yhi - (xhi - xlo)
    if (extra > 0) { xlo <- xlo - extra / 2; xhi <- xhi + extra / 2 }
  }

  p +
    coord_fixed(xlim = c(xlo, xhi), ylim = c(0, yhi), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

pA <- plate(list(list(path = file.path(ASSETS, "Fig1A_IF_B220_CD4_DCIR.png"), h = 1)))
pB <- plate(list(list(path = file.path(ASSETS, "Fig1B_IF_F480_Laminin.png"), h = 1)))
pC <- plate(list(
  list(path = file.path(ASSETS, "Fig1C_flow_gating.png"),  h = 1.00),
  list(path = file.path(ASSETS, "Fig1C_flow_gating2.png"), h = 0.60)
), gap = 0.12, pad_l = 0.02, pad_t = 0.09,
   fill_aspect = FIG_W_MM / ROW2_MM)

# -----------------------------------------------------------------------------
# Panel d : UMAP with RNA-rule CD4 T / CD8 T identities
# T subsets occupy the same island, so colour + legend is used rather than
# overlapping median labels. Every cell is drawn; none is subsampled.
# -----------------------------------------------------------------------------
cat("Panel d (UMAP) ...\n")
df_d <- readRDS(file.path(SRC, "Fig1D_umap_reannotated.rds"))
n_total <- nrow(df_d)
stopifnot(n_total == 15271)
df_d$lineage <- factor(as.character(df_d$lineage), levels = names(CT_COL))
df_d$IV_lab <- factor(as.character(df_d$IV_lab), levels = unname(IV_LAB))

n_by_group <- df_d %>% count(IV_lab, name = "n")
strip_lab  <- setNames(
  sprintf("%s\nn = %s cells", n_by_group$IV_lab, comma(n_by_group$n)),
  as.character(n_by_group$IV_lab)
)
df_d$facet <- factor(strip_lab[as.character(df_d$IV_lab)],
                     levels = unname(strip_lab))

pD <- ggplot(df_d, aes(umap_1, umap_2, colour = lineage)) +
  geom_point(size = 0.22, stroke = 0, shape = 16, alpha = 0.70) +
  scale_colour_manual(values = CT_COL, name = NULL, drop = FALSE) +
  facet_wrap(~ facet, nrow = 1) +
  coord_equal(clip = "off") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1),
                               ncol = 4, byrow = TRUE)) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(0.8, "mm"),
    legend.key.size = unit(2.4, "mm"),
    legend.text = element_text(size = BASE - 1.2),
    legend.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = BASE - 1.0, lineheight = 1.05),
    plot.margin = margin(1, 1, 1, 1)
  )

# -----------------------------------------------------------------------------
# Panel e : composition on the same identities and colours
# -----------------------------------------------------------------------------
cat("Panel e (composition) ...\n")
prop <- read.csv(file.path(SRC, "Fig1E_proportions.csv"), stringsAsFactors = FALSE)
ord <- read.csv(file.path(SRC, "lineage_order.csv"), stringsAsFactors = FALSE)
ct_order <- as.character(ord$lineage)
prop$lineage <- factor(as.character(prop$lineage), levels = rev(ct_order))
prop$IV_lab <- factor(as.character(prop$IV_lab), levels = unname(IV_LAB))

pE <- ggplot(prop, aes(Percentage, lineage, fill = lineage)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", Percentage)), hjust = -0.12,
            size = 1.85, family = FAM, colour = "black") +
  facet_wrap(~ IV_lab, nrow = 1) +
  scale_fill_manual(values = CT_COL, guide = "none") +
  scale_x_continuous(name = "Cells (%)", limits = c(0, 50),
                     breaks = c(0, 20, 40),
                     expand = expansion(mult = c(0, 0.28))) +
  labs(y = NULL) +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = BASE - 1.0),
    strip.text = element_text(size = BASE - 1.0),
    plot.margin = margin(1, 2, 1, 1)
  )

# -----------------------------------------------------------------------------
# Assemble (no marker row; that panel is now Fig. S1a)
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
tag_it <- function(p, t) {
  p + labs(tag = t) +
    theme(plot.tag = element_text(size = 8, face = "bold", family = FAM,
                                  hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}
seal <- function(x) wrap_elements(full = x)

row1 <- seal(tag_it(pA, "a") + tag_it(pB, "b") +
               plot_layout(widths = c(0.962, 1.111)))
row2 <- seal(tag_it(pC, "c"))
row3 <- seal(tag_it(pD, "d") + tag_it(pE, "e") +
               plot_layout(widths = c(1.20, 0.90)))

fig <- row1 / row2 / row3 +
  plot_layout(heights = grid::unit.c(grid::unit(ROW1_MM, "mm"),
                                     grid::unit(ROW2_MM, "mm"),
                                     grid::unit(1, "null"))) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "Figure1")

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

SD <- file.path(OUTDIR, "source_data_plotted")
dir.create(SD, showWarnings = FALSE)
write.csv(df_d[, c("barcode", "umap_1", "umap_2", "lineage", "IV_lab")],
          file.path(SD, "SourceData_Fig1d_umap.csv"), row.names = FALSE)
write.csv(prop[, c("IV_lab", "lineage", "Freq", "total", "Percentage")],
          file.path(SD, "SourceData_Fig1e_composition.csv"), row.names = FALSE)

dpi_report <- function(path, printed_h_mm) {
  d <- dim(png::readPNG(path))
  data.frame(image = basename(path), px_w = d[2], px_h = d[1],
             printed_h_mm = round(printed_h_mm, 1),
             printed_w_mm = round(printed_h_mm * d[2] / d[1], 1),
             dpi = round(d[1] / printed_h_mm * 25.4))
}
row1_img_h <- ROW1_MM / 1.075
row2_img_h <- ROW2_MM / 1.09
dpi_tab <- rbind(
  dpi_report(file.path(ASSETS, "Fig1A_IF_B220_CD4_DCIR.png"), row1_img_h),
  dpi_report(file.path(ASSETS, "Fig1B_IF_F480_Laminin.png"),  row1_img_h),
  dpi_report(file.path(ASSETS, "Fig1C_flow_gating.png"),      row2_img_h),
  dpi_report(file.path(ASSETS, "Fig1C_flow_gating2.png"),     row2_img_h * 0.60)
)
cat("\nEffective print resolution of embedded images:\n")
print(dpi_tab, row.names = FALSE)
write.csv(dpi_tab, file.path(QA, "embedded_image_dpi.csv"), row.names = FALSE)

cat("\nCells plotted in panel d:", n_total, "\n")
cat("DONE ->", OUTDIR, "\n")
