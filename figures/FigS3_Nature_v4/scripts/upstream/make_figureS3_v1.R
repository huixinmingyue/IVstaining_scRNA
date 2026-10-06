#!/usr/bin/env Rscript
# =============================================================================
# Fig. S3 v1 - Nature-style lineage-marker support for Fig. 3
#
# Core conclusion
#   Canonical lineage genes recover the hard RNA gate used in Fig. 3:
#   Cd4, Cd8b1 and Klrb1c peak in true CD4, true CD8 and NKT; Ms4a1 peaks
#   in the cluster-2 doublet.
#
# Evidence chain
#   a-d  Cd4 / Cd8b1 / Klrb1c / Ms4a1 (former Fig. 3 v2 panels d-g)
#
# Archetype : quantitative grid
# Backend   : R (ggplot2 / patchwork / svglite / ragg)
# Size      : 183 x 130 mm
# New folder
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
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

set.seed(42)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS3_Nature_v1")
}
SRC <- file.path(OUTDIR, "source_data")
QA  <- file.path(OUTDIR, "qa")
SDP <- file.path(OUTDIR, "source_data_plotted")
dir.create(QA, showWarnings = FALSE)
dir.create(SDP, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 130
width_mm <- FIG_W_MM
height_mm <- FIG_H_MM
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

CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
CLASS_COL <- c(
  true_CD4 = "#2F6F8F", true_CD8 = "#5B8F4E", NKT = "#B2473E",
  DP = "#C4A35A", unassigned = "#9A9A9A", doublet_c2 = "#5A5A5A"
)
CLASS_SHORT <- c(true_CD4 = "CD4", true_CD8 = "CD8", NKT = "NKT",
                 DP = "DP", unassigned = "Unass.", doublet_c2 = "Dbl.")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.0, "pt"),
      axis.title = element_text(size = BASE, colour = "black"),
      axis.text = element_text(size = BASE - 0.5, colour = "black"),
      legend.title = element_text(size = BASE - 0.2),
      legend.text = element_text(size = BASE - 0.6),
      legend.key.size = unit(2.4, "mm"),
      legend.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE - 0.7, hjust = 0.5, colour = "#555555"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      panel.grid = element_blank()
    )
)

fmt_p_sci1 <- function(p) {
  if (is.na(p)) return("NA")
  if (p >= 0.001) return(sprintf("%.3g", p))
  sprintf("%.1e", p)
}
tag_it <- function(p, t) {
  p + labs(tag = t) +
    theme(plot.tag = element_text(size = 8, face = "bold", family = FAM,
                                  hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}

cat("Load cells ...\n")
cells <- readRDS(file.path(SRC, "Fig3_T_cells.rds"))
stopifnot(nrow(cells) == 6098)
cells$t_class <- factor(cells$t_class, levels = CLASS_LV)
cells$t_short <- factor(unname(CLASS_SHORT[as.character(cells$t_class)]),
                        levels = unname(CLASS_SHORT))

feat_vln <- function(df, gene, pal) {
  zmax <- max(as.numeric(quantile(df$expr, 0.99, names = FALSE)), 0.15)
  peak <- names(which.max(tapply(df$expr, df$t_short, mean)))
  wt <- suppressWarnings(wilcox.test(df$expr[as.character(df$t_short) == peak],
                                     df$expr[as.character(df$t_short) != peak]))
  p_f <- ggplot(df, aes(UMAP1, UMAP2, colour = expr)) +
    geom_point(size = 0.12, alpha = 0.78, shape = 16) +
    scale_colour_gradientn(
      colours = c(SEQ_LO, SEQ_HI), limits = c(0, zmax), oob = squish,
      name = NULL, breaks = c(0, zmax),
      labels = c("0", sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = unit(1.8, "mm"),
                             barheight = unit(11, "mm"),
                             ticks.colour = "black", frame.colour = "black",
                             frame.linewidth = 0.2)
    ) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = gene) +
    theme(plot.title = element_text(face = "italic", size = BASE + 0.2, hjust = 0.5),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
          legend.position = "right",
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0.6, 0.4, 0.4, 0.4, "mm"))
  p_v <- ggplot(df, aes(t_short, expr, fill = t_class)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Expression",
         subtitle = sprintf("highest in %s", peak)) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5,
                                       colour = "#555555"),
          axis.text.x = element_text(size = 5.6, face = "bold", angle = 45,
                                     hjust = 1, vjust = 1),
          plot.margin = margin(0.6, 0.6, 0.4, 0.6, "mm"))
  list(plot = wrap_elements(full = p_f + p_v + plot_layout(widths = c(1.15, 0.85))),
       peak = peak, peak_p = wt$p.value)
}

cat("Panels a-d ...\n")
mk <- list(a = "Cd4", b = "Cd8b1", c = "Klrb1c", d = "Ms4a1")
built <- list()
for (nm in names(mk)) {
  d <- cells
  d$expr <- cells[[mk[[nm]]]]
  built[[nm]] <- feat_vln(d, mk[[nm]], CLASS_COL)
  cat(mk[[nm]], "peak", built[[nm]]$peak,
      "P", fmt_p_sci1(built[[nm]]$peak_p), "\n")
}

fig <- (tag_it(built$a$plot, "a") | tag_it(built$b$plot, "b")) /
  (tag_it(built$c$plot, "c") | tag_it(built$d$plot, "d")) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(0.6, 0.6, 0.6, 0.6, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS3")

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

write.csv(cells[, c("cell", "UMAP1", "UMAP2", "t_class", "IV_group",
                    "Cd4", "Cd8b1", "Klrb1c", "Ms4a1")],
          file.path(SDP, "SourceData_FigS3_markers.csv"), row.names = FALSE)

cat("DONE ->", OUTDIR, "\n")
