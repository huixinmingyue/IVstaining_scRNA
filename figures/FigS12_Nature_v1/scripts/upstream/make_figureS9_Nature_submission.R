#!/usr/bin/env Rscript
# =============================================================================
# Fig. S9 — Nature Communications SI submission (R / ggplot2 / patchwork)
# Archetype: image plate + quant
# Core conclusion: Independent mouse Visium sections recapitulate n-/p-
# microglia subtype maps and POS-spot enrichment of Lgals3/Prdx1/Mki67.
# Image integrity: title-band crop only; no contrast/local edits.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
  library(svglite)
  library(ragg)
  library(systemfonts)
})

set.seed(42)
options(warn = 1)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS9_Nature_submission")
ASSET  <- file.path(OUTDIR, "assets")
SRC    <- file.path(OUTDIR, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 180
FIG_H_MM <- 280
BASE <- 6.5
FAM  <- "Arial"
fm <- systemfonts::match_fonts(FAM)
cat("font resolve:", FAM, "->", fm$path, "\n")
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM
    do.call(reg, args)
  }
}
PT <- function(x) x / .pt

theme_set(
  theme_void(base_size = BASE, base_family = FAM) +
    theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM))
)

# -----------------------------------------------------------------------------
# Raster helpers
# -----------------------------------------------------------------------------
load_rast <- function(path) {
  stopifnot(file.exists(path))
  png::readPNG(path)
}

raster_panel <- function(path, title = NULL, title_face = "bold",
                         title_size = BASE) {
  rast <- load_rast(path)
  g <- rasterGrob(rast, interpolate = TRUE)
  ar <- nrow(rast) / ncol(rast)
  p <- ggplot() +
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void(base_size = BASE, base_family = FAM) +
    theme(
      plot.margin = margin(1, 1, 1, 1),
      aspect.ratio = ar
    )
  if (!is.null(title) && nzchar(title)) {
    p <- p +
      labs(title = title) +
      theme(plot.title = element_text(
        size = title_size, face = title_face, hjust = 0.5, family = FAM,
        margin = margin(b = 1, t = 0)
      ))
  }
  p
}

hdr_plot <- function(txt, face = "bold", size = BASE) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = txt,
      size = PT(size), fontface = face, family = FAM
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

ap <- function(...) file.path(ASSET, ...)

# -----------------------------------------------------------------------------
# Panel a — six n-/p- spatial maps (two sections × clusters 5/3/0)
# -----------------------------------------------------------------------------
# Row 1: mg_52 / mg_32 / mg_02 ; Row 2: mg_45 / mg_43 / mg_40
a_top <- wrap_plots(
  raster_panel(ap("brain_mouse_Spatial_mg_52.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_32.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_02.png")),
  nrow = 1
)
a_bot <- wrap_plots(
  raster_panel(ap("brain_mouse_Spatial_mg_45.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_43.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_40.png")),
  nrow = 1
)
a_col_hdr <- wrap_plots(
  hdr_plot("Cluster 5"),
  hdr_plot("Cluster 3"),
  hdr_plot("Cluster 0"),
  nrow = 1
)
pA <- (a_col_hdr / a_top / a_bot) + plot_layout(heights = c(0.08, 1, 1))

# -----------------------------------------------------------------------------
# Panel b — gene spatial stacks + boxplots (section set 1)
# -----------------------------------------------------------------------------
make_gene_quant_col <- function(gene, gene_png, box_png, cluster_lab) {
  (hdr_plot(gene, face = "bold.italic", size = BASE + 0.5) /
     raster_panel(gene_png, title = paste0("POS / NEG · ", cluster_lab),
                  title_face = "plain", title_size = BASE - 0.8) /
     raster_panel(box_png)) +
    plot_layout(heights = c(0.07, 1.15, 0.95))
}

pB <- wrap_plots(
  make_gene_quant_col(
    "Lgals3",
    ap("mouse_5_Lgals3_clean.png"),
    ap("n5_p5_Lgals3_mouse_clean.png"),
    "cluster 5"
  ),
  make_gene_quant_col(
    "Prdx1",
    ap("mouse_5_Prdx14new_clean.png"),
    ap("n5_p5_Prdx1_mouse_clean.png"),
    "cluster 5"
  ),
  make_gene_quant_col(
    "Mki67",
    ap("1_3_Mki67n_clean.png"),
    ap("n3_p3_Mki67_mouse2new_clean.png"),
    "cluster 3"
  ),
  nrow = 1
)

# -----------------------------------------------------------------------------
# Panel c — second-section spatial maps (clusters 5/3/0)
# -----------------------------------------------------------------------------
c_maps <- wrap_plots(
  raster_panel(ap("brain_mouse_Spatial_mg_45.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_43.png")),
  raster_panel(ap("brain_mouse_Spatial_mg_40.png")),
  nrow = 1
)
c_hdr <- wrap_plots(
  hdr_plot("Cluster 5"),
  hdr_plot("Cluster 3"),
  hdr_plot("Cluster 0"),
  nrow = 1
)
pC <- (c_hdr / c_maps) + plot_layout(heights = c(0.1, 1))

# -----------------------------------------------------------------------------
# Panel d — gene spatial stacks + boxplots (section set 2)
# -----------------------------------------------------------------------------
pD <- wrap_plots(
  make_gene_quant_col(
    "Lgals3",
    ap("mouse_5_Lgals34_clean.png"),
    ap("n5_p5_Lgals3_mouse2_clean.png"),
    "cluster 5"
  ),
  make_gene_quant_col(
    "Prdx1",
    ap("mouse_5_Prdx14_clean.png"),
    ap("n5_p5_Prdx1_mouse2_clean.png"),
    "cluster 5"
  ),
  make_gene_quant_col(
    "Mki67",
    ap("MCAO_43_Mki67_clean.png"),
    ap("n3_p3_Mki67_mouse2_clean.png"),
    "cluster 3"
  ),
  nrow = 1
)

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
fig <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  wrap_elements(full = pC) /
  wrap_elements(full = pD) +
  plot_layout(heights = c(1.15, 1.55, 0.65, 1.55)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.004, 0.995)
  )

# -----------------------------------------------------------------------------
# Provenance / legend / image map
# -----------------------------------------------------------------------------
img_map <- data.frame(
  panel = c(
    rep("a", 6), rep("b", 6), rep("c", 3), rep("d", 6)
  ),
  role = c(
    "n5_p5_section1", "n3_p3_section1", "n0_p0_section1",
    "n5_p5_section2", "n3_p3_section2", "n0_p0_section2",
    "Lgals3_spatial", "Prdx1_spatial", "Mki67_spatial",
    "Lgals3_boxplot", "Prdx1_boxplot", "Mki67_boxplot",
    "n5_p5_section2", "n3_p3_section2", "n0_p0_section2",
    "Lgals3_spatial", "Prdx1_spatial", "Mki67_spatial",
    "Lgals3_boxplot", "Prdx1_boxplot", "Mki67_boxplot"
  ),
  source_tiff = c(
    "brain_mouse_Spatial_mg_52.tiff", "brain_mouse_Spatial_mg_32.tiff",
    "brain_mouse_Spatial_mg_02.tiff", "brain_mouse_Spatial_mg_45.tiff",
    "brain_mouse_Spatial_mg_43.tiff", "brain_mouse_Spatial_mg_40.tiff",
    "mouse_5_Lgals3.tiff", "mouse_5_Prdx14new.tiff", "1_3_Mki67n.tiff",
    "n5_p5_Lgals3_mouse.tiff", "n5_p5_Prdx1_mouse.tiff",
    "n3_p3_Mki67_mouse2new.tiff",
    "brain_mouse_Spatial_mg_45.tiff", "brain_mouse_Spatial_mg_43.tiff",
    "brain_mouse_Spatial_mg_40.tiff",
    "mouse_5_Lgals34.tiff", "mouse_5_Prdx14.tiff", "MCAO_43_Mki67.tiff",
    "n5_p5_Lgals3_mouse2.tiff", "n5_p5_Prdx1_mouse2.tiff",
    "n3_p3_Mki67_mouse2.tiff"
  ),
  processing = c(
    rep("resize_only", 6),
    rep("title_band_crop+resize", 6),
    rep("resize_only", 3),
    rep("title_band_crop+resize", 6)
  ),
  stringsAsFactors = FALSE
)
write.csv(img_map, file.path(SRC, "FigS9_image_map.csv"), row.names = FALSE)

writeLines(c(
  "Fig. S9 | Spatial microglia subtype maps and effector-gene enrichment across mouse Visium sections.",
  "a, Spatial n- and p- prediction scores for clusters 5, 3 and 0 on two independent mouse Visium sections (top row: mg_52/32/02; bottom row: mg_45/43/40).",
  "b, Spatial expression of Lgals3 and Prdx1 (cluster-5 POS vs NEG spots) and Mki67 (cluster-3 POS vs NEG spots), with matching boxplots (Wilcoxon; ****P < 1×10⁻⁴ as displayed on source plates).",
  "c, Spatial n-/p- maps for clusters 5, 3 and 0 on the second section (mg_45/43/40).",
  "d, Independent spatial gene maps and boxplots for Lgals3, Prdx1 and Mki67 (second plate set).",
  "Image integrity: global resize; Chinese/title bands cropped from gene and boxplot plates; no local contrast edits.",
  "Source data / image map are provided as a Source Data file."
), file.path(OUTDIR, "FigureS9_legend.txt"))

writeLines(c(
  "Figure S9 | Nature Communications SI image-plate assembly",
  "Backend: R (ggplot2/patchwork/png/grid)",
  "Primary deliverables: FigureS9.svg, FigureS9.png, FigureS9.pdf",
  paste0("Size: ", FIG_W_MM, " x ", FIG_H_MM, " mm; font: ", FAM, " ", BASE, " pt"),
  paste0("Output: ", OUTDIR)
), file.path(OUTDIR, "FigureS9_stats.txt"))

writeLines(c(
  "Core conclusion: Independent mouse Visium sections recapitulate n-/p- microglia subtype spatial maps and POS-spot enrichment of Lgals3/Prdx1/Mki67.",
  "Archetype: image plate + quant",
  "Backend: R exclusive",
  "Final size: 180 x 280 mm",
  "Panel a: 2x3 spatial n/p maps; Cluster 5/3/0 headers",
  "Panel b: italic gene titles; POS/NEG spatial + boxplots",
  "Panel c: second-section Cluster 5/3/0 maps",
  "Panel d: second gene/boxplot plate set",
  "Integrity: title-band crop only on gene/boxplot plates; spatial maps unmodified except resize"
), file.path(OUTDIR, "FigureS9_QA.txt"))

writeLines(c(
  "Figure S9 provenance",
  "Source directory: author MCAO-spatial TIFF plates",
  "Assets converted to PNG for R embedding; cleaned gene/box plates stored as *_clean.png",
  "Does not overwrite prior FigS* folders."
), file.path(OUTDIR, "FigureS9_provenance.txt"))

# -----------------------------------------------------------------------------
# Export SVG + PNG
# -----------------------------------------------------------------------------
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS9")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig)
dev.off()

cat("Writing PNG (600 dpi) ...\n")
ragg::agg_png(
  paste0(base, ".png"), width = w_in, height = h_in,
  units = "in", res = 600, background = "white"
)
print(fig)
dev.off()

cat("Writing PDF ...\n")
pdf_path <- paste0(base, ".pdf")
pdf_alt  <- paste0(base, "_native.pdf")
write_pdf <- function(path) {
  grDevices::cairo_pdf(path, width = w_in, height = h_in, family = FAM, fallback_resolution = 600)
  print(fig)
  grDevices::dev.off()
}
ok <- tryCatch({
  if (file.exists(pdf_path) && !isTRUE(file.remove(pdf_path))) stop("locked")
  write_pdf(pdf_path)
  TRUE
}, error = function(e) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  write_pdf(pdf_alt)
  FALSE
})
cat("PDF:", if (ok) pdf_path else pdf_alt, "\n")

cat("DONE Nature submission Figure S9 ->", OUTDIR, "\n")
