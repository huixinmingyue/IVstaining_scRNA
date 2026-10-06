#!/usr/bin/env Rscript
# =============================================================================
# Fig. S6 — Nature SI native ggplot vector (R / ggplot2 / patchwork)
# Archetype: asymmetric mixed-modality figure
# Core conclusion: Microglia traverse Monocle states with vascular-contact
# enrichment; pseudotime-ordered programs culminate in inflammatory/Cxcl10
# modules that map to perivascular EAE lesions on ST.
# =============================================================================
options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(ggnewscale)
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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS6_Nature_vector")
SRC <- file.path(OUTDIR, "source_data")
FIG_W_MM <- 183
FIG_H_MM <- 255
FAM <- "Arial"
BASE <- 6.5
PT <- function(x) x / ggplot2::.pt

has_scattermore <- requireNamespace("scattermore", quietly = TRUE)

traj <- readRDS(file.path(SRC, "traj_slim.rds"))
hm   <- readRDS(file.path(SRC, "heatmap_slim.rds"))
sp_d <- readRDS(file.path(SRC, "spatial_d_slim.rds"))
sp_e <- readRDS(file.path(SRC, "spatial_e_slim.rds"))

state_cols <- c("1" = "#E64B35", "2" = "#00A087", "3" = "#3C5488")
group_cols <- c("CNS_IV" = "#E64B35", "CNS_IVC" = "#00A087")
grp_cols <- c(homeostatic = "#4C72B0", inflammatory = "#C44E52", proliferative = "#DD8452")
pal_spatial <- c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF8C00", "#FF0000")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = BASE),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.3),
    legend.text = element_text(size = BASE - 0.7),
    legend.key.size = unit(3.2, "mm"),
    plot.title = element_text(size = BASE + 0.5, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 9, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid = element_blank()
  )

# ---- Panel a / b: Monocle trajectory (native ggplot) ----
traj_plot <- function(color_col, pal, title, letter, legend_title) {
  df <- traj$cells
  df$.col <- factor(as.character(df[[color_col]]), levels = names(pal))
  p <- ggplot(df, aes(comp1, comp2, colour = .col)) +
    geom_segment(
      data = traj$edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, linewidth = 0.35, colour = "grey20"
    ) +
    geom_point(size = 0.70, stroke = 0, alpha = 0.88) +
    geom_point(
      data = traj$branch_points,
      aes(comp1, comp2), inherit.aes = FALSE,
      shape = 21, size = 2.4, fill = "black", colour = "white", stroke = 0.4
    ) +
    geom_text(
      data = traj$branch_points,
      aes(comp1, comp2, label = label), inherit.aes = FALSE,
      colour = "white", size = PT(5.5), fontface = "bold", family = FAM
    ) +
    scale_colour_manual(values = pal, name = legend_title, drop = FALSE) +
    labs(x = "Component 1", y = "Component 2", title = title, tag = letter) +
    coord_equal() +
    theme_pub +
    theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))
  p
}

pA <- traj_plot("State", state_cols, NULL, "a", "State")
pB <- traj_plot("group", group_cols, NULL, "b", "group")

# ---- Panel c: pseudotime heatmap (single ggplot; vector geom_tile) ----
mat <- hm$mat
gmeta <- hm$gene_meta
gene_levels <- rev(as.character(gmeta$gene))
gmeta$gene <- factor(as.character(gmeta$gene), levels = gene_levels)
df_hm <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
colnames(df_hm) <- c("gene", "bin", "z")
df_hm$gene <- factor(as.character(df_hm$gene), levels = gene_levels)
df_hm$x <- match(df_hm$bin, colnames(mat))
side <- data.frame(
  gene = factor(as.character(gmeta$gene), levels = gene_levels),
  group = factor(gmeta$group, levels = names(grp_cols)),
  x = -0.35
)
vmax <- max(0.5, as.numeric(quantile(abs(df_hm$z), 0.98, names = FALSE)))

pC <- ggplot() +
  geom_tile(data = side, aes(x = x, y = gene, fill = group), width = 0.55, height = 1, colour = NA) +
  scale_fill_manual(values = grp_cols, name = "functional group", guide = guide_legend(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = df_hm, aes(x = x, y = gene, fill = z), width = 1, height = 1, colour = NA) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-vmax, vmax), oob = squish,
    name = "z-score (bin mean)", guide = guide_colorbar(order = 2)
  ) +
  scale_x_continuous(
    breaks = seq_len(ncol(mat)),
    labels = colnames(mat),
    expand = expansion(mult = c(0.02, 0.02)),
    limits = c(-0.75, ncol(mat) + 0.5)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Pseudotime bins",
    y = NULL,
    title = "Pseudotime-ordered gene programs (all microglia)",
    subtitle = "early homeostatic \u2192 late inflammatory/proliferative",
    tag = "c"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 5.2, face = "italic", family = FAM),
    axis.text.x = element_text(size = 5.2),
    axis.title.x = element_text(size = 5.8),
    plot.subtitle = element_text(size = 5.2, hjust = 0.5, family = FAM, margin = margin(b = 2)),
    legend.key.height = unit(7, "mm"),
    legend.key.width = unit(2.4, "mm"),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3),
    plot.margin = margin(2, 4, 6, 2)
  )

# ---- H&E helpers (raster background only; spots are vector) ----
fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}
he_ann <- function(he, alpha = 0.90) {
  arr <- fade_he(he$image, alpha = alpha)
  nr <- dim(arr)[1]; nc <- dim(arr)[2]
  arr <- arr[rev(seq_len(nr)), , , drop = FALSE]
  list(
    raster = grDevices::as.raster(arr),
    xmin = 0.5, xmax = nc + 0.5,
    ymin = 0.5, ymax = nr + 0.5
  )
}

spot_spatial <- function(coords, values, cells = NULL, limits, title = NULL,
                         letter = NULL, pt = 0.35, he = NULL, he_alpha = 0.90,
                         show_legend = TRUE, legend_name = NULL,
                         crop = NULL, use_scatter = FALSE) {
  if (is.null(cells)) cells <- rownames(coords)
  cells <- intersect(cells, names(values))
  cells <- intersect(cells, rownames(coords))
  df <- data.frame(
    x = coords[cells, "x"],
    y = coords[cells, "y"],
    val = pmin(pmax(as.numeric(values[cells]), limits[1]), limits[2])
  )
  df <- df[order(df$val), , drop = FALSE]
  # Coord convention: x=imagecol, y=+imagerow after flip for display
  df$y <- -df$y
  p <- ggplot(df, aes(x, y, colour = val))
  if (!is.null(he)) {
    ha <- he_ann(he, alpha = he_alpha)
    # Flip HE to match y = -imagerow
    # he_ann already flipped once for +imagerow; for -imagerow flip again = original order
    arr <- fade_he(he$image, he_alpha)
    # y plotted as -imagerow: raster row1 at ymax should be small imagerow
    # array row 1 = top of image = small imagerow typically → place at ymax (-small)
    # With y = -imagerow: ymin = -nr-0.5, ymax = -0.5
    rast <- grDevices::as.raster(arr)
    ymin <- -(dim(arr)[1] + 0.5)
    ymax <- -0.5
    xmin <- 0.5; xmax <- dim(arr)[2] + 0.5
    p <- p + annotation_raster(rast, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    xlim <- c(xmin, xmax); ylim <- c(ymin, ymax)
  } else {
    xlim <- range(df$x); ylim <- range(df$y)
  }
  if (!is.null(crop)) {
    # crop in original imagerow/imagecol before y flip
    xlim <- c(crop$xmin, crop$xmax)
    ylim <- c(-crop$ymax, -crop$ymin)
  }
  if (use_scatter && has_scattermore && nrow(df) > 20000) {
    p <- p + scattermore::geom_scattermore(pointsize = 1.1, pixels = c(900, 900), interpolate = TRUE)
  } else {
    p <- p + geom_point(size = pt, stroke = 0)
  }
  p <- p +
    scale_colour_gradientn(colours = pal_spatial, limits = limits, name = legend_name, oob = squish) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, tag = letter) +
    theme_void(base_family = FAM) +
    theme(
      plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
      plot.tag = element_text(size = 9, face = "bold", family = FAM),
      plot.tag.position = c(0.01, 0.99),
      legend.title = element_text(size = BASE - 0.5),
      legend.text = element_text(size = BASE - 0.8),
      legend.key.height = unit(6, "mm"),
      legend.key.width = unit(2.2, "mm"),
      plot.margin = margin(1, 1, 1, 1)
    )
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

# Lesion cell sets from barcode row/col
bc <- sp_d$barcode_rc
lesion_cells <- function(reg) {
  hit <- rownames(bc)[bc$row >= reg$row[1] & bc$row <= reg$row[2] &
                        bc$col >= reg$col[1] & bc$col <= reg$col[2]]
  intersect(hit, rownames(sp_d$coords))
}

# Crop window in HE pixel space from selected spots
crop_from_cells <- function(cells, pad = 8) {
  xy <- sp_d$coords[cells, , drop = FALSE]
  list(
    xmin = max(0.5, min(xy$x) - pad),
    xmax = min(sp_d$he$dim[2] + 0.5, max(xy$x) + pad),
    ymin = max(0.5, min(xy$y) - pad),
    ymax = min(sp_d$he$dim[1] + 0.5, max(xy$y) + pad)
  )
}

make_lesion_row <- function(reg_name, letter = NULL, show_legend = FALSE) {
  cells <- lesion_cells(sp_d$regions[[reg_name]])
  crop <- crop_from_cells(cells)
  cat(reg_name, "spots:", length(cells), "\n")
  p_pos <- spot_spatial(
    sp_d$coords, sp_d$MgIV, cells = cells, limits = c(0, 0.3),
    title = "MgIV", letter = letter, pt = 0.9, he = sp_d$he,
    show_legend = FALSE, legend_name = "MgIV", crop = crop
  )
  p_neg <- spot_spatial(
    sp_d$coords, sp_d$MgIVC, cells = cells, limits = c(0, 0.3),
    title = "MgIVC", letter = NULL, pt = 0.9, he = sp_d$he,
    show_legend = show_legend, legend_name = "MgIVC", crop = crop
  )
  p_pos | p_neg
}

cat("Building panel d lesions ...\n")
pD1 <- make_lesion_row("Mg1", letter = "d", show_legend = FALSE)
pD2 <- make_lesion_row("Mg2", letter = NULL, show_legend = FALSE)
pD3 <- make_lesion_row("Mg3", letter = NULL, show_legend = TRUE)
pD <- (pD1 / pD2 / pD3) + plot_layout(heights = c(1, 1, 1))

cat("Building panel e: Cxcl10 | nCount_Spatial ...\n")
stopifnot(!is.null(sp_e$nCount_Spatial))
# Author QC spatial used viridis-like UMI map; keep shared spatial palette, limits ~0–6500
pE_cx <- spot_spatial(
  sp_e$coords, sp_e$Cxcl10, limits = c(0, 2.5),
  title = "Cxcl10", letter = "e", pt = 0.12, he = sp_e$he,
  show_legend = TRUE, legend_name = "Cxcl10",
  use_scatter = TRUE
)
pE_nc <- spot_spatial(
  sp_e$coords, sp_e$nCount_Spatial, limits = c(0, 6500),
  title = "nCount_Spatial", letter = NULL, pt = 0.12, he = sp_e$he,
  show_legend = TRUE, legend_name = "nCount",
  use_scatter = TRUE
)
pE <- (pE_cx | pE_nc) + plot_layout(widths = c(1, 1))

# ---- Layout ----
left_top <- (pA | pB) + plot_layout(widths = c(1, 1))
left <- (left_top / pD) + plot_layout(heights = c(0.95, 2.6))
top <- (left | pC) + plot_layout(widths = c(1.45, 1.05))
# Panel e is now two full-section maps side-by-side — give it a bit more height
fig <- (top / pE) + plot_layout(heights = c(2.95, 1.25))
fig <- fig + plot_annotation(
  theme = theme(plot.margin = margin(3, 3, 3, 3), plot.background = element_rect(fill = "white", colour = NA))
)

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS6")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig)
dev.off()

cat("Writing PDF ...\n")
pdf_path <- paste0(base, ".pdf")
pdf_alt <- paste0(base, "_native.pdf")
write_pdf <- function(path) {
  grDevices::cairo_pdf(path, width = w_in, height = h_in, family = FAM, fallback_resolution = 600)
  print(fig)
  grDevices::dev.off()
}
ok <- tryCatch({
  if (file.exists(pdf_path) && !isTRUE(file.remove(pdf_path))) stop("locked")
  write_pdf(pdf_path); TRUE
}, error = function(e) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  write_pdf(pdf_alt); FALSE
})
cat("PDF:", if (ok) pdf_path else pdf_alt, "\n")

cat("Writing PNG 600 dpi ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in", res = 600, background = "white")
print(fig)
dev.off()

writeLines(c(
  "Fig. S6 | Region-specific microglia subtypes and pseudotime differentiation.",
  "a,b, Native ggplot Monocle trajectories of EAE microglia coloured by State (a) or vascular-contact group (CNS_IV / CNS_IVC; b).",
  "c, Pseudotime-ordered gene-program heatmap (ggplot geom_tile) across all microglia; genes ranked by peak along pseudotime; sidebar marks homeostatic, inflammatory and proliferative programs; colour, z-scored bin-mean expression.",
  "d, Spatial ggplot maps of MgIV / MgIVC prediction scores in three perivascular lesions (rows), overlaid on Visium H&E.",
  "e, Spatial ggplot maps on the EAE brain section (Visium H&E): Cxcl10 expression (left) and nCount_Spatial UMI density (right; from the QC SpatialFeaturePlot companion to the nCount violin). Dense spots rasterized for performance.",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS6_legend.txt"), useBytes = TRUE)

writeLines(c(
  "Figure S6 | Nature SI native ggplot vector",
  "Backend: R exclusive (ggplot2/patchwork)",
  "Reuse: build anew from monocle slim + s6c export + Visium slim",
  sprintf("Size: %d x %d mm", FIG_W_MM, FIG_H_MM),
  "a/b: geom_point + geom_segment trajectory (vector)",
  "c: geom_tile heatmap (vector)",
  "d/e: geom_point/scattermore on H&E annotation_raster",
  "Panel e: Cxcl10 | nCount_Spatial side-by-side",
  "No author PDF/TIFF plate embeds; prior FigS6_Nature_submission untouched."
), file.path(OUTDIR, "FigureS6_QA.txt"), useBytes = TRUE)

writeLines(c(
  "Figure S6 | Nature vector submission",
  "Rendered in R (native ggplot2).",
  "Deliverables: FigureS6.svg, FigureS6.pdf, FigureS6.png",
  paste0("Output: ", OUTDIR)
), file.path(OUTDIR, "FigureS6_provenance.txt"), useBytes = TRUE)

cat("DONE ->", OUTDIR, "\n")