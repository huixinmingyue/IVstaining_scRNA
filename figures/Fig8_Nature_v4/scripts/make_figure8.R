#!/usr/bin/env Rscript
# =============================================================================
# Fig. 8 v4 — Nature main figure (R / ggplot2 exclusive)
#
# Core conclusion
#   Six perivascular windows mark vessel-associated niches; T-cell and
#   microglia POS versus NEG scores differ among all immune spots and
#   inside L1–L3 versus L4–L6; microglia POS/NEG geography is resolved
#   inside each window.
#
# Evidence chain
#   a    Vector schematic: EAE section, slice, Visium ST
#   b    Sagittal Laminin/DAPI
#   c    Coronal IF with L1–L6 labels
#   d    L1–L3 IF zooms plus L4–L6 crops
#   e    H&E and clusters (no L1–L6 zooms)
#   f    T-cell maps; all-spot and window boxplots
#   g    Microglia maps; all-spot and window boxplots
#   h    Microglia POS/NEG window crops L1–L6; same boxplots as f–g
#
# Archetype : asymmetric mixed-modality
# Size      : 183 x 470 mm
# New folder — does not overwrite Fig8_Nature_v1/v2/v3
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
  library(systemfonts)
  library(grid)
  library(png)
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "Fig8_Nature_v4")
}
SRC     <- file.path(OUTDIR, "source_data")
ASSET   <- file.path(OUTDIR, "assets")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 470
width_mm = 183
height_mm = 470
FAM      <- "Arial"
BASE     <- 6.5
PT5      <- 5.2
PT       <- function(x) x / ggplot2::.pt

pos_neg_cols <- c("EAE_POS" = "#9B59B6", "EAE_NEG" = "#E67E22")
cluster_cols <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3",
  "4" = "#FF7F00", "5" = "#A65628", "6" = "#F781BF", "7" = "#999999",
  "8" = "#66C2A5"
)
cluster_annot <- c(
  "0" = "Neurons", "1" = "Neurons", "2" = "Hippocampal neurons",
  "3" = "Oligodendrocytes", "4" = "Purkinje cells", "5" = "Granule cells",
  "6" = "Meninges", "7" = "Vascular-associated", "8" = "Plasma cells"
)
cluster_legend_lab <- setNames(
  paste0(names(cluster_annot), ": ", unname(cluster_annot)),
  names(cluster_annot)
)
score_pal <- c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725")

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

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.title = element_text(size = 5.8, colour = "black"),
    axis.text = element_text(size = 5.2, colour = "black"),
    legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.2),
    legend.key.size = unit(2.2, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid = element_blank()
  )
theme_set(theme_pub)
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = FAM))

star_lab <- function(p) {
  if (is.na(p)) return("")
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 0.01) "**" else
    if (p < 0.05) "*" else "ns"
}

in_rect <- function(x, y, r) {
  x >= r$xmin & x <= r$xmax & y >= r$ymin & y <= r$ymax
}

embed_png <- function(path, letter = NULL, ink = "black") {
  img <- png::readPNG(path)
  p <- embed_array(img, title = NULL)
  if (!is.null(letter)) {
    p <- p + labs(tag = letter) +
      theme(
        plot.tag = element_text(size = 8, face = "bold", family = FAM, colour = ink),
        plot.tag.position = c(0.045, 0.90),
        plot.margin = margin(12, 6, 2, 8)
      )
  }
  p
}

embed_array <- function(arr, title = NULL) {
  rast <- grDevices::as.raster(pmin(pmax(arr, 0), 1))
  nr <- dim(arr)[1]; nc <- dim(arr)[2]
  p <- ggplot() +
    annotation_raster(rast, xmin = 0.5, xmax = nc + 0.5,
                      ymin = -(nr + 0.5), ymax = -0.5) +
    coord_fixed(xlim = c(0.5, nc + 0.5), ylim = c(-(nr + 0.5), -0.5),
                expand = FALSE) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(1, 1, 1, 1))
  if (!is.null(title)) {
    p <- p + labs(title = title) +
      theme(plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5,
                                      family = FAM))
  }
  p
}

flow_arrow <- function(vertical = FALSE) {
  if (vertical) {
    p <- ggplot() +
      annotate(
        "segment", x = 0.5, xend = 0.5, y = 0.92, yend = 0.08,
        arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
        linewidth = 0.5, colour = "#2C5AA0"
      )
  } else {
    p <- ggplot() +
      annotate(
        "segment", x = 0.06, xend = 0.94, y = 0.5, yend = 0.5,
        arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
        linewidth = 0.5, colour = "#2C5AA0"
      )
  }
  p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

panel_sagittal <- function() {
  outline <- data.frame(
    x = c(0.06, 0.16, 0.32, 0.50, 0.66, 0.78, 0.88, 0.92, 0.86, 0.76,
          0.66, 0.52, 0.36, 0.20, 0.10, 0.06),
    y = c(0.46, 0.62, 0.72, 0.74, 0.66, 0.54, 0.44, 0.32, 0.24, 0.22,
          0.28, 0.32, 0.34, 0.36, 0.40, 0.46)
  )
  cereb <- data.frame(
    x = c(0.70, 0.80, 0.90, 0.91, 0.84, 0.74, 0.70),
    y = c(0.50, 0.60, 0.50, 0.36, 0.30, 0.36, 0.50)
  )
  ggplot() +
    geom_polygon(data = outline, aes(x, y),
                 fill = "#F3C4B3", colour = "#B56B5C", linewidth = 0.45) +
    geom_polygon(data = cereb, aes(x, y),
                 fill = "#E39A88", colour = "#B56B5C", linewidth = 0.35) +
    annotate(
      "segment", x = 0.20, xend = 0.74, y = 0.44, yend = 0.44,
      linetype = "22", linewidth = 0.55, colour = "#1F3A5F",
      arrow = arrow(length = unit(1.5, "mm"), type = "closed")
    ) +
    annotate("text", x = 0.48, y = 0.90, label = "EAE brain",
             size = PT(5.8), fontface = "bold", family = FAM) +
    annotate("text", x = 0.48, y = 0.10, label = "Coronal cut",
             size = PT(5.2), family = FAM, colour = "#1F3A5F") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(1, 1, 1, 1))
}

panel_coronal <- function() {
  th <- seq(0, 2 * pi, length.out = 90)
  outline <- data.frame(x = 0.50 + 0.34 * cos(th), y = 0.50 + 0.40 * sin(th))
  ggplot() +
    geom_polygon(data = outline, aes(x, y),
                 fill = "#F3C4B3", colour = "#B56B5C", linewidth = 0.45) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.12, yend = 0.88,
             colour = "#B56B5C", linewidth = 0.25) +
    annotate("rect", xmin = 0.30, xmax = 0.70, ymin = 0.16, ymax = 0.56,
             fill = NA, colour = "#1F3A5F", linewidth = 0.5, linetype = "22") +
    annotate("text", x = 0.50, y = 0.92, label = "Coronal slice",
             size = PT(5.8), fontface = "bold", family = FAM) +
    annotate("text", x = 0.50, y = 0.06, label = "Sampling ROI",
             size = PT(5.2), family = FAM, colour = "#1F3A5F") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(1, 1, 1, 1))
}

panel_slide <- function(he) {
  rast <- grDevices::as.raster(pmin(pmax(he, 0), 1))
  ggplot() +
    annotate("rect", xmin = 0.06, xmax = 0.94, ymin = 0.24, ymax = 0.76,
             fill = "#EEF3F8", colour = "#5D6D7E", linewidth = 0.5) +
    annotate("rect", xmin = 0.72, xmax = 0.94, ymin = 0.24, ymax = 0.76,
             fill = "#F4F4F4", colour = NA) +
    annotation_raster(rast, xmin = 0.14, xmax = 0.66, ymin = 0.32, ymax = 0.68) +
    annotate("text", x = 0.50, y = 0.90, label = "Visium slide",
             size = PT(5.8), fontface = "bold", family = FAM) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(1, 1, 1, 1))
}

spat_cluster <- function(df, he, nr, nc, title, pt = 0.25, show_lab = TRUE,
                         rects = NULL, show_legend = FALSE) {
  d <- df
  d$cluster <- factor(as.character(d$cluster),
                      levels = as.character(sort(as.integer(unique(as.character(d$cluster))))))
  labs_df <- d %>%
    group_by(cluster) %>%
    summarise(x_img = median(x_img), y_img = median(y_img), .groups = "drop")
  lev <- levels(d$cluster)
  cols <- cluster_cols[lev]
  names(cols) <- lev
  leg_lab <- cluster_legend_lab[lev]
  leg_lab[is.na(leg_lab)] <- lev[is.na(leg_lab)]
  p <- ggplot(d, aes(x_img, -y_img, colour = cluster)) +
    annotation_raster(he, xmin = 0, xmax = nc, ymin = -nr, ymax = 0) +
    scattermore::geom_scattermore(pointsize = max(1.0, pt * 4), pixels = c(900, 900))
  if (show_legend) {
    p <- p + scale_colour_manual(
      values = cols, labels = unname(leg_lab), name = "Cluster",
      guide = guide_legend(override.aes = list(size = 2.2, alpha = 1),
                           keywidth = unit(2.2, "mm"), keyheight = unit(2.2, "mm"))
    )
  } else {
    p <- p + scale_colour_manual(values = cols, guide = "none")
  }
  p <- p +
    coord_fixed(xlim = c(0, nc), ylim = c(-nr, 0), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = BASE, base_family = FAM) +
    theme(plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5,
                                    family = FAM))
  if (show_legend) {
    p <- p + theme(
      legend.position = "right",
      legend.title = element_text(size = 5.4, face = "bold", family = FAM),
      legend.text = element_text(size = 5.2, family = FAM)
    )
  }
  if (show_lab) {
    p <- p + geom_text(
      data = labs_df, aes(x = x_img, y = -y_img, label = cluster),
      colour = "black", size = PT(5.2), fontface = "bold", family = FAM,
      inherit.aes = FALSE
    )
  }
  if (!is.null(rects) && nrow(rects) > 0) {
    p <- p + geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
      inherit.aes = FALSE, fill = NA, colour = "yellow",
      linewidth = 0.45, linetype = "dashed"
    )
    p <- p + geom_text(
      data = rects,
      aes(x = (xmin + xmax) / 2, y = -ymin + 10, label = lesion),
      inherit.aes = FALSE, colour = "yellow", size = PT(5.2),
      fontface = "bold", family = FAM
    )
  }
  p
}

spat_score <- function(df, score_col, title, max_cut, he, nr, nc,
                       pt = 0.35, rects = NULL) {
  d <- df
  d$val <- pmin(d[[score_col]], max_cut)
  p <- ggplot(d, aes(x_img, -y_img, colour = val)) +
    annotation_raster(he, xmin = 0, xmax = nc, ymin = -nr, ymax = 0) +
    scattermore::geom_scattermore(pointsize = max(1.2, pt * 4), pixels = c(800, 800)) +
    scale_colour_gradientn(
      colours = score_pal, limits = c(0, max_cut), name = NULL, oob = squish
    ) +
    coord_fixed(xlim = c(0, nc), ylim = c(-nr, 0), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = BASE, base_family = FAM) +
    theme(
      plot.title = element_text(size = 5.4, face = "bold", hjust = 0.5, family = FAM),
      legend.key.height = unit(8, "mm"),
      legend.key.width = unit(1.8, "mm"),
      legend.text = element_text(size = 5.2)
    )
  if (!is.null(rects) && nrow(rects) > 0) {
    p <- p + geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
      inherit.aes = FALSE, fill = NA, colour = "yellow",
      linewidth = 0.4, linetype = "dashed"
    )
  }
  p
}

make_box <- function(df, title) {
  d <- df
  d$group <- factor(d$group, levels = c("EAE_POS", "EAE_NEG"))
  wt <- wilcox.test(score ~ group, data = d)
  ymax <- max(d$score, na.rm = TRUE)
  ggplot(d, aes(group, score, fill = group)) +
    geom_boxplot(width = 0.55, outlier.size = 0.22, outlier.alpha = 0.3,
                 linewidth = 0.3, alpha = 0.85) +
    geom_text(data = data.frame(x = 1.5, y = ymax * 1.08, lab = star_lab(wt$p.value)),
              aes(x, y, label = lab), inherit.aes = FALSE,
              size = PT(5.2), family = FAM) +
    scale_fill_manual(values = pos_neg_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    labs(x = NULL, y = "Prediction score", title = title) +
    theme(
      axis.text.x = element_text(size = 5.2),
      plot.title = element_text(size = 5.4, face = "bold", family = FAM,
                                margin = margin(b = 1)),
      plot.margin = margin(1, 2, 8, 2)
    )
}

make_box_split <- function(df, title) {
  d <- df
  d$group <- factor(d$group, levels = c("EAE_POS", "EAE_NEG"))
  d$set <- factor(d$set, levels = c("L1-L3", "L4-L6"))
  stats <- d %>%
    group_by(set) %>%
    summarise(
      p = wilcox.test(score ~ group)$p.value,
      y = max(score, na.rm = TRUE) * 1.08,
      .groups = "drop"
    ) %>%
    mutate(lab = vapply(p, star_lab, character(1)), x = as.numeric(set))
  ggplot(d, aes(set, score, fill = group)) +
    geom_boxplot(width = 0.7, outlier.size = 0.18, outlier.alpha = 0.28,
                 linewidth = 0.3, alpha = 0.85,
                 position = position_dodge(width = 0.75)) +
    geom_text(data = stats, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, size = PT(5.2), family = FAM) +
    scale_fill_manual(values = pos_neg_cols, name = NULL) +
    scale_x_discrete(labels = c("L1-L3" = "L1–L3", "L4-L6" = "L4–L6")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.24))) +
    labs(x = NULL, y = "Prediction score", title = title) +
    theme(
      legend.position = "top",
      legend.margin = margin(0, 0, 0, 0),
      legend.key.size = grid::unit(2.4, "mm"),
      legend.text = element_text(size = 5.2),
      axis.text.x = element_text(size = 5.2),
      plot.title = element_text(size = 5.4, face = "bold", family = FAM,
                                margin = margin(t = 2, b = 1)),
      plot.margin = margin(2, 2, 16, 2)
    )
}

crop_view <- function(df, he, pad = 12) {
  nr0 <- dim(he)[1]; nc0 <- dim(he)[2]
  xmin <- max(0, min(df$x_img) - pad); xmax <- min(nc0, max(df$x_img) + pad)
  ymin <- max(0, min(df$y_img) - pad); ymax <- min(nr0, max(df$y_img) + pad)
  r1 <- max(1, floor(ymin) + 1); r2 <- min(nr0, ceiling(ymax))
  c1 <- max(1, floor(xmin) + 1); c2 <- min(nc0, ceiling(xmax))
  he_c <- he[r1:r2, c1:c2, , drop = FALSE]
  d2 <- df
  d2$x_img <- d2$x_img - (c1 - 1)
  d2$y_img <- d2$y_img - (r1 - 1)
  list(df = d2, he = he_c, nr = dim(he_c)[1], nc = dim(he_c)[2])
}

mk_map_box <- function(map_l, map_r, box) {
  maps <- wrap_elements(panel = map_l | map_r)
  wrap_elements(panel = (maps | box) + plot_layout(widths = c(2.15, 0.85)))
}

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
cat("Load ...\n")
rects <- read.csv(file.path(SRC, "lesion_rects_expanded.csv"), stringsAsFactors = FALSE)
rects$lesion <- factor(rects$lesion, levels = paste0("L", 1:6))
df_ef <- readRDS(file.path(SRC, "Fig5EF_clusters.rds"))
he_e  <- readRDS(file.path(SRC, "Fig5E_HE_lowres.rds"))
he_g  <- readRDS(file.path(SRC, "HE_Brain_lowres.rds"))
nr <- dim(he_e)[1]; nc <- dim(he_e)[2]
nrg <- dim(he_g)[1]; ncg <- dim(he_g)[2]
gj <- readRDS(file.path(SRC, "Fig5GJ_spatial_scores.rds"))
box_long <- readRDS(file.path(SRC, "Fig5GJ_boxplot_long.rds"))

assign_win <- function(d) {
  d$lesion <- NA_character_
  for (i in seq_len(nrow(rects))) {
    r <- rects[i, ]
    hit <- in_rect(d$x_img, d$y_img, r) & is.na(d$lesion)
    d$lesion[hit] <- as.character(r$lesion)
  }
  d$set <- ifelse(d$lesion %in% c("L1", "L2", "L3"), "L1-L3",
                  ifelse(d$lesion %in% c("L4", "L5", "L6"), "L4-L6", NA))
  d
}
df_ef <- assign_win(df_ef)
gj <- assign_win(gj)
gj_les <- gj[!is.na(gj$set), , drop = FALSE]
cat("ST spots", nrow(df_ef), " immune", nrow(gj), " in windows", nrow(gj_les), "\n")
print(table(gj$lesion, useNA = "ifany"))

feat_keep <- c("Tcell EAE_POS", "Tcell EAE_NEG", "Mg EAE_POS", "Mg EAE_NEG")
box_les <- gj_les %>%
  dplyr::select(cell, set, lesion, dplyr::all_of(feat_keep)) %>%
  tidyr::pivot_longer(dplyr::all_of(feat_keep), names_to = "feature", values_to = "score") %>%
  dplyr::mutate(
    celltype = ifelse(grepl("^Tcell", feature), "Tcell", "Mg"),
    group = ifelse(grepl("EAE_POS", feature), "EAE_POS", "EAE_NEG")
  )

# Visium -> Fig5C IF pixels (right-panel tissue bbox)
img_c <- png::readPNG(file.path(ASSET, "Fig5C_IF.png"))
nr_c <- dim(img_c)[1]; nc_c <- dim(img_c)[2]
ix0 <- 542; ix1 <- 1073; iy0 <- 19; iy1 <- 506
sx <- (ix1 - ix0) / 600; sy <- (iy1 - iy0) / 595
if_rects <- transform(
  rects,
  xmin = ix0 + xmin * sx,
  xmax = ix0 + xmax * sx,
  ymin = iy0 + ymin * sy,
  ymax = iy0 + ymax * sy
)
if_new <- if_rects[if_rects$set == "added", ]
write.csv(if_rects, file.path(PLOTTED, "SourceData_Fig8c_if_rects.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# a–c
# -----------------------------------------------------------------------------
cat("Panels a-c ...\n")
pA_st <- spat_cluster(df_ef, he_e, nr, nc, "Spatial transcriptome",
                      pt = 0.12, show_lab = FALSE, rects = NULL,
                      show_legend = FALSE)
pA <- wrap_elements(
  panel = (
    (panel_sagittal() | flow_arrow() | panel_coronal()) /
      (panel_slide(he_e) | flow_arrow() | pA_st)
  ) + plot_layout(widths = c(1.15, 0.22, 1.15), heights = c(1, 1.05))
) + labs(tag = "a") + tag_theme
pB <- embed_png(file.path(ASSET, "Fig5B_IF.png"), "b", ink = "white")

pC <- ggplot() +
  annotation_raster(grDevices::as.raster(pmin(pmax(img_c, 0), 1)),
                    xmin = 0.5, xmax = nc_c + 0.5,
                    ymin = -(nr_c + 0.5), ymax = -0.5) +
  geom_rect(
    data = if_new,
    aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
    fill = NA, colour = "yellow", linewidth = 0.5, linetype = "dashed"
  ) +
  geom_text(
    data = if_new,
    aes(x = (xmin + xmax) / 2, y = -ymin + 12, label = lesion),
    colour = "yellow", size = PT(5.2), fontface = "bold", family = FAM
  ) +
  geom_text(
    data = data.frame(
      lesion = c("L1", "L2", "L3"),
      x = c(655, 766, 930),
      y = c(92, 30, 80)
    ),
    aes(x = x, y = -y, label = lesion),
    inherit.aes = FALSE, colour = "white", size = PT(5.2),
    fontface = "bold", family = FAM
  ) +
  coord_fixed(xlim = c(0.5, nc_c + 0.5), ylim = c(-(nr_c + 0.5), -0.5),
              expand = FALSE) +
  theme_void(base_family = FAM) +
  labs(tag = "c") + tag_theme +
  theme(plot.margin = margin(2, 2, 2, 2))

# -----------------------------------------------------------------------------
# d  original L1–L3 IF composite + L4–L6 crops from the high-res plate
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
pD_old <- embed_png(file.path(ASSET, "Fig5D_IF.png")) +
  labs(title = "L1–L3") +
  theme(plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5,
                                  family = FAM))
pic <- png::readPNG(file.path(ASSET, "Fig8d_L456_source.png"))
nr_p <- dim(pic)[1]; nc_p <- dim(pic)[2]
added <- rects[rects$set == "added", ]
pic_rects <- data.frame(
  lesion = added$lesion,
  xmin = added$xmin / 600 * nc_p,
  xmax = added$xmax / 600 * nc_p,
  ymin = added$ymin / 595 * nr_p,
  ymax = added$ymax / 595 * nr_p,
  stringsAsFactors = FALSE
)
write.csv(pic_rects, file.path(PLOTTED, "SourceData_Fig8d_pic1_crops.csv"),
          row.names = FALSE)
crop_pic <- function(r, pad = 18) {
  x1 <- max(1L, floor(r$xmin) - pad); x2 <- min(nc_p, ceiling(r$xmax) + pad)
  y1 <- max(1L, floor(r$ymin) - pad); y2 <- min(nr_p, ceiling(r$ymax) + pad)
  pic[y1:y2, x1:x2, 1:3, drop = FALSE]
}
d_new <- lapply(seq_len(nrow(pic_rects)), function(i) {
  r <- pic_rects[i, ]
  embed_array(crop_pic(r), title = as.character(r$lesion))
})
right_d <- d_new[[1]] | d_new[[2]] | d_new[[3]]
pD <- wrap_elements(
  panel = (pD_old | right_d) + plot_layout(widths = c(1.05, 1.45))
) + labs(tag = "d") + tag_theme

# -----------------------------------------------------------------------------
# e–f
# -----------------------------------------------------------------------------
cat("Panels e-f ...\n")
pE_he <- ggplot() +
  annotation_raster(he_e, xmin = 0, xmax = nc, ymin = -nr, ymax = 0) +
  geom_rect(
    data = rects, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin),
    fill = NA, colour = "yellow", linewidth = 0.4, linetype = "dashed"
  ) +
  geom_text(
    data = rects, aes(x = (xmin + xmax) / 2, y = -ymin + 10, label = lesion),
    colour = "yellow", size = PT(5.2), fontface = "bold", family = FAM
  ) +
  coord_fixed(xlim = c(0, nc), ylim = c(-nr, 0), expand = FALSE) +
  labs(title = "H&E") +
  theme_void(base_family = FAM) +
  theme(plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5, family = FAM))

pE_cl <- spat_cluster(df_ef, he_e, nr, nc, "Clusters", pt = 0.2,
                      rects = rects, show_legend = TRUE, show_lab = FALSE)
pE <- wrap_elements(
  panel = (pE_he | pE_cl) + plot_layout(widths = c(1, 1.38))
) + labs(tag = "e") + tag_theme

# -----------------------------------------------------------------------------
# g–j
# -----------------------------------------------------------------------------
cat("Panels f-g ...\n")
# Whole-section maps only. Window-restricted maps repeat the same
# geography, so the window contrast is kept in the boxplots.
score_panel <- function(map_l, map_r, box_all, box_win) {
  maps <- wrap_elements(panel = map_l | map_r)
  wrap_elements(panel = (maps | (box_all / box_win)) +
                  plot_layout(widths = c(2.15, 0.95), guides = "keep"))
}
pF <- score_panel(
  spat_score(gj, "Tcell EAE_POS", "T cell POS", 0.3, he_g, nrg, ncg, 0.22, rects),
  spat_score(gj, "Tcell EAE_NEG", "T cell NEG", 0.3, he_g, nrg, ncg, 0.22, rects),
  make_box(box_long %>% filter(celltype == "Tcell"), "All immune spots"),
  make_box_split(box_les %>% filter(celltype == "Tcell"), "Windows")
) + labs(tag = "f") + tag_theme
pG <- score_panel(
  spat_score(gj, "Mg EAE_POS", "Microglia POS", 0.45, he_g, nrg, ncg, 0.22, rects),
  spat_score(gj, "Mg EAE_NEG", "Microglia NEG", 0.45, he_g, nrg, ncg, 0.22, rects),
  make_box(box_long %>% filter(celltype == "Mg"), "All immune spots"),
  make_box_split(box_les %>% filter(celltype == "Mg"), "Windows")
) + labs(tag = "g") + tag_theme

# -----------------------------------------------------------------------------
# h  Microglia POS/NEG window crops (Fig. S6 d style) + f–g boxplots
# -----------------------------------------------------------------------------
cat("Panel h ...\n")
mg_cut <- 0.45
les_map <- function(nm, score_col, title) {
  sub <- gj[!is.na(gj$lesion) & gj$lesion == as.character(nm), , drop = FALSE]
  ht <- max(sub$y_img, na.rm = TRUE) - min(sub$y_img, na.rm = TRUE)
  pad_use <- if (is.finite(ht) && ht < 50) 18 else 12
  cv <- crop_view(sub, he_g, pad = pad_use)
  spat_score(cv$df, score_col, title, mg_cut, cv$he, cv$nr, cv$nc,
             pt = 0.9, rects = NULL) +
    theme(legend.position = "none")
}
score_legend <- function(max_cut) {
  dummy <- ggplot(data.frame(x = 1, y = 1, z = c(0, max_cut)),
                  aes(x, y, colour = z)) +
    geom_point(size = 0.01, alpha = 0) +
    scale_colour_gradientn(
      colours = score_pal, limits = c(0, max_cut), name = NULL, oob = squish
    ) +
    guides(colour = guide_colourbar(
      barwidth = unit(1.8, "mm"),
      barheight = unit(28, "mm"),
      ticks.colour = "black",
      frame.colour = "black",
      frame.linewidth = 0.25
    )) +
    theme_void(base_family = FAM) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 5.2, family = FAM, colour = "black"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
  g <- ggplotGrob(dummy)
  idx <- which(grepl("guide-box", g$layout$name))
  if (!length(idx)) stop("score_legend: no guide-box grob")
  wrap_elements(full = g$grobs[[idx[[1]]]])
}
maps_h <- wrap_elements(panel = wrap_plots(
  les_map("L1", "Mg EAE_POS", "L1 POS"),
  les_map("L1", "Mg EAE_NEG", "L1 NEG"),
  les_map("L4", "Mg EAE_POS", "L4 POS"),
  les_map("L4", "Mg EAE_NEG", "L4 NEG"),
  les_map("L2", "Mg EAE_POS", "L2 POS"),
  les_map("L2", "Mg EAE_NEG", "L2 NEG"),
  les_map("L5", "Mg EAE_POS", "L5 POS"),
  les_map("L5", "Mg EAE_NEG", "L5 NEG"),
  les_map("L3", "Mg EAE_POS", "L3 POS"),
  les_map("L3", "Mg EAE_NEG", "L3 NEG"),
  les_map("L6", "Mg EAE_POS", "L6 POS"),
  les_map("L6", "Mg EAE_NEG", "L6 NEG"),
  ncol = 4, nrow = 3
) + plot_layout(heights = c(1, 1, 1)))
pH <- wrap_elements(
  panel = (maps_h | score_legend(mg_cut) | (
    make_box(box_long %>% filter(celltype == "Mg"), "All immune spots") /
      make_box_split(box_les %>% filter(celltype == "Mg"), "Windows")
  )) + plot_layout(widths = c(2.00, 0.20, 0.88), guides = "keep")
) + labs(tag = "h") + tag_theme

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
fig <- ((pA | pB) + plot_layout(widths = c(1.12, 1.20))) / pC / pD / pE / pF / pG / pH +
  plot_layout(heights = grid::unit.c(
    grid::unit(58, "mm"),
    grid::unit(54, "mm"),
    grid::unit(68, "mm"),
    grid::unit(54, "mm"),
    grid::unit(72, "mm"),
    grid::unit(72, "mm"),
    grid::unit(92, "mm")
  ), guides = "keep")
# 58+54+68+54+72+72+92 = 470

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "Figure8")

cat("Writing SVG ...\n")
svglite::svglite(paste0(stem, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in,
                     family = FAM, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PNG ...\n")
ragg::agg_png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())

cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(stem, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); invisible(dev.off())

write.csv(rects, file.path(PLOTTED, "SourceData_Fig8_lesion_rects.csv"), row.names = FALSE)
n_win <- as.data.frame(table(gj$lesion, useNA = "ifany"), stringsAsFactors = FALSE)
names(n_win) <- c("lesion", "n_immune")
n_st <- as.data.frame(table(df_ef$lesion, useNA = "ifany"), stringsAsFactors = FALSE)
names(n_st) <- c("lesion", "n_st")
write.csv(merge(n_st, n_win, by = "lesion", all = TRUE),
          file.path(PLOTTED, "SourceData_Fig8_window_n.csv"), row.names = FALSE)
write.csv(box_les, file.path(PLOTTED, "SourceData_Fig8hj_window_scores.csv"),
          row.names = FALSE)

wt_rows <- list()
for (ct in c("Tcell", "Mg")) {
  d0 <- box_long[box_long$celltype == ct, ]
  wt_rows[[length(wt_rows) + 1]] <- data.frame(
    celltype = ct, scope = "all", set = "all",
    n = nrow(d0), p = wilcox.test(score ~ group, data = d0)$p.value
  )
  for (st in c("L1-L3", "L4-L6")) {
    d1 <- box_les[box_les$celltype == ct & box_les$set == st, ]
    wt_rows[[length(wt_rows) + 1]] <- data.frame(
      celltype = ct, scope = "windows", set = st,
      n = nrow(d1), p = wilcox.test(score ~ group, data = d1)$p.value
    )
  }
}
write.csv(bind_rows(wt_rows), file.path(PLOTTED, "SourceData_Fig8_wilcoxon.csv"),
          row.names = FALSE)
print(bind_rows(wt_rows))

n_imm <- nrow(gj)
n_l13 <- sum(gj$set == "L1-L3", na.rm = TRUE)
n_l46 <- sum(gj$set == "L4-L6", na.rm = TRUE)
writeLines(c(
  "Fig. 8 | Perivascular windows and vessel-contact scores in the EAE brain",
  "",
  "a, Schematic of EAE brain sectioning, slice collection and Visium spatial transcriptomics. b, Laminin/DAPI sagittal overview. c, DAPI and Laminin/CD45. Dashed yellow boxes mark six windows: original L1–L3 and added L4–L6. d, Laminin/CD45 zooms. L1–L3 are the original three fields; L4–L6 are the added windows cropped from the high-resolution coronal plate. e, Visium H&E and cluster map with the six windows. Cluster identity is in the legend (0/1 neurons, 2 hippocampal neurons, 3 oligodendrocytes, 4 Purkinje cells, 5 granule cells, 6 meninges, 7 vascular-associated, 8 plasma cells). Cluster numbers are not printed on the maps. f, T-cell POS versus NEG prediction scores on the whole section (immune spots, n = 32,811). Boxplots compare all immune spots and, separately, spots inside the six windows, split as L1–L3 versus L4–L6. g, The same layout for microglia. h, Microglia POS versus NEG scores inside each window (L1–L3 left, L4–L6 right), with the same all-spot and window boxplots as g. Boxes show the median and interquartile range; whiskers extend to 1.5\u00d7IQR. Stars are two-sided Wilcoxon rank-sum tests of POS versus NEG within each set. Window sets are tested separately because L4–L6 need not repeat the L1–L3 direction.",
  "",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "Figure8_legend.txt"))
writeLines(c(
  "Figure 8 statistics / integrity (Fig8_Nature_v4)",
  sprintf("Immune Visium spots: %d.", n_imm),
  sprintf("Windows: L1-L3 n = %d; L4-L6 n = %d.", n_l13, n_l46),
  "Panel a: vector schematic (not the v3 raster plate).",
  "Panel e: L1-L6 cluster zooms omitted.",
  "Panel h: microglia window crops + Wilcoxon as in g.",
  "Image integrity: IF author rasters; H&E unmodified; no contrast stretch.",
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "Figure8_stats.txt"))
writeLines(c(
  "Figure 8 | new folder Fig8_Nature_v4 (does not overwrite v1-v3).",
  "Reuse: structural adaptation of v3; h from Fig. S6 d crop layout.",
  "Panel a redrawn as vector sectioning/Visium schematic.",
  "Backend: R exclusive."
), file.path(OUTDIR, "Figure8_provenance.txt"))
cat("DONE ->", OUTDIR, "\n")
