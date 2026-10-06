#!/usr/bin/env Rscript
# =============================================================================
# Fig. 8 v2 — Nature main figure (R / ggplot2 exclusive)
#
# Core conclusion
#   Six perivascular windows (original L1–L3 plus added L4–L6) mark
#   vessel-associated niches; T-cell and microglia POS versus NEG
#   prediction scores are shown for the whole section and, in the same
#   panel, inside Original versus Added windows.
#
# Evidence chain
#   a    Sampling schematic
#   b    Sagittal Laminin/DAPI
#   c    Coronal IF with L1–L6 labels
#   d    L1–L3 IF zooms plus L4–L6 crops from the high-res plate
#   e    H&E and clusters (external legend; no on-plot cluster numbers) plus zooms
#   f    T-cell maps once; all-spot and window boxplots
#   g    Microglia maps once; all-spot and window boxplots
#
# Archetype : image plate + quant
# Size      : 183 x 450 mm
# New folder — does not overwrite Fig8_Nature_v1
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
  file.path(IV_FIG_ROOT, "Fig8_Nature_v3")
}
SRC     <- file.path(OUTDIR, "source_data")
ASSET   <- file.path(OUTDIR, "assets")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183  # fig_width_mm = 183
FIG_H_MM <- 450
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
pA <- embed_png(file.path(ASSET, "Fig5A_IF.png"), "a")
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
les_plots <- lapply(paste0("L", 1:6), function(nm) {
  sub <- df_ef[!is.na(df_ef$lesion) & df_ef$lesion == nm, , drop = FALSE]
  cv <- crop_view(sub, he_e, pad = 8)
  spat_cluster(cv$df, cv$he, cv$nr, cv$nc, nm, pt = 0.55,
               show_lab = FALSE, rects = NULL, show_legend = FALSE)
})
pE <- wrap_elements(
  panel = ((pE_he | pE_cl) + plot_layout(widths = c(1, 1.38))) /
    ((les_plots[[1]] | les_plots[[2]] | les_plots[[3]]) /
       (les_plots[[4]] | les_plots[[5]] | les_plots[[6]])) +
    plot_layout(heights = c(1.05, 1.15))
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
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
fig <- (pA | pB) / pC / pD / pE / pF / pG +
  plot_layout(heights = grid::unit.c(
    grid::unit(56, "mm"),
    grid::unit(56, "mm"),
    grid::unit(72, "mm"),
    grid::unit(96, "mm"),
    grid::unit(78, "mm"),
    grid::unit(92, "mm")
  ), guides = "keep")
# 56+56+72+96+78+92 = 450

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
cat("DONE ->", OUTDIR, "\n")
