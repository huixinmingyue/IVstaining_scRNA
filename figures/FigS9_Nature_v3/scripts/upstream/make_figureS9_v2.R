#!/usr/bin/env Rscript
# =============================================================================
# Fig. S9 v2 — Nature multipanel (FigS9_Nature_v2)
# New folder; does not overwrite FigS9_Nature_v1.
#
# Changed relative to v1
#   former h (Mg maps) removed
#   former i is now h, without Mg
#   c-g titles: POS / NEG (Fig. 8 f-g style), not IV / IVC
#   h  All immune spots and windows L1-L3 vs L4-L6 (Fig. 8 f-g stats)
#
# Core conclusion
#   EAE Visium resolves cluster marker programmes and maps POS versus NEG
#   immune-lineage prediction scores onto the section, with quantified
#   differences among all immune spots and inside L1-L3 versus L4-L6.
#
# Archetype : image plate + quant
# Backend   : R exclusive
# Size      : 183 x 248 mm
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(grid)
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
  file.path(IV_FIG_ROOT, "FigS9_Nature_v2")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
dir.create(PLOTTED, showWarnings = FALSE, recursive = TRUE)
dir.create(QA, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 248
width_mm = 183
height_mm = 248
FAM  <- "Arial"
BASE <- 6.5
PT_SIZE <- 1.00
PT <- function(x) x / ggplot2::.pt
has_sm <- requireNamespace("scattermore", quietly = TRUE)

pal_feat <- c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF8C00", "#FF0000")
clust_lvls <- as.character(0:8)
clust_pal <- setNames(
  c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
    "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"),
  clust_lvls
)
# Fig. 8 f-g fill: POS purple, NEG orange
pos_neg_cols <- c(POS = "#9B59B6", NEG = "#E67E22")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = BASE),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.4),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(3.0, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid = element_blank()
  )
theme_set(theme_pub)
tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.01, 0.99)
)

fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}

fix_images <- function(sce) {
  for (nm in names(sce@images)) {
    img <- sce@images[[nm]]
    ok <- tryCatch({ slot(img, "misc"); TRUE }, error = function(e) FALSE)
    if (!ok) {
      slot(img, "misc", check = FALSE) <- list()
      sce@images[[nm]] <- img
    }
  }
  sce
}

fmtp <- function(p) {
  if (!is.finite(p) || p < 2.2e-16) "<2.2e-16" else sprintf("%.3g", p)
}

star_lab <- function(p) {
  if (is.na(p)) return("")
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 0.01) "**" else
    if (p < 0.05) "*" else "ns"
}

in_rect <- function(x, y, r) {
  x >= r$xmin & x <= r$xmax & y >= r$ymin & y <= r$ymax
}

hm <- readRDS(file.path(SRC, "heatmap_slim.rds"))
full <- readRDS(file.path(SRC, "spatial_full_slim.rds"))
imm_slim <- readRDS(file.path(SRC, "spatial_immune_slim.rds"))
rects <- read.csv(file.path(SRC, "lesion_rects_expanded.csv"), stringsAsFactors = FALSE)
rects$lesion <- factor(rects$lesion, levels = paste0("L", 1:6))

# -----------------------------------------------------------------------------
# a — cluster marker heatmap
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
set.seed(1)
mk <- hm$markers
mk$cluster <- as.character(mk$cluster)
gene_keep <- unique(unlist(lapply(split(mk, mk$cluster), function(d) {
  head(unique(as.character(d$gene)), 5)
}), use.names = FALSE))
gene_ord2 <- intersect(hm$gene_order, gene_keep)
mat <- hm$mat[gene_ord2, , drop = FALSE]
cl <- hm$cell_cluster
n_a_before <- length(cl)
keep <- unlist(lapply(split(seq_along(cl), cl), function(ix) {
  if (length(ix) <= 60) ix else sample(ix, 60)
}), use.names = FALSE)
mat <- mat[, keep, drop = FALSE]
cl <- cl[keep]
n_a_after <- length(cl)
mat_c <- pmax(pmin(mat, 2.5), -2.5)

write.csv(
  mk[, intersect(c("gene", "cluster", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"),
                 names(mk)), drop = FALSE],
  file.path(PLOTTED, "SourceData_FigS9a_markers.csv"),
  row.names = FALSE
)

df_a <- as.data.frame(as.table(mat_c), stringsAsFactors = FALSE)
colnames(df_a) <- c("gene", "cell", "z")
df_a$gene <- factor(df_a$gene, levels = rev(rownames(mat)))
df_a$cell <- factor(df_a$cell, levels = colnames(mat_c))
bar <- data.frame(
  cell = factor(colnames(mat_c), levels = colnames(mat_c)),
  cluster = factor(cl, levels = clust_lvls),
  y = 1
)
pA_bar <- ggplot(bar, aes(cell, y, fill = cluster)) +
  geom_tile(height = 1, width = 1, colour = NA) +
  scale_fill_manual(values = clust_pal, name = "Cluster", drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_family = FAM) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = BASE - 0.5),
    legend.text = element_text(size = BASE - 1),
    legend.key.size = unit(2.5, "mm"),
    plot.margin = margin(0, 0, 0, 0)
  )
pA_hm <- ggplot(df_a, aes(cell, gene, fill = z)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradientn(
    colours = c("#FF00FF", "#000000", "#FFFF00"),
    limits = c(-2.5, 2.5), oob = squish, name = "Scaled\nexpr."
  ) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = "a") + theme_pub +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 5.2, face = "italic", family = FAM),
    legend.key.height = unit(6, "mm"), legend.key.width = unit(2.2, "mm"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.3)
  )
pA <- (pA_bar / pA_hm) + plot_layout(heights = c(0.06, 1), guides = "collect")

# -----------------------------------------------------------------------------
# b — gene spatial maps
# -----------------------------------------------------------------------------
spot_gene <- function(gene, letter = NULL, limits = NULL, pt = 0.12, show_legend = TRUE) {
  vals <- as.numeric(full$expr[gene, ]); names(vals) <- colnames(full$expr)
  if (is.null(limits)) {
    hi <- as.numeric(quantile(vals[vals > 0], 0.995, names = FALSE, na.rm = TRUE))
    if (!is.finite(hi) || hi <= 0) hi <- max(vals, na.rm = TRUE)
    limits <- c(0, max(hi, 0.1))
  }
  df <- data.frame(
    x = full$coords[, "x"], y = -full$coords[, "y"],
    val = pmin(pmax(vals[rownames(full$coords)], limits[1]), limits[2])
  )
  df <- df[order(df$val), , drop = FALSE]
  arr <- fade_he(full$he$image, 0.90)
  rast <- grDevices::as.raster(arr)
  ymin <- -(dim(arr)[1] + 0.5); ymax <- -0.5
  xmin <- 0.5; xmax <- dim(arr)[2] + 0.5
  p <- ggplot(df, aes(x, y, colour = val)) +
    annotation_raster(rast, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  if (has_sm) {
    p <- p + scattermore::geom_scattermore(pointsize = 1.05, pixels = c(800, 800),
                                           interpolate = TRUE)
  } else {
    p <- p + geom_point(size = pt, stroke = 0)
  }
  p <- p +
    scale_colour_gradientn(colours = pal_feat, limits = limits, oob = squish,
                           name = "Expr.") +
    coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    labs(title = gene, tag = letter) +
    theme_void(base_family = FAM) +
    theme(
      plot.title = element_text(size = BASE, face = "bold.italic", hjust = 0.5,
                                family = FAM),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.tag.position = c(0.01, 0.99),
      legend.title = element_text(size = BASE - 0.5),
      legend.text = element_text(size = BASE - 0.8),
      legend.key.height = unit(5, "mm"), legend.key.width = unit(2.0, "mm"),
      plot.margin = margin(1, 1, 1, 1)
    )
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

cat("Panel b ...\n")
lim_mki <- c(0, max(2, as.numeric(quantile(full$expr["Mki67", ], 0.995))))
lim_lga <- c(0, max(2, as.numeric(quantile(full$expr["Lgals3", ], 0.995))))
lim_ccl <- c(0, max(2, as.numeric(quantile(full$expr["Ccl2", ], 0.995))))
lim_hex <- c(0, max(2, as.numeric(quantile(full$expr["Hexb", ], 0.995))))
pB <- (spot_gene("Mki67", "b", lim_mki) | spot_gene("Lgals3", NULL, lim_lga)) /
  (spot_gene("Ccl2", NULL, lim_ccl) | spot_gene("Hexb", NULL, lim_hex))

# -----------------------------------------------------------------------------
# c-g — SpatialFeaturePlot; POS / NEG titles (Fig. 8 f-g style)
# -----------------------------------------------------------------------------
cat("Loading predictions object for c-g ...\n")
spe_i <- fix_images(readRDS(file.path(SRC, "spe_immune_predictions.rds")))
DefaultAssay(spe_i) <- "predictions"
n_immune <- ncol(spe_i)
cat("spe_immune dim:", paste(dim(spe_i), collapse = "x"), "\n")

theme_sf <- theme_void(base_family = FAM) +
  theme(
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.01, 0.99),
    legend.title = element_text(size = BASE - 0.5, family = FAM),
    legend.text = element_text(size = BASE - 0.8, family = FAM),
    legend.key.height = unit(5, "mm"),
    legend.key.width = unit(2, "mm"),
    plot.margin = margin(1, 1, 1, 1)
  )

sf_one <- function(feat, max.cutoff, title, letter = NULL, show_legend = TRUE) {
  args <- list(
    object = spe_i, features = feat, pt.size.factor = PT_SIZE,
    min.cutoff = 0, crop = TRUE, stroke = 0,
    image.alpha = 0.9, alpha = c(0.1, 1)
  )
  if (!is.null(max.cutoff)) args$max.cutoff <- max.cutoff
  p <- do.call(SpatialFeaturePlot, args)
  lim <- if (is.null(max.cutoff)) {
    v <- as.numeric(GetAssayData(spe_i, assay = "predictions", layer = "data")[feat, ])
    c(0, max(0.1, as.numeric(quantile(v, 0.995, names = FALSE, na.rm = TRUE))))
  } else c(0, max.cutoff)
  p <- p +
    scale_fill_gradientn(colours = pal_feat, limits = lim, name = "Score") +
    ggtitle(title) + labs(tag = letter) + theme_sf
  if (!show_legend) p <- p + guides(fill = "none")
  p
}

sf_pair <- function(f_pos, f_neg, max.cutoff, letter, lineage) {
  cat("  SpatialFeaturePlot", letter, lineage, "\n")
  left <- sf_one(f_pos, max.cutoff, paste(lineage, "POS"), letter, FALSE)
  right <- sf_one(f_neg, max.cutoff, paste(lineage, "NEG"), NULL, TRUE)
  wrap_elements(full = (left | right))
}

q995 <- function(feat) {
  as.numeric(quantile(
    GetAssayData(spe_i, assay = "predictions", layer = "data")[feat, ],
    0.995, names = FALSE, na.rm = TRUE
  ))
}

cat("Building panels c-g ...\n")
pC <- sf_pair("DCIV", "DCIVC", 0.65, "c", "DC")
pD <- sf_pair("BcellIV", "BcellIVC", 0.15, "d", "B cell")
pE <- sf_pair("Neu1IV", "Neu1IVC", 0.08, "e", "Neu1")
pc_hi <- max(q995("PCIV"), q995("PCIVC"), 0.1)
pF <- sf_pair("PCIV", "PCIVC", pc_hi, "f", "PC")
pG <- sf_pair("MacIV", "MacIVC", 0.45, "g", "Mac")
rm(spe_i); gc(verbose = FALSE)

# -----------------------------------------------------------------------------
# h — All immune spots + windows L1-L3 / L4-L6 (Fig. 8 f-g; no Mg)
# -----------------------------------------------------------------------------
cat("Panel h ...\n")
spots <- data.frame(
  cell = colnames(imm_slim$pred),
  x_img = imm_slim$coords[colnames(imm_slim$pred), "x"],
  y_img = imm_slim$coords[colnames(imm_slim$pred), "y"],
  stringsAsFactors = FALSE
)
spots$lesion <- NA_character_
for (i in seq_len(nrow(rects))) {
  r <- rects[i, ]
  hit <- in_rect(spots$x_img, spots$y_img, r) & is.na(spots$lesion)
  spots$lesion[hit] <- as.character(r$lesion)
}
spots$set <- ifelse(spots$lesion %in% c("L1", "L2", "L3"), "L1-L3",
                    ifelse(spots$lesion %in% c("L4", "L5", "L6"), "L4-L6", NA))
n_win <- as.data.frame(table(spots$lesion, useNA = "ifany"), stringsAsFactors = FALSE)
names(n_win) <- c("lesion", "n_immune")
print(n_win)
n_l13 <- sum(spots$set == "L1-L3", na.rm = TRUE)
n_l46 <- sum(spots$set == "L4-L6", na.rm = TRUE)
cat("windows L1-L3", n_l13, " L4-L6", n_l46, "\n")
write.csv(rects, file.path(PLOTTED, "SourceData_FigS9_lesion_rects.csv"),
          row.names = FALSE)
write.csv(n_win, file.path(PLOTTED, "SourceData_FigS9h_window_n.csv"),
          row.names = FALSE)
write.csv(spots[, c("cell", "lesion", "set")],
          file.path(PLOTTED, "SourceData_FigS9h_spot_windows.csv"),
          row.names = FALSE)

pairs <- list(
  DC = c("DCIV", "DCIVC"), B = c("BcellIV", "BcellIVC"),
  Neu1 = c("Neu1IV", "Neu1IVC"), PC = c("PCIV", "PCIVC"),
  Mac = c("MacIV", "MacIVC")
)
labs_h <- c(DC = "DC", B = "B cell", Neu1 = "Neu1", PC = "PC", Mac = "Mac")

long_all <- lapply(names(pairs), function(nm) {
  f_pos <- pairs[[nm]][1]; f_neg <- pairs[[nm]][2]
  data.frame(
    cell = rep(colnames(imm_slim$pred), 2),
    lineage = nm,
    group = factor(rep(c("POS", "NEG"), each = ncol(imm_slim$pred)),
                   levels = c("POS", "NEG")),
    score = c(as.numeric(imm_slim$pred[f_pos, ]),
              as.numeric(imm_slim$pred[f_neg, ])),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
long_all <- long_all[is.finite(long_all$score), , drop = FALSE]
long_win <- merge(long_all, spots[, c("cell", "set", "lesion")], by = "cell")
long_win <- long_win[!is.na(long_win$set), , drop = FALSE]
long_win$set <- factor(long_win$set, levels = c("L1-L3", "L4-L6"))

make_box_all <- function(df, title, letter = NULL, show_y = TRUE, show_legend = FALSE) {
  d <- df
  wt <- wilcox.test(score ~ group, data = d)
  ymax <- max(d$score, na.rm = TRUE)
  p <- ggplot(d, aes(group, score, fill = group)) +
    geom_boxplot(width = 0.55, outlier.size = 0.22, outlier.alpha = 0.3,
                 linewidth = 0.3, alpha = 0.85) +
    geom_text(data = data.frame(x = 1.5, y = ymax * 1.08, lab = star_lab(wt$p.value)),
              aes(x, y, label = lab), inherit.aes = FALSE,
              size = PT(5.2), family = FAM) +
    scale_fill_manual(values = pos_neg_cols, guide = if (show_legend) "legend" else "none",
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    labs(x = title, y = if (show_y) "Prediction score" else NULL, tag = letter,
         title = NULL) +
    theme_pub +
    theme(
      axis.text.x = element_text(size = 5.2),
      plot.margin = margin(2, 2, 4, 2)
    )
  list(plot = p, p = wt$p.value, n = nrow(d),
       n_POS = sum(d$group == "POS"), n_NEG = sum(d$group == "NEG"))
}

make_box_split <- function(df, title, show_y = TRUE, show_legend = FALSE) {
  d <- df
  stats <- d %>%
    group_by(set) %>%
    summarise(
      p = wilcox.test(score ~ group)$p.value,
      y = max(score, na.rm = TRUE) * 1.08,
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(lab = vapply(p, star_lab, character(1)), x = as.numeric(set))
  p <- ggplot(d, aes(set, score, fill = group)) +
    geom_boxplot(width = 0.7, outlier.size = 0.18, outlier.alpha = 0.28,
                 linewidth = 0.3, alpha = 0.85,
                 position = position_dodge(width = 0.75)) +
    geom_text(data = stats, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, size = PT(5.2), family = FAM) +
    scale_fill_manual(values = pos_neg_cols, name = NULL,
                      guide = if (show_legend) "legend" else "none") +
    scale_x_discrete(labels = c("L1-L3" = "L1\u2013L3", "L4-L6" = "L4\u2013L6")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.24))) +
    labs(x = title, y = if (show_y) "Prediction score" else NULL) +
    theme_pub +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.margin = margin(0, 0, 0, 0),
      legend.key.size = unit(2.4, "mm"),
      legend.text = element_text(size = 5.2),
      axis.text.x = element_text(size = 5.2),
      plot.margin = margin(2, 2, 4, 2)
    )
  list(plot = p, stats = stats)
}

all_plots <- list(); win_plots <- list()
wt_rows <- list()
i <- 0
for (nm in names(pairs)) {
  i <- i + 1
  d_all <- long_all[long_all$lineage == nm, ]
  d_win <- long_win[long_win$lineage == nm, ]
  res_a <- make_box_all(
    d_all, labs_h[[nm]],
    letter = NULL,
    show_y = i == 1
  )
  res_w <- make_box_split(
    d_win, labs_h[[nm]],
    show_y = i == 1, show_legend = i == 1
  )
  all_plots[[nm]] <- res_a$plot
  win_plots[[nm]] <- res_w$plot
  wt_rows[[length(wt_rows) + 1]] <- data.frame(
    lineage = nm, scope = "all", set = "all",
    n = res_a$n, n_spots = res_a$n_POS, wilcox_p = res_a$p,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(nrow(res_w$stats))) {
    st <- res_w$stats[j, ]
    wt_rows[[length(wt_rows) + 1]] <- data.frame(
      lineage = nm, scope = "windows", set = as.character(st$set),
      n = st$n, n_spots = st$n / 2, wilcox_p = st$p,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("  %s all p=%s n=%d | L1-L3 p=%s | L4-L6 p=%s\n",
              nm, fmtp(res_a$p), res_a$n_POS,
              fmtp(res_w$stats$p[res_w$stats$set == "L1-L3"]),
              fmtp(res_w$stats$p[res_w$stats$set == "L4-L6"])))
}
stats_h <- bind_rows(wt_rows)
write.csv(stats_h, file.path(PLOTTED, "SourceData_FigS9h_wilcoxon.csv"),
          row.names = FALSE)
write.csv(
  long_all[, c("cell", "lineage", "group", "score")],
  file.path(PLOTTED, "SourceData_FigS9h_all_scores.csv"),
  row.names = FALSE
)
write.csv(
  long_win[, c("cell", "lineage", "group", "score", "set", "lesion")],
  file.path(PLOTTED, "SourceData_FigS9h_window_scores.csv"),
  row.names = FALSE
)

hdr <- function(txt, letter = NULL) {
  p <- ggplot() +
    annotate(
      "text", x = 0.045, y = 0.42, label = txt,
      fontface = "bold", size = PT(6.2), family = FAM, hjust = 0
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(0, 0, 0, 0))
  if (!is.null(letter)) {
    p <- p + labs(tag = letter) +
      theme(
        plot.tag = element_text(size = 8, face = "bold", family = FAM),
        plot.tag.position = c(0.008, 0.5)
      )
  }
  p
}

all_row <- wrap_plots(all_plots, nrow = 1)
win_row <- wrap_plots(win_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "top")
hdr_all <- hdr("All immune spots", "h")
hdr_win <- hdr("Windows")

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- (pA | pB) + plot_layout(widths = c(1.0, 1.15))
row2 <- (pC | pD | pE)
row3 <- (pF | pG)
# mm-locked header rows at the top-level patchwork (same device as Fig. 8).
fig <- (row1 / row2 / row3 / hdr_all / all_row / hdr_win / win_row) +
  plot_layout(heights = grid::unit.c(
    grid::unit(1.22, "null"),
    grid::unit(0.78, "null"),
    grid::unit(0.78, "null"),
    grid::unit(8, "mm"),
    grid::unit(0.58, "null"),
    grid::unit(8, "mm"),
    grid::unit(0.70, "null")
  )) +
  plot_annotation(theme = theme(
    plot.margin = margin(3, 3, 3, 3),
    plot.background = element_rect(fill = "white", colour = NA)
  ))

p_all_txt <- paste(
  sprintf("%s %s", stats_h$lineage[stats_h$scope == "all"],
          vapply(stats_h$wilcox_p[stats_h$scope == "all"], fmtp, "")),
  collapse = "; "
)
p_w_txt <- paste(
  apply(stats_h[stats_h$scope == "windows", ], 1, function(r) {
    sprintf("%s %s %s", r[["lineage"]], r[["set"]], fmtp(as.numeric(r[["wilcox_p"]])))
  }),
  collapse = "; "
)

stats <- c(
  "Figure S9 statistics / integrity (FigS9_Nature_v2)",
  sprintf("Immune Visium spots: %d.", n_immune),
  sprintf("Windows: L1-L3 n = %d; L4-L6 n = %d (same rects as Fig. 8).", n_l13, n_l46),
  sprintf("Panel a heatmap display cells: %d -> %d (<=60 spots/cluster)", n_a_before, n_a_after),
  "Panels c-g: POS/NEG titles; pt.size.factor = 1.00; Mg maps omitted.",
  sprintf("Panel h all-spot Wilcoxon: %s", p_all_txt),
  sprintf("Panel h window Wilcoxon: %s", p_w_txt),
  "Image integrity: H&E global alpha 0.9; no contrast stretch.",
  sprintf("Figure size: %d x %d mm", FIG_W_MM, FIG_H_MM)
)
writeLines(stats, file.path(OUTDIR, "FigureS9_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. S9 | Spatial cluster markers and POS versus NEG immune-lineage mapping in the EAE Visium section.",
  "",
  "a, Heatmap of top marker genes per Visium cluster (0-8); expression scaled and clipped to [-2.5, 2.5].",
  "b, Spatial maps of Mki67, Lgals3, Ccl2 and Hexb on Visium H&E (image alpha = 0.9).",
  "c-g, Prediction-score maps for DC (c), B cells (d), Neu1 (e), plasma cells (f) and macrophages (g); left, POS; right, NEG; image alpha = 0.9.",
  paste0(
    "h, Prediction scores for the same lineages. Top, all immune spots (n = ",
    format(n_immune, big.mark = ","),
    "). Bottom, spots inside perivascular windows, split as L1\u2013L3 (n = ",
    format(n_l13, big.mark = ","),
    ") versus L4\u2013L6 (n = ",
    format(n_l46, big.mark = ","),
    "). Boxes show median and IQR; whiskers 1.5\u00d7IQR. Two-sided Wilcoxon rank-sum tests of POS versus NEG within each set. Window sets are tested separately because L4\u2013L6 need not repeat the L1\u2013L3 direction."
  ),
  "",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS9_legend.txt"))

writeLines(c(
  "Figure S9 | new folder FigS9_Nature_v2 (does not overwrite FigS9_Nature_v1).",
  "Reuse: structural adaptation of v1; window stats follow Fig. 8 f-g.",
  "Dropped: Mg spatial maps and Mg boxplots.",
  "c-g titles relabelled IV/IVC -> POS/NEG.",
  "Panel h (former i): all immune spots and L1-L3 vs L4-L6.",
  "Backend: R exclusive."
), file.path(OUTDIR, "FigureS9_provenance.txt"))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
base <- file.path(OUTDIR, "FigureS9")

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
  "Build notes (full visual QA is in FigureS9_QA.txt)",
  "- Backend: R exclusive",
  "- c-g: POS/NEG titles; no Mg maps",
  "- h: All immune spots and L1-L3 vs L4-L6; no Mg",
  sprintf("- Exports: PDF, SVG, PNG, TIFF at %d x %d mm, raster 600 dpi",
          FIG_W_MM, FIG_H_MM)
), file.path(OUTDIR, "FigureS9_build_notes.txt"))

cat("DONE\n")
