#!/usr/bin/env Rscript
# =============================================================================
# Figure 10 | Vascular-contact microglia after MCAO: states, territory, modules
#
#   a  UMAP of all sorted immune cells, microglial clusters 0/3/5 highlighted
#   b  Hexb / Siglech / P2ry12 / Mki67 on the same embedding
#   c  microglial cluster composition, MCAO_NEG versus MCAO_POS, per mouse
#   d  marker dot plot, 9 genes x 3 clusters
#   e  volcano, cluster 5 versus cluster 0
#   f  transferred MCAO_POS fraction plus Lgals3 / Mki67 / P2ry12 / Prdx1 maps
#   g  detection rate by POS dominance (Wilson 95% CI, Fisher)
#   h  Spearman rho against the continuous MCAO_POS fraction
#   i  mean z-scored module scores by cluster x contact
#
# New folder. Does not overwrite Fig9 or any earlier package.
# Backend: R only.
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(scattermore)
  library(cowplot)
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
cowplot::set_null_device("agg")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "Fig10_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- page geometry ----------------------------------------------------------
FIG_W_MM <- 183
width_mm = 183
R1 <- 42; R2 <- 49; R3 <- 54; R4 <- 66; PAD_B <- 3
FIG_H_MM <- R1 + R2 + R3 + R4 + PAD_B               # 214 mm
height_mm = 214

FAM  <- "Arial"
BASE <- 6.0
TAG  <- 8
ANN_PT <- 5.2
ann <- ANN_PT / .pt
MATH_PT <- 7.2

if (!FAM %in% systemfonts::system_fonts()$family) {
  stop("Font '", FAM, "' is not registered with systemfonts; refusing to render ",
       "with a substituted font.")
}

TISSUE_GREY <- "#DCDCDC"
OTHER_GREY  <- "#D9D9D9"
seq_cols    <- c("#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#7F0000")
zero_ramp   <- c("#E8E8E8", seq_cols)
CL_COLS     <- c(`0` = "#0072B2", `3` = "#009E73", `5` = "#A6301F")
ARM_COLS    <- c(MCAO_NEG = "#E67E22", MCAO_POS = "#9B59B6")
DOM_COLS    <- c(`NEG-dominant` = "#E67E22", `POS-dominant` = "#9B59B6")
DIV_LOW     <- "#0072B2"; DIV_MID <- "#F5F5F5"; DIV_HIGH <- "#A6301F"

ARM_LAB  <- c(n = "MCAO_NEG", p = "MCAO_POS")
CL_LEVELS <- c("0", "3", "5")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(0.9, "pt"),
    axis.title   = element_text(size = BASE, colour = "black"),
    axis.text    = element_text(size = BASE - 0.6, colour = "black"),
    legend.title = element_text(size = BASE - 0.4, colour = "black"),
    legend.text  = element_text(size = BASE - 0.6, colour = "black"),
    legend.key.size = unit(2.2, "mm"),
    legend.margin   = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    plot.title    = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5, colour = "grey25"),
    plot.margin   = margin(1.2, 1.2, 1.2, 1.2),
    panel.grid    = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = BASE, face = "bold")
  )

theme_map <- theme_void(base_family = FAM) +
  theme(
    plot.title    = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5, colour = "grey25"),
    legend.title  = element_text(size = BASE - 0.4),
    legend.text   = element_text(size = BASE - 0.6),
    plot.margin   = margin(1.2, 1.2, 1.2, 1.2)
  )

cbar_v <- function(title = NULL, height = 9) {
  guides(colour = guide_colourbar(
    title = title, title.position = "top",
    theme = theme(legend.key.width  = unit(1.3, "mm"),
                  legend.key.height = unit(height, "mm")),
    ticks.colour = "white", ticks.linewidth = 0.3, frame.colour = NA))
}

inside_legend <- function(x, y, just_x = x, just_y = y) {
  theme(legend.position = "inside",
        legend.position.inside = c(x, y),
        legend.justification.inside = c(just_x, just_y))
}

cap_at <- function(v, q = 0.99) {
  cp <- as.numeric(stats::quantile(v[is.finite(v)], q, na.rm = TRUE))
  if (!is.finite(cp) || cp <= 0) cp <- max(v, na.rm = TRUE)
  cp
}

wilson_ci <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  d <- 1 + z^2 / n
  c(centre = p,
    lo = max(0, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d),
    hi = min(1, (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d))
}

rank_biserial <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  U <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))$statistic
  as.numeric(2 * U / (n1 * n2) - 1)
}

fmt_q <- function(q) {
  ifelse(q < 1e-300, "q < 1e-300",
         ifelse(q < 0.001, sprintf("q = %.0e", q), sprintf("q = %.3f", q)))
}

cat("== Figure 10 ==\n")

# =============================================================================
# a | UMAP of all sorted immune cells
# =============================================================================
cat("a | UMAP\n")
umap <- readRDS(file.path(SRC, "panel_a_umap.rds"))
umap$cluster <- factor(as.character(umap$cluster), levels = as.character(0:9))
umap$arm <- unname(ARM_LAB[substr(umap$orig.ident, 1, 1)])
stopifnot(!anyNA(umap$arm))
N_CELLS <- nrow(umap)

umap$state <- ifelse(as.character(umap$cluster) %in% CL_LEVELS,
                     as.character(umap$cluster), "other")
umap$state <- factor(umap$state, levels = c(CL_LEVELS, "other"))
umap_ord <- umap[order(umap$state == "other", decreasing = TRUE), ]

lab_a <- umap %>%
  group_by(cluster) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2),
            mg = first(as.character(cluster)) %in% CL_LEVELS, .groups = "drop")

state_cols <- c(CL_COLS, other = OTHER_GREY)

pA <- ggplot(umap_ord, aes(umap_1, umap_2, colour = state)) +
  geom_scattermore(pointsize = 1.6, pixels = c(1800, 1800), alpha = 0.9) +
  geom_text(data = lab_a, aes(umap_1, umap_2, label = cluster,
                              colour = NULL, fontface = ifelse(mg, "bold", "plain")),
            inherit.aes = FALSE, size = 2.05, family = FAM,
            colour = ifelse(lab_a$mg, "black", "grey35")) +
  scale_colour_manual(values = state_cols, breaks = CL_LEVELS,
                      labels = paste("Cluster", CL_LEVELS),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.1, alpha = 1))) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "Sorted brain immune cells",
       subtitle = sprintf("%s cells, 2 mice per arm; microglial clusters in colour",
                          comma(N_CELLS))) +
  theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.height = unit(2.4, "mm")) +
  inside_legend(0.99, 0.02)

# =============================================================================
# b | marker feature plots
# =============================================================================
cat("b | feature plots\n")
feat <- readRDS(file.path(SRC, "panel_b_features.rds"))
GENES_B <- c("Hexb", "Siglech", "P2ry12", "Mki67")
stopifnot(all(GENES_B %in% names(feat)))

panels_b <- lapply(GENES_B, function(g) {
  d <- feat[order(feat[[g]]), , drop = FALSE]
  cp <- cap_at(d[[g]], 0.99)
  ggplot(d, aes(umap_1, umap_2, colour = pmin(.data[[g]], cp))) +
    geom_scattermore(pointsize = 1.4, pixels = c(1300, 1300)) +
    scale_colour_gradientn(colours = zero_ramp, limits = c(0, cp),
                           oob = squish, name = NULL,
                           breaks = c(0, floor(cp * 10) / 10)) +
    coord_equal() +
    labs(title = g, x = NULL, y = NULL) +
    cbar_v(height = 6) +
    theme_map +
    theme(plot.title = element_text(face = "bold.italic", size = BASE),
          legend.position = "right",
          legend.box.spacing = unit(0.3, "mm"),
          legend.margin = margin(0, 0, 0, 0))
})
pB <- wrap_plots(panels_b, nrow = 1) +
  plot_annotation(
    title = "Microglial identity and proliferation markers",
    theme = theme(plot.title = element_text(size = BASE, face = "bold",
                                            hjust = 0.5, family = FAM),
                  plot.margin = margin(1, 1, 1, 1)))

# =============================================================================
# c | microglial composition per sort arm, both mice shown
# =============================================================================
cat("c | composition\n")
mg <- umap %>% filter(as.character(cluster) %in% CL_LEVELS)
mg$cluster <- factor(as.character(mg$cluster), levels = CL_LEVELS)

per_mouse <- mg %>%
  count(orig.ident, arm, cluster, name = "cells") %>%
  group_by(orig.ident) %>%
  mutate(pct = 100 * cells / sum(cells), mouse_total = sum(cells)) %>%
  ungroup()

arm_summary <- per_mouse %>%
  group_by(arm, cluster) %>%
  summarise(pct_mean = mean(pct), n_mice = n(), .groups = "drop")

arm_cells <- mg %>% count(arm, name = "cells")
N_MG <- nrow(mg)

write.csv(per_mouse, file.path(PLOTTED, "SourceData_Fig10c_microglial_composition_per_mouse.csv"),
          row.names = FALSE)

xpos <- function(cl, arm) match(cl, CL_LEVELS) + ifelse(arm == "MCAO_NEG", -0.19, 0.19)
arm_summary$x   <- xpos(as.character(arm_summary$cluster), arm_summary$arm)
per_mouse$x     <- xpos(as.character(per_mouse$cluster), per_mouse$arm)
cl_tick <- data.frame(cluster = CL_LEVELS, x = seq_along(CL_LEVELS))

pC <- ggplot() +
  geom_col(data = arm_summary, aes(x, pct_mean, fill = arm),
           width = 0.34, colour = "black", linewidth = 0.22) +
  geom_point(data = per_mouse, aes(x, pct), shape = 21, size = 0.8, stroke = 0.22,
             fill = "white", colour = "black",
             position = position_jitter(width = 0.07, height = 0, seed = 42)) +
  geom_segment(data = cl_tick, aes(x = x - 0.37, xend = x + 0.37, y = -4.5, yend = -4.5),
               colour = CL_COLS[CL_LEVELS], linewidth = 1.1, lineend = "butt") +
  scale_fill_manual(values = ARM_COLS, name = NULL) +
  scale_x_continuous(breaks = cl_tick$x, labels = paste("Cluster", CL_LEVELS),
                     expand = expansion(add = 0.45)) +
  scale_y_continuous(limits = c(-6, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Microglia in the arm (%)",
       title = "Microglial state composition",
       subtitle = sprintf("%s microglia; points, individual mice (2 per arm)",
                          comma(N_MG))) +
  theme_pub +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.height = unit(2.2, "mm"), legend.key.width = unit(2.6, "mm"),
        legend.box.spacing = unit(0.4, "mm"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank())

# =============================================================================
# d | marker dot plot
# =============================================================================
cat("d | dot plot\n")
vln <- readRDS(file.path(SRC, "panel_f_violin.rds"))
GENE_ORDER <- c("Hexb", "P2ry12", "Klf4", "Mki67",
                "Prdx1", "Fth1", "Spp1", "Lgals3", "Cxcl2")
stopifnot(all(GENE_ORDER %in% unique(vln$gene)))

dot <- vln %>%
  filter(gene %in% GENE_ORDER, as.character(cluster) %in% CL_LEVELS) %>%
  group_by(gene, cluster) %>%
  summarise(n_cells = n(), pct_expr = 100 * mean(expression > 0),
            mean_expr = mean(expression), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(z = as.numeric(scale(mean_expr))) %>%
  ungroup()
dot$gene    <- factor(dot$gene, levels = rev(GENE_ORDER))
dot$cluster <- factor(as.character(dot$cluster), levels = CL_LEVELS)

write.csv(dot, file.path(PLOTTED, "SourceData_Fig10d_marker_summary.csv"), row.names = FALSE)

pD <- ggplot(dot, aes(cluster, gene)) +
  geom_point(aes(size = pct_expr, fill = z), shape = 21, colour = "grey35",
             stroke = 0.18) +
  scale_size_area(max_size = 2.8, name = "Cells expressing (%)",
                  breaks = c(5, 50, 100), limits = c(0, 100)) +
  scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                       limits = c(-1.2, 1.2), oob = squish,
                       name = "Mean expression (z across clusters)",
                       breaks = c(-1, 0, 1),
                       guide = guide_colourbar(
                         title.position = "top",
                         theme = theme(legend.key.width  = unit(14, "mm"),
                                       legend.key.height = unit(1.3, "mm")),
                         ticks.colour = "white", ticks.linewidth = 0.3)) +
  scale_x_discrete(labels = paste("C", CL_LEVELS, sep = "")) +
  labs(x = "Microglial cluster", y = NULL, title = "State-defining markers",
       subtitle = "log-normalised expression") +
  theme_pub +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.direction = "horizontal",
        legend.spacing.y = unit(0.4, "mm"),
        legend.box.spacing = unit(0.8, "mm"),
        legend.key.height = unit(2.1, "mm"),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.18),
        axis.line = element_line(linewidth = 0.3))

# =============================================================================
# e | volcano, cluster 5 versus cluster 0
# =============================================================================
cat("e | volcano\n")
deg <- read.csv(file.path(SRC, "panel_e_DEG_C5_vs_C0.csv"), stringsAsFactors = FALSE)
P_FLOOR  <- min(deg$p_val_adj[deg$p_val_adj > 0], na.rm = TRUE)
P_CAP    <- -log10(P_FLOOR)
deg$censored  <- deg$p_val_adj <= 0
deg$neglog10q <- ifelse(deg$censored, P_CAP, -log10(deg$p_val_adj))
LFC_T <- 1; Q_T <- 0.05
deg$call <- with(deg, ifelse(p_val_adj < Q_T & avg_log2FC >  LFC_T, "Up in cluster 5",
                      ifelse(p_val_adj < Q_T & avg_log2FC < -LFC_T, "Up in cluster 0",
                             "Not significant")))
deg$call <- factor(deg$call,
                   levels = c("Up in cluster 0", "Not significant", "Up in cluster 5"))
N_UP   <- sum(deg$call == "Up in cluster 5")
N_DOWN <- sum(deg$call == "Up in cluster 0")
N_CENS <- sum(deg$censored)

write.csv(deg[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj",
                  "censored", "call")],
          file.path(PLOTTED, "SourceData_Fig10e_DEG_cluster5_vs_cluster0.csv"),
          row.names = FALSE)

lab_top <- bind_rows(
  deg %>% filter(call == "Up in cluster 5") %>% arrange(desc(avg_log2FC)) %>% head(6),
  deg %>% filter(call == "Up in cluster 0") %>% arrange(avg_log2FC) %>% head(4),
  deg %>% filter(gene %in% GENE_ORDER, call != "Not significant")
) %>% distinct(gene, .keep_all = TRUE)

vol_cols <- c(`Up in cluster 0` = CL_COLS[["0"]], `Not significant` = "grey78",
              `Up in cluster 5` = CL_COLS[["5"]])

pE <- ggplot(deg, aes(avg_log2FC, neglog10q, colour = call)) +
  geom_vline(xintercept = c(-LFC_T, LFC_T), linetype = "22", linewidth = 0.2,
             colour = "grey55") +
  geom_hline(yintercept = -log10(Q_T), linetype = "22", linewidth = 0.2,
             colour = "grey55") +
  geom_point(data = ~ subset(.x, !censored), size = 0.32, stroke = 0, shape = 16,
             alpha = 0.65) +
  geom_point(data = ~ subset(.x, censored), size = 0.55, stroke = 0.2, shape = 24,
             fill = NA) +
  geom_text_repel(data = lab_top, aes(label = gene), size = ann, family = FAM,
                  fontface = "italic", colour = "black", max.overlaps = 60,
                  box.padding = 0.14, point.padding = 0.05, segment.size = 0.15,
                  segment.colour = "grey50", min.segment.length = 0, seed = 42,
                  ylim = c(NA, P_CAP * 1.30), force = 2.5, max.iter = 20000) +
  scale_colour_manual(values = vol_cols, name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.1, alpha = 1,
                                                               shape = 16))) +
  scale_y_continuous(limits = c(0, P_CAP * 1.33), breaks = c(0, 100, 200, 300),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = expression(log[2]~"fold change (cluster 5 / cluster 0)"),
       y = expression(-log[10]~"adjusted "*italic(P)),
       title = "Cluster 5 versus cluster 0",
       subtitle = sprintf(
         "%s genes tested; %s up in cluster 5, %s up in cluster 0\n\u25b3, %d genes at the numeric floor of the test",
         comma(nrow(deg)), N_UP, N_DOWN, N_CENS)) +
  theme_pub +
  theme(legend.direction = "horizontal",
        legend.key.width = unit(1.6, "mm"),
        axis.title = element_text(size = MATH_PT),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA)) +
  inside_legend(0.5, 0.02)

# =============================================================================
# f | MCAO_POS fraction + state-marker maps on the registered section
# =============================================================================
cat("f | spatial maps\n")
spat <- readRDS(file.path(SRC, "panel_d_spatial_pred.rds"))
gspat <- readRDS(file.path(SRC, "panel_g_spatial_genes.rds"))

FEATS <- c("n-0", "p-0", "n-3", "p-3", "n-5", "p-5")
stopifnot(all(FEATS %in% names(spat)))
wsum <- rowSums(spat[, FEATS])
stopifnot(max(abs(wsum - 1)) < 1e-6)

spat$pos_frac <- rowSums(spat[, c("p-0", "p-3", "p-5")])
spat$y_plot   <- -spat$y
N_SPOT <- nrow(spat)

PADPX <- 12
XLIM <- c(min(spat$x) - PADPX, max(spat$x) + PADPX)
YLIM <- c(min(spat$y_plot) - PADPX, max(spat$y_plot) + PADPX)

map_base <- function() {
  list(coord_equal(xlim = XLIM, ylim = YLIM, expand = FALSE),
       theme_map)
}

sil <- spat[, c("x", "y_plot")]
SPOT_SZ <- 0.38

spot_map <- function(df, value, title, subtitle = NULL, cap = NULL,
                     bar_title = NULL, bar_labels = waiver(), italic = TRUE) {
  d <- df
  d$val <- as.numeric(value)
  if (is.null(cap)) cap <- cap_at(d$val, 0.99)
  d <- d[order(d$val), , drop = FALSE]
  ggplot() +
    geom_point(data = sil, aes(x, y_plot), colour = TISSUE_GREY, size = SPOT_SZ,
               stroke = 0, shape = 16) +
    geom_point(data = subset(d, val > 0), aes(x, y_plot, colour = pmin(val, cap)),
               size = SPOT_SZ, stroke = 0, shape = 16) +
    scale_colour_gradientn(colours = seq_cols, limits = c(0, cap), oob = squish,
                           name = NULL, labels = bar_labels) +
    labs(title = title, subtitle = subtitle) +
    cbar_v(title = bar_title, height = 7) +
    map_base() +
    theme(legend.position = "right", legend.margin = margin(0, 0, 0, -1.2),
          plot.title = element_text(face = if (italic) "bold.italic" else "bold",
                                    size = BASE))
}

SPLIT <- 0.5
spat$dominance <- factor(ifelse(spat$pos_frac >= SPLIT, "POS-dominant", "NEG-dominant"),
                         levels = c("NEG-dominant", "POS-dominant"))
n_dom <- table(spat$dominance)

wt <- data.frame(feature = FEATS, mean_weight = colMeans(spat[, FEATS]),
                 spots_gt0 = sapply(spat[, FEATS], function(v) sum(v > 0)),
                 row.names = NULL)
wt$arm     <- unname(ARM_LAB[substr(wt$feature, 1, 1)])
wt$cluster <- factor(sub("^._", "", sub("-", "_", wt$feature)), levels = CL_LEVELS)
write.csv(wt[, c("feature", "arm", "cluster", "mean_weight", "spots_gt0")],
          file.path(PLOTTED, "SourceData_Fig10_reference_weight_composition.csv"),
          row.names = FALSE)

pF_pos <- ggplot() +
  geom_point(data = spat[order(spat$pos_frac), ],
             aes(x, y_plot, colour = pos_frac), size = SPOT_SZ, stroke = 0,
             shape = 16) +
  scale_colour_gradientn(colours = seq_cols, limits = c(0, 1),
                         breaks = c(0, 0.5, 1), name = NULL) +
  labs(title = "MCAO_POS fraction",
       subtitle = sprintf("%s spots; split at 0.5 for g", comma(N_SPOT))) +
  cbar_v(height = 7) +
  map_base() +
  theme(legend.position = "right", legend.margin = margin(0, 0, 0, -1.2),
        plot.title = element_text(face = "bold", size = BASE))

GENES_F <- c("Lgals3", "Mki67", "P2ry12", "Prdx1")
stopifnot(all(GENES_F %in% names(gspat)))
stopifnot(identical(gspat$cell, spat$cell))
gspat$y_plot <- -gspat$y

g_caps <- sapply(GENES_F, function(g) {
  v <- gspat[[g]]
  max(1, ceiling(stats::quantile(v[v > 0], 0.99)))
})
panels_genes <- lapply(GENES_F, function(g) {
  spot_map(gspat, gspat[[g]], title = g,
           subtitle = sprintf("detected in %.0f%%; max %d",
                              100 * mean(gspat[[g]] > 0), max(gspat[[g]])),
           cap = g_caps[[g]], bar_title = sprintf("counts\n(cap %d)", g_caps[[g]]),
           italic = TRUE)
})
pF <- wrap_plots(c(list(pF_pos), panels_genes), nrow = 1)

# =============================================================================
# g-h | detection and Spearman against the POS fraction
# =============================================================================
cat("g-h | quantification\n")
gspat$dominance <- spat$dominance
gspat$pos_frac  <- spat$pos_frac

det <- lapply(GENES_F, function(g) {
  tab <- table(factor(gspat[[g]] > 0, levels = c(FALSE, TRUE)), gspat$dominance)
  ft  <- stats::fisher.test(tab)
  do.call(rbind, lapply(levels(gspat$dominance), function(dm) {
    k <- sum(gspat[[g]][gspat$dominance == dm] > 0)
    n <- sum(gspat$dominance == dm)
    ci <- wilson_ci(k, n)
    data.frame(gene = g, dominance = dm, k = k, n = n,
               pct = 100 * ci[["centre"]], lo = 100 * ci[["lo"]], hi = 100 * ci[["hi"]],
               fisher_p = ft$p.value, odds_ratio = unname(ft$estimate),
               or_lo = ft$conf.int[1], or_hi = ft$conf.int[2])
  }))
}) %>% bind_rows()
det$fisher_q <- p.adjust(det$fisher_p, method = "BH")

rho_tab <- lapply(GENES_F, function(g) {
  ct <- suppressWarnings(stats::cor.test(gspat[[g]], gspat$pos_frac, method = "spearman"))
  b <- replicate(2000, {
    i <- sample.int(N_SPOT, N_SPOT, replace = TRUE)
    suppressWarnings(stats::cor(gspat[[g]][i], gspat$pos_frac[i], method = "spearman"))
  })
  x <- gspat[[g]][gspat$dominance == "POS-dominant"]
  y <- gspat[[g]][gspat$dominance == "NEG-dominant"]
  ws <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
  data.frame(gene = g, rho = unname(ct$estimate),
             rho_lo = unname(stats::quantile(b, 0.025, na.rm = TRUE)),
             rho_hi = unname(stats::quantile(b, 0.975, na.rm = TRUE)),
             rho_p = ct$p.value,
             ranksum_p = ws$p.value, r_rb = rank_biserial(x, y),
             median_pos = median(x), median_neg = median(y),
             mean_pos = mean(x), mean_neg = mean(y))
}) %>% bind_rows()
rho_tab$rho_q      <- p.adjust(rho_tab$rho_p, method = "BH")
rho_tab$ranksum_q  <- p.adjust(rho_tab$ranksum_p, method = "BH")

write.csv(det, file.path(PLOTTED, "SourceData_Fig10g_detection_by_dominance.csv"),
          row.names = FALSE)
write.csv(rho_tab, file.path(PLOTTED, "SourceData_Fig10h_spearman_and_ranksum.csv"),
          row.names = FALSE)
write.csv(data.frame(cell = spat$cell, x = spat$x, y = spat$y,
                     spat[, FEATS], pos_frac = spat$pos_frac,
                     dominance = as.character(spat$dominance),
                     gspat[, GENES_F], check.names = FALSE),
          file.path(PLOTTED, "SourceData_Fig10f_spot_level.csv"), row.names = FALSE)

det$gene <- factor(det$gene, levels = GENES_F)
det$dominance <- factor(det$dominance, levels = names(DOM_COLS))
q_lab <- det %>% group_by(gene) %>%
  summarise(hi = max(hi), q = first(fisher_q), or = first(odds_ratio), .groups = "drop") %>%
  mutate(
    lab = sprintf("%s\nOR %.2f", fmt_q(q), or),
    y = ifelse(hi >= 85, 103, hi + 3)
  )

pG <- ggplot(det, aes(gene, pct, fill = dominance)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.6,
           colour = "black", linewidth = 0.22) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.72),
                width = 0.16, linewidth = 0.28) +
  geom_text(data = q_lab, aes(gene, y, label = lab), inherit.aes = FALSE,
            size = ann, family = FAM, colour = "grey20", lineheight = 0.95,
            vjust = 0) +
  scale_fill_manual(values = DOM_COLS, name = NULL) +
  scale_y_continuous(limits = c(0, 130), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Spots with \u2265 1 count (%)",
       title = "Detection by POS dominance",
       subtitle = "Wilson 95% CI; Fisher, BH over 4 genes") +
  theme_pub +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.height = unit(2.2, "mm"), legend.key.width = unit(2.6, "mm"),
        legend.box.spacing = unit(0.4, "mm"))

rho_tab$gene <- factor(rho_tab$gene, levels = rev(GENES_F))
rho_tab$lab  <- sprintf("%+.2f\n%s", rho_tab$rho, fmt_q(rho_tab$rho_q))

pH <- ggplot(rho_tab, aes(rho, gene)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = rho_lo, xmax = rho_hi), orientation = "y", width = 0.14,
                linewidth = 0.32, colour = "grey25") +
  geom_point(size = 1.15, colour = "black") +
  geom_text(aes(x = 0.72, label = lab), hjust = 0, vjust = 0.5,
            size = ann, family = FAM, colour = "grey20", lineheight = 0.92) +
  scale_x_continuous(limits = c(-0.28, 1.12), breaks = c(-0.2, 0, 0.2, 0.4, 0.6),
                     expand = expansion(mult = 0.01)) +
  labs(x = expression("Spearman "*rho*" vs POS fraction"), y = NULL,
       title = "Association with POS fraction",
       subtitle = sprintf("all %s spots; bootstrap 95%% CI", comma(N_SPOT))) +
  theme_pub +
  theme(axis.text.y = element_text(face = "italic"))

# =============================================================================
# i | module scores by cluster x contact
# =============================================================================
cat("i | modules\n")
mod <- read.csv(file.path(SRC, "module_cluster_by_contact.csv"), stringsAsFactors = FALSE)
mod <- mod %>%
  filter(set %in% c("DAM", "DAM_core", "ICAM", "SAMC", "IFN", "G2M")) %>%
  mutate(
    cluster = factor(as.character(cluster), CL_LEVELS),
    contact = recode(contact, NEG = "MCAO_NEG", POS = "MCAO_POS"),
    contact = factor(contact, c("MCAO_NEG", "MCAO_POS")),
    set = recode(set, DAM = "DAM-like", DAM_core = "DAM core"),
    set = factor(set, c("DAM-like", "DAM core", "ICAM", "SAMC", "IFN", "G2M")),
    cluster_lab = factor(paste0("C", cluster),
                         levels = paste0("C", CL_LEVELS)),
    lab_col = ifelse(mean > 0.90, "white", "black")
  )
write.csv(mod[, c("set", "cluster", "contact", "n", "mean", "median")],
          file.path(PLOTTED, "SourceData_Fig10i_module_cluster_by_contact.csv"),
          row.names = FALSE)

mod_genes <- data.frame(
  set = c("DAM-like", "DAM core", "ICAM", "SAMC", "IFN", "G2M"),
  n_used = c(8, 3, 10, 9, 9, 8),
  genes = c(
    "Clec7a, Itgax, Cst7, Trem2, Apoe, Lpl, Tyrobp, Axl",
    "Itgax, Cst7, Trem2 (Clec7a not detected)",
    "Lgals3, Plau, Srxn1, Ankrd33b, Gas2l3, Nes, Edn1, Fam129b, Vat1, Fmn1",
    "Spp1, Fabp5, Gpnmb, Lgals3, Lpl, Fth1, Cd63, Ctsb, Mmp12",
    "Ifit1, Ifit2, Ifit3, Isg15, Irf7, Rsad2, Oasl1, Oasl2, Stat1",
    "Mki67, Top2a, Birc5, Cdk1, Pcna, Tyms, Pclaf, Ube2c"
  ),
  stringsAsFactors = FALSE
)
write.csv(mod_genes, file.path(PLOTTED, "SourceData_Fig10i_module_genes.csv"),
          row.names = FALSE)

pI <- ggplot(mod, aes(contact, set, fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", mean), colour = lab_col),
            size = ann, family = FAM, show.legend = FALSE) +
  scale_colour_identity() +
  facet_grid(. ~ cluster_lab) +
  scale_fill_gradient2(
    low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
    limits = c(-0.5, 1.6), oob = squish, name = "Mean z",
    breaks = c(0, 0.5, 1, 1.5),
    guide = guide_colourbar(
      title.position = "top",
      theme = theme(legend.key.width  = unit(1.6, "mm"),
                    legend.key.height = unit(10, "mm")),
      ticks.colour = "white", ticks.linewidth = 0.3)
  ) +
  scale_x_discrete(expand = c(0, 0), labels = c("NEG", "POS")) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Modules by cluster \u00d7 contact",
       subtitle = "Mean z; DAM core = Itgax, Cst7, Trem2") +
  theme_pub +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(size = BASE - 0.4),
    strip.text = element_text(face = "bold", size = BASE - 0.2),
    panel.spacing.x = unit(1.4, "mm"),
    plot.margin = margin(1.2, 1.2, 1.2, 1.2)
  )

# =============================================================================
# assemble
# =============================================================================
cat("assemble\n")
row1 <- plot_grid(pA, pB, nrow = 1, rel_widths = c(0.52, 1.31))
row2 <- plot_grid(pC, pD, pE, nrow = 1, rel_widths = c(0.52, 0.62, 0.69))
row3 <- pF
row4 <- plot_grid(pG, pH, pI, nrow = 1, rel_widths = c(1.08, 0.88, 1.12),
                  align = "h", axis = "tb")

page <- plot_grid(row1, row2, row3, row4, ncol = 1,
                  rel_heights = c(R1, R2, R3, R4 + PAD_B))

tag_at <- function(lab, x_mm, y_mm) {
  draw_plot_label(lab, x = x_mm / FIG_W_MM, y = 1 - y_mm / FIG_H_MM,
                  size = TAG, fontface = "bold", family = FAM,
                  hjust = 0, vjust = 1)
}

fig <- ggdraw(page) +
  tag_at("a", 0.8, 0.8) +
  tag_at("b", 52.0, 0.8) +
  tag_at("c", 0.8, R1 + 0.8) +
  tag_at("d", 52.0, R1 + 0.8) +
  tag_at("e", 114.0, R1 + 0.8) +
  tag_at("f", 0.8, R1 + R2 + 0.8) +
  tag_at("g", 0.8, R1 + R2 + R3 + 0.8) +
  tag_at("h", 65.0, R1 + R2 + R3 + 0.8) +
  tag_at("i", 117.0, R1 + R2 + R3 + 0.8) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# =============================================================================
# export
# =============================================================================
save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  cat("  SVG\n")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); grDevices::dev.off()
  cat("  PDF\n")
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = FAM)
  print(plot); grDevices::dev.off()
  cat("  PNG\n")
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in",
                res = dpi, background = "white")
  print(plot); grDevices::dev.off()
  cat("  TIFF\n")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = dpi, compression = "lzw", background = "white")
  print(plot); grDevices::dev.off()
}

save_pub(fig, file.path(OUTDIR, "Figure10"), FIG_W_MM, FIG_H_MM)

saveRDS(list(per_mouse = per_mouse, arm_summary = arm_summary, arm_cells = arm_cells,
             dot = dot, det = det, rho_tab = rho_tab, wt = wt, mod = mod,
             n_cells = N_CELLS, n_mg = N_MG, n_spot = N_SPOT, n_dom = n_dom,
             p_floor = P_FLOOR, p_cap = P_CAP, n_up = N_UP, n_down = N_DOWN,
             n_cens = N_CENS, g_caps = g_caps,
             fig_w = FIG_W_MM, fig_h = FIG_H_MM),
        file.path(QA, "figure10_numbers.rds"))

cat("DONE ->", OUTDIR, "\n")
cat(sprintf("final size %.0f x %.0f mm\n", FIG_W_MM, FIG_H_MM))
