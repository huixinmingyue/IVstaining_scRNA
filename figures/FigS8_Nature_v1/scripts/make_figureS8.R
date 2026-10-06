#!/usr/bin/env Rscript
# =============================================================================
# Fig. S8 v1 — Nature SI multipanel (R / ggplot2 exclusive)
#
# Core conclusion
#   EAE microglia traverse a branched Monocle trajectory in which
#   vessel-contacting cells occupy States 1–2; cluster-wise POS versus NEG
#   expression and a 50-gene pseudotime heatmap show a shift from
#   homeostatic to inflammatory and proliferative programmes.
#
# Evidence chain
#   a  DDRTree coloured by State
#   b  DDRTree coloured by CNS_IV / CNS_IVC
#   c  Six-gene violins by cluster and EAE_POS / EAE_NEG
#   d  Pseudotime-ordered gene-programme heatmap
#
# Archetype : asymmetric mixed-modality (hero = d)
# Size      : 183 x 270 mm
# New folder — does not overwrite FigS6 or Fig6 packages
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS8_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183  # fig_width_mm = 183
FIG_H_MM <- 418
FAM      <- "Arial"
BASE     <- 6.5
PT5      <- 5.2
PT       <- function(x) x / ggplot2::.pt
LABEL_PT <- PT(PT5)

state_cols <- c("1" = "#E64B35", "2" = "#00A087", "3" = "#3C5488")
group_cols <- c("CNS_IV" = "#E64B35", "CNS_IVC" = "#00A087")
grp_cols   <- c("EAE_POS" = "#E64B35", "EAE_NEG" = "#00A087")
prog_cols  <- c(
  homeostatic   = "#4C72B0",
  inflammatory  = "#C44E52",
  proliferative = "#DD8452"
)
gene_lv <- c("P2ry12", "Hexb", "Il1b", "Cxcl2", "Cxcl10", "Id2")

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
    axis.ticks.length = unit(1.0, "pt"),
    axis.title = element_text(size = BASE - 0.2, colour = "black"),
    axis.text = element_text(size = BASE - 0.6, colour = "black"),
    legend.title = element_text(size = BASE - 0.4),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 5.2, hjust = 0.5, colour = "grey30"),
    plot.tag = element_text(size = 8, face = "bold", family = FAM,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(6, 4, 2, 8),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = BASE - 0.6, face = "bold.italic", family = FAM)
  )
theme_set(theme_pub)
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = FAM))

# -----------------------------------------------------------------------------
# a, b  Monocle DDRTree
# -----------------------------------------------------------------------------
cat("Panels a-b ...\n")
traj <- readRDS(file.path(SRC, "traj_slim.rds"))
cells <- traj$cells
stopifnot(nrow(cells) == 1807L)
n_iv  <- sum(cells$group == "CNS_IV")
n_ivc <- sum(cells$group == "CNS_IVC")
cat("  n cells", nrow(cells), " CNS_IV", n_iv, " CNS_IVC", n_ivc, "\n")
print(table(cells$State, cells$group))

traj_plot <- function(color_col, pal, legend_title) {
  df <- cells
  df$.col <- factor(as.character(df[[color_col]]), levels = names(pal))
  ggplot(df, aes(comp1, comp2, colour = .col)) +
    geom_segment(
      data = traj$edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, linewidth = 0.35, colour = "grey25"
    ) +
    geom_point(size = 0.55, stroke = 0, alpha = 0.82) +
    geom_point(
      data = traj$branch_points,
      aes(comp1, comp2), inherit.aes = FALSE,
      shape = 21, size = 2.3, fill = "black", colour = "white", stroke = 0.35
    ) +
    geom_text(
      data = traj$branch_points,
      aes(comp1, comp2, label = label), inherit.aes = FALSE,
      colour = "white", size = PT(5.2), fontface = "bold", family = FAM
    ) +
    scale_colour_manual(values = pal, name = legend_title, drop = FALSE) +
    labs(x = "Component 1", y = "Component 2") +
    coord_equal(clip = "off") +
    theme(
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-2, 0, 0, 0),
      axis.text = element_text(size = 5.2),
      axis.title = element_text(size = 5.8),
      plot.title = element_text(size = 8, face = "bold", hjust = 0, family = FAM,
                                margin = margin(0, 0, 1, 0)),
      plot.title.position = "plot",
      plot.margin = margin(2, 6, 2, 4)
    )
}

pA <- traj_plot("State", state_cols, "State") + labs(title = "a")
pB <- traj_plot("group", group_cols, NULL) + labs(title = "b")

# -----------------------------------------------------------------------------
# c  Violins
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
df_f <- read.csv(file.path(SRC, "Fig6f_violin.csv"), stringsAsFactors = FALSE)
stat_f <- read.csv(file.path(SRC, "Fig6f_wilcoxon_BH.csv"), stringsAsFactors = FALSE)
df_f$gene            <- factor(df_f$gene, levels = gene_lv)
df_f$seurat_clusters <- factor(as.character(df_f$seurat_clusters), levels = as.character(0:3))
df_f$group_label     <- factor(df_f$group_label, levels = c("EAE_POS", "EAE_NEG"))
stopifnot(!anyNA(df_f$gene), !anyNA(df_f$group_label), !anyNA(df_f$seurat_clusters))

n_cells_v <- n_distinct(df_f$cell)
stopifnot(n_cells_v == 1807L)

stat_f$gene            <- factor(stat_f$gene, levels = gene_lv)
stat_f$seurat_clusters <- factor(as.character(stat_f$seurat_clusters),
                                 levels = as.character(0:3))
facet_max <- df_f %>% group_by(gene) %>% summarise(gmax = max(expression), .groups = "drop")
stat_f <- stat_f %>%
  left_join(facet_max, by = "gene") %>%
  mutate(
    star = ifelse(is.na(star), "", star),
    y = gmax * 1.10
  )
stat_star <- stat_f %>% filter(nzchar(star))
cat("  Wilcoxon BH q < 0.05:", sum(stat_f$q < 0.05, na.rm = TRUE), "/ 24\n")
cat("  cluster 3 EAE_POS n =", unique(stat_f$n_pos[stat_f$seurat_clusters == "3"]), "\n")

pC <- ggplot(df_f, aes(seurat_clusters, expression, fill = group_label)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA,
              position = position_dodge(width = 0.78), width = 0.72) +
  stat_summary(fun = median, geom = "point", size = 0.45, colour = "white",
               position = position_dodge(width = 0.78)) +
  geom_text(data = stat_star, aes(seurat_clusters, y, label = star),
            inherit.aes = FALSE, size = LABEL_PT, family = FAM, vjust = 0.15) +
  scale_fill_manual(
    values = grp_cols, name = NULL,
    labels = c("EAE_POS (CNS_IV)", "EAE_NEG (CNS_IVC)")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  labs(x = "Cluster", y = "Normalized expression", title = "c") +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.margin = margin(0, 0, -2, 0),
    axis.text = element_text(size = 5.2),
    axis.title = element_text(size = 5.8),
    plot.title = element_text(size = 8, face = "bold", hjust = 0, family = FAM,
                              margin = margin(0, 0, 1, 0)),
    plot.title.position = "plot",
    plot.margin = margin(2, 6, 4, 4),
    strip.text = element_text(size = 5.8, face = "bold.italic", family = FAM,
                              margin = margin(1, 0, 2, 0))
  )

# -----------------------------------------------------------------------------
# d  Heatmap
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
hm    <- readRDS(file.path(SRC, "heatmap_slim.rds"))
mat   <- hm$mat
gmeta <- hm$gene_meta
stopifnot(nrow(mat) == 50L, ncol(mat) == 8L, nrow(gmeta) == 50L)
stopifnot(identical(rownames(mat), as.character(gmeta$gene)))

gene_levels <- rev(as.character(gmeta$gene))
gmeta$gene <- factor(as.character(gmeta$gene), levels = gene_levels)
df_hm <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
colnames(df_hm) <- c("gene", "bin", "z")
df_hm$gene <- factor(as.character(df_hm$gene), levels = gene_levels)
df_hm$x <- match(df_hm$bin, colnames(mat))
side <- data.frame(
  gene  = factor(as.character(gmeta$gene), levels = gene_levels),
  group = factor(gmeta$group, levels = names(prog_cols)),
  x     = 0.18
)
vmax <- max(0.5, as.numeric(quantile(abs(df_hm$z), 0.98, names = FALSE)))

pD <- ggplot() +
  geom_tile(data = side, aes(x = x, y = gene, fill = group),
            width = 0.42, height = 1, colour = NA) +
  scale_fill_manual(
    values = prog_cols, name = "Programme",
    breaks = names(prog_cols),
    labels = c("Homeostatic", "Inflammatory", "Proliferative"),
    guide = guide_legend(order = 1, override.aes = list(colour = NA))
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = df_hm, aes(x = x, y = gene, fill = z),
            width = 1, height = 1, colour = NA) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-vmax, vmax), oob = squish,
    name = "z-score",
    guide = guide_colorbar(
      order = 2,
      barheight = unit(28, "mm"),
      barwidth = unit(2.2, "mm")
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(ncol(mat)),
    labels = colnames(mat),
    expand = expansion(mult = c(0.01, 0.02)),
    limits = c(-0.12, ncol(mat) + 0.5)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Pseudotime (relative bin)", y = NULL, title = "d") +
  theme(
    axis.text.y = element_text(size = 5.2, face = "italic", family = FAM,
                               colour = "black", lineheight = 0.92,
                               margin = margin(t = 0, r = 1.2, b = 0, l = 0)),
    axis.text.x = element_text(size = 5.2, margin = margin(t = 1.2)),
    axis.title.x = element_text(size = 5.8, margin = margin(t = 2)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.2),
    legend.key.height = unit(3.0, "mm"),
    legend.key.width = unit(3.0, "mm"),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3),
    plot.title = element_text(size = 8, face = "bold", hjust = 0, family = FAM,
                              margin = margin(0, 0, 2, 0)),
    plot.title.position = "plot",
    plot.margin = margin(2, 4, 6, 4)
  )

# -----------------------------------------------------------------------------
# Assemble and export
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
fig <- (pA | pB) / pC / pD +
  plot_layout(heights = grid::unit.c(
    grid::unit(68, "mm"),
    grid::unit(110, "mm"),
    grid::unit(194, "mm")
  ), guides = "keep")
# 68 + 110 + 194 = 372

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "FigureS8")

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

write.csv(cells[, c("cell", "comp1", "comp2", "State", "group", "Pseudotime")],
          file.path(PLOTTED, "SourceData_FigS8ab_trajectory.csv"), row.names = FALSE)
write.csv(traj$edges,
          file.path(PLOTTED, "SourceData_FigS8ab_mst_edges.csv"), row.names = FALSE)
write.csv(df_f,
          file.path(PLOTTED, "SourceData_FigS8c_violin.csv"), row.names = FALSE)
write.csv(stat_f[, c("gene", "seurat_clusters", "n_pos", "n_neg", "p", "q", "star")],
          file.path(PLOTTED, "SourceData_FigS8c_wilcoxon_BH.csv"), row.names = FALSE)
write.csv(df_hm,
          file.path(PLOTTED, "SourceData_FigS8d_heatmap.csv"), row.names = FALSE)
write.csv(gmeta,
          file.path(PLOTTED, "SourceData_FigS8d_gene_meta.csv"), row.names = FALSE)

cat("DONE ->", OUTDIR, "\n")
