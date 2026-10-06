#!/usr/bin/env Rscript
# =============================================================================
# Fig. 2 v5 - Nature-style multipanel figure (R exclusive)
#
# Core conclusion
#   Intravascular-enriched B-cell subcluster 1 is MZ-like: it is Cr2-high /
#   Fcer2a-low, the S1a MZ∩Intra set maps to ImmGen MZ versus Fo, and MZ
#   module scores are highest in subcluster 1 and higher in IV+ cells.
#
# Evidence chain
#   a  Do IV- and IV+ B cells occupy the same space, and does cluster
#      composition shift?                                            (UMAP + %)
#   b  Is subcluster 1 Cr2-high / Fcer2a-low?                       (feature+violin)
#   c  Do independent MZ∩Intra genes peak at ImmGen MZ?             (heatmap)
#   d  Are MZ module scores highest in subcluster 1 and in IV+?     (violins)
#
# Archetype : quantitative grid (panel c S1a mapping as hero)
# Backend   : R (ggplot2 / patchwork / svglite / ragg)
# Size      : 183 x 196 mm
# New folder — does not overwrite Fig2_Nature_v4
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

OUTDIR <- file.path(IV_FIG_ROOT, "Fig2_Nature_v5")
SRC    <- file.path(OUTDIR, "source_data")
QA     <- file.path(OUTDIR, "qa")
dir.create(QA, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 178
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

CL_COL <- c("0" = "#4C6E8A", "1" = "#D1495B", "2" = "#E9A03B")
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"
DIV_LO <- "#3B6FA0"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#C0392B"
SUBSET_ORDER <- c("T1", "T2", "T3", "Fo", "MZ", "Mem", "GC-CB", "GC-CC")

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
      legend.key.size = unit(2.5, "mm"),
      legend.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = BASE + 0.3, face = "bold", hjust = 0.5),
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
seal <- function(x) wrap_elements(full = x)

# -----------------------------------------------------------------------------
# Panel a : B-cell UMAP + subcluster composition (Fig. 1e style)
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
cells <- readRDS(file.path(SRC, "Fig2ab_cells.rds"))
n_total <- nrow(cells)
stopifnot(n_total == 5975)
n_iv <- table(cells$IV_group)
lab_pos <- cells %>%
  group_by(cluster) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop")
xr <- range(cells$UMAP1); yr <- range(cells$UMAP2)
pad <- 0.05 * diff(xr)
iv_n <- as.integer(n_iv[names(IV_LAB)])
facet_lab <- setNames(
  sprintf("%s\nn = %s", unname(IV_LAB), format(iv_n, big.mark = ",")),
  names(IV_LAB)
)
cells$facet <- factor(facet_lab[as.character(cells$IV_group)],
                      levels = unname(facet_lab))

pA_umap <- ggplot(cells, aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 0.18, alpha = 0.72, shape = 16) +
  geom_text(data = lab_pos, aes(x, y, label = cluster), inherit.aes = FALSE,
            family = FAM, fontface = "bold", size = PT(6.4), colour = "black") +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = CL_COL, guide = "none") +
  coord_fixed(xlim = c(xr[1] - pad, xr[2] + pad),
              ylim = c(yr[1] - pad, yr[2] + pad), expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
        axis.line = element_blank(),
        strip.text = element_text(size = BASE - 0.2, face = "bold", lineheight = 0.95),
        panel.spacing = unit(1.4, "mm"),
        plot.margin = margin(0.4, 1, 0.2, 1, "mm"))

prop <- cells %>%
  count(IV_group, cluster, name = "n") %>%
  group_by(IV_group) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(
    cluster = factor(as.character(cluster), levels = rev(c("0", "1", "2"))),
    IV_lab = factor(unname(IV_LAB[as.character(IV_group)]), levels = unname(IV_LAB))
  )
stopifnot(sum(prop$n) == n_total)

pA_bar <- ggplot(prop, aes(pct, cluster, fill = cluster)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.12,
            size = 1.85, family = FAM, colour = "black") +
  facet_wrap(~ IV_lab, nrow = 1) +
  scale_fill_manual(values = CL_COL, guide = "none") +
  scale_x_continuous(name = "Cells (%)", limits = c(0, 100),
                     breaks = c(0, 40, 80),
                     expand = expansion(mult = c(0, 0.28))) +
  labs(y = NULL) +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = BASE - 0.4, face = "bold",
                               margin = margin(r = 1.8, unit = "mm")),
    strip.text = element_text(size = BASE - 0.8, face = "bold"),
    panel.spacing = unit(1.2, "mm"),
    plot.margin = margin(1, 1.5, 1, 0.8, "mm")
  )

# -----------------------------------------------------------------------------
# Panel b : Cr2 and Fcer2a (former Fig. S2 a-b)
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
feat_vln <- function(df, gene, y_lab = "Expression") {
  zmax <- max(as.numeric(quantile(df$expr, 0.99, names = FALSE)), 0.15)
  peak <- names(which.max(tapply(df$expr, df$cluster, mean)))
  wt <- suppressWarnings(wilcox.test(df$expr[df$cluster == peak],
                                     df$expr[df$cluster != peak]))
  p_f <- ggplot(df, aes(UMAP1, UMAP2, colour = expr)) +
    geom_point(size = 0.13, alpha = 0.78, shape = 16) +
    scale_colour_gradientn(
      colours = c(SEQ_LO, SEQ_HI), limits = c(0, zmax), oob = squish,
      name = NULL, breaks = c(0, zmax),
      labels = c("0", sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = unit(2.0, "mm"),
                             barheight = unit(12, "mm"),
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
          plot.margin = margin(1, 0.5, 0.5, 0.5, "mm"))
  p_v <- ggplot(df, aes(cluster, expr, fill = cluster)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_fill_manual(values = CL_COL, guide = "none") +
    labs(x = "Cluster", y = y_lab,
         subtitle = sprintf("highest in %s", peak)) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5,
                                       colour = "#555555"),
          axis.text.x = element_text(size = BASE - 0.6, face = "bold"),
          plot.margin = margin(1, 1, 0.5, 1, "mm"))
  list(feat = p_f, vln = p_v, peak = peak, peak_p = wt$p.value)
}

d_cr2 <- cells; d_cr2$expr <- cells$Cr2
d_fc  <- cells; d_fc$expr  <- cells$Fcer2a
g_cr2 <- feat_vln(d_cr2, "Cr2")
g_fc  <- feat_vln(d_fc, "Fcer2a")
pB <- wrap_elements(full = (
  (g_cr2$feat | g_cr2$vln) / (g_fc$feat | g_fc$vln) +
    plot_layout(widths = c(1.20, 0.80), heights = c(1, 1))
))

# -----------------------------------------------------------------------------
# Panel c : S1a MZ intersect Intra vs ImmGen (former Fig. 2d / P0-2)
# -----------------------------------------------------------------------------
cat("Heatmap c ...\n")
read_z <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  mat <- as.matrix(d[, SUBSET_ORDER, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- d$gene
  list(mat = mat, memb = d$gene_set)
}
read_delta <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$subset <- factor(d$subset, levels = SUBSET_ORDER)
  d[order(d$subset), ]
}

z_c <- read_z(file.path(SRC, "SourceData_P0-2_S1a_zscore.csv"))
z_c$memb <- factor(
  ifelse(grepl("^MZ", z_c$memb), "MZ", "Fo"),
  levels = c("Fo", "MZ")
)
dlt_c <- read_delta(file.path(SRC, "SourceData_P0-2_S1a_delta.csv"))

order_genes <- function(mat, memb) {
  memb <- droplevels(as.factor(memb))
  unlist(lapply(levels(memb), function(lv) {
    idx <- which(as.character(memb) == lv)
    if (length(idx) < 2) return(rownames(mat)[idx])
    hc <- hclust(dist(mat[idx, , drop = FALSE]), method = "complete")
    rownames(mat)[idx][hc$order]
  }), use.names = FALSE)
}

gg_heatmap <- function(z, delta, show_names = TRUE, zlim = NULL) {
  mat <- z$mat
  if (is.null(zlim)) {
    zlim <- max(2.5, ceiling(max(abs(mat), na.rm = TRUE) * 2) / 2)
  }
  gene_ord <- order_genes(mat, z$memb)
  dd <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(dd) <- c("gene", "subset", "z")
  dd$gene <- factor(as.character(dd$gene), levels = rev(gene_ord))
  dd$subset <- factor(as.character(dd$subset), levels = SUBSET_ORDER)
  dd$memb <- z$memb[match(as.character(dd$gene), rownames(mat))]
  dlt <- delta
  dlt$subset <- factor(as.character(dlt$subset), levels = SUBSET_ORDER)
  dmax <- max(0.5, ceiling(max(abs(dlt$delta_mean_z), na.rm = TRUE) * 10) / 10)
  p_dlt <- ggplot(dlt, aes(subset, delta_mean_z, fill = delta_mean_z > 0)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#888888") +
    geom_col(width = 0.78, colour = NA) +
    scale_fill_manual(values = c("TRUE" = unname(CL_COL[["1"]]),
                                 "FALSE" = unname(CL_COL[["0"]])),
                      guide = "none") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(-dmax, dmax),
                       breaks = c(-dmax, 0, dmax),
                       expand = expansion(mult = c(0.05, 0.08))) +
    labs(x = NULL, y = "\u0394 mean z") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = BASE - 0.8),
          axis.text.y = element_text(size = 5.2),
          plot.margin = margin(0.4, 1, 0, 2.2, "mm"))
  p_hm <- ggplot(dd, aes(subset, gene, fill = z)) +
    geom_tile() +
    facet_grid(memb ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(position = if (show_names) "right" else "left") +
    scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI,
                         midpoint = 0, limits = c(-zlim, zlim),
                         oob = squish, name = "z",
                         breaks = c(-zlim, 0, zlim),
                         guide = guide_colorbar(barwidth = unit(2.2, "mm"),
                                                barheight = unit(16, "mm"),
                                                ticks.colour = "black",
                                                frame.colour = "black",
                                                frame.linewidth = 0.2)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = BASE - 0.5),
          axis.text.y = if (show_names) {
            element_text(face = "italic", size = 5.3)
          } else {
            element_blank()
          },
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(size = 5.6, face = "bold", angle = 90),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0, 1, 1, 2.2, "mm"))
  p_dlt / p_hm + plot_layout(heights = c(0.26, 1.00))
}

as_wrapped <- function(p) wrap_elements(full = patchwork::patchworkGrob(p))

pC <- as_wrapped(gg_heatmap(z_c, dlt_c, show_names = TRUE, zlim = 4))

# -----------------------------------------------------------------------------
# Panel d : ImmGen module scores (former Fig. 2e / P0-4)
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
sc <- read.csv(file.path(SRC, "SourceData_P0-4_module_scores_per_cell.csv"),
               stringsAsFactors = FALSE)
sc <- sc %>%
  filter(module %in% c("MZ", "Fo", "Mem", "GC")) %>%
  mutate(
    module = factor(module, levels = c("MZ", "Fo", "Mem", "GC"),
                    labels = c("MZ signature", "Fo signature",
                               "Mem signature", "GC signature")),
    cluster = factor(as.character(cluster), levels = c("0", "1", "2")),
    IV = recode(IV_group, !!!IV_LAB)
  )
stopifnot(nrow(sc) == 5975 * 4)

pE <- ggplot(sc, aes(cluster, score, fill = IV)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, width = 0.86,
              position = position_dodge(width = 0.85)) +
  geom_boxplot(width = 0.14, outlier.size = 0.12, outlier.alpha = 0.20,
               position = position_dodge(width = 0.85), linewidth = 0.25) +
  facet_wrap(~ module, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("IV\u2212 (extravascular)" = IV_COL[["IntraNeg"]],
                               "IV+ (intravascular)" = IV_COL[["IntraPos"]]),
                    name = NULL) +
  labs(x = "B-cell subcluster", y = "Module score") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        strip.text = element_text(size = BASE - 0.2, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- seal(
  tag_it(pA_umap, "a") +
    (plot_spacer() / pA_bar / plot_spacer() + plot_layout(heights = c(0.38, 1.00, 0.38))) +
    tag_it(pB, "b") +
    plot_layout(widths = c(1.18, 0.95, 1.00))
)
row2 <- seal(tag_it(pC, "c") + tag_it(pE, "d") + plot_layout(widths = c(0.95, 1.15)))

fig <- row1 / row2 +
  plot_layout(heights = c(0.95, 1.12)) &
  theme(plot.margin = margin(0.6, 0.6, 0.6, 0.6, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "Figure2")

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
write.csv(cells[, c("cell", "UMAP1", "UMAP2", "cluster", "IV_group", "Cr2", "Fcer2a")],
          file.path(SD, "SourceData_Fig2ab_cells.csv"), row.names = FALSE)
write.csv(prop[, c("IV_group", "cluster", "n", "total", "pct")],
          file.path(SD, "SourceData_Fig2a_composition.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(z_c$mat), gene_set = as.character(z_c$memb),
                     z_c$mat, check.names = FALSE),
          file.path(SD, "SourceData_Fig2c_S1a_zscore.csv"), row.names = FALSE)
write.csv(dlt_c, file.path(SD, "SourceData_Fig2c_S1a_delta.csv"),
          row.names = FALSE)
write.csv(sc[, c("cell", "cluster", "IV_group", "module", "score")],
          file.path(SD, "SourceData_Fig2d_module_scores.csv"), row.names = FALSE)

cat("\nCr2 peak cluster:", g_cr2$peak, " P vs rest", fmt_p_sci1(g_cr2$peak_p), "\n")
cat("Fcer2a peak cluster:", g_fc$peak, " P vs rest", fmt_p_sci1(g_fc$peak_p), "\n")
cat("Cells:", n_total, " IV-:", n_iv[["IntraNeg"]], " IV+:", n_iv[["IntraPos"]], "\n")
cat("DONE ->", OUTDIR, "\n")
