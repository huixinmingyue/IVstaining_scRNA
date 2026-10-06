#!/usr/bin/env Rscript
# =============================================================================
# Fig. S12 v1 — Nature SI multipanel (FigS12_Nature_v1)
# New folder; does not overwrite FigS9 packages, FigS10, or any prior folder.
#
# Core conclusion
#   Two independent mouse Visium sections recapitulate n-/p- microglia
#   subtype maps and POS-spot enrichment of Lgals3. Prdx1 is not
#   significant on either section.
#
# Evidence chain
#   a  Do cluster 5/3/0 NEG vs POS scores occupy distinct territories
#      on section 1?
#   b  Are Lgals3 / Prdx1 / Mki67 enriched on the top 100 POS spots
#      of section 1?
#   c  Does an independent section repeat the score geography (same
#      0-0.5 scale as a)?
#   d  Do the same genes repeat on section 2?
#
# Archetype : image plate + quant
# Backend   : R exclusive
# Size      : 183 x 245 mm
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
  file.path(IV_FIG_ROOT, "FigS12_Nature_v1")
}
SRC     <- file.path(OUTDIR, "source_data")
PLOTTED <- file.path(OUTDIR, "source_data_plotted")
QA      <- file.path(OUTDIR, "qa")
for (d in c(PLOTTED, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 245
width_mm = 183
height_mm = 245
FAM  <- "Arial"
BASE <- 6.5
PT5  <- 5.2
PT <- function(x) x / ggplot2::.pt
N_TOP <- 100L
LIM_PRED <- c(0, 0.5)
LIM_GENE <- c(0, 2)
PAL_PRED <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
BOX_COLS <- c(NEG = "#D0D0D0", POS = "#4D4D4D")

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

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
)

scale_spot <- function(limits = NULL, colours = PAL_PRED) {
  scale_colour_gradientn(
    colours = colours,
    limits = limits, oob = squish,
    breaks = if (!is.null(limits)) seq(limits[1], limits[2], length.out = 6) else waiver(),
    guide = guide_colourbar(
      barwidth = unit(1.5, "mm"),
      barheight = unit(9, "mm"),
      ticks.colour = "black",
      frame.colour = "black",
      nbin = 100
    )
  )
}

pt_for_section <- function(slim, ref_slim, base_pt = 0.40) {
  area_per <- function(s) {
    xr <- diff(range(s$coords$x, na.rm = TRUE))
    yr <- diff(range(s$coords$y, na.rm = TRUE))
    (xr * yr) / max(nrow(s$coords), 1)
  }
  base_pt * sqrt(area_per(slim) / area_per(ref_slim))
}

theme_spot <- function() {
  theme_void(base_size = BASE, base_family = FAM) +
    theme(
      plot.title = element_text(
        size = BASE - 0.2, face = "bold", hjust = 0.5, family = FAM,
        margin = margin(b = 1, t = 0)
      ),
      plot.subtitle = element_text(
        size = BASE - 0.8, face = "plain", hjust = 0.5, family = FAM,
        colour = "grey30", margin = margin(b = 1, t = 0)
      ),
      legend.text = element_text(size = PT5, family = FAM),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1),
      legend.key.height = unit(6, "mm"),
      legend.key.width = unit(1.4, "mm")
    )
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

section_lab <- function(txt) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = txt, angle = 90,
      size = PT(BASE), fontface = "bold", family = FAM
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void(base_family = FAM) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# -----------------------------------------------------------------------------
# Load slim payloads + H&E
# -----------------------------------------------------------------------------
s2 <- readRDS(file.path(SRC, "slim_immune2.rds"))
s4 <- readRDS(file.path(SRC, "slim_immune4.rds"))
he2 <- readRDS(file.path(SRC, "he_immune2.rds"))
he4 <- readRDS(file.path(SRC, "he_immune4.rds"))
n_s2 <- as.integer(s2$n_spots)
n_s4 <- as.integer(s4$n_spots)
cat("section1 spots:", n_s2, "; section2 spots:", n_s4, "\n")
cat("HE dims: s1", paste(he2$dim, collapse = "x"), "; s2", paste(he4$dim, collapse = "x"), "\n")
PT_S2 <- pt_for_section(s2, s2, base_pt = 0.40)
PT_S4 <- pt_for_section(s4, s2, base_pt = 0.40)
cat(sprintf("spot pt: section1=%.3f; section2=%.3f\n", PT_S2, PT_S4))

fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}

he_annotation <- function(he, alpha = 0.55) {
  arr <- fade_he(he$image, alpha = alpha)
  nr <- dim(arr)[1]
  nc <- dim(arr)[2]
  arr <- arr[rev(seq_len(nr)), , , drop = FALSE]
  list(
    raster = grDevices::as.raster(arr),
    xmin = 0.5,
    xmax = nc + 0.5,
    ymin = 0.5,
    ymax = nr + 0.5
  )
}

spot_plot <- function(slim, values, cells = NULL, title = NULL, subtitle = NULL,
                      limits, pt = 0.45, show_legend = TRUE,
                      colours = PAL_PRED, he = NULL, he_alpha = 0.90) {
  coords <- slim$coords
  if (is.null(cells)) cells <- rownames(coords)
  cells <- intersect(cells, rownames(coords))
  v <- as.numeric(values[cells])
  if (!is.null(limits)) v <- pmin(pmax(v, limits[1]), limits[2])
  names(v) <- cells
  df <- data.frame(
    x = coords[cells, "x"],
    y = coords[cells, "y"],
    val = v
  )
  df <- df[order(df$val), , drop = FALSE]

  p <- ggplot(df, aes(x, y, colour = val))
  if (!is.null(he)) {
    ha <- he_annotation(he, alpha = he_alpha)
    p <- p +
      annotation_raster(
        ha$raster, xmin = ha$xmin, xmax = ha$xmax,
        ymin = ha$ymin, ymax = ha$ymax
      ) +
      geom_point(size = pt, stroke = 0) +
      coord_fixed(
        xlim = c(ha$xmin, ha$xmax),
        ylim = c(ha$ymin, ha$ymax),
        expand = FALSE
      )
  } else {
    p <- p +
      geom_point(size = pt, stroke = 0) +
      coord_fixed()
  }

  p <- p +
    scale_spot(limits = limits, colours = colours) +
    labs(title = title, subtitle = subtitle) +
    theme_spot()
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

top_cells <- function(slim, feat, n = N_TOP) {
  names(sort(slim$pred[feat, ], decreasing = TRUE)[seq_len(n)])
}

feat_neg <- function(cluster) paste0("n-", cluster)
feat_pos <- function(cluster) paste0("p-", cluster)

make_df <- function(slim, gene, cluster, n = N_TOP) {
  n_feat <- feat_neg(cluster)
  p_feat <- feat_pos(cluster)
  cells_n <- top_cells(slim, n_feat, n)
  cells_p <- top_cells(slim, p_feat, n)
  data.frame(
    gene = gene,
    cluster = cluster,
    group = factor(c(rep("NEG", length(cells_n)), rep("POS", length(cells_p))),
                   levels = c("NEG", "POS")),
    expression = c(slim$expr[gene, cells_n], slim$expr[gene, cells_p]),
    barcode = c(cells_n, cells_p),
    n_before = slim$n_spots,
    n_after = n,
    stringsAsFactors = FALSE
  )
}

wilcox_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 0.05) return("*")
  "ns"
}

make_box <- function(df) {
  wt <- wilcox.test(expression ~ group, data = df, exact = FALSE)
  star <- wilcox_stars(wt$p.value)
  y_max <- max(df$expression, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  y_star <- max(y_max * 1.08, y_max + 0.35)
  ggplot(df, aes(group, expression, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(
      width = 0.18, height = 0, size = 0.45, stroke = 0,
      colour = "grey20", alpha = 0.55
    ) +
    scale_fill_manual(values = BOX_COLS, guide = "none") +
    annotate(
      "text", x = 1.5, y = y_star, label = star,
      size = PT(6.5), fontface = "bold", family = FAM
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.22), clip = "off") +
    labs(x = NULL, y = "Expression") +
    theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.title.y = element_text(size = 6.1),
      axis.text = element_text(size = 5.2),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid = element_blank(),
      plot.margin = margin(2, 3, 2, 3)
    )
}

# =============================================================================
# a / c — prediction maps
# =============================================================================
make_np_pair <- function(slim, cluster, pt, show_legend = FALSE) {
  nf <- feat_neg(cluster)
  pf <- feat_pos(cluster)
  (spot_plot(slim, slim$pred[nf, ], title = "NEG", limits = LIM_PRED,
             pt = pt, show_legend = FALSE, colours = PAL_PRED) |
   spot_plot(slim, slim$pred[pf, ], title = "POS", limits = LIM_PRED,
             pt = pt, show_legend = show_legend, colours = PAL_PRED))
}

make_section_pred_row <- function(slim, section_txt, pt) {
  clusters <- wrap_plots(
    (hdr_plot("Cluster 5") / make_np_pair(slim, 5, pt, show_legend = FALSE)) +
      plot_layout(heights = c(0.12, 1)),
    (hdr_plot("Cluster 3") / make_np_pair(slim, 3, pt, show_legend = FALSE)) +
      plot_layout(heights = c(0.12, 1)),
    (hdr_plot("Cluster 0") / make_np_pair(slim, 0, pt, show_legend = TRUE)) +
      plot_layout(heights = c(0.12, 1)),
    nrow = 1
  )
  (section_lab(section_txt) | clusters) + plot_layout(widths = c(0.055, 1))
}

pA <- make_section_pred_row(s2, "Section 1", PT_S2)
pC <- make_section_pred_row(s4, "Section 2", PT_S4)

# =============================================================================
# b / d — gene spatial (top 100) + boxplots
# =============================================================================
make_gene_block <- function(slim, gene, cluster, cluster_lab, pt_gene = 1.25,
                            show_legend = FALSE, he = NULL) {
  n_feat <- feat_neg(cluster)
  p_feat <- feat_pos(cluster)
  cells_pos <- top_cells(slim, p_feat, N_TOP)
  cells_neg <- top_cells(slim, n_feat, N_TOP)
  stopifnot(length(cells_pos) == N_TOP, length(cells_neg) == N_TOP)
  v <- slim$expr[gene, ]
  p_pos <- spot_plot(
    slim, v, cells_pos, title = gene,
    subtitle = paste0("POS · ", cluster_lab),
    limits = LIM_GENE, pt = pt_gene, show_legend = FALSE, he = he
  ) +
    theme(plot.title = element_text(
      size = BASE + 0.4, face = "bold.italic", hjust = 0.5, family = FAM
    ))
  p_neg <- spot_plot(
    slim, v, cells_neg, title = "NEG",
    limits = LIM_GENE, pt = pt_gene, show_legend = show_legend, he = he
  )
  df <- make_df(slim, gene, cluster, N_TOP)
  p_box <- make_box(df)
  (p_pos / p_neg / p_box) + plot_layout(heights = c(1, 1, 0.85))
}

pB <- wrap_plots(
  make_gene_block(s2, "Lgals3", 5, "cluster 5",
                  pt_gene = 1.25 * PT_S2 / 0.40, he = he2),
  make_gene_block(s2, "Prdx1",  5, "cluster 5",
                  pt_gene = 1.25 * PT_S2 / 0.40, he = he2),
  make_gene_block(s2, "Mki67",  3, "cluster 3",
                  pt_gene = 1.25 * PT_S2 / 0.40, show_legend = TRUE, he = he2),
  nrow = 1
)

pD <- wrap_plots(
  make_gene_block(s4, "Lgals3", 5, "cluster 5",
                  pt_gene = 1.25 * PT_S4 / 0.40, he = he4),
  make_gene_block(s4, "Prdx1",  5, "cluster 5",
                  pt_gene = 1.25 * PT_S4 / 0.40, he = he4),
  make_gene_block(s4, "Mki67",  3, "cluster 3",
                  pt_gene = 1.25 * PT_S4 / 0.40, show_legend = TRUE, he = he4),
  nrow = 1
)

# =============================================================================
# Assemble
# =============================================================================
tag_theme <- theme(
  plot.tag = element_text(size = 8, face = "bold", family = FAM),
  plot.tag.position = c(0.004, 0.995)
)
fig <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  wrap_elements(full = pC) /
  wrap_elements(full = pD) +
  plot_layout(heights = c(0.78, 1.50, 0.78, 1.50)) +
  plot_annotation(tag_levels = "a") &
  tag_theme

# -----------------------------------------------------------------------------
# Source data
# -----------------------------------------------------------------------------
pred_long <- function(slim, section) {
  feats <- rownames(slim$pred)
  coords <- slim$coords
  bind_rows(lapply(feats, function(ft) {
    cl <- sub("^[np]-", "", ft)
    grp <- ifelse(startsWith(ft, "p-"), "POS", "NEG")
    data.frame(
      section = section,
      cluster = cl,
      group = grp,
      feature = ft,
      barcode = colnames(slim$pred),
      x = coords[colnames(slim$pred), "x"],
      y = coords[colnames(slim$pred), "y"],
      score = as.numeric(slim$pred[ft, ]),
      stringsAsFactors = FALSE
    )
  }))
}

dfs_b <- list(
  Lgals3 = make_df(s2, "Lgals3", 5),
  Prdx1  = make_df(s2, "Prdx1",  5),
  Mki67  = make_df(s2, "Mki67",  3)
)
dfs_d <- list(
  Lgals3 = make_df(s4, "Lgals3", 5),
  Prdx1  = make_df(s4, "Prdx1",  5),
  Mki67  = make_df(s4, "Mki67",  3)
)

pval_row <- function(df, gene, panel) {
  data.frame(
    panel = panel,
    gene = gene,
    cluster = unique(df$cluster),
    wilcox_p = wilcox.test(expression ~ group, data = df, exact = FALSE)$p.value,
    n_per_group = N_TOP,
    n_before = unique(df$n_before),
    stringsAsFactors = FALSE
  )
}
pvals <- rbind(
  pval_row(dfs_b$Lgals3, "Lgals3", "b"),
  pval_row(dfs_b$Prdx1,  "Prdx1",  "b"),
  pval_row(dfs_b$Mki67,  "Mki67",  "b"),
  pval_row(dfs_d$Lgals3, "Lgals3", "d"),
  pval_row(dfs_d$Prdx1,  "Prdx1",  "d"),
  pval_row(dfs_d$Mki67,  "Mki67",  "d")
)
pvals$star <- vapply(pvals$wilcox_p, wilcox_stars, character(1))
write.csv(pvals, file.path(PLOTTED, "SourceData_FigS12_wilcoxon.csv"), row.names = FALSE)
write.csv(pvals, file.path(SRC, "FigS12_wilcoxon.csv"), row.names = FALSE)

for (nm in names(dfs_b)) {
  write.csv(dfs_b[[nm]], file.path(PLOTTED, paste0("SourceData_FigS12b_", nm, ".csv")),
            row.names = FALSE)
  write.csv(dfs_d[[nm]], file.path(PLOTTED, paste0("SourceData_FigS12d_", nm, ".csv")),
            row.names = FALSE)
}

write.csv(pred_long(s2, "section1"),
          file.path(PLOTTED, "SourceData_FigS12a_prediction.csv"), row.names = FALSE)
write.csv(pred_long(s4, "section2"),
          file.path(PLOTTED, "SourceData_FigS12c_prediction.csv"), row.names = FALSE)

fmt_p <- function(panel) {
  sub <- pvals[pvals$panel == panel, ]
  paste(sprintf("%s=%.3g", sub$gene, sub$wilcox_p), collapse = "; ")
}

# -----------------------------------------------------------------------------
# Legend / QA / stats
# -----------------------------------------------------------------------------
writeLines(c(
  "Fig. S12 | Spatial microglia subtype maps and effector-gene enrichment on two independent mouse Visium sections.",
  "a, Spatial NEG and POS prediction scores for microglia clusters 5, 3 and 0 on Visium section 1. Colour scale 0-0.5; values outside the range are clipped.",
  sprintf("b, Spatial expression of Lgals3 and Prdx1 on the top %d cluster-5 POS or NEG spots, and of Mki67 on the top %d cluster-3 POS or NEG spots (section 1), overlaid on Visium H&E (image alpha = 0.90), with matching boxplots. Centre line, median; box, IQR; whiskers, 1.5xIQR; points, individual spots (n = %d/group). Two-sided Wilcoxon rank-sum tests. Exact P: Lgals3, %.3g; Prdx1, %.3g (ns); Mki67, %.3g.",
          N_TOP, N_TOP, N_TOP, pvals$wilcox_p[pvals$panel == "b" & pvals$gene == "Lgals3"],
          pvals$wilcox_p[pvals$panel == "b" & pvals$gene == "Prdx1"],
          pvals$wilcox_p[pvals$panel == "b" & pvals$gene == "Mki67"]),
  "c, Spatial NEG and POS prediction scores for clusters 5, 3 and 0 on Visium section 2, on the same colour scale as a.",
  sprintf("d, As in b for section 2. Exact P: Lgals3, %.3g; Prdx1, %.3g (ns); Mki67, %.3g (ns).",
          pvals$wilcox_p[pvals$panel == "d" & pvals$gene == "Lgals3"],
          pvals$wilcox_p[pvals$panel == "d" & pvals$gene == "Prdx1"],
          pvals$wilcox_p[pvals$panel == "d" & pvals$gene == "Mki67"]),
  "Prediction scores describe transcriptional resemblance to a sorted reference population. n is spatial spots on two sections, not biological mouse replicates. Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS12_legend.txt"))

writeLines(c(
  "Figure S12 | Nature SI — native ggplot Visium maps",
  sprintf("a: Section 1 prediction maps; spots = %d; limits 0-0.5", n_s2),
  sprintf("b: Section 1 top-%d POS/NEG gene maps on H&E + Wilcoxon boxplots", N_TOP),
  sprintf("c: Section 2 prediction maps; spots = %d; same scale as a", n_s4),
  sprintf("d: Section 2 top-%d gene maps + Wilcoxon", N_TOP),
  paste0("Wilcoxon b: ", fmt_p("b")),
  paste0("Wilcoxon d: ", fmt_p("d")),
  sprintf("Size: %d x %d mm; font: %s %.1f pt", FIG_W_MM, FIG_H_MM, FAM, BASE),
  "Deliverables: SVG + PNG + PDF + TIFF; no TIFF plate embeds"
), file.path(OUTDIR, "FigureS12_stats.txt"))

writeLines(c(
  "Core conclusion: Two independent Visium sections recapitulate n-/p- microglia subtype maps and POS-spot enrichment of Lgals3.",
  "Archetype: image plate + quant",
  "Backend: R ggplot2/patchwork exclusive",
  sprintf("Final size: %d x %d mm", FIG_W_MM, FIG_H_MM),
  "Panel logic: a=section1 maps; b=section1 genes+quant; c=section2 maps; d=section2 genes+quant",
  "Colour: a and c share viridis, limits 0-0.5, hard-clip, matched spot size",
  sprintf("Statistics: n=%d spots/group; two-sided Wilcoxon; honest ns retained", N_TOP),
  paste0("Wilcoxon b: ", fmt_p("b")),
  paste0("Wilcoxon d: ", fmt_p("d")),
  sprintf("n_before section1=%d section2=%d; n_after gene maps=%d/group", n_s2, n_s4, N_TOP),
  "H&E: global alpha 0.90; no local contrast",
  "Display labels NEG/POS replace n-/p- titles; tests use the same groups"
), file.path(OUTDIR, "FigureS12_QA.txt"))

writeLines(c(
  "Figure S12 provenance",
  "Reuse level: structural adaptation of the native ggplot two-section Visium figure",
  "The named raster-plate SVG is stored under source_figures and was not redrawn as embeds",
  "This folder is new and does not overwrite any previous FigS9 or FigS1-S10 package"
), file.path(OUTDIR, "FigureS12_provenance.txt"))

writeLines(c(
  "Changed relative to the source Visium figure",
  "- Renumbered as Fig. S12",
  "- NEG/POS display labels instead of n-5/p-5 axis titles",
  "- Section identity as a left-side row label",
  "- TIFF export added; plotted CSVs written to source_data_plotted",
  "Unchanged: score limits, top-100 rule, Wilcoxon P, H&E alpha, gene set"
), file.path(OUTDIR, "FigureS12_build_notes.txt"))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "FigureS12")

cat("Writing SVG ...\n")
svglite::svglite(paste0(stem, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in,
                     family = FAM, bg = "white", fallback_resolution = 600)
print(fig); invisible(dev.off())

cat("Writing PNG ...\n")
ragg::agg_png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())

cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(stem, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); invisible(dev.off())

cat("DONE ->", OUTDIR, "\n")
cat("Wilcoxon b:", fmt_p("b"), "\n")
cat("Wilcoxon d:", fmt_p("d"), "\n")
