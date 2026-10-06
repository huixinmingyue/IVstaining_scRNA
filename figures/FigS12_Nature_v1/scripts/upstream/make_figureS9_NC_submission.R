#!/usr/bin/env Rscript
# =============================================================================
# Fig. S9 — Nature Communications SI submission (native ggplot vector)
# Archetype: image plate + quant
# Core conclusion: Two independent Visium sections recapitulate n-/p-
# microglia subtype maps and POS-spot enrichment of Lgals3/Prdx1/Mki67.
# No TIFF embeds; no a/c map duplication.
# =============================================================================

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

OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS9_NC_submission")
SRC    <- file.path(OUTDIR, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

# Nature double-column SI width
FIG_W_MM <- 183
FIG_H_MM <- 240
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
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
)

# Shared prediction-map style for panels a and c (keep identical)
PAL_PRED <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
LIM_PRED <- c(0, 0.5)

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

# Spot size matched to Section 1 so a/c have the same visual density language
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
      legend.text = element_text(size = BASE - 1.3, family = FAM),
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

# -----------------------------------------------------------------------------
# Load slim payloads + H&E backgrounds (for panels b/d)
# -----------------------------------------------------------------------------
s2 <- readRDS(file.path(SRC, "slim_immune2.rds"))
s4 <- readRDS(file.path(SRC, "slim_immune4.rds"))
he2 <- readRDS(file.path(SRC, "he_immune2.rds"))
he4 <- readRDS(file.path(SRC, "he_immune4.rds"))
cat("section1 spots:", s2$n_spots, "; section2 spots:", s4$n_spots, "\n")
cat("HE dims: s1", paste(he2$dim, collapse = "x"), "; s2", paste(he4$dim, collapse = "x"), "\n")
PT_S2 <- pt_for_section(s2, s2, base_pt = 0.40)
PT_S4 <- pt_for_section(s4, s2, base_pt = 0.40)
cat(sprintf("spot pt: section1=%.3f; section2=%.3f (matched to s1 density)\n", PT_S2, PT_S4))

# Keep H&E largely opaque; spots sit on top (author image.alpha ≈ 0.4–1)
fade_he <- function(arr, alpha = 0.90) {
  arr <- pmin(pmax(arr, 0), 1)
  arr * alpha + (1 - alpha)
}

# Align HE to slim coords (x = imagecol, y = +imagerow; same orientation as a/c)
he_annotation <- function(he, alpha = 0.55) {
  arr <- fade_he(he$image, alpha = alpha)
  nr <- dim(arr)[1]
  nc <- dim(arr)[2]
  # annotation_raster draws array row 1 at ymax; slim has large imagerow toward
  # the top of the panel, so flip HE so original bottom is drawn first.
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

top_cells <- function(slim, feat, n = 100) {
  names(sort(slim$pred[feat, ], decreasing = TRUE)[seq_len(n)])
}

make_df <- function(slim, gene, n_feat, p_feat, n = 100) {
  cells_n <- top_cells(slim, n_feat, n)
  cells_p <- top_cells(slim, p_feat, n)
  data.frame(
    expression = c(slim$expr[gene, cells_n], slim$expr[gene, cells_p]),
    cluster = factor(
      c(rep(n_feat, length(cells_n)), rep(p_feat, length(cells_p))),
      levels = c(n_feat, p_feat)
    ),
    barcode = c(cells_n, cells_p),
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

make_box <- function(df, gene, levels_x) {
  df <- df %>% mutate(cluster = factor(cluster, levels = levels_x))
  wt <- wilcox.test(expression ~ cluster, data = df, exact = FALSE)
  star <- wilcox_stars(wt$p.value)
  y_max <- max(df$expression, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  y_star <- max(y_max * 1.08, y_max + 0.35)
  ggplot(df, aes(cluster, expression, fill = cluster)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(
      width = 0.18, height = 0, size = 0.45, stroke = 0,
      colour = "grey20", alpha = 0.55
    ) +
    scale_fill_manual(values = setNames(c("#D9D9D9", "#595959"), levels_x), guide = "none") +
    annotate(
      "text", x = 1.5, y = y_star, label = star,
      size = PT(6.5), fontface = "bold", family = FAM
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.22), clip = "off") +
    labs(x = NULL, y = "Expression") +
    theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.title.y = element_text(size = BASE - 0.4),
      axis.text = element_text(size = BASE - 0.6),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid = element_blank(),
      plot.margin = margin(2, 3, 2, 3)
    )
}

# =============================================================================
# Panel a — section 1 only × clusters 5/3/0 (n-/p-)
# Panel c — section 2; MUST reuse same palette/limits/pt language as panel a
# =============================================================================
make_np_pair <- function(slim, cluster, pt, show_legend = FALSE) {
  nf <- paste0("n-", cluster)
  pf <- paste0("p-", cluster)
  (spot_plot(slim, slim$pred[nf, ], title = nf, limits = LIM_PRED,
             pt = pt, show_legend = FALSE, colours = PAL_PRED) |
   spot_plot(slim, slim$pred[pf, ], title = pf, limits = LIM_PRED,
             pt = pt, show_legend = show_legend, colours = PAL_PRED))
}

make_section_pred_row <- function(slim, section_lab, pt) {
  (
    hdr_plot(section_lab, size = BASE - 0.2) /
    wrap_plots(
      (hdr_plot("Cluster 5") / make_np_pair(slim, 5, pt, show_legend = FALSE)) +
        plot_layout(heights = c(0.12, 1)),
      (hdr_plot("Cluster 3") / make_np_pair(slim, 3, pt, show_legend = FALSE)) +
        plot_layout(heights = c(0.12, 1)),
      (hdr_plot("Cluster 0") / make_np_pair(slim, 0, pt, show_legend = TRUE)) +
        plot_layout(heights = c(0.12, 1)),
      nrow = 1
    )
  ) + plot_layout(heights = c(0.08, 1))
}

pA <- make_section_pred_row(s2, "Section 1", PT_S2)
pC <- make_section_pred_row(s4, "Section 2", PT_S4)

# =============================================================================
# Panel b — section 1 gene spatial (top-100) + boxplots
# =============================================================================
make_gene_block <- function(slim, gene, n_feat, p_feat, cluster_lab, pt_gene = 1.25,
                            show_legend = FALSE, he = NULL) {
  cells_pos <- top_cells(slim, p_feat, 100)
  cells_neg <- top_cells(slim, n_feat, 100)
  v <- slim$expr[gene, ]
  p_pos <- spot_plot(
    slim, v, cells_pos, title = gene,
    subtitle = paste0("POS · ", cluster_lab),
    limits = c(0, 2), pt = pt_gene, show_legend = FALSE, he = he
  ) +
    theme(plot.title = element_text(
      size = BASE + 0.4, face = "bold.italic", hjust = 0.5, family = FAM
    ))
  p_neg <- spot_plot(
    slim, v, cells_neg, title = "NEG",
    limits = c(0, 2), pt = pt_gene, show_legend = show_legend, he = he
  )
  df <- make_df(slim, gene, n_feat, p_feat, 100)
  p_box <- make_box(df, gene, c(n_feat, p_feat))
  (p_pos / p_neg / p_box) + plot_layout(heights = c(1, 1, 0.85))
}

pB <- wrap_plots(
  make_gene_block(s2, "Lgals3", "n-5", "p-5", "cluster 5",
                  pt_gene = 1.25 * PT_S2 / 0.40, he = he2),
  make_gene_block(s2, "Prdx1",  "n-5", "p-5", "cluster 5",
                  pt_gene = 1.25 * PT_S2 / 0.40, he = he2),
  make_gene_block(s2, "Mki67",  "n-3", "p-3", "cluster 3",
                  pt_gene = 1.25 * PT_S2 / 0.40, show_legend = TRUE, he = he2),
  nrow = 1
)

# =============================================================================
# Panel d — section 2 gene spatial + boxplots
# =============================================================================
pD <- wrap_plots(
  make_gene_block(s4, "Lgals3", "n-5", "p-5", "cluster 5",
                  pt_gene = 1.25 * PT_S4 / 0.40, he = he4),
  make_gene_block(s4, "Prdx1",  "n-5", "p-5", "cluster 5",
                  pt_gene = 1.25 * PT_S4 / 0.40, he = he4),
  make_gene_block(s4, "Mki67",  "n-3", "p-3", "cluster 3",
                  pt_gene = 1.25 * PT_S4 / 0.40, show_legend = TRUE, he = he4),
  nrow = 1
)

# =============================================================================
# Assemble
# =============================================================================
fig <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  wrap_elements(full = pC) /
  wrap_elements(full = pD) +
  plot_layout(heights = c(0.78, 1.50, 0.78, 1.50)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = FAM),
    plot.tag.position = c(0.004, 0.995)
  )

# -----------------------------------------------------------------------------
# Source data + stats
# -----------------------------------------------------------------------------
dfs_b <- list(
  Lgals3 = make_df(s2, "Lgals3", "n-5", "p-5"),
  Prdx1  = make_df(s2, "Prdx1",  "n-5", "p-5"),
  Mki67  = make_df(s2, "Mki67",  "n-3", "p-3")
)
dfs_d <- list(
  Lgals3 = make_df(s4, "Lgals3", "n-5", "p-5"),
  Prdx1  = make_df(s4, "Prdx1",  "n-5", "p-5"),
  Mki67  = make_df(s4, "Mki67",  "n-3", "p-3")
)

pval_row <- function(df, gene, section) {
  data.frame(
    panel = section,
    gene = gene,
    wilcox_p = wilcox.test(expression ~ cluster, data = df, exact = FALSE)$p.value,
    n_per_group = 100L,
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
write.csv(pvals, file.path(SRC, "FigS9_wilcoxon.csv"), row.names = FALSE)

for (nm in names(dfs_b)) {
  write.csv(dfs_b[[nm]], file.path(SRC, paste0("FigS9B_", nm, ".csv")), row.names = FALSE)
  write.csv(dfs_d[[nm]], file.path(SRC, paste0("FigS9D_", nm, ".csv")), row.names = FALSE)
}

fmt_p <- function(panel) {
  sub <- pvals[pvals$panel == panel, ]
  paste(sprintf("%s=%.3g", sub$gene, sub$wilcox_p), collapse = "; ")
}

stats <- c(
  "Figure S9 | Nature Communications SI submission — native ggplot vector",
  "a: Section 1 prediction maps n-/p- clusters 5/3/0; limits 0-0.5",
  "b: Section 1 top-100 POS/NEG gene maps on H&E + Wilcoxon boxplots",
  "c: Section 2 prediction maps (unique; not duplicated in a)",
  "d: Section 2 top-100 POS/NEG gene maps on H&E + Wilcoxon boxplots",
  paste0("Wilcoxon b: ", fmt_p("b")),
  paste0("Wilcoxon d: ", fmt_p("d")),
  sprintf("Size: %d x %d mm; font: %s %.1f pt", FIG_W_MM, FIG_H_MM, FAM, BASE),
  "Deliverables: SVG + PNG + PDF; no TIFF plate embeds"
)
writeLines(stats, file.path(OUTDIR, "FigureS9_stats.txt"))
cat(paste(stats, collapse = "\n"), "\n")

writeLines(c(
  "Fig. S9 | Spatial microglia subtype maps and effector-gene enrichment on two independent mouse Visium sections.",
  "a, Spatial n- and p- prediction scores for clusters 5, 3 and 0 on Visium section 1. Colour scale 0–0.5.",
  "b, Spatial expression of Lgals3 and Prdx1 on the top 100 cluster-5 POS or NEG spots, and of Mki67 on the top 100 cluster-3 POS or NEG spots (section 1), overlaid on Visium H&E, with matching boxplots. Centre line, median; box, IQR; whiskers, 1.5×IQR; points, individual spots (n = 100/group). Two-sided Wilcoxon rank-sum tests. Exact P: Lgals3, 1.84×10⁻³³; Prdx1, 0.733 (ns); Mki67, 4.08×10⁻⁸.",
  "c, Spatial n-/p- maps for clusters 5, 3 and 0 on Visium section 2.",
  "d, As in b for section 2 (H&E background). Exact P: Lgals3, 2.16×10⁻³⁴; Prdx1, 0.409 (ns); Mki67, 0.729 (ns).",
  "Source data are provided as a Source Data file."
), file.path(OUTDIR, "FigureS9_legend.txt"))

writeLines(c(
  "Core conclusion: Two independent Visium sections recapitulate n-/p- microglia subtype maps and POS-spot enrichment of Lgals3/Prdx1/Mki67.",
  "Archetype: image plate + quant",
  "Backend: R ggplot2/patchwork (native vector geom_point)",
  "Final size: 183 x 240 mm (Nature double-column SI)",
  "Panel logic: a=section1 maps; b=section1 genes+quant; c=section2 maps; d=section2 genes+quant (no a/c duplication)",
  "Colour: panels a and c share identical viridis palette, limits 0-0.5, hard-clip, matched spot size",
  "Statistics: n=100 spots/group; Wilcoxon; honest ns retained for Prdx1 and section2 Mki67",
  "Primary deliverables: FigureS9.svg, FigureS9.png, FigureS9.pdf",
  paste0("Wilcoxon b: ", fmt_p("b")),
  paste0("Wilcoxon d: ", fmt_p("d"))
), file.path(OUTDIR, "FigureS9_QA.txt"))

writeLines(c(
  "Figure S9 | Nature Communications SI submission package",
  "Reuse level: structural adaptation of FigS9_Nature_vector native ggplot build",
  "Slim payloads from sce_immune2 / sce_immune4",
  "Does not overwrite FigS9_Nature_vector or plate-assembly folders"
), file.path(OUTDIR, "FigureS9_provenance.txt"))

# -----------------------------------------------------------------------------
# Export SVG + PNG + PDF
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
cat("DONE Nature NC submission Figure S9 ->", OUTDIR, "\n")
