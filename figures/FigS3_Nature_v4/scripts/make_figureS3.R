#!/usr/bin/env Rscript
# =============================================================================
# Fig. S3 v4 - Nature-style support for Fig. 3 (R exclusive)
#
# Core conclusion
#   The hard RNA gate is recovered by lineage detection and by Cd4 / Cd8b1 /
#   Klrb1c / Ms4a1. The Fig. 3 functional labels are recovered by the
#   table marker set (dot plot) and by Ccr7 / S100a4 / Eomes / Ccl5.
#
# Evidence chain
#   a      RNA detection after the hard gate
#   b-e    Cd4 / Cd8b1 / Klrb1c / Ms4a1
#   f      Dot plot of table markers on the new identities
#   g-j    Ccr7 / S100a4 / Eomes / Ccl5 (same format as b-e)
#
# Archetype : quantitative grid
# Backend   : R (ggplot2 / patchwork / svglite / ragg)
# Size      : 183 x 305 mm
# New folder — does not overwrite FigS3_Nature_v3
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS3_Nature_v4")
}
SRC <- file.path(OUTDIR, "source_data")
QA  <- file.path(OUTDIR, "qa")
SDP <- file.path(OUTDIR, "source_data_plotted")
dir.create(QA, showWarnings = FALSE)
dir.create(SDP, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 305
BASE <- 6.5
FAM  <- "Arial"
PT <- function(x) x / .pt
zlim <- 2.2

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

CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
CLASS_COL <- c(
  true_CD4 = "#2F6F8F", true_CD8 = "#5B8F4E", NKT = "#B2473E",
  DP = "#C4A35A", unassigned = "#9A9A9A", doublet_c2 = "#5A5A5A"
)
CLASS_LAB <- c(true_CD4 = "True CD4", true_CD8 = "True CD8", NKT = "NKT",
               DP = "DP", unassigned = "Unassigned", doublet_c2 = "Doublet")
CLASS_SHORT <- c(true_CD4 = "CD4", true_CD8 = "CD8", NKT = "NKT",
                 DP = "DP", unassigned = "Unass.", doublet_c2 = "Dbl.")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"

FUN_LV <- c("CD4-1b", "CD4-1a", "CD8-0", "CD8-3a", "CD8-3b", "NKT", "Doublet")
FUN_LAB <- c(
  "CD4-1b"  = "Naive CD4 T",
  "CD4-1a"  = "Memory/activated CD4 T",
  "CD8-0"   = "Naive CD8 T",
  "CD8-3a"  = "Memory-like CD8 T",
  "CD8-3b"  = "Effector/cytotoxic CD8 T",
  "NKT"     = "NKT-like T",
  "Doublet" = "B\u2013T doublet"
)
FUN_SHORT <- c(
  "CD4-1b"  = "Naive CD4",
  "CD4-1a"  = "Mem. CD4",
  "CD8-0"   = "Naive CD8",
  "CD8-3a"  = "Mem. CD8",
  "CD8-3b"  = "Eff. CD8",
  "NKT"     = "NKT-like",
  "Doublet" = "B\u2013T dbl."
)
FUN_COL <- c(
  "CD4-1b"  = "#F2B6A0",
  "CD4-1a"  = "#E07A5F",
  "CD8-0"   = "#3D5A80",
  "CD8-3a"  = "#2A9D8F",
  "CD8-3b"  = "#7BC9BE",
  "NKT"     = "#B2473E",
  "Doublet" = "#9AA0A6"
)

GENE_BLOCKS <- list(
  "Naive"      = c("Cd4", "Cd3d", "Cd3e", "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r"),
  "Mem. CD4"   = c("S100a4", "Cxcr3", "Cxcr6"),
  "Treg"       = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"),
  "CD8"        = c("Cd8b1", "Cd44", "Eomes", "Il2rb"),
  "Cytotoxic"  = c("Nkg7", "Gzma", "Ccl5", "Prf1", "Gzmb"),
  "NKT"        = c("Trac", "Klrb1c"),
  "gdT"        = c("Trdc", "Trgc1", "Trgc2"),
  "Doublet"    = c("Ms4a1", "Cd79a", "Lyz2", "S100a8")
)

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
      legend.key.size = unit(2.4, "mm"),
      legend.background = element_blank(),
      strip.text = element_text(size = BASE, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = BASE + 0.2, face = "bold", hjust = 0.5),
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

# -----------------------------------------------------------------------------
# a  RNA detection after the hard gate
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
det <- read.csv(file.path(SRC, "SourceData_lineage_detection.csv"),
                stringsAsFactors = FALSE)
stopifnot(nrow(det) == 42)
det$t_class <- factor(det$t_class, levels = CLASS_LV)
det$gene <- factor(det$gene, levels = c("Cd4", "Cd8a", "Cd8b1", "Klrb1c",
                                        "Cd3e", "Ncr1", "Ms4a1"))
stopifnot(!anyNA(det$t_class), !anyNA(det$gene))

pA <- ggplot(det, aes(gene, t_class, fill = pct)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f", pct)), size = PT(5.26),
            family = FAM, colour = "grey15") +
  scale_fill_gradient(low = "white", high = "#2F6F8F", name = "% detected",
                      limits = c(0, 100),
                      guide = guide_colorbar(barwidth = unit(2.4, "mm"),
                                             barheight = unit(18, "mm"),
                                             ticks.colour = "black",
                                             frame.colour = "black",
                                             frame.linewidth = 0.2)) +
  scale_y_discrete(labels = CLASS_LAB) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(face = "italic", size = BASE,
                                   margin = margin(t = 1.2, unit = "mm")),
        axis.text.y = element_text(size = BASE),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.margin = margin(0, 0, 0, 2, "mm"),
        plot.margin = margin(4, 3, 4, 2, "mm"))

# -----------------------------------------------------------------------------
# b-e  lineage markers
# -----------------------------------------------------------------------------
cat("Panels b-e ...\n")
cells <- readRDS(file.path(SRC, "Fig3_T_cells.rds"))
stopifnot(nrow(cells) == 6098)
cells$t_class <- factor(cells$t_class, levels = CLASS_LV)
cells$t_short <- factor(unname(CLASS_SHORT[as.character(cells$t_class)]),
                        levels = unname(CLASS_SHORT))

feat_vln <- function(df, gene, pal, xvar, xlev, xlab) {
  zmax <- max(as.numeric(quantile(df$expr, 0.99, names = FALSE)), 0.15)
  peak <- names(which.max(tapply(df$expr, df[[xvar]], mean)))
  wt <- suppressWarnings(wilcox.test(df$expr[as.character(df[[xvar]]) == peak],
                                     df$expr[as.character(df[[xvar]]) != peak]))
  p_f <- ggplot(df, aes(UMAP1, UMAP2, colour = expr)) +
    geom_point(size = 0.12, alpha = 0.78, shape = 16) +
    scale_colour_gradientn(
      colours = c(SEQ_LO, SEQ_HI), limits = c(0, zmax), oob = squish,
      name = NULL, breaks = c(0, zmax),
      labels = c("0", sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = unit(1.8, "mm"),
                             barheight = unit(11, "mm"),
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
          plot.margin = margin(0.6, 0.4, 0.4, 0.4, "mm"))
  p_v <- ggplot(df, aes(.data[[xvar]], expr, fill = .data[[xvar]])) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_x_discrete(labels = xlab) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Expression",
         subtitle = sprintf("highest in %s", unname(xlab[peak]))) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5,
                                       colour = "#555555"),
          axis.text.x = element_text(size = 5.4, face = "bold", angle = 45,
                                     hjust = 1, vjust = 1),
          plot.margin = margin(0.6, 0.6, 0.4, 0.6, "mm"))
  list(plot = wrap_elements(full = p_f + p_v + plot_layout(widths = c(1.15, 0.85))),
       peak = peak, peak_p = wt$p.value)
}

mk_be <- list(b = "Cd4", c = "Cd8b1", d = "Klrb1c", e = "Ms4a1")
built_be <- list()
for (nm in names(mk_be)) {
  d <- cells
  d$expr <- cells[[mk_be[[nm]]]]
  built_be[[nm]] <- feat_vln(d, mk_be[[nm]], CLASS_COL, "t_class", CLASS_LV, CLASS_SHORT)
  cat(mk_be[[nm]], "peak", built_be[[nm]]$peak,
      "P", fmt_p_sci1(built_be[[nm]]$peak_p), "\n")
}

# -----------------------------------------------------------------------------
# f  Dot plot of table markers on new identities
# -----------------------------------------------------------------------------
cat("Panel f (dot plot) ...\n")
ann <- readRDS(file.path(SRC, "FigS3_annot_expr.rds"))
stopifnot(nrow(ann) == 6098)
ann$fun <- factor(as.character(ann$subtype), levels = FUN_LV)
ann_f <- ann %>% filter(!is.na(fun))
stopifnot(nlevels(droplevels(ann_f$fun)) == 7)

present_genes <- unlist(GENE_BLOCKS, use.names = FALSE)
present_genes <- present_genes[present_genes %in% colnames(ann)]
gene_blocks <- lapply(GENE_BLOCKS, function(g) g[g %in% present_genes])
gene_blocks <- gene_blocks[lengths(gene_blocks) > 0]
gene_order <- unlist(gene_blocks, use.names = FALSE)

dot <- ann_f %>%
  select(fun, all_of(gene_order)) %>%
  pivot_longer(all_of(gene_order), names_to = "gene", values_to = "expr") %>%
  group_by(fun, gene) %>%
  summarise(avg = mean(expr), pct = 100 * mean(expr > 0), n = n(),
            .groups = "drop") %>%
  group_by(gene) %>%
  mutate(avg_scaled = as.numeric(scale(avg))) %>%
  ungroup()
dot$avg_plot <- pmax(pmin(dot$avg_scaled, zlim), -zlim)
dot$gene <- factor(dot$gene, levels = gene_order)
dot$fun <- factor(dot$fun, levels = rev(FUN_LV))
block_edges <- head(cumsum(vapply(gene_blocks, length, integer(1))), -1) + 0.5
stopifnot(nrow(dot) == length(gene_order) * 7)

pF <- ggplot(dot, aes(gene, fun)) +
  geom_vline(xintercept = block_edges, colour = "grey88", linewidth = 0.25) +
  geom_point(aes(size = pct, colour = avg_plot), shape = 16) +
  scale_size_area(max_size = 3.4, breaks = c(10, 30, 50, 75),
                  limits = c(0, 100), name = "Expressing\ncells (%)") +
  scale_colour_gradientn(
    colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(-zlim, zlim), breaks = c(-2, 0, 2),
    name = "Scaled mean\nexpression",
    guide = guide_colorbar(barwidth = unit(2.2, "mm"), barheight = unit(16, "mm"),
                           ticks.colour = "black", frame.colour = "black",
                           frame.linewidth = 0.2, order = 1)
  ) +
  scale_y_discrete(labels = FUN_LAB) +
  guides(size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", size = 5.4),
    axis.text.y = element_text(size = 5.6, colour = "black"),
    legend.position = "right",
    legend.box.spacing = unit(1.0, "mm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.25),
    plot.margin = margin(2, 2, 1, 2, "mm")
  )

# -----------------------------------------------------------------------------
# g-j  diagnostic genes, same format as b-e, violins on new identities
# -----------------------------------------------------------------------------
cat("Panels g-j ...\n")
ann$fun <- factor(as.character(ann$subtype), levels = FUN_LV)
mk_gj <- list(g = "Ccr7", h = "S100a4", i = "Eomes", j = "Ccl5")
built_gj <- list()
for (nm in names(mk_gj)) {
  gene <- mk_gj[[nm]]
  stopifnot(gene %in% colnames(ann))
  d <- ann %>% filter(!is.na(fun))
  d$expr <- d[[gene]]
  # Feature plot uses all cells; violin uses annotated cells only.
  d_all <- ann
  d_all$expr <- d_all[[gene]]
  d_all$fun_plot <- d_all$fun
  zmax <- max(as.numeric(quantile(d_all$expr, 0.99, names = FALSE)), 0.15)
  peak <- names(which.max(tapply(d$expr, d$fun, mean)))
  wt <- suppressWarnings(wilcox.test(d$expr[as.character(d$fun) == peak],
                                     d$expr[as.character(d$fun) != peak]))
  p_f <- ggplot(d_all, aes(UMAP1, UMAP2, colour = expr)) +
    geom_point(size = 0.12, alpha = 0.78, shape = 16) +
    scale_colour_gradientn(
      colours = c(SEQ_LO, SEQ_HI), limits = c(0, zmax), oob = squish,
      name = NULL, breaks = c(0, zmax),
      labels = c("0", sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = unit(1.8, "mm"),
                             barheight = unit(11, "mm"),
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
          plot.margin = margin(0.6, 0.4, 0.4, 0.4, "mm"))
  p_v <- ggplot(d, aes(fun, expr, fill = fun)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_x_discrete(labels = FUN_SHORT) +
    scale_fill_manual(values = FUN_COL, guide = "none") +
    labs(x = NULL, y = "Expression",
         subtitle = sprintf("highest in %s", unname(FUN_SHORT[[peak]]))) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5,
                                       colour = "#555555"),
          axis.text.x = element_text(size = 5.3, face = "bold", angle = 45,
                                     hjust = 1, vjust = 1),
          plot.margin = margin(0.6, 0.6, 0.4, 0.6, "mm"))
  built_gj[[nm]] <- list(
    plot = wrap_elements(full = p_f + p_v + plot_layout(widths = c(1.15, 0.85))),
    peak = peak, peak_p = wt$p.value
  )
  cat(gene, "peak", peak, "P", fmt_p_sci1(wt$p.value), "\n")
}

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row_a <- tag_it(pA, "a")
row_be <- (tag_it(built_be$b$plot, "b") | tag_it(built_be$c$plot, "c")) /
  (tag_it(built_be$d$plot, "d") | tag_it(built_be$e$plot, "e"))
row_f <- tag_it(pF, "f")
row_gj <- (tag_it(built_gj$g$plot, "g") | tag_it(built_gj$h$plot, "h")) /
  (tag_it(built_gj$i$plot, "i") | tag_it(built_gj$j$plot, "j"))

fig <- row_a / row_be / row_f / row_gj +
  plot_layout(heights = c(0.88, 1.32, 1.18, 1.32)) &
  theme(plot.margin = margin(1.2, 1.2, 1.2, 1.2, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS3")

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

write.csv(det, file.path(SDP, "SourceData_FigS3a_detection.csv"), row.names = FALSE)
write.csv(cells[, c("cell", "UMAP1", "UMAP2", "t_class", "IV_group",
                    "Cd4", "Cd8b1", "Klrb1c", "Ms4a1")],
          file.path(SDP, "SourceData_FigS3be_markers.csv"), row.names = FALSE)
write.csv(dot, file.path(SDP, "SourceData_FigS3f_dotplot.csv"), row.names = FALSE)
write.csv(ann[, c("cell", "UMAP1", "UMAP2", "subtype", "IV_group",
                  intersect(c("Ccr7", "S100a4", "Eomes", "Ccl5"), colnames(ann)))],
          file.path(SDP, "SourceData_FigS3gj_markers.csv"), row.names = FALSE)

cat("DONE ->", OUTDIR, "\n")
