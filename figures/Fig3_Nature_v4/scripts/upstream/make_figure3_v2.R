#!/usr/bin/env Rscript
# =============================================================================
# Fig. 3 v2 - Nature-style multipanel figure (R exclusive)
#
# Core conclusion
#   A hard RNA gate splits the original CD4_T object into bona fide CD4, CD8
#   and NKT. Those identities occupy both IV compartments but with a strong
#   composition shift (true CD4 in IV-, NKT/CD8 in IV+). Canonical lineage
#   genes recover the gate, and ImmGen naive-CD4 / CD8 / NKT / NK modules
#   land on the matching gated classes. Residual true-CD4 Intra DEGs map to
#   ImmGen NKT/NK versus naive CD4/Treg (heatmap rows are genes, not cells).
#
# Evidence chain
#   a    Hard-gated identity on one UMAP
#   b    Same identities, split by IV
#   c    Class composition by IV (hero)
#   d-g  Cd4 / Cd8b1 / Klrb1c / Ms4a1
#   h    ImmGen z of true-CD4 Intra DEGs (gene-to-reference mapping)
#   i    ImmGen identity modules on the gated classes (Fig. 2d method)
#
# Archetype : quantitative grid (panel c hero)
# Backend   : R (ggplot2 / patchwork / svglite / ragg)
# Size      : 183 x 236 mm
# New folder — does not overwrite Fig3_Nature_v1
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

OUTDIR <- file.path(IV_FIG_ROOT, "Fig3_Nature_v2")
SRC    <- file.path(OUTDIR, "source_data")
QA     <- file.path(OUTDIR, "qa")
dir.create(QA, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 236
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

CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
CLASS_COL <- c(
  true_CD4 = "#2F6F8F", true_CD8 = "#5B8F4E", NKT = "#B2473E",
  DP = "#C4A35A", unassigned = "#9A9A9A", doublet_c2 = "#5A5A5A"
)
CLASS_LAB <- c(true_CD4 = "True CD4", true_CD8 = "True CD8", NKT = "NKT",
               DP = "DP", unassigned = "Unassigned", doublet_c2 = "Doublet")
CLASS_SHORT <- c(true_CD4 = "CD4", true_CD8 = "CD8", NKT = "NKT",
                 DP = "DP", unassigned = "Unass.", doublet_c2 = "Dbl.")
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"
DIV_LO <- "#2F6F8F"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#B2473E"

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
pstar1 <- function(p) {
  if (is.na(p)) "n.s." else if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < 0.05) "*" else "n.s."
}
pstar <- function(p) vapply(p, pstar1, character(1))
tag_it <- function(p, t) {
  p + labs(tag = t) +
    theme(plot.tag = element_text(size = 8, face = "bold", family = FAM,
                                  hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}
seal <- function(x) wrap_elements(full = x)

cat("Load cells ...\n")
cells <- readRDS(file.path(SRC, "Fig3_T_cells.rds"))
stopifnot(nrow(cells) == 6098)
cells$t_class <- factor(cells$t_class, levels = CLASS_LV)
cells$IV_group <- factor(cells$IV_group, levels = c("IntraNeg", "IntraPos"))
n_iv <- table(cells$IV_group)
cells$t_short <- factor(unname(CLASS_SHORT[as.character(cells$t_class)]),
                        levels = unname(CLASS_SHORT))

# -----------------------------------------------------------------------------
# a  identity UMAP
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
pA <- ggplot(cells, aes(UMAP1, UMAP2, colour = t_class)) +
  geom_point(size = 0.16, stroke = 0, alpha = 0.72, shape = 16) +
  scale_colour_manual(values = CLASS_COL, labels = CLASS_LAB, name = NULL) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = BASE - 0.8),
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), nrow = 2))

# -----------------------------------------------------------------------------
# b  same identities, split by IV
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
iv_n <- as.integer(n_iv[names(IV_LAB)])
facet_lab <- setNames(
  sprintf("%s\nn = %s", unname(IV_LAB), format(iv_n, big.mark = ",")),
  names(IV_LAB)
)
cells$facet <- factor(facet_lab[as.character(cells$IV_group)],
                      levels = unname(facet_lab))
xr <- range(cells$UMAP1); yr <- range(cells$UMAP2)
pad <- 0.05 * diff(xr)

pB <- ggplot(cells, aes(UMAP1, UMAP2, colour = t_class)) +
  geom_point(size = 0.14, alpha = 0.70, shape = 16) +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  coord_fixed(xlim = c(xr[1] - pad, xr[2] + pad),
              ylim = c(yr[1] - pad, yr[2] + pad), expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
        axis.line = element_blank(),
        strip.text = element_text(size = BASE - 0.3, face = "bold", lineheight = 0.95),
        panel.spacing = unit(1.2, "mm"))

# -----------------------------------------------------------------------------
# c  composition of gated classes (hero)
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
comp <- read.csv(file.path(SRC, "n_by_class_IV.csv"), stringsAsFactors = FALSE)
fish <- read.csv(file.path(SRC, "fisher_class_IV.csv"), stringsAsFactors = FALSE)
comp$t_class <- factor(comp$t_class, levels = CLASS_LV)
comp$IV_group <- factor(comp$IV_group, levels = c("IntraNeg", "IntraPos"))
fish$t_class <- factor(fish$t_class, levels = CLASS_LV)
comp$t_short <- factor(unname(CLASS_SHORT[as.character(comp$t_class)]),
                       levels = unname(CLASS_SHORT))
fish$t_short <- factor(unname(CLASS_SHORT[as.character(fish$t_class)]),
                       levels = unname(CLASS_SHORT))
chi <- chisq.test(xtabs(n ~ t_class + IV_group, data = comp))

pC <- ggplot(comp, aes(t_short, pct, fill = IV_group)) +
  geom_col(aes(colour = t_class), position = position_dodge(width = 0.72),
           width = 0.64, linewidth = 0.40) +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  geom_errorbar(aes(ymin = pmax(pct - se, 0), ymax = pct + se),
                position = position_dodge(width = 0.72), width = 0.18,
                linewidth = 0.28) +
  geom_text(data = fish, aes(x = t_short, y = 88, label = pstar(p_holm)),
            inherit.aes = FALSE, family = FAM, size = PT(5.6), vjust = 0.35) +
  geom_text(data = fish, aes(x = t_short, y = 95,
                             label = sprintf("OR %.2f", OR)),
            inherit.aes = FALSE, family = FAM, size = PT(5.2), colour = "#444444") +
  scale_fill_manual(values = IV_COL, name = NULL, labels = unname(IV_LAB)) +
  scale_y_continuous(limits = c(0, 102), breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Hard-gated class", y = "% of cells") +
  theme(axis.text.x = element_text(face = "bold", size = BASE - 0.3),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.3, "mm"),
        legend.margin = margin(t = -1.5, b = 0, unit = "mm"))

# -----------------------------------------------------------------------------
# d-g  lineage markers of the gated classes
# -----------------------------------------------------------------------------
cat("Panels d-g ...\n")
feat_vln <- function(df, gene, pal) {
  zmax <- max(as.numeric(quantile(df$expr, 0.99, names = FALSE)), 0.15)
  peak <- names(which.max(tapply(df$expr, df$t_short, mean)))
  wt <- suppressWarnings(wilcox.test(df$expr[as.character(df$t_short) == peak],
                                     df$expr[as.character(df$t_short) != peak]))
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
  p_v <- ggplot(df, aes(t_short, expr, fill = t_class)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Expression",
         subtitle = sprintf("highest in %s", peak)) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5,
                                       colour = "#555555"),
          axis.text.x = element_text(size = 5.6, face = "bold", angle = 45,
                                     hjust = 1, vjust = 1),
          plot.margin = margin(0.6, 0.6, 0.4, 0.6, "mm"))
  list(plot = wrap_elements(full = p_f + p_v + plot_layout(widths = c(1.15, 0.85))),
       peak = peak, peak_p = wt$p.value)
}

mk <- list(d = "Cd4", e = "Cd8b1", f = "Klrb1c", g = "Ms4a1")
built <- list()
for (nm in names(mk)) {
  d <- cells
  d$expr <- cells[[mk[[nm]]]]
  built[[nm]] <- feat_vln(d, mk[[nm]], CLASS_COL)
  cat(mk[[nm]], "peak", built[[nm]]$peak,
      "P", fmt_p_sci1(built[[nm]]$peak_p), "\n")
}

# -----------------------------------------------------------------------------
# h  ImmGen heatmap of true-CD4 Intra DEGs (rows = genes, not cells)
# -----------------------------------------------------------------------------
cat("Panel h ...\n")
zlong <- read.csv(file.path(SRC, "SourceData_ImmGen_trueCD4.csv"),
                  stringsAsFactors = FALSE)
sub_lv <- c("Naive CD4", "Total CD4", "Treg", "Naive CD8", "NKT",
            "gdT", "NK immature", "NK int", "NK mature")
zlong$subset <- factor(zlong$subset, levels = sub_lv)
zlong$direction <- recode(zlong$direction,
                          `IV- up` = "IV\u2212 up",
                          `IV+ up` = "IV+ up")
zlong$direction <- factor(zlong$direction, levels = c("IV\u2212 up", "IV+ up"))
src_ord <- unique(zlong$gene)
zlong$gene <- factor(zlong$gene, levels = rev(src_ord))
n_gene <- n_distinct(zlong$gene)
stopifnot(n_gene == 54)

pH <- ggplot(zlong, aes(subset, gene, fill = z)) +
  geom_tile() +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI,
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                       name = "z", breaks = c(-2.5, 0, 2.5),
                       guide = guide_colorbar(barwidth = unit(2.2, "mm"),
                                              barheight = unit(16, "mm"),
                                              ticks.colour = "black",
                                              frame.colour = "black",
                                              frame.linewidth = 0.2)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = BASE - 0.5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 5.6, face = "bold", angle = 90),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(1, 1, 1, 2, "mm"))

# -----------------------------------------------------------------------------
# i  ImmGen modules on gated classes (Fig. 2d method)
# -----------------------------------------------------------------------------
cat("Panel i ...\n")
sc <- read.csv(file.path(SRC, "SourceData_Fig3_module_scores_per_cell.csv"),
               stringsAsFactors = FALSE)
sc$t_class <- factor(sc$t_class, levels = CLASS_LV)
sc$t_short <- factor(unname(CLASS_SHORT[as.character(sc$t_class)]),
                     levels = unname(CLASS_SHORT))
sc$module <- factor(sc$module,
                    levels = c("Naive_CD4", "Naive_CD8", "NKT", "NK"),
                    labels = c("Naive CD4 signature", "Naive CD8 signature",
                               "NKT signature", "NK signature"))
sc$IV <- recode(sc$IV_group, !!!IV_LAB)
stopifnot(nrow(sc) == 6098 * 4)

pI <- ggplot(sc, aes(t_short, score, fill = IV)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, width = 0.86,
              position = position_dodge(width = 0.85)) +
  geom_boxplot(width = 0.12, outlier.size = 0.10, outlier.alpha = 0.18,
               position = position_dodge(width = 0.85), linewidth = 0.25) +
  facet_wrap(~ module, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("IV\u2212 (extravascular)" = IV_COL[["IntraNeg"]],
                               "IV+ (intravascular)" = IV_COL[["IntraPos"]]),
                    name = NULL) +
  labs(x = "Hard-gated class", y = "Module score") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        axis.text.x = element_text(size = 5.6, face = "bold", angle = 45,
                                   hjust = 1, vjust = 1),
        strip.text = element_text(size = BASE - 0.2, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- seal(tag_it(pA, "a") + tag_it(pB, "b") + plot_layout(widths = c(0.92, 1.08)))
row2 <- tag_it(pC, "c")
row3 <- seal(
  (tag_it(built$d$plot, "d") | tag_it(built$e$plot, "e")) /
    (tag_it(built$f$plot, "f") | tag_it(built$g$plot, "g"))
)
row4 <- seal(tag_it(pH, "h") + tag_it(pI, "i") + plot_layout(widths = c(1.05, 0.95)))

fig <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(0.95, 0.88, 1.32, 1.20)) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "Figure3")

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
write.csv(cells[, c("cell", "UMAP1", "UMAP2", "t_class", "IV_group",
                    "Cd4", "Cd8b1", "Klrb1c", "Ms4a1")],
          file.path(SD, "SourceData_Fig3_cells.csv"), row.names = FALSE)
write.csv(comp, file.path(SD, "SourceData_Fig3c_composition.csv"), row.names = FALSE)
write.csv(fish, file.path(SD, "SourceData_Fig3c_fisher.csv"), row.names = FALSE)
write.csv(zlong, file.path(SD, "SourceData_Fig3h_ImmGen.csv"), row.names = FALSE)
write.csv(sc[, c("cell", "t_class", "IV_group", "module", "score")],
          file.path(SD, "SourceData_Fig3i_module_scores.csv"), row.names = FALSE)

cat("Cells:", nrow(cells), " IV-:", n_iv[["IntraNeg"]], " IV+:", n_iv[["IntraPos"]], "\n")
cat("chi-square panel c P", fmt_p_sci1(chi$p.value), "\n")
cat("DONE ->", OUTDIR, "\n")
