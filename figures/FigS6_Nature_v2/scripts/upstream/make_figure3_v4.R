#!/usr/bin/env Rscript
# =============================================================================
# Fig. 3 v4 - Nature-style multipanel figure (R exclusive)
#
# Core conclusion
#   A hard RNA gate splits the original CD4_T object into bona fide CD4, CD8
#   and NKT. CD4 1a/1b and CD8 0/3a/3b match naive versus memory/effector
#   programmes. Identities occupy both IV compartments but shift in
#   composition. Residual true-CD4 Intra DEGs map to ImmGen NKT/NK versus
#   naive CD4/Treg, and ImmGen modules land on the matching gated classes.
#
# Evidence chain
#   a    IV-split UMAP; CD4/CD8 islands labelled by function
#   b    Subtype composition by IV (Fig. 1e style; hero)
#   c    ImmGen z of true-CD4 Intra DEGs
#   d    Nine ImmGen identity modules matching panel c columns
#
# Archetype : quantitative grid (panel b hero)
# Backend   : R (ggplot2 / patchwork / svglite / ragg)
# Size      : 183 x 216 mm
# New folder — does not overwrite Fig3_Nature_v3
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
  file.path(IV_FIG_ROOT, "Fig3_Nature_v4")
}
SRC <- file.path(OUTDIR, "source_data")
QA  <- file.path(OUTDIR, "qa")
SDP <- file.path(OUTDIR, "source_data_plotted")
dir.create(QA, showWarnings = FALSE)
dir.create(SDP, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 216
width_mm <- FIG_W_MM
height_mm <- FIG_H_MM
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
CLASS_SHORT <- c(true_CD4 = "CD4", true_CD8 = "CD8", NKT = "NKT",
                 DP = "DP", unassigned = "Unass.", doublet_c2 = "Dbl.")
SUB_LV <- c(
  "CD8-0", "CD4-1a", "CD4-1b", "Doublet", "CD8-3a", "CD8-3b",
  "CD4-other", "CD8-other", "NKT", "DP", "Unassigned"
)
# Cluster colours from Fig. S2 panel c for 0/1a/1b/2/3a/3b; extras for residual classes.
SUB_COL <- c(
  "CD8-0"      = "#3D5A80",
  "CD4-1a"     = "#E07A5F",
  "CD4-1b"     = "#F2B6A0",
  "Doublet"    = "#9AA0A6",
  "CD8-3a"     = "#2A9D8F",
  "CD8-3b"     = "#7BC9BE",
  "CD4-other"  = "#A65D45",
  "CD8-other"  = "#7A9BB8",
  "NKT"        = "#B2473E",
  "DP"         = "#C4A35A",
  "Unassigned" = "#6E6E6E"
)
SUB_LAB <- c(
  "CD8-0"      = "Naive CD8 T",
  "CD4-1a"     = "Memory/activated CD4 T",
  "CD4-1b"     = "Naive CD4 T",
  "Doublet"    = "B\u2013T doublet",
  "CD8-3a"     = "Memory-like CD8 T",
  "CD8-3b"     = "Effector/cytotoxic CD8 T",
  "CD4-other"  = "CD4 other",
  "CD8-other"  = "CD8 other",
  "NKT"        = "NKT-like T",
  "DP"         = "DP",
  "Unassigned" = "Unassigned"
)
UMAP_CL_LAB <- c(
  "0"  = "Naive CD8",
  "1a" = "Mem. CD4",
  "1b" = "Naive CD4",
  "2"  = "Doublet",
  "3a" = "Mem. CD8",
  "3b" = "Eff. CD8"
)
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")
DIV_LO <- "#2F6F8F"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#B2473E"

MOD_LV <- c("Naive_CD4", "Total_CD4", "Treg", "Naive_CD8", "NKT",
            "gdT", "NK_immature", "NK_int", "NK_mature")
MOD_LAB <- c(
  Naive_CD4 = "Naive CD4", Total_CD4 = "Total CD4", Treg = "Treg",
  Naive_CD8 = "Naive CD8", NKT = "NKT", gdT = "gdT",
  NK_immature = "NK immature", NK_int = "NK int", NK_mature = "NK mature"
)
HEAT_LV <- c("Naive CD4", "Total CD4", "Treg", "Naive CD8", "NKT",
             "gdT", "NK immature", "NK int", "NK mature")

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
seal <- function(x) wrap_elements(full = x)

cat("Load cells ...\n")
cells <- readRDS(file.path(SRC, "Fig3_T_cells.rds"))
stopifnot(nrow(cells) == 6098)
cells$t_class <- factor(cells$t_class, levels = CLASS_LV)
cells$subtype <- factor(as.character(cells$subtype), levels = SUB_LV)
cells$IV_group <- factor(cells$IV_group, levels = c("IntraNeg", "IntraPos"))
n_iv <- table(cells$IV_group)
cells$t_short <- factor(unname(CLASS_SHORT[as.character(cells$t_class)]),
                        levels = unname(CLASS_SHORT))
stopifnot(!anyNA(cells$subtype))

# -----------------------------------------------------------------------------
# a  IV-split UMAP with CD4/CD8 subtypes (legend from former panel a)
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
iv_n <- as.integer(n_iv[names(IV_LAB)])
facet_lab <- setNames(
  sprintf("%s\nn = %s", unname(IV_LAB), format(iv_n, big.mark = ",")),
  names(IV_LAB)
)
cells$facet <- factor(facet_lab[as.character(cells$IV_group)],
                      levels = unname(facet_lab))
xr <- range(cells$UMAP1); yr <- range(cells$UMAP2)
pad <- 0.06 * diff(xr)
lab_pos <- cells %>%
  group_by(cluster_ref) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop") %>%
  mutate(lab = unname(UMAP_CL_LAB[as.character(cluster_ref)]))
lab_pos$lab[is.na(lab_pos$lab)] <- as.character(lab_pos$cluster_ref[is.na(lab_pos$lab)])

pA <- ggplot(cells, aes(UMAP1, UMAP2, colour = subtype)) +
  geom_point(size = 0.14, alpha = 0.70, shape = 16) +
  geom_text(data = lab_pos, aes(x, y, label = lab), inherit.aes = FALSE,
            family = FAM, fontface = "bold", size = PT(5.4), colour = "black") +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = SUB_COL, labels = SUB_LAB, name = NULL,
                      breaks = SUB_LV) +
  coord_fixed(xlim = c(xr[1] - pad, xr[2] + pad),
              ylim = c(yr[1] - pad, yr[2] + pad * 1.15), expand = FALSE) +
  labs(x = NULL, y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
        axis.line = element_blank(),
        strip.text = element_text(size = BASE - 0.3, face = "bold", lineheight = 0.95),
        panel.spacing = unit(1.2, "mm"),
        legend.position = "bottom",
        legend.text = element_text(size = 5.4),
        legend.margin = margin(t = 1.5, b = 0.5, unit = "mm"),
        legend.box.spacing = unit(1.2, "mm"),
        plot.margin = margin(1, 1, 2.5, 1, "mm")) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1),
                               nrow = 3, byrow = TRUE))

# -----------------------------------------------------------------------------
# b  composition, Fig. 1e style (horizontal %, IV facets, no error bars)
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
prop <- cells %>%
  count(subtype, IV_group, name = "Freq") %>%
  group_by(IV_group) %>%
  mutate(total = sum(Freq), Percentage = 100 * Freq / total) %>%
  ungroup()
# Keep every subtype in both facets (zeros if absent).
prop <- tidyr::complete(prop, subtype = factor(SUB_LV, levels = SUB_LV),
                        IV_group, fill = list(Freq = 0, total = 0, Percentage = 0))
prop$total <- ifelse(prop$IV_group == "IntraNeg", n_iv[["IntraNeg"]], n_iv[["IntraPos"]])
prop$Percentage <- 100 * prop$Freq / prop$total
mean_pct <- prop %>%
  group_by(subtype) %>%
  summarise(m = mean(Percentage), .groups = "drop") %>%
  arrange(desc(m))
ct_order <- as.character(mean_pct$subtype)
prop$lineage <- factor(as.character(prop$subtype), levels = rev(ct_order))
prop$IV_lab <- factor(unname(IV_LAB[as.character(prop$IV_group)]),
                      levels = unname(IV_LAB))
xmax <- max(ceiling(max(prop$Percentage) / 10) * 10, 40)
x_breaks <- pretty(c(0, xmax), n = 4)
x_breaks <- x_breaks[x_breaks <= xmax]

pB <- ggplot(prop, aes(Percentage, lineage, fill = lineage)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", Percentage)), hjust = -0.12,
            size = PT(5.26), family = FAM, colour = "black") +
  facet_wrap(~ IV_lab, nrow = 1) +
  scale_fill_manual(values = SUB_COL, guide = "none") +
  scale_y_discrete(labels = SUB_LAB) +
  scale_x_continuous(name = "Cells (%)", limits = c(0, xmax),
                     breaks = x_breaks,
                     expand = expansion(mult = c(0, 0.28))) +
  labs(y = NULL) +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.4),
    strip.text = element_text(size = BASE - 0.5),
    plot.margin = margin(1, 2, 1, 1, "mm")
  )

chi <- chisq.test(xtabs(Freq ~ subtype + IV_group, data = prop))

# -----------------------------------------------------------------------------
# c  ImmGen heatmap of true-CD4 Intra DEGs (rows = genes, not cells)
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
zlong <- read.csv(file.path(SRC, "SourceData_ImmGen_trueCD4.csv"),
                  stringsAsFactors = FALSE)
zlong$subset <- factor(zlong$subset, levels = HEAT_LV)
zlong$direction <- recode(zlong$direction,
                          `IV- up` = "IV\u2212 up",
                          `IV+ up` = "IV+ up")
zlong$direction <- factor(zlong$direction, levels = c("IV\u2212 up", "IV+ up"))
src_ord <- unique(zlong$gene)
zlong$gene <- factor(zlong$gene, levels = rev(src_ord))
n_gene <- n_distinct(zlong$gene)
stopifnot(n_gene == 54)
stopifnot(all(HEAT_LV %in% unique(as.character(zlong$subset))))

pC <- ggplot(zlong, aes(subset, gene, fill = z)) +
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
# d  ImmGen modules for every panel-c cell type (Fig. 2d method)
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
sc <- read.csv(file.path(SRC, "SourceData_Fig3_module_scores_per_cell.csv"),
               stringsAsFactors = FALSE)
sc$t_class <- factor(sc$t_class, levels = CLASS_LV)
sc$t_short <- factor(unname(CLASS_SHORT[as.character(sc$t_class)]),
                     levels = unname(CLASS_SHORT))
sc$module <- factor(sc$module, levels = MOD_LV, labels = unname(MOD_LAB[MOD_LV]))
sc$IV <- recode(sc$IV_group, !!!IV_LAB)
stopifnot(nrow(sc) == 6098 * 9)
stopifnot(nlevels(sc$module) == 9)
stopifnot(identical(levels(sc$module), unname(MOD_LAB[MOD_LV])))

pD <- ggplot(sc, aes(t_short, score, fill = IV)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, width = 0.86,
              position = position_dodge(width = 0.85)) +
  geom_boxplot(width = 0.12, outlier.size = 0.10, outlier.alpha = 0.18,
               position = position_dodge(width = 0.85), linewidth = 0.25) +
  facet_wrap(~ module, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = c("IV\u2212 (extravascular)" = IV_COL[["IntraNeg"]],
                               "IV+ (intravascular)" = IV_COL[["IntraPos"]]),
                    name = NULL) +
  labs(x = "Hard-gated class", y = "Module score") +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        axis.text.x = element_text(size = 5.6, face = "bold", angle = 45,
                                   hjust = 1, vjust = 1),
        strip.text = element_text(size = BASE - 0.4, face = "bold"),
        panel.spacing = unit(1.0, "mm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row1 <- tag_it(pA, "a")
row2 <- tag_it(pB, "b")
row3 <- seal(tag_it(pC, "c") + tag_it(pD, "d") + plot_layout(widths = c(0.95, 1.05)))

fig <- row1 / row2 / row3 +
  plot_layout(heights = c(1.18, 1.08, 1.38)) &
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

write.csv(cells[, c("cell", "UMAP1", "UMAP2", "t_class", "subtype",
                    "IV_group", "cluster_ref")],
          file.path(SDP, "SourceData_Fig3a_cells.csv"), row.names = FALSE)
write.csv(prop, file.path(SDP, "SourceData_Fig3b_composition.csv"), row.names = FALSE)
write.csv(zlong, file.path(SDP, "SourceData_Fig3c_ImmGen.csv"), row.names = FALSE)
write.csv(sc[, c("cell", "t_class", "subtype", "IV_group", "module", "score")],
          file.path(SDP, "SourceData_Fig3d_module_scores.csv"), row.names = FALSE)

cat("Cells:", nrow(cells), " IV-:", n_iv[["IntraNeg"]], " IV+:", n_iv[["IntraPos"]], "\n")
cat("chi-square panel b P", fmt_p_sci1(chi$p.value), "\n")
cat("DONE ->", OUTDIR, "\n")
