#!/usr/bin/env Rscript
# =============================================================================
# Fig. S6 v2 - Nature rebuild in final figures (R / ggplot2 exclusive)
#
# Core conclusion
#   In the same oligo-CD45 CITE-seq experiment as Fig. 5, splenic T-cell
#   subtypes (Fig. 3 v4 taxonomy) occupy the labelling gradient more weakly
#   than B cells: memory-like CD8 and memory/activated CD4 shares rise along
#   oligo-CD45 copies, and no gene within memory/activated CD4 T cells reaches
#   BH-adjusted significance at high versus low tails.
#
# Evidence chain
#   a  Reclustered T-cell UMAP after QC drop of B-T doublet and myeloid
#   b  Oligo-CD45 copies by subtype
#   c  Subtype share along the continuous labelling gradient
#   d  Same reclustered embedding, subset at high (>=6) vs low (<=1) copies
#   e  Composition of those two tails
#   f  Memory/activated CD4 T differential expression at the same thresholds
#
# Archetype : quantitative grid
# Size      : 183 x 220 mm
# New folder — does not overwrite FigS6_Nature_v1
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
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
  file.path(IV_FIG_ROOT, "FigS6_Nature_v2")
}
SRC <- file.path(OUTDIR, "source_data")
SD  <- file.path(OUTDIR, "source_data_plotted")
QA  <- file.path(OUTDIR, "qa")
for (d in c(SD, QA)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 183
FIG_H_MM <- 220
FAM      <- "Arial"
BASE     <- 6.5
FC_CUT   <- 0.25
PCT_SHOW <- 0.20

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

# Fig. 3 v4 display-class colours, remapped onto CITE T cluster labels.
t_lv <- c("Naive CD4 T", "Naive CD8 T", "Memory/activated CD4 T",
          "Memory-like CD8 T", "Unassigned")
t_cols <- c(
  "Naive CD4 T"            = "#F2B6A0",
  "Naive CD8 T"            = "#3D5A80",
  "Memory/activated CD4 T" = "#E07A5F",
  "Memory-like CD8 T"      = "#2A9D8F",
  "Unassigned"             = "#6E6E6E"
)
short_t <- c(
  "Naive CD4 T"            = "Naive CD4",
  "Naive CD8 T"            = "Naive CD8",
  "Memory/activated CD4 T" = "Mem./act. CD4",
  "Memory-like CD8 T"      = "Mem.-like CD8",
  "Unassigned"             = "Unassigned"
)
OLIGO <- "Oligo-CD45 copies"
# CD4 identity/activation genes only (DEG class is Mem./act. CD4).
# Label if present in the plotted table; never inject missing genes.
label_genes <- c("Ccr7", "Sell", "S100a4", "Cxcr3", "Foxp3", "Lef1", "Tcf7",
                 "Il7r", "Cd44", "Ikzf2", "Ctla4", "Il2ra", "Klf2", "S1pr1")
pf_levels <- c("Low (PF<=1)", "High (PF>=6)")
pf_lab <- c("Low (PF<=1)" = "Low (\u2264 1)",
            "High (PF>=6)" = "High (\u2265 6)")

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line    = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title   = element_text(size = BASE - 0.3, colour = "black"),
    axis.text    = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.5),
    legend.text  = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title    = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5,
                                 colour = "grey30"),
    plot.tag      = element_text(size = 8, face = "bold", family = FAM,
                                 hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    plot.margin  = margin(1, 2, 1, 1),
    panel.grid   = element_blank(),
    strip.background = element_blank(),
    strip.text   = element_text(size = BASE - 0.6, face = "bold", family = FAM)
  )

fmtp <- function(p) {
  vapply(p, function(x) {
    if (!is.finite(x)) return("NA")
    if (x < 1e-16) return("< 1e-16")
    sprintf("%.2e", x)
  }, character(1))
}

STATS <- list()

# -----------------------------------------------------------------------------
# a  UMAP
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
df_ab <- readRDS(file.path(SRC, "panel_a_umap.rds"))
df_ab[["Tcelltype"]] <- factor(as.character(df_ab[["Tcelltype"]]), levels = t_lv)
stopifnot(!anyNA(df_ab$Tcelltype))
STATS$n_cells <- nrow(df_ab)
STATS$n_by_tt <- table(df_ab[["Tcelltype"]])

lab_a <- df_ab %>%
  group_by(Tcelltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2),
            .groups = "drop") %>%
  mutate(
    lab = unname(short_t[as.character(Tcelltype)]),
    umap_1 = umap_1 + dplyr::recode(as.character(Tcelltype),
      "Naive CD4 T" = 1.9, "Naive CD8 T" = -1.9,
      "Memory/activated CD4 T" = 0.15, "Memory-like CD8 T" = -0.2,
      "Unassigned" = 1.6, .default = 0),
    umap_2 = umap_2 + dplyr::recode(as.character(Tcelltype),
      "Naive CD4 T" = -1.1, "Naive CD8 T" = -0.9,
      "Memory/activated CD4 T" = 1.15, "Memory-like CD8 T" = -1.2,
      "Unassigned" = 0.35, .default = 0)
  )

pA <- ggplot(df_ab, aes(umap_1, umap_2, colour = Tcelltype)) +
  geom_point(size = 0.18, stroke = 0, shape = 16, alpha = 0.75) +
  geom_text_repel(
    data = lab_a, aes(umap_1, umap_2, label = lab, colour = Tcelltype),
    size = 1.95, family = FAM, fontface = "bold", inherit.aes = FALSE,
    bg.colour = "white", bg.r = 0.12, max.overlaps = Inf,
    box.padding = 0.05, point.padding = 0, force = 0.2,
    min.segment.length = Inf, seed = 42, show.legend = FALSE
  ) +
  scale_colour_manual(values = t_cols, guide = "none") +
  coord_equal(expand = TRUE, clip = "off") +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "T-cell subtypes (reclustered)",
       subtitle = sprintf("n = %s cells", comma(STATS$n_cells)), tag = "a") +
  theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# -----------------------------------------------------------------------------
# b  oligo-CD45 by subtype
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
df_b <- df_ab
kw <- kruskal.test(Proteinfeatures ~ Tcelltype, data = df_b)
STATS$kruskal <- kw
df_b[["tt_lab"]] <- factor(
  unname(short_t[as.character(df_b[["Tcelltype"]])]),
  levels = unname(short_t[t_lv])
)
n_b <- df_b %>% count(tt_lab, name = "n") %>%
  mutate(lab = sprintf("n = %s", comma(n)))
y_top_b <- max(df_b[["Proteinfeatures"]], na.rm = TRUE)

pB <- ggplot(df_b, aes(tt_lab, Proteinfeatures, fill = Tcelltype)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.13, outlier.size = 0.12, outlier.stroke = 0,
               fill = "white", alpha = 0.9, linewidth = 0.25) +
  geom_text(data = n_b, aes(tt_lab, y_top_b * 1.9, label = lab),
            inherit.aes = FALSE, size = 1.8, family = FAM, colour = "grey25") +
  scale_fill_manual(values = t_cols, guide = "none") +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = NULL, y = paste0(OLIGO, " (log1p scale)"),
       title = paste0(OLIGO, " by subtype"),
       subtitle = paste0("Kruskal-Wallis P = ", fmtp(kw[["p.value"]])),
       tag = "b") +
  theme_pub +
  theme(axis.text.x = element_text(size = BASE - 1.0, angle = 32,
                                   hjust = 1, vjust = 1))

# -----------------------------------------------------------------------------
# c  subtype share along the gradient
# -----------------------------------------------------------------------------
cat("Panel c ...\n")
prop_c <- read.csv(file.path(SRC, "panel_c_bin_prop.csv"),
                   stringsAsFactors = FALSE)
cor_c  <- read.csv(file.path(SRC, "panel_c_spearman.csv"),
                   stringsAsFactors = FALSE)
prop_c[["group"]] <- factor(as.character(prop_c[["group"]]), levels = t_lv)
prop_c[["facet"]] <- factor(unname(short_t[as.character(prop_c[["group"]])]),
                            levels = unname(short_t[t_lv]))
cor_c[["facet"]]  <- factor(unname(short_t[as.character(cor_c[["group"]])]),
                            levels = levels(prop_c[["facet"]]))
STATS$cor_c <- cor_c
ann_c <- cor_c %>%
  mutate(txt = sprintf("rho = %.2f\nP = %s", spearman_rho, fmtp(spearman_p)),
         x = 1, y = Inf)

pC <- ggplot(prop_c, aes(bin_mid, pct)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.35,
              colour = "grey25", alpha = 0.12, formula = y ~ x) +
  geom_line(aes(colour = group, group = group), linewidth = 0.45) +
  geom_point(aes(colour = group, size = n), alpha = 0.9, stroke = 0) +
  geom_text(data = ann_c, aes(x = x, y = y, label = txt), inherit.aes = FALSE,
            vjust = 1.15, hjust = 0, size = 1.8, family = FAM, lineheight = 0.95,
            colour = "grey20") +
  facet_wrap(~ facet, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = t_cols, guide = "none") +
  scale_x_continuous(trans = "log1p", breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_size_continuous(name = "Cells per bin", range = c(0.9, 2.8)) +
  labs(x = paste0(OLIGO, " (bin median, log1p scale)"),
       y = "Share of T cells in bin (%)",
       title = "Subtype share along the labelling gradient", tag = "c") +
  theme_pub +
  theme(legend.position = "inside",
        legend.position.inside = c(0.84, 0.18),
        legend.direction = "vertical",
        legend.title = element_text(size = BASE - 0.8),
        legend.text = element_text(size = BASE - 1.1),
        legend.key.size = unit(2.2, "mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.spacing = unit(1.6, "mm"))

# -----------------------------------------------------------------------------
# d  split UMAP, high >=6 vs low <=1
# -----------------------------------------------------------------------------
cat("Panel d ...\n")
df_d <- readRDS(file.path(SRC, "panel_d_umap_split.rds"))
df_d[["Tcelltype"]] <- factor(as.character(df_d[["Tcelltype"]]), levels = t_lv)
df_d[["PF_class"]] <- factor(as.character(df_d[["PF_class"]]), levels = pf_levels)
stopifnot(!anyNA(df_d$PF_class), !anyNA(df_d$Tcelltype))
n_d <- as.integer(table(df_d$PF_class))
STATS$n_d <- n_d

pD <- ggplot(df_d, aes(umap_1, umap_2, colour = Tcelltype)) +
  geom_point(size = 0.22, stroke = 0, shape = 16, alpha = 0.85) +
  facet_wrap(~ PF_class, nrow = 1, labeller = as_labeller(pf_lab)) +
  scale_colour_manual(values = t_cols, name = NULL,
                      labels = unname(short_t[t_lv])) +
  coord_equal(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "Split by oligo-CD45 tails",
       subtitle = sprintf("Low \u2264 1 (n = %s) vs high \u2265 6 (n = %s)",
                          comma(n_d[1]), comma(n_d[2])),
       tag = "d") +
  theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = BASE - 1.2),
        legend.key.size = unit(2.0, "mm")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6, alpha = 1)))

# -----------------------------------------------------------------------------
# e  composition of the same tails
# -----------------------------------------------------------------------------
cat("Panel e ...\n")
prop_e <- readRDS(file.path(SRC, "panel_e_proportion.rds"))
prop_e[["Tcelltype"]] <- factor(as.character(prop_e[["Tcelltype"]]),
                                levels = t_lv)
prop_e[["PF_class"]] <- factor(as.character(prop_e[["PF_class"]]),
                               levels = pf_levels)
stopifnot(nrow(prop_e) == length(t_lv) * length(pf_levels))
xlab_e <- paste0(unname(pf_lab[pf_levels]), "\n(n = ",
                 comma(n_d), ")")

pE <- ggplot(prop_e, aes(PF_class, pct, fill = Tcelltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = t_cols, guide = "none") +
  scale_x_discrete(labels = xlab_e) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (%)",
       title = "Subtype composition",
       subtitle = "Same thresholds as d", tag = "e") +
  theme_pub +
  theme(axis.text.x = element_text(size = BASE - 0.8, lineheight = 0.95))

# -----------------------------------------------------------------------------
# f  Mem./act. CD4 DEG, PF >=6 vs <=1
# -----------------------------------------------------------------------------
cat("Panel f ...\n")
deg_class <- read.csv(file.path(SRC, "panel_f_deg_class.csv"),
                      stringsAsFactors = FALSE)
deg <- read.csv(file.path(SRC, "panel_f_DEG.csv"), stringsAsFactors = FALSE)
STATS$n_deg_export <- nrow(deg)
STATS$deg_class <- deg_class$deg_class[1]
STATS$deg_n_low <- deg_class$n_low[1]
STATS$deg_n_high <- deg_class$n_high[1]

# Display filter: drop sparse-gene Wilcoxon noise. Recorded in QA.
deg <- deg[is.finite(deg$avg_log2FC) &
             abs(deg$avg_log2FC) >= FC_CUT &
             (deg[["pct.1"]] >= PCT_SHOW | deg[["pct.2"]] >= PCT_SHOW), ]
# y uses unadjusted P because BH-adjusted P is 1 for every gene in this table.
deg[["neglog10p"]] <- -log10(pmax(deg[["p_val"]], 1e-300))
deg[["sig"]] <- ifelse(
  !is.na(deg[["p_val_adj"]]) & deg[["p_val_adj"]] < 0.05 &
    abs(deg[["avg_log2FC"]]) >= FC_CUT, "sig", "ns")
STATS$n_deg <- nrow(deg)
STATS$n_sig <- sum(deg[["sig"]] == "sig", na.rm = TRUE)
STATS$n_nom <- sum(deg[["p_val"]] < 0.05, na.rm = TRUE)

lab_f <- deg %>% filter(gene %in% label_genes)
if (!nrow(lab_f)) {
  lab_f <- deg %>% arrange(p_val) %>% slice_head(n = 6)
  STATS$label_fallback <- TRUE
} else {
  STATS$label_fallback <- FALSE
}
cat("  labelling:", paste(sort(lab_f[["gene"]]), collapse = ", "),
    " fallback=", STATS$label_fallback, " sig=", STATS$n_sig, "\n", sep = "")

pF <- ggplot(deg, aes(avg_log2FC, neglog10p)) +
  geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = 2,
             colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2,
             colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = sig), size = 0.7, alpha = 0.7, stroke = 0) +
  geom_point(data = lab_f, shape = 21, size = 1.7, stroke = 0.5,
             colour = "#08306B", fill = NA) +
  scale_colour_manual(values = c(ns = "grey72", sig = "#B2182B"),
                      guide = "none") +
  geom_text_repel(
    data = lab_f, aes(label = gene), size = 2.1, family = FAM,
    colour = "black", fontface = "italic", bg.colour = "white", bg.r = 0.12,
    max.overlaps = Inf, box.padding = 0.45, point.padding = 0.25,
    segment.size = 0.25, segment.colour = "grey40",
    min.segment.length = 0, seed = 42, force = 5
  ) +
  scale_x_continuous(expand = expansion(mult = 0.09)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  labs(x = "avg log2 fold change (high vs low oligo-CD45)",
       y = "-log10 unadjusted P",
       title = "Mem./act. CD4: high vs low",
       subtitle = sprintf(">= 6 vs <= 1 copies; BH-adjusted P < 0.05: %d / %d genes",
                          STATS$n_sig, STATS$n_deg),
       tag = "f") +
  theme_pub +
  theme(plot.subtitle = element_text(size = BASE - 1.4, hjust = 0.5,
                                     colour = "grey30", lineheight = 1.05))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
seal <- function(x) wrap_elements(full = x)
row1 <- seal(pA + pB + plot_layout(widths = c(1, 1.02)))
row2 <- seal(pC)
row3 <- seal((pD | pE | pF) + plot_layout(widths = c(1.48, 0.68, 1.18)))

fig <- row1 / row2 / row3 +
  plot_layout(heights = grid::unit.c(grid::unit(1, "null"),
                                     grid::unit(1.35, "null"),
                                     grid::unit(1.15, "null"))) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
stem <- file.path(OUTDIR, "FigureS6")

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

write.csv(df_ab[, c("cell", "umap_1", "umap_2", "cluster", "Tcelltype",
                    "Proteinfeatures")],
          file.path(SD, "SourceData_FigS6ab_umap_oligo.csv"), row.names = FALSE)
write.csv(prop_c, file.path(SD, "SourceData_FigS6c_bins.csv"), row.names = FALSE)
write.csv(cor_c, file.path(SD, "SourceData_FigS6c_spearman.csv"),
          row.names = FALSE)
write.csv(df_d, file.path(SD, "SourceData_FigS6d_umap_split.csv"),
          row.names = FALSE)
write.csv(prop_e, file.path(SD, "SourceData_FigS6e_proportion.csv"),
          row.names = FALSE)
write.csv(deg, file.path(SD, "SourceData_FigS6f_DEG.csv"), row.names = FALSE)

cnt <- as.data.frame(STATS$n_by_tt)
colnames(cnt) <- c("Tcelltype", "n")
write.csv(cnt, file.path(QA, "subtype_counts.csv"), row.names = FALSE)
print(cnt, row.names = FALSE)
cat("T cells plotted:", STATS$n_cells, "\n")
cat("d n low/high:", paste(n_d, collapse = "/"), "\n")
cat("DEG class:", STATS$deg_class, STATS$deg_n_low, "/", STATS$deg_n_high, "\n")
cat("DEG export:", STATS$n_deg_export, " plotted:", STATS$n_deg,
    " BH-sig:", STATS$n_sig, " nominal p<0.05:", STATS$n_nom, "\n")
cat("Kruskal-Wallis P =", fmtp(kw[["p.value"]]), "\n")
cat("DONE ->", OUTDIR, "\n")
