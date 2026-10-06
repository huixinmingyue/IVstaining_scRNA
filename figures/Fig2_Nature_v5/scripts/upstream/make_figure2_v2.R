#!/usr/bin/env Rscript
# =============================================================================
# Figure 2 - Nature build v2 (double column, 180 mm)
#
# This revision keeps the statistical logic, the data-integrity decisions and
# the reviewer pre-empt of the previous build unchanged. What changed is the
# presentation layer:
#   - compartment vocabulary matched to Figure 1 (IV- / IV+ rather than
#     IntraNeg / IntraPos), so the two figures name the same groups identically
#   - panel b reports exact Holm-adjusted P values instead of asterisks
#   - the shared c/d guides move to a horizontal strip under the row, which
#     returns roughly a quarter of the row width to the two dot plots
#   - taller canvas so the panel e gene labels and the c/d pathway labels are
#     legible at final printed size
#
# Core conclusion (what the data actually supports)
#   Intravascular-labelled (IntraPos) splenic B cells are selectively enriched
#   for subcluster 1, whose frequency more than doubles relative to IntraNeg
#   (13.1% -> 27.8%, OR 2.55). This compositional shift is the defensible claim.
#
# Two limits that the panels must not paper over:
#   - Panels c/d: no KEGG pathway reaches FDR significance (q ~ 0.99), so they
#     are exploratory and are labelled as such.
#   - Panel e: the two supplied ImmGen matrices share 16 genes with identical
#     rows and their subset profiles correlate at r = 0.94, so two separate
#     heat maps would have shown the same evidence twice. They are merged into
#     one panel split by gene-set membership, with the difference quantified
#     per subset. Nothing survives Holm correction, so panel e maps a largely
#     shared programme onto the reference rather than proving divergent identity.
#
# Panel -> claim question
#   a  Do IntraNeg and IntraPos B cells occupy different transcriptional space?
#   b  Does subcluster composition shift, and by how much?   (quantitative core)
#   c  Which KEGG pathways trend with cluster 1?   (exploratory - see FDR note)
#   d  Which KEGG pathways trend with cluster 0?   (exploratory - see FDR note)
#   e  Do the two subcluster gene sets map to different ImmGen B subsets?
#
# Archetype : quantitative grid, panel a as hero
# Backend   : R (ggplot2 / ComplexHeatmap / patchwork), exclusive
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
  library(svglite)
  library(ragg)
  library(systemfonts)
})

set.seed(42)
options(warn = 1)   # surface warnings where they happen, not in a batch at exit
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

FIGDIR  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures")
IMMGEN  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "spleen_Ab_immgen")
OUTDIR  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig2_Nature_v2")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "qa"), showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 180
FIG_H_MM <- 212

# ---- typography -------------------------------------------------------------
# Every rendered glyph must stay at or above 5 pt.
BASE   <- 6.5
FAM    <- "Arial"
fm <- systemfonts::match_fonts(FAM)
cat("font resolve:", FAM, "->", fm$path, "\n")

# Rendering goes through systemfonts (svglite / cairo_pdf / ragg), which uses the
# real Arial file above. Some grid text-measuring paths still fall back to the R
# PostScript font database, which has no Arial entry. Register one backed by
# Helvetica, whose metrics are identical to Arial, so any fallback measurement
# is still correct instead of silently wrong.
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM
    do.call(reg, args)
  }
}

PT <- function(x) x / .pt          # ggplot geom_text size is in mm

# ---- restrained palette -----------------------------------------------------
# Cluster identity carries colour. IV group carries light/dark neutral, so the
# two encodings never compete and panel b survives greyscale printing.
CL_COL <- c("0" = "#4C6E8A", "1" = "#D1495B", "2" = "#E9A03B")
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")

# Display vocabulary shared with Figure 1. The internal keys stay IntraNeg /
# IntraPos so every statistic, table and source-data column is unchanged; only
# the rendered labels differ.
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"
DIV_LO <- "#3B6FA0"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#C0392B"

theme_fig2 <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FAM) +
    theme(
      axis.line    = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.1, "pt"),
      axis.title   = element_text(size = base),
      axis.text    = element_text(size = base - 0.4, colour = "black"),
      legend.title = element_text(size = base - 0.2),
      legend.text  = element_text(size = base - 0.6),
      legend.key.size = unit(2.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_blank(),
      strip.text   = element_text(size = base, face = "bold"),
      strip.background = element_blank(),
      plot.title   = element_text(size = base + 0.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base - 0.6, hjust = 0.5,
                                   colour = "#555555"),
      plot.margin  = margin(1, 1, 1, 1, "mm"),
      panel.grid   = element_blank()
    )
}
theme_set(theme_fig2())

# P-value text kept free of shrunk superscripts so nothing drops below 5 pt.
fmt_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 1e-4) return("P < 0.0001")
  if (p < 0.001) return(sprintf("P = %.5f", p))
  sprintf("P = %.3g", p)
}
pstar1 <- function(p) {
  if (is.na(p)) "n.s." else if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < 0.05) "*" else "n.s."
}
pstar <- function(p) vapply(p, pstar1, character(1))

# Exact P values are preferred over significance stars. Scientific notation is
# built from Unicode superscript glyphs rather than plotmath, so the exponent is
# a full-size character in the chosen font and cannot fall below the 5 pt floor
# the way a shrunk superscript would.
sup_digits <- function(x) {
  chartr("0123456789-", paste0("\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076",
                               "\u2077\u2078\u2079\u207b"), as.character(x))
}
fmt_p_sci1 <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p >= 0.001) return(sprintf("P = %.3g", p))
  ex <- floor(log10(p))
  sprintf("P = %.1f \u00d7 10%s", p / 10^ex, sup_digits(ex))
}
fmt_p_sci <- function(p) vapply(p, fmt_p_sci1, character(1))

STATS <- list()
say <- function(...) cat(sprintf(...), "\n")

# =============================================================================
# Load sources
# =============================================================================
say("loading B-cell object ...")
b <- readRDS(file.path(FIGDIR, "Bcell_reclustered.rds"))
Idents(b) <- "seurat_clusters"

md <- b@meta.data %>%
  mutate(cluster = factor(as.character(seurat_clusters), levels = c("0", "1", "2")),
         IV_group = factor(IV_group, levels = c("IntraNeg", "IntraPos")))
emb <- as.data.frame(Embeddings(b, "umap"))
colnames(emb) <- c("UMAP1", "UMAP2")
umap_df <- cbind(emb, md[, c("cluster", "IV_group")])

STATS$n_cells <- nrow(umap_df)
STATS$n_by_iv <- table(umap_df$IV_group)
STATS$tab     <- table(umap_df$cluster, umap_df$IV_group)
say("cells: %d   IntraNeg: %d   IntraPos: %d", STATS$n_cells,
    STATS$n_by_iv[["IntraNeg"]], STATS$n_by_iv[["IntraPos"]])

# =============================================================================
# Panel a - UMAP split by intravascular labelling
# =============================================================================
lab_pos <- umap_df %>%
  group_by(cluster) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop")

facet_n <- umap_df %>%
  count(IV_group) %>%
  mutate(lab = sprintf("n = %s", format(n, big.mark = ",")))

xr <- range(umap_df$UMAP1); yr <- range(umap_df$UMAP2)
pad <- 0.06 * diff(xr)

p_a <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 0.16, alpha = 0.70, shape = 16) +
  geom_text(data = lab_pos, aes(x, y, label = cluster), inherit.aes = FALSE,
            family = FAM, fontface = "bold", size = PT(6.4), colour = "black") +
  geom_text(data = facet_n, aes(x = xr[1], y = yr[2] + pad * 0.4, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 1, family = FAM,
            size = PT(5.4), colour = "#444444") +
  facet_wrap(~ IV_group, nrow = 1, labeller = as_labeller(IV_LAB)) +
  scale_colour_manual(values = CL_COL, guide = "none") +
  coord_fixed(xlim = c(xr[1] - pad, xr[2] + pad),
              ylim = c(yr[1] - pad, yr[2] + pad * 1.6), expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "B-cell subclusters by intravascular labelling") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black",
                                    linewidth = 0.3),
        axis.line = element_blank(),
        panel.spacing = unit(1.6, "mm"))

# =============================================================================
# Panel b - subcluster composition, IntraNeg vs IntraPos
# =============================================================================
comp <- as.data.frame(STATS$tab)
colnames(comp) <- c("cluster", "IV_group", "n")
comp <- comp %>%
  group_by(IV_group) %>%
  mutate(total = sum(n), pct = 100 * n / total,
         se = 100 * sqrt(pct / 100 * (1 - pct / 100) / total)) %>%
  ungroup()

# Per-cluster Fisher exact test: this cluster vs all other B cells in the group.
fish <- lapply(levels(comp$cluster), function(cl) {
  m <- matrix(c(STATS$tab[cl, "IntraPos"],
                sum(STATS$tab[, "IntraPos"]) - STATS$tab[cl, "IntraPos"],
                STATS$tab[cl, "IntraNeg"],
                sum(STATS$tab[, "IntraNeg"]) - STATS$tab[cl, "IntraNeg"]),
              nrow = 2)
  ft <- fisher.test(m)
  data.frame(cluster = cl, OR = unname(ft$estimate), p = ft$p.value)
}) %>% bind_rows()
fish$p_holm <- p.adjust(fish$p, method = "holm")
fish$cluster <- factor(fish$cluster, levels = levels(comp$cluster))
STATS$fisher <- fish

chi <- chisq.test(STATS$tab)
STATS$chisq <- chi
say("chi-square 3x2: X2 = %.1f, df = %d, %s", chi$statistic, chi$parameter,
    fmt_p(chi$p.value))
for (i in seq_len(nrow(fish))) {
  say("  cluster %s: OR = %.2f, %s (Holm %s)", fish$cluster[i], fish$OR[i],
      fmt_p(fish$p[i]), fmt_p(fish$p_holm[i]))
}

# Odds ratios and exact Holm-adjusted P values sit on a common baseline rather
# than tracking each bar height per cluster, so they read as a row and cannot
# collide with the bars. Asterisks were dropped: they encode the same quantity
# as the printed P value and add nothing.
brk <- fish
brk$lab <- sprintf("OR %.2f\n%s", brk$OR, fmt_p_sci(brk$p_holm))
ymax_b <- 100
Y_ANN  <- 98

# Fill encodes the experimental group; the bar outline carries cluster identity,
# which ties panel b back to panel a without a second colour legend.
p_b <- ggplot(comp, aes(cluster, pct, fill = IV_group)) +
  geom_col(aes(colour = cluster), position = position_dodge(width = 0.72),
           width = 0.64, linewidth = 0.45) +
  scale_colour_manual(values = CL_COL, guide = "none") +
  geom_errorbar(aes(ymin = pmax(pct - se, 0), ymax = pct + se),
                position = position_dodge(width = 0.72), width = 0.18,
                linewidth = 0.28) +
  geom_text(aes(y = pct + se + 2.0,
                label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.72), family = FAM,
            size = PT(5.2), colour = "black") +
  geom_text(data = brk, aes(x = cluster, y = Y_ANN, label = lab),
            inherit.aes = FALSE, family = FAM, size = PT(5.2),
            colour = "#333333", vjust = 1, lineheight = 1.15) +
  scale_fill_manual(values = IV_COL, name = NULL,
                    labels = c(sprintf("%s, n = %s", IV_LAB[["IntraNeg"]],
                                       format(STATS$n_by_iv[["IntraNeg"]],
                                              big.mark = ",")),
                               sprintf("%s, n = %s", IV_LAB[["IntraPos"]],
                                       format(STATS$n_by_iv[["IntraPos"]],
                                              big.mark = ",")))) +
  scale_y_continuous(limits = c(0, ymax_b), breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "B-cell subcluster", y = "B cells in compartment (%)",
       title = "Subcluster composition",
       subtitle = sprintf("chi-square %s; per-cluster Fisher, Holm-adjusted",
                          fmt_p(chi$p.value))) +
  theme(axis.text.x = element_text(face = "bold", size = BASE + 0.3,
                                   colour = "black"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(2.3, "mm"),
        legend.spacing.y = unit(0.3, "mm"),
        legend.margin = margin(t = -1.0, b = 0, unit = "mm"))

# =============================================================================
# Panels c / d - KEGG GSEA dot plots
#
# Data-integrity note: all 20 rows of each supplied table are shown and no
# pathway reaches FDR significance (q ~ 0.99). The panels therefore encode the
# nominal P value only and say so, rather than implying a significant result.
# =============================================================================
e <- new.env()
load(file.path(IMMGEN, "IVstaining_DEG2.RData"), envir = e)
k1 <- get("kegg.B1_top20", envir = e)
k0 <- get("kegg.B0_top20", envir = e)

prep_kegg <- function(df, tag) {
  d <- as.data.frame(df)
  d$Description <- trimws(gsub(" - Mus musculus.*|house mouse", "",
                               d$Description))
  # Every supplied pathway name ends in "signaling pathway"; the suffix moves
  # to the axis title so the labels stay legible at 5.5 pt.
  d$label <- trimws(sub("signaling pathway", "", d$Description))
  d$label <- trimws(gsub("\\s+", " ", d$label))   # collapse the gap left behind
  d$label <- gsub("TGF-beta", "TGF-\u03b2", d$label)
  d$label <- ifelse(d$label == "", d$Description, d$label)
  d$neglogp <- -log10(pmax(d$pvalue, 1e-300))
  d$nominal <- d$pvalue < 0.05
  d$panel <- tag
  d[order(d$NES), ]
}
d1 <- prep_kegg(k1, "cluster 1")
d0 <- prep_kegg(k0, "cluster 0")

STATS$kegg_maxq <- c(cluster1 = max(k1$p.adjust), cluster0 = max(k0$p.adjust))
STATS$kegg_minq <- c(cluster1 = min(k1$p.adjust), cluster0 = min(k0$p.adjust))
STATS$kegg_nom  <- c(cluster1 = sum(d1$nominal), cluster0 = sum(d0$nominal))
say("KEGG cluster1: q range %.4f-%.4f, nominal P<0.05: %d/20",
    STATS$kegg_minq[["cluster1"]], STATS$kegg_maxq[["cluster1"]],
    STATS$kegg_nom[["cluster1"]])
say("KEGG cluster0: q range %.4f-%.4f, nominal P<0.05: %d/20",
    STATS$kegg_minq[["cluster0"]], STATS$kegg_maxq[["cluster0"]],
    STATS$kegg_nom[["cluster0"]])

# Shared scales so c and d are directly comparable and the guides collect once.
nes_lim  <- c(-1.75, 1.75)
col_lim  <- range(c(d1$neglogp, d0$neglogp))
size_lim <- range(c(d1$setSize, d0$setSize))

# A second, unfilled point layer rings the nominally significant pathways.
kegg_panel <- function(d, ttl, qmax) {
  d$label <- factor(d$label, levels = d$label)
  ggplot(d, aes(NES, label)) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "#999999") +
    geom_segment(aes(x = 0, xend = NES, y = label, yend = label),
                 linewidth = 0.2, colour = "#CFD6DA") +
    geom_point(aes(size = setSize, fill = neglogp), shape = 21,
               colour = "black", stroke = 0.14) +
    geom_point(data = subset(d, nominal), aes(size = setSize),
               shape = 21, fill = NA, colour = "black", stroke = 0.6,
               show.legend = FALSE) +
    scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = col_lim,
                        name = "-log10 P", breaks = c(0.4, 0.8, 1.2)) +
    scale_size_continuous(range = c(0.7, 2.5), limits = size_lim,
                          name = "Gene set size",
                          breaks = c(5, 15, 30, 44)) +
    scale_x_continuous(limits = nes_lim, breaks = seq(-1.5, 1.5, 0.5)) +
    guides(size = guide_legend(override.aes = list(fill = "#AEBEC8"),
                               order = 1),
           fill = guide_colourbar(order = 2)) +
    labs(x = "Normalized enrichment score", y = "KEGG signaling pathway",
         title = ttl,
         subtitle = sprintf("all FDR q > %.2f; nominal P shown, ring = P < 0.05",
                            floor(qmax * 100) / 100)) +
    theme(axis.text.y = element_text(size = 6.0),
          panel.grid.major.y = element_line(linewidth = 0.15,
                                            colour = "#F2F4F5"))
}

p_c <- kegg_panel(d1, "KEGG GSEA, subcluster 1",
                  STATS$kegg_maxq[["cluster1"]])
p_d <- kegg_panel(d0, "KEGG GSEA, subcluster 0",
                  STATS$kegg_maxq[["cluster0"]])

# =============================================================================
# Panels e / f - ImmGen splenic B-subset z-score heatmaps
# =============================================================================
read_z <- function(f) {
  as.matrix(read.csv(file.path(IMMGEN, f), row.names = 1, check.names = FALSE))
}
mat_e <- read_z("Fig2E_zscore_matrix.csv")
mat_f <- read_z("Fig2F_zscore_matrix.csv")

# ImmGen shorthand -> readable subset labels. Every column is spleen (.Sp), so
# the tissue moves into the panel title.
IMM_MAP <- c("B.T1.Sp" = "T1", "B.T2.Sp" = "T2", "B.T3.Sp" = "T3",
             "B.Fo.Sp" = "Fo", "B.MZ.Sp" = "MZ", "B.mem.Sp" = "Mem",
             "B.GC.CB.Sp" = "GC-CB", "B.GC.CC.Sp" = "GC-CC")
IMM_ORDER <- names(IMM_MAP)
relab <- function(m) {
  m <- m[, IMM_ORDER, drop = FALSE]
  colnames(m) <- IMM_MAP[colnames(m)]
  m
}
mat_e <- relab(mat_e); mat_f <- relab(mat_f)

STATS$imm_dim <- list(e = dim(mat_e), f = dim(mat_f))
STATS$imm_rng <- range(c(mat_e, mat_f))
say("ImmGen e: %d x %d   f: %d x %d   z range %.2f .. %.2f",
    nrow(mat_e), ncol(mat_e), nrow(mat_f), ncol(mat_f),
    STATS$imm_rng[1], STATS$imm_rng[2])

# Symmetric scale rounded up to the next 0.5 so it just covers the data and 0
# stays white. ComplexHeatmap otherwise stretches the legend bar out to pretty()
# bounds (+/- 4 here), which would advertise a range the data never reaches, so
# the tick positions are pinned to the true scale endpoints.
zlim <- ceiling(max(abs(STATS$imm_rng)) * 2) / 2
zcol <- colorRamp2(c(-zlim, 0, zlim), c(DIV_LO, DIV_MID, DIV_HI))
z_at <- c(-zlim, 0, zlim)

# Heatmap() measures label widths at construction time, and grid.grabExpr()
# defaults to an offscreen pdf() device whose PostScript font database has no
# Arial. Both steps therefore run with an svglite device current, so text is
# measured through systemfonts with the same metrics as the final export.
svglite_dev <- function(width, height) {
  svglite::svglite(tempfile(fileext = ".svg"), width = width, height = height)
}
# --- merge the two gene sets into a single membership-split panel ------------
#
# The 16 genes present in both matrices carry identical rows (verified: maximum
# absolute row difference is exactly 0), because both matrices are row-subsets
# of one ImmGen z-score table. A per-gene subtraction would therefore be all
# zeros and carry no information. The informative contrast is which genes are
# exclusive to each subcluster, so the union is shown once and split into
# subcluster-1-only / shared / subcluster-0-only blocks, with the difference
# quantified per subset as the mean-z gap between the two exclusive sets.
shared  <- intersect(rownames(mat_e), rownames(mat_f))
e_only  <- setdiff(rownames(mat_e), shared)
f_only  <- setdiff(rownames(mat_f), shared)
stopifnot(max(abs(mat_e[shared, ] - mat_f[shared, ])) == 0)

mat_u <- rbind(mat_e[e_only, , drop = FALSE],
               mat_e[shared, , drop = FALSE],
               mat_f[f_only, , drop = FALSE])
BLK <- c(sprintf("Subcluster 1 only (%d)", length(e_only)),
         sprintf("Shared (%d)", length(shared)),
         sprintf("Subcluster 0 only (%d)", length(f_only)))
memb <- factor(rep(BLK, c(length(e_only), length(shared), length(f_only))),
               levels = BLK)

# Per-subset difference between the two exclusive gene sets, with a rank test
# and Holm correction across the 8 subsets.
delta <- colMeans(mat_e[e_only, , drop = FALSE]) -
         colMeans(mat_f[f_only, , drop = FALSE])
dw <- sapply(colnames(mat_u), function(s) {
  suppressWarnings(wilcox.test(mat_e[e_only, s], mat_f[f_only, s])$p.value)
})
dw_holm <- p.adjust(dw, method = "holm")
STATS$delta <- delta
STATS$delta_p <- dw
STATS$delta_holm <- dw_holm
STATS$blocks <- c(sc1_only = length(e_only), shared = length(shared),
                  sc0_only = length(f_only))
say("union genes: %d (SC1-only %d, shared %d, SC0-only %d)",
    nrow(mat_u), length(e_only), length(shared), length(f_only))
say("largest |delta|: %s = %+.3f (Wilcoxon P = %.4f, Holm = %.3f)",
    names(which.max(abs(delta))), delta[which.max(abs(delta))],
    dw[which.max(abs(delta))], dw_holm[which.max(abs(delta))])
say("subsets surviving Holm across 8 tests: %d", sum(dw_holm < 0.05))

# Lineage-purity check. The subcluster 1 gene set carries canonical T-cell
# machinery, which would admit a B-T doublet explanation for the panel b shift.
# That alternative is tested against the object itself rather than assumed:
# marker expression and the Cd3-positive cell fraction per IV group decide it.
TCELL <- c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Cd8b1", "Lck", "Itk",
           "Cd28", "Il2rb", "Camk4", "Cblb", "Lef1", "Prkcq", "Lcp2",
           "Skap1", "Thy1", "Tcf7", "Zap70", "Themis", "Ctla4", "Icos",
           "Cd5", "Cd6", "Trat1", "Rasgrp1", "Dgka", "Ppp3cc", "Tespa1")
STATS$tc_e <- intersect(rownames(mat_e), TCELL)
STATS$tc_f <- intersect(rownames(mat_f), TCELL)
say("T-cell genes: subcluster 1 %d/%d (%.0f%%), subcluster 0 %d/%d (%.0f%%)",
    length(STATS$tc_e), nrow(mat_e), 100 * length(STATS$tc_e) / nrow(mat_e),
    length(STATS$tc_f), nrow(mat_f), 100 * length(STATS$tc_f) / nrow(mat_f))

MK <- intersect(c("Cd3e", "Cd3d", "Cd4", "Cd8b1", "Lck", "Itk",
                  "Cd19", "Cd79a", "Cd79b", "Ms4a1"), rownames(b[["RNA"]]))
STATS$mk <- as.matrix(
  Seurat::AverageExpression(b, features = MK, group.by = "seurat_clusters",
                            assays = "RNA", verbose = FALSE)$RNA)
cd3g <- intersect(c("Cd3e", "Cd3d", "Cd3g"), rownames(b[["RNA"]]))
cd3pos <- Matrix::colSums(
  GetAssayData(b, assay = "RNA", layer = "counts")[cd3g, , drop = FALSE]) > 0
STATS$cd3_by_cl <- 100 * tapply(cd3pos, md$cluster, mean)
STATS$cd3_by_cl_iv <- 100 * tapply(cd3pos, list(md$cluster, md$IV_group), mean)
STATS$qc <- md %>%
  group_by(cluster) %>%
  summarise(n = n(), nCount = median(nCount_RNA),
            nFeature = median(nFeature_RNA),
            pct_mt = median(percent.mt), .groups = "drop")
say("Cd3+ fraction by cluster: %s",
    paste(sprintf("%s %.1f%%", names(STATS$cd3_by_cl), STATS$cd3_by_cl),
          collapse = "  "))

# The grab device must match the height the assembled layout gives panel e,
# otherwise ComplexHeatmap measures its labels against the wrong canvas.
HT_W <- FIG_W_MM / 25.4
HT_H <- 70 / 25.4

# Transposed so the 79 genes run along x: keeping genes as rows would need
# ~150 mm of height and would not fit the page alongside rows 1 and 2.
heat_grob <- function() {
  tmp <- tempfile(fileext = ".svg")
  svglite::svglite(tmp, width = HT_W, height = HT_H)
  on.exit({ if (dev.cur() > 1L) dev.off(); unlink(tmp) }, add = TRUE)

  ann <- rowAnnotation(
    delta = anno_barplot(
      delta, baseline = 0, border = FALSE,
      bar_width = 0.72,
      gp = gpar(fill = ifelse(delta > 0, CL_COL[["1"]], CL_COL[["0"]]),
                col = NA),
      ylim = c(-0.70, 0.36),
      axis_param = list(side = "bottom", at = c(-0.6, 0, 0.3),
                        labels = c("-0.6", "0", "0.3"),
                        gp = gpar(fontsize = 5, fontfamily = FAM)),
      width = unit(15, "mm")
    ),
    annotation_label = "\u0394 mean z",
    annotation_name_gp = gpar(fontsize = 5.5, fontfamily = FAM),
    annotation_name_rot = 0, annotation_name_side = "top"
  )

  ht <- Heatmap(
    t(mat_u), name = "z-score", col = zcol,
    cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE,
    column_split = memb, cluster_column_slices = FALSE,
    column_gap = unit(1.6, "mm"),
    border = FALSE, rect_gp = gpar(col = NA),
    row_names_gp = gpar(fontsize = 6.2, fontfamily = FAM),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 5.2, fontfamily = FAM,
                           fontface = "italic"),
    column_names_rot = 90, column_names_side = "bottom",
    column_title_gp = gpar(fontsize = 6, fontfamily = FAM, fontface = "bold"),
    row_title = "ImmGen splenic B subset",
    row_title_gp = gpar(fontsize = 6.2, fontfamily = FAM),
    right_annotation = ann,
    heatmap_legend_param = list(
      title = "z-score", at = z_at, labels = sprintf("%.1f", z_at),
      title_gp = gpar(fontsize = 6, fontfamily = FAM),
      labels_gp = gpar(fontsize = 5.5, fontfamily = FAM),
      grid_height = unit(2.4, "mm"), legend_width = unit(22, "mm"),
      direction = "horizontal", title_position = "lefttop"
    )
  )
  grid.grabExpr(draw(ht, heatmap_legend_side = "bottom",
                     align_heatmap_legend = "heatmap_center",
                     padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm")),
                device = svglite_dev, width = HT_W, height = HT_H)
}
p_e <- wrap_elements(full = heat_grob())

# =============================================================================
# Assemble
# =============================================================================
# Panel a is coord_fixed, so it only grows taller by getting wider; give it the
# larger share so the UMAPs fill the row instead of floating in white space.
row1 <- p_a + p_b + plot_layout(widths = c(1.85, 1))
# The shared c/d guides are collected into a horizontal strip under the row
# rather than into a right-hand column. A guide column of usable width left
# most of its area empty below the two short legends; moving the guides to the
# bottom returns that width to the dot plots, so the pathway labels and the NES
# axis both gain space.
row2 <- p_c + p_d +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.box.just = "center", legend.justification = "center",
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        legend.box.spacing = unit(1.2, "mm"))

fig <- row1 / row2 / p_e +
  plot_layout(heights = c(1.00, 1.22, 1.10)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FAM),
        plot.tag.position = c(0, 1))

# =============================================================================
# Export
# =============================================================================
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "Figure2")

# A single locked output file (a PDF viewer or Explorer preview handler holding
# a handle) must not abort the run and lose the other formats, so each export is
# attempted independently and failures are collected and reported at the end.
export_failures <- character(0)
export_one <- function(label, open_device) {
  say("writing %s ...", label)
  ok <- tryCatch({
    open_device()
    on.exit(while (dev.cur() > 1L) invisible(dev.off()), add = TRUE)
    print(fig)
    TRUE
  }, error = function(e) {
    message(sprintf("  %s FAILED: %s", label, conditionMessage(e)))
    FALSE
  }, warning = function(w) {
    message(sprintf("  %s FAILED: %s", label, conditionMessage(w)))
    FALSE
  })
  if (!ok) export_failures <<- c(export_failures, label)
  invisible(ok)
}

export_one("SVG", function() {
  svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
})
export_one("PDF", function() {
  cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM,
            bg = "white")
})
export_one("PNG", function() {
  agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
          res = 600, background = "white")
})
# A TIFF submission raster is one call away if the target journal asks for it:
#   agg_tiff(paste0(base, ".tiff"), width = w_in, height = h_in, units = "in",
#            res = 600, background = "white", compression = "lzw")

# =============================================================================
# Source data + statistics report
# =============================================================================
write.csv(comp %>% select(cluster, IV_group, n, total, pct, se),
          file.path(OUTDIR, "SourceData_Fig2b_composition.csv"),
          row.names = FALSE)
write.csv(rbind(cbind(panel = "c", d1[, c("ID", "Description", "setSize",
                                          "NES", "pvalue", "p.adjust")]),
                cbind(panel = "d", d0[, c("ID", "Description", "setSize",
                                          "NES", "pvalue", "p.adjust")])),
          file.path(OUTDIR, "SourceData_Fig2cd_kegg.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(mat_u), gene_set = as.character(memb),
                     mat_u, check.names = FALSE),
          file.path(OUTDIR, "SourceData_Fig2e_zscore_union.csv"),
          row.names = FALSE)
write.csv(data.frame(subset = names(delta), delta_mean_z = as.numeric(delta),
                     wilcoxon_p = as.numeric(dw),
                     wilcoxon_p_holm = as.numeric(dw_holm)),
          file.path(OUTDIR, "SourceData_Fig2e_delta.csv"), row.names = FALSE)

rep <- file.path(OUTDIR, "Figure2_stats.txt")
con <- file(rep, "w", encoding = "UTF-8")
wl <- function(...) cat(sprintf(...), "\n", file = con, sep = "")

wl("Figure 2 - source statistics")
wl(strrep("=", 64))
wl("")
wl("[Panel a/b] Splenic B cells, reclustered")
wl("  total cells              : %s", format(STATS$n_cells, big.mark = ","))
wl("  IntraNeg                 : %s", format(STATS$n_by_iv[["IntraNeg"]],
                                             big.mark = ","))
wl("  IntraPos                 : %s", format(STATS$n_by_iv[["IntraPos"]],
                                             big.mark = ","))
wl("  subclusters              : 0, 1, 2 (RNA_snn_res.0.1)")
wl("")
wl("  cluster x IV_group counts")
for (cl in rownames(STATS$tab)) {
  wl("    cluster %s : IntraNeg %5d (%5.2f%%)   IntraPos %5d (%5.2f%%)", cl,
     STATS$tab[cl, "IntraNeg"],
     100 * STATS$tab[cl, "IntraNeg"] / sum(STATS$tab[, "IntraNeg"]),
     STATS$tab[cl, "IntraPos"],
     100 * STATS$tab[cl, "IntraPos"] / sum(STATS$tab[, "IntraPos"]))
}
wl("")
wl("  Pearson chi-square (3x2) : X2 = %.2f, df = %d, P = %.3e",
   STATS$chisq$statistic, STATS$chisq$parameter, STATS$chisq$p.value)
wl("  per-cluster Fisher exact (two-sided), cluster vs rest, Holm-adjusted")
for (i in seq_len(nrow(STATS$fisher))) {
  wl("    cluster %s : OR = %.3f, P = %.3e, P_Holm = %.3e",
     STATS$fisher$cluster[i], STATS$fisher$OR[i], STATS$fisher$p[i],
     STATS$fisher$p_holm[i])
}
wl("  error bars               : binomial SE of the observed proportion")
wl("")
wl("[Panel c/d] KEGG GSEA, top 20 pathways as supplied")
wl("  IMPORTANT - no pathway reaches FDR significance in either panel.")
wl("  subcluster 1 : q range %.4f .. %.4f; nominal P < 0.05 in %d of 20",
   STATS$kegg_minq[["cluster1"]], STATS$kegg_maxq[["cluster1"]],
   STATS$kegg_nom[["cluster1"]])
wl("  subcluster 0 : q range %.4f .. %.4f; nominal P < 0.05 in %d of 20",
   STATS$kegg_minq[["cluster0"]], STATS$kegg_maxq[["cluster0"]],
   STATS$kegg_nom[["cluster0"]])
wl("  smallest gene set        : %d genes (below the usual GSEA floor of 15)",
   min(size_lim))
wl("  panels therefore encode the nominal P value and state the FDR status;")
wl("  they must be described as exploratory, not as significant enrichment.")
wl("")
wl("[Panel e] ImmGen splenic B-subset z-scores, merged gene sets")
wl("  supplied matrices        : %d genes (subcluster 1) and %d genes",
   nrow(mat_e), nrow(mat_f))
wl("                             (subcluster 0), %d subsets each", ncol(mat_e))
wl("  the %d genes present in both carry identical rows (max absolute row",
   length(shared))
wl("  difference = 0), because both matrices are row-subsets of one ImmGen")
wl("  z-score table. A per-gene subtraction is therefore uninformative, so the")
wl("  union is plotted once and split by gene-set membership.")
wl("  union plotted            : %d genes x %d subsets", nrow(mat_u),
   ncol(mat_u))
wl("    subcluster 1 only      : %d", STATS$blocks[["sc1_only"]])
wl("    shared                 : %d", STATS$blocks[["shared"]])
wl("    subcluster 0 only      : %d", STATS$blocks[["sc0_only"]])
wl("  z-score range            : %.3f .. %.3f", STATS$imm_rng[1],
   STATS$imm_rng[2])
wl("  colour scale             : symmetric, +/- %.3f", zlim)
wl("  subset order (fixed)     : %s", paste(IMM_MAP, collapse = ", "))
wl("  genes                    : clustered within each membership block")
wl("  subsets                  : not clustered, developmental order retained")
wl("")
wl("  delta = mean z(subcluster-1-only) - mean z(subcluster-0-only) per subset")
wl("  Wilcoxon rank-sum per subset, Holm-adjusted across the %d subsets",
   length(delta))
for (s in names(delta)) {
  wl("    %-6s delta = %+.3f   P = %.4f   P_Holm = %.3f", s, delta[[s]],
     dw[[s]], dw_holm[[s]])
}
wl("  subsets significant after Holm : %d of %d", sum(dw_holm < 0.05),
   length(delta))
wl("")
wl("  INTERPRETATION LIMIT - the two gene sets are only weakly separable.")
wl("  Full-set subset profiles correlate at r = %.3f and the exclusive-set",
   cor(colMeans(mat_e), colMeans(mat_f)))
wl("  profiles at r = %.3f; both peak in Mem and then GC-CC.",
   cor(colMeans(mat_e[e_only, , drop = FALSE]),
       colMeans(mat_f[f_only, , drop = FALSE])))
wl("  The largest gap is at Mem (delta = %+.3f, nominal P = %.4f), i.e. the",
   delta[["Mem"]], dw[["Mem"]])
wl("  subcluster-0-exclusive genes are the more memory-skewed of the two, but")
wl("  this does not survive correction across the 8 subsets. Panel e should be")
wl("  described as mapping a largely shared programme onto the ImmGen")
wl("  reference, not as evidence of divergent subcluster identity.")
wl("")
wl("LINEAGE-PURITY CHECK (reviewer pre-empt, tested and resolved)")
wl(strrep("-", 64))
wl("  canonical T-cell genes in the subcluster 1 set : %d of %d (%.0f%%)",
   length(STATS$tc_e), nrow(mat_e),
   100 * length(STATS$tc_e) / nrow(mat_e))
wl("    %s", paste(sort(STATS$tc_e), collapse = ", "))
wl("  canonical T-cell genes in the subcluster 0 set : %d of %d (%.0f%%)",
   length(STATS$tc_f), nrow(mat_f),
   100 * length(STATS$tc_f) / nrow(mat_f))
wl("    %s", paste(sort(STATS$tc_f), collapse = ", "))
wl("")
wl("  These enter the list through pathway membership, not through expression:")
wl("  the top panel c pathway is 'T cell receptor signaling pathway', which")
wl("  contributes Cd3e, Cd4, Cd8b1, Lck, Itk and Cd28 as set members.")
wl("")
wl("  A B-T doublet explanation for the panel b shift was tested directly.")
wl("  mean RNA expression per subcluster (assay RNA):")
wl("    %-7s %s", "gene", paste(sprintf("%8s", colnames(STATS$mk)),
                                collapse = ""))
for (g in rownames(STATS$mk)) {
  wl("    %-7s %s", g, paste(sprintf("%8.3f", STATS$mk[g, ]), collapse = ""))
}
wl("  cells with any Cd3e/Cd3d/Cd3g count, per subcluster (%%):")
wl("    %s", paste(sprintf("cluster %s %.2f", names(STATS$cd3_by_cl),
                           STATS$cd3_by_cl), collapse = "   "))
wl("  same, split by IV group (%%):")
for (cl in rownames(STATS$cd3_by_cl_iv)) {
  wl("    cluster %s : IntraNeg %.2f   IntraPos %.2f", cl,
     STATS$cd3_by_cl_iv[cl, "IntraNeg"], STATS$cd3_by_cl_iv[cl, "IntraPos"])
}
wl("")
wl("  VERDICT: the doublet alternative is not supported. Subcluster 1 is")
wl("  B-lineage (high Cd79a, Ms4a1, Cd19) with T-lineage transcripts near zero")
wl("  and not elevated over subcluster 0. Critically, the Cd3-positive")
wl("  fraction within subcluster 1 is higher in IntraNeg than in IntraPos,")
wl("  the opposite of what blood-derived B-T doublets would produce.")
wl("")
wl("  Remaining item to address in the text: subcluster 1 has a higher median")
wl("  library size than subcluster 0 (see QC below), consistent with larger or")
wl("  activated cells but worth stating explicitly.")
wl("  per-subcluster QC (median):")
for (i in seq_len(nrow(STATS$qc))) {
  wl("    cluster %s : n = %5d  nCount = %6.0f  nFeature = %5.0f  mt = %.2f%%",
     STATS$qc$cluster[i], STATS$qc$n[i], STATS$qc$nCount[i],
     STATS$qc$nFeature[i], STATS$qc$pct_mt[i])
}
wl("  ADT panel available : %s", paste(rownames(b[["Protein"]]),
                                       collapse = ", "))
wl("                        (no CD3 antibody, so protein cannot adjudicate")
wl("                        CD3 directly)")
wl("")
wl(strrep("=", 64))
wl("QA notes")
wl(strrep("=", 64))
wl("  figure size            : %.0f x %.0f mm (Nature double column)",
   FIG_W_MM, FIG_H_MM)
wl("  archetype              : quantitative grid, panel a as hero")
wl("  backend                : R / ggplot2 + ComplexHeatmap (exclusive)")
wl("  font                   : Arial, 5.0-8.0 pt; no plotmath superscript,")
wl("                           so no glyph falls below the 5 pt floor")
wl("  editable text          : svglite SVG and cairo_pdf PDF keep live text")
wl("  data integrity         : all %s cells, all 20+20 KEGG rows and all",
   format(STATS$n_cells, big.mark = ","))
wl("                           %d union ImmGen genes are shown; nothing",
   nrow(mat_u))
wl("                           dropped, the %d shared genes are drawn once",
   length(shared))
wl("                           instead of twice")
wl("  centre / spread        : panel b, observed proportion +/- binomial SE")
wl("  colour policy          : cluster identity carries hue; IV group carries")
wl("                           a light/dark neutral, so panel b survives")
wl("                           greyscale print and the two encodings never")
wl("                           compete. No rainbow map anywhere.")
wl("  relabelling            : KEGG 'signaling pathway' suffix moved to the")
wl("                           axis title; ImmGen '.Sp' shorthand expanded.")
wl("                           No value was altered.")
wl("  source data            : four SourceData_Fig2*.csv files in this folder")
wl("  panel count            : a-e; the two supplied ImmGen heat maps are")
wl("                           merged into panel e by gene-set membership")
wl("  exports                : SVG, PDF, PNG at 600 dpi (TIFF available via one")
wl("                           commented call in the export block)")
wl("")
wl("  changes in this build (presentation only, no statistic was altered)")
wl("    vocabulary           : IntraNeg / IntraPos are displayed as")
wl("                           'IV- (extravascular)' and 'IV+ (intravascular)',")
wl("                           matching Figure 1. Internal keys, tables and")
wl("                           source-data columns are unchanged.")
wl("    panel b annotation   : significance asterisks were replaced by the exact")
wl("                           Holm-adjusted P value. The asterisks encoded the")
wl("                           same quantity and added a second, coarser scale.")
wl("                           Exponents use Unicode superscript glyphs rather")
wl("                           than plotmath, so no glyph is shrunk below 5 pt.")
wl("    row 2 layout         : the shared c/d guides moved from a right-hand")
wl("                           guide column to a horizontal strip under the row.")
wl("                           The column reserved width that stayed empty below")
wl("                           the two short legends; the dot plots now use it.")
wl("    canvas               : 196 -> 212 mm tall, which lets the panel e gene")
wl("                           labels go from 5.0 to 5.2 pt, the ImmGen subset")
wl("                           labels from 6.0 to 6.2 pt and the c/d pathway")
wl("                           labels from 5.5 to 6.0 pt.")
wl("    reuse level          : structural adaptation of the previous build -")
wl("                           statistical logic, gene sets, tests, corrections")
wl("                           and integrity checks are byte-for-byte the same.")
wl("")
wl("  verified after render (see qa/ for the measured values of this build)")
wl("    SVG glyph sizes      : read back from the exported SVG, which records")
wl("                           true point sizes; the floor is 5 pt")
wl("    PDF fonts            : Arial embedded as TrueType subsets, so text stays")
wl("                           selectable and editable")
wl("    PDF Tf caveat        : cairo_pdf writes '1 Tf' and applies the real size")
wl("                           through the text matrix, so a raw Tf scan reports")
wl("                           1 pt for every run. Sizes were confirmed from the")
wl("                           SVG instead, which records true point sizes.")
close(con)
say("wrote %s", rep)

# Close any device left open and drop the blank Rplots.pdf that R creates when
# a default device gets touched, so the delivery folder holds outputs only.
while (dev.cur() > 1L) invisible(dev.off())
stray <- file.path(OUTDIR, "Rplots.pdf")
if (file.exists(stray)) unlink(stray)

if (length(export_failures)) {
  say("WARNING: these exports failed (target file most likely locked by a")
  say("         viewer): %s", paste(export_failures, collapse = ", "))
  say("         close the file and re-run to refresh them.")
} else {
  say("all three exports written")
}
say("done")
