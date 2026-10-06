#!/usr/bin/env Rscript
# =============================================================================
# Figure S2 — Nature quantitative grid (R exclusive, submission-ready)
#
# Source: FigS2_NC_submission logic; NEW T-cell resolution that splits former
#   clusters 1 and 3 into two populations each (panels c–h).
# Panels a–b: unchanged B-cell Cr2 / Fcer2a from Bcell_reclustered.rds
#
# Core conclusion
#   IntraPos T cells (CD4_T annotation) show composition shifts among refined
#   subclusters after splitting prior clusters 1 (Maf-high) and 3 (Cd160-high)
#   into two groups each. B-cell IntraPos-enriched subcluster 1 remains
#   Cr2-high / Fcer2a-low (a–b).
#
# Panel map
#   a  Cr2 marks B-cell subcluster 1
#   b  Fcer2a marks B-cell subcluster 0
#   c  T UMAP by IV group with refined clusters (0, 1a, 1b, 2, 3a, 3b)
#   d  Composition bars + Fisher/Holm (hero metric)
#   e–h  Tcf7 / Bcl6 / Maf / Cd160 feature+violin on refined clusters
#
# Archetype : quantitative grid (panel d hero)
# Backend   : R native ggplot2 vector
# Size      : 180 x 230 mm
# Exports   : PDF / SVG / PNG (600 dpi)
# New folder — does not overwrite FigS2_NC_submission
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
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
options(future.globals.maxSize = 8 * 1024^3)

FIGDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures")
OUTDIR <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS2_Nature_split_c1c3")
dir.create(file.path(OUTDIR, "scripts"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "source_data"), showWarnings = FALSE, recursive = TRUE)

FIG_W_MM <- 180
FIG_H_MM <- 230
BASE <- 6.2
FAM <- "Arial"
fm <- systemfonts::match_fonts(FAM)
cat("font:", FAM, "->", fm$path, "\n")

for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM
    do.call(reg, args)
  }
}
PT <- function(x) x / .pt

B_COL <- c("0" = "#4C6E8A", "1" = "#D1495B", "2" = "#E9A03B")
# Refined T palette: 1a/1b warm (ex-1); 3a/3b teal (ex-3)
T_COL <- c(
  "0"  = "#3D5A80",
  "1a" = "#E07A5F",
  "1b" = "#F2B6A0",
  "2"  = "#9AA0A6",
  "3a" = "#2A9D8F",
  "3b" = "#7BC9BE"
)
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")
SEQ_LO <- "#EAEFF2"; SEQ_HI <- "#26506B"
T_LEVELS <- c("0", "1a", "1b", "2", "3a", "3b")

theme_s2 <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.1, "pt"),
      axis.title = element_text(size = base),
      axis.text = element_text(size = base - 0.4, colour = "black"),
      legend.title = element_text(size = base - 0.2),
      legend.text = element_text(size = base - 0.6),
      legend.key.size = unit(2.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_blank(),
      strip.text = element_text(size = base, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = base + 0.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base - 0.6, hjust = 0.5, colour = "#555555"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      panel.grid = element_blank()
    )
}
theme_set(theme_s2())

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
STATS <- list()
say <- function(...) cat(sprintf(...), "\n")

# =============================================================================
# Split former clusters 1 and 3 into two groups each
# =============================================================================
bipartition_cluster <- function(obj, cl, gene_rank, res = 0.1, prefix) {
  Idents(obj) <- "seurat_clusters_old"
  sub <- subset(obj, idents = cl)
  DefaultAssay(sub) <- "RNA"
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
  dims <- 1:min(15, ncol(Embeddings(sub, "pca")))
  sub <- FindNeighbors(sub, dims = dims, verbose = FALSE)
  sub <- FindClusters(sub, resolution = res, verbose = FALSE)
  tab <- sort(table(Idents(sub)), decreasing = TRUE)
  say("  old%s @ res=%.2f raw: %s", cl, res,
      paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "))
  # Keep two largest; merge remainder into nearest large by PCA centroid
  keep <- names(tab)[seq_len(min(2L, length(tab)))]
  if (length(tab) > 2L) {
    emb <- Embeddings(sub, "pca")[, dims, drop = FALSE]
    cents <- lapply(keep, function(k) colMeans(emb[Idents(sub) == k, , drop = FALSE]))
    names(cents) <- keep
    new_id <- as.character(Idents(sub))
    for (k in setdiff(names(tab), keep)) {
      ck <- colMeans(emb[Idents(sub) == k, , drop = FALSE])
      d1 <- sum((ck - cents[[keep[1]]])^2)
      d2 <- sum((ck - cents[[keep[2]]])^2)
      new_id[Idents(sub) == k] <- if (d1 <= d2) keep[1] else keep[2]
    }
    Idents(sub) <- factor(new_id)
  }
  # Order by ranking gene mean: higher -> prefix+"a", lower -> prefix+"b"
  means <- sapply(levels(Idents(sub)), function(k) {
    mean(FetchData(sub, vars = gene_rank)[Idents(sub) == k, 1])
  })
  ord <- names(sort(means, decreasing = TRUE))
  map <- setNames(paste0(prefix, c("a", "b")), ord[1:2])
  out <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
  lab <- as.character(Idents(sub))
  out[colnames(sub)] <- unname(map[lab])
  say("  old%s -> %s (by %s: %s)", cl, paste(unique(out[!is.na(out)]), collapse = "/"),
      gene_rank,
      paste(sprintf("%s mean=%.3f", names(means), means), collapse = "; "))
  out
}

say("loading objects ...")
b  <- readRDS(file.path(FIGDIR, "Bcell_reclustered.rds"))
t4 <- readRDS(file.path(FIGDIR, "CD4T_reclustered.rds"))
t4$seurat_clusters_old <- as.character(t4$seurat_clusters)
Idents(t4) <- "seurat_clusters_old"
say("old T clusters: %s", paste(sprintf("%s=%d", names(table(t4$seurat_clusters_old)),
                                         as.integer(table(t4$seurat_clusters_old))), collapse = ", "))

say("splitting old cluster 1 (rank by Maf) ...")
lab1 <- bipartition_cluster(t4, "1", "Maf", res = 0.1, prefix = "1")
say("splitting old cluster 3 (rank by Cd160) ...")
lab3 <- bipartition_cluster(t4, "3", "Cd160", res = 0.1, prefix = "3")

new_lab <- t4$seurat_clusters_old
new_lab[!is.na(lab1)] <- lab1[!is.na(lab1)]
new_lab[!is.na(lab3)] <- lab3[!is.na(lab3)]
# old 0 and 2 stay "0" and "2"
t4$seurat_clusters <- factor(new_lab, levels = T_LEVELS)
Idents(t4) <- "seurat_clusters"
say("new T clusters: %s", paste(sprintf("%s=%d", names(table(t4$seurat_clusters)),
                                         as.integer(table(t4$seurat_clusters))), collapse = ", "))
saveRDS(t4, file.path(OUTDIR, "source_data", "CD4T_reclustered_split_c1c3.rds"))
write.csv(
  data.frame(cell = colnames(t4),
             cluster_old = t4$seurat_clusters_old,
             cluster_new = as.character(t4$seurat_clusters),
             IV_group = as.character(t4$IV_group)),
  file.path(OUTDIR, "source_data", "CD4T_cluster_map_old_to_new.csv"),
  row.names = FALSE
)

Idents(b) <- "seurat_clusters"
prep_cells <- function(obj, cl_levels, cluster_col = "seurat_clusters") {
  md <- obj@meta.data
  md$cluster <- factor(as.character(md[[cluster_col]]), levels = cl_levels)
  md$IV_group <- factor(md$IV_group, levels = c("IntraNeg", "IntraPos"))
  emb <- as.data.frame(Embeddings(obj, "umap"))
  colnames(emb) <- c("UMAP1", "UMAP2")
  cbind(emb, md[, c("cluster", "IV_group")])
}

b_df <- prep_cells(b, c("0", "1", "2"))
t_df <- prep_cells(t4, T_LEVELS)
STATS$n_b <- nrow(b_df)
STATS$n_t <- nrow(t_df)
STATS$n_t_iv <- table(t_df$IV_group)
STATS$tab_t <- table(t_df$cluster, t_df$IV_group)
STATS$n_b_iv <- table(b_df$IV_group)

expr_vec <- function(obj, gene) {
  as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[gene, ])
}
pct_pos <- function(obj, gene, cluster_col = "seurat_clusters") {
  cnt <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "counts")[gene, ])
  100 * tapply(cnt > 0, factor(as.character(obj@meta.data[[cluster_col]]), levels = T_LEVELS), mean)
}

LIN <- c("Cd3e", "Cd4", "Cd8b1", "Ms4a1", "Klrb1c")
STATS$lin_pct <- sapply(LIN, function(g) pct_pos(t4, g))
STATS$n_t_cl <- as.numeric(table(t_df$cluster))
names(STATS$n_t_cl) <- levels(t_df$cluster)
say("T lineage %% positive:")
print(round(STATS$lin_pct, 1))

# =============================================================================
# Feature + violin
# =============================================================================
marker_stats <- function(df, gene) {
  kw <- kruskal.test(expr ~ cluster, df)
  means <- tapply(df$expr, df$cluster, mean)
  peak <- names(which.max(means))
  wt <- suppressWarnings(
    wilcox.test(df$expr[df$cluster == peak], df$expr[df$cluster != peak])
  )
  list(gene = gene, kw_p = kw$p.value, peak = peak,
       peak_mean = unname(means[[peak]]), rest_mean = mean(df$expr[df$cluster != peak]),
       peak_p = wt$p.value, n = nrow(df),
       n_cl = as.numeric(table(df$cluster)),
       mean_cl = as.numeric(means),
       pct_pos = 100 * tapply(df$expr > 0, df$cluster, mean))
}

feat_vln <- function(df, gene, pal, y_lab = "Expression", vln_sub = NULL) {
  zmax <- max(as.numeric(quantile(df$expr, 0.99, names = FALSE)), 0.15)
  st <- marker_stats(df, gene)
  if (is.null(vln_sub)) {
    vln_sub <- sprintf("highest in %s, vs rest %s", st$peak, fmt_p(st$peak_p))
  }
  p_f <- ggplot(df, aes(UMAP1, UMAP2, colour = expr)) +
    geom_point(size = 0.13, alpha = 0.78, shape = 16) +
    scale_colour_gradientn(
      colours = c(SEQ_LO, SEQ_HI), limits = c(0, zmax), oob = squish,
      name = NULL, breaks = c(0, zmax), labels = c("0", sprintf("%.1f", zmax)),
      guide = guide_colorbar(barwidth = unit(2.2, "mm"), barheight = unit(16, "mm"),
                             ticks.colour = "black", frame.colour = "black",
                             frame.linewidth = 0.2)
    ) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", title = gene) +
    theme(plot.title = element_text(face = "italic", size = BASE + 0.6),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
          legend.position = "right",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -1, "mm"))
  p_v <- ggplot(df, aes(cluster, expr, fill = cluster)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "black", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = "Cluster", y = y_lab, subtitle = vln_sub) +
    theme(plot.subtitle = element_text(size = BASE - 0.8, hjust = 0.5, colour = "#555555"),
          axis.text.x = element_text(size = BASE - 0.8))
  list(plot = wrap_elements(full = p_f + p_v + plot_layout(widths = c(1.28, 0.72))),
       stats = st)
}

# =============================================================================
# Panels a / b
# =============================================================================
b_cr2 <- b_df; b_cr2$expr <- expr_vec(b, "Cr2")
b_fc  <- b_df; b_fc$expr  <- expr_vec(b, "Fcer2a")
pa <- feat_vln(b_cr2, "Cr2", B_COL)
pb <- feat_vln(b_fc, "Fcer2a", B_COL)
STATS$a <- pa$stats; STATS$b <- pb$stats

# =============================================================================
# Panel c
# =============================================================================
lab_pos <- t_df %>% group_by(cluster) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = "drop")
facet_n <- t_df %>% count(IV_group) %>%
  mutate(lab = sprintf("n = %s", format(n, big.mark = ",")))
xr <- range(t_df$UMAP1); yr <- range(t_df$UMAP2)
pad <- 0.06 * diff(xr)

p_c <- ggplot(t_df, aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 0.15, alpha = 0.70, shape = 16) +
  geom_text(data = lab_pos, aes(x, y, label = cluster), inherit.aes = FALSE,
            family = FAM, fontface = "bold", size = PT(5.8), colour = "black") +
  geom_text(data = facet_n, aes(x = xr[1], y = yr[2] + pad * 0.35, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 1, family = FAM,
            size = PT(5.4), colour = "#444444") +
  facet_wrap(~ IV_group, nrow = 1) +
  scale_colour_manual(values = T_COL, guide = "none") +
  coord_fixed(xlim = c(xr[1] - pad, xr[2] + pad),
              ylim = c(yr[1] - pad, yr[2] + pad * 1.55), expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "T-cell subclusters by intravascular labelling",
       subtitle = "CD4_T annotation; former clusters 1 and 3 each split into a/b. RNA lineage note on integrity.") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
        axis.line = element_blank(),
        panel.spacing = unit(1.6, "mm"),
        plot.subtitle = element_text(size = 5.2, hjust = 0.5, colour = "#555555"))

# =============================================================================
# Panel d
# =============================================================================
comp <- as.data.frame(STATS$tab_t)
colnames(comp) <- c("cluster", "IV_group", "n")
comp <- comp %>%
  group_by(IV_group) %>%
  mutate(total = sum(n), pct = 100 * n / total,
         se = 100 * sqrt(pct / 100 * (1 - pct / 100) / total)) %>%
  ungroup()

fish <- lapply(levels(comp$cluster), function(cl) {
  m <- matrix(c(STATS$tab_t[cl, "IntraPos"],
                sum(STATS$tab_t[, "IntraPos"]) - STATS$tab_t[cl, "IntraPos"],
                STATS$tab_t[cl, "IntraNeg"],
                sum(STATS$tab_t[, "IntraNeg"]) - STATS$tab_t[cl, "IntraNeg"]),
              nrow = 2)
  ft <- fisher.test(m)
  data.frame(cluster = cl, OR = unname(ft$estimate), p = ft$p.value,
             n_neg = STATS$tab_t[cl, "IntraNeg"],
             n_pos = STATS$tab_t[cl, "IntraPos"])
}) %>% bind_rows()
fish$p_holm <- p.adjust(fish$p, method = "holm")
fish$cluster <- factor(fish$cluster, levels = levels(comp$cluster))
STATS$fisher <- fish
chi <- chisq.test(STATS$tab_t)
STATS$chisq <- chi
say("chi-square: X2=%.1f df=%d %s", unname(chi$statistic), chi$parameter, fmt_p(chi$p.value))
print(fish)

Y_STAR <- 90; Y_OR <- 96; ymax_d <- 102
p_d <- ggplot(comp, aes(cluster, pct, fill = IV_group)) +
  geom_col(aes(colour = cluster), position = position_dodge(width = 0.72),
           width = 0.64, linewidth = 0.45) +
  scale_colour_manual(values = T_COL, guide = "none") +
  geom_errorbar(aes(ymin = pmax(pct - se, 0), ymax = pct + se),
                position = position_dodge(width = 0.72), width = 0.18, linewidth = 0.28) +
  geom_text(aes(y = pct + se + 1.8, label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.72), family = FAM,
            size = PT(4.8), colour = "black") +
  geom_text(data = fish, aes(x = cluster, y = Y_STAR, label = pstar(p_holm)),
            inherit.aes = FALSE, family = FAM, size = PT(5.6), vjust = 0.35) +
  geom_text(data = fish, aes(x = cluster, y = Y_OR, label = sprintf("OR %.2f", OR)),
            inherit.aes = FALSE, family = FAM, size = PT(4.6), colour = "#444444") +
  scale_fill_manual(values = IV_COL, name = NULL,
                    labels = c(sprintf("IntraNeg (n = %s)", format(STATS$n_t_iv[["IntraNeg"]], big.mark = ",")),
                               sprintf("IntraPos (n = %s)", format(STATS$n_t_iv[["IntraPos"]], big.mark = ",")))) +
  scale_y_continuous(limits = c(0, ymax_d), breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "T-cell subcluster", y = "Cells in group (%)",
       title = "Subcluster composition",
       subtitle = sprintf("chi-square %s; per-cluster Fisher, Holm-adjusted", fmt_p(chi$p.value))) +
  theme(axis.text.x = element_text(face = "bold", size = BASE - 0.2),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.3, "mm"),
        legend.margin = margin(t = -1.5, b = 0, unit = "mm"))

# =============================================================================
# Panels e–h
# =============================================================================
mk_t <- list(
  e = list(gene = "Tcf7"),
  f = list(gene = "Bcl6"),
  g = list(gene = "Maf"),
  h = list(gene = "Cd160")
)
built <- list()
for (nm in names(mk_t)) {
  d <- t_df
  d$expr <- expr_vec(t4, mk_t[[nm]]$gene)
  vln_sub <- if (mk_t[[nm]]$gene == "Bcl6") {
    "weakly higher in 1a/1b; modest restriction"
  } else NULL
  built[[nm]] <- feat_vln(d, mk_t[[nm]]$gene, T_COL, vln_sub = vln_sub)
  STATS[[nm]] <- built[[nm]]$stats
  say("%s peak %s %s", mk_t[[nm]]$gene, STATS[[nm]]$peak, fmt_p(STATS[[nm]]$peak_p))
}

# =============================================================================
# Assemble + export
# =============================================================================
fig <- (pa$plot | pb$plot) /
  p_c /
  p_d /
  (built$e$plot | built$f$plot) /
  (built$g$plot | built$h$plot) +
  plot_layout(heights = c(1.00, 1.12, 0.95, 1.00, 1.00)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FAM),
        plot.tag.position = c(0.01, 0.98))

base <- file.path(OUTDIR, "FigureS2")
w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4

svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM, bg = "white")
print(fig); invisible(dev.off())
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())

# Source data
write.csv(comp, file.path(OUTDIR, "source_data", "SourceData_FigS2d_composition.csv"), row.names = FALSE)
write.csv(comp, file.path(OUTDIR, "SourceData_FigS2d_composition.csv"), row.names = FALSE)
write.csv(fish, file.path(OUTDIR, "source_data", "SourceData_FigS2d_fisher.csv"), row.names = FALSE)

summ_row <- function(st, panel, cell) {
  mc <- st$mean_cl
  names(mc) <- levels(t_df$cluster)[seq_along(mc)]
  data.frame(
    panel = panel, cell = cell, gene = st$gene,
    peak_cluster = st$peak, peak_vs_rest_p = st$peak_p, kruskal_p = st$kw_p,
    t(mc), stringsAsFactors = FALSE, check.names = FALSE
  )
}
# B panels use 3 clusters
summ_b <- function(st, panel) {
  data.frame(panel = panel, cell = "B", gene = st$gene,
             peak_cluster = st$peak, peak_vs_rest_p = st$peak_p, kruskal_p = st$kw_p,
             mean_cl0 = st$mean_cl[1], mean_cl1 = st$mean_cl[2], mean_cl2 = st$mean_cl[3],
             stringsAsFactors = FALSE)
}
marker_summ <- dplyr::bind_rows(
  summ_b(STATS$a, "a"),
  summ_b(STATS$b, "b"),
  summ_row(STATS$e, "e", "T"),
  summ_row(STATS$f, "f", "T"),
  summ_row(STATS$g, "g", "T"),
  summ_row(STATS$h, "h", "T")
)
write.csv(marker_summ, file.path(OUTDIR, "source_data", "SourceData_FigS2_markers.csv"), row.names = FALSE)
write.csv(marker_summ, file.path(OUTDIR, "SourceData_FigS2_markers.csv"), row.names = FALSE)

lin_df <- as.data.frame(STATS$lin_pct)
lin_df$cluster <- rownames(lin_df)
write.csv(lin_df, file.path(OUTDIR, "source_data", "SourceData_FigS2c_lineage_pct.csv"), row.names = FALSE)

writeLines(c(
  "Fig. S2 | B- and T-cell subcluster markers and IntraPos/IntraNeg composition.",
  "",
  "a,b, Feature and violin plots for Cr2 and Fcer2a in reclustered B cells.",
  "c, UMAP of the CD4_T annotation split by intravascular labelling, after increasing",
  "resolution so that former clusters 1 and 3 each form two groups (1a/1b; 3a/3b).",
  "Clusters 0 and 2 are retained. 1a ranked by higher Maf; 3a by higher Cd160.",
  "d, Subcluster composition (% of cells within each IV group) with binomial SE;",
  "stars, per-cluster Fisher exact tests (Holm-adjusted); OR shown.",
  "e–h, Tcf7, Bcl6, Maf and Cd160 feature and violin plots on the refined clusters.",
  "Source data: SourceData_FigS2*.csv."
), file.path(OUTDIR, "FigureS2_legend.txt"))

writeLines(c(
  "Core conclusion: IntraPos T cells shift among refined subclusters after splitting",
  "former clusters 1 and 3; B-cell panels a–b unchanged.",
  "Archetype: quantitative grid (panel d hero)",
  "Backend: R ggplot2 + patchwork, native vector",
  sprintf("Size: %d x %d mm", FIG_W_MM, FIG_H_MM),
  "T clusters: 0, 1a, 1b, 2, 3a, 3b (old 1 and 3 bipartitioned at subset res=0.1)",
  paste("New sizes:", paste(sprintf("%s=%d", names(STATS$n_t_cl), STATS$n_t_cl), collapse = ", "))
), file.path(OUTDIR, "figure_contract.txt"))

# Stats report
rep <- file.path(OUTDIR, "FigureS2_stats.txt")
con <- file(rep, "wt")
wl <- function(...) cat(sprintf(...), "\n", file = con)
wl("Figure S2 statistics (split c1/c3)")
wl("T new clusters: %s", paste(sprintf("%s=%d", names(STATS$n_t_cl), STATS$n_t_cl), collapse = "  "))
wl("chi-square %s", fmt_p(chi$p.value))
for (i in seq_len(nrow(fish))) {
  wl("  %s: OR=%.2f %s (Holm %s) n_neg=%d n_pos=%d",
     fish$cluster[i], fish$OR[i], fmt_p(fish$p[i]), fmt_p(fish$p_holm[i]),
     fish$n_neg[i], fish$n_pos[i])
}
close(con)

say("DONE -> %s", OUTDIR)
