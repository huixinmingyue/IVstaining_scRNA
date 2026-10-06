#!/usr/bin/env Rscript
# =============================================================================
# Fig. 2e P0 reanalysis — identity genes, not KEGG core, mapped to ImmGen
#
# Core conclusion
#   Cluster-identity and MZ/Fo marker genes separate ImmGen MZ from Fo;
#   the original KEGG-filtered gene set does not. MZ module scores are
#   higher in subcluster 1 and in IV+ B cells.
#
# Panel map
#   Figure_P0_identity
#     a  Do cluster 1 vs 0 marker genes map to different ImmGen B subsets?
#     b  Do the 21 S1a MZ ∩ Intra genes separate MZ from Fo?   (positive control)
#     c  Does a textbook MZ/Fo panel separate the same columns?
#     d  Does dropping T-lineage genes from the original KEGG set rescue MZ?
#   Figure_P0_intra
#     a  Do unfiltered within-cluster Intra DEGs map to different ImmGen subsets?
#   Figure_P0_modules
#     a–d  ImmGen MZ / Fo / Mem / GC module scores by subcluster and IV group
#
# Archetype : quantitative grid
# Backend   : R (ggplot2 / ComplexHeatmap / patchwork), exclusive
# Output    : this folder only; does not overwrite Fig2_Nature_v2
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
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

# ---- paths (new folder only) ------------------------------------------------
FIGDIR  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures")
IMMGEN  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "spleen_Ab_immgen")
S1DIR   <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS1_Nature_b_allDEG/source_data")
FIG2V2  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig2_Nature_v2")
OUTDIR  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig2e_P0_ImmGen_reanalysis")
SRC     <- file.path(OUTDIR, "source_data")
QA      <- file.path(OUTDIR, "qa")
for (d in c(OUTDIR, SRC, QA, file.path(OUTDIR, "scripts"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

FIG_W_MM <- 180
FAM  <- "Arial"
BASE <- 6.5
fm <- systemfonts::match_fonts(FAM)
cat("font resolve:", FAM, "->", fm$path, "\n")
for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
  if (!FAM %in% names(reg())) {
    args <- list(reg()$Helvetica); names(args) <- FAM
    do.call(reg, args)
  }
}

CL_COL <- c("0" = "#4C6E8A", "1" = "#D1495B", "2" = "#E9A03B")
IV_COL <- c("IntraNeg" = "#C9D1D6", "IntraPos" = "#4A5A66")
IV_LAB <- c("IntraNeg" = "IV\u2212 (extravascular)",
            "IntraPos" = "IV+ (intravascular)")
DIV_LO <- "#3B6FA0"; DIV_MID <- "#FFFFFF"; DIV_HI <- "#C0392B"
SUBSET_ORDER <- c("T1", "T2", "T3", "Fo", "MZ", "Mem", "GC-CB", "GC-CC")
TARGET_PREFIX <- c("B.T1.Sp", "B.T2.Sp", "B.T3.Sp",
                   "B.Fo.Sp", "B.MZ.Sp", "B.mem.Sp",
                   "B.GC.CB.Sp", "B.GC.CC.Sp")
DROP_PREFIX <- c("B.PB.Sp", "B.PC.Sp", "B.Sp", "B.Fem.Sp")
TCELL <- c("Camk4", "Cblb", "Cd28", "Cd3e", "Cd4", "Cd8b1", "Il2rb",
           "Itk", "Lck", "Lcp2", "Lef1", "Prkcq", "Skap1", "Cd3d",
           "Cd3g", "Cd8a", "Zap70", "Thy1", "Tcf7")
# The 12 genes named in the P0-3 contract (subcluster-1 KEGG set).
TCELL_P03 <- c("Camk4", "Cblb", "Cd28", "Cd3e", "Cd4", "Cd8b1",
               "Il2rb", "Itk", "Lck", "Lcp2", "Lef1", "Prkcq")

theme_p0 <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FAM) +
    theme(
      axis.line  = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.title = element_text(size = base),
      axis.text  = element_text(size = base - 0.4, colour = "black"),
      legend.title = element_text(size = base - 0.2),
      legend.text  = element_text(size = base - 0.6),
      legend.key.size = unit(2.6, "mm"),
      strip.text = element_text(size = base, face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(size = base + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "#555555", hjust = 0),
      plot.margin = margin(1.5, 2, 1.5, 2, "mm"),
      panel.grid = element_blank()
    )
}
theme_set(theme_p0())

fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 1e-4) return(formatC(x, format = "e", digits = 1))
    sprintf("%.3f", x)
  }, character(1))
}

STATS <- list()
log_lines <- character()
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  log_lines <<- c(log_lines, msg)
}

# ---- ImmGen helpers ---------------------------------------------------------
imm_raw <- read.csv(file.path(IMMGEN, "B.Sp_immgen.csv"),
                    row.names = 1, check.names = FALSE)

merge_immgen <- function(expr) {
  pref0 <- sub("(.*)\\.\\d+$", "\\1", colnames(expr))
  keep <- !pref0 %in% DROP_PREFIX
  expr <- expr[, keep, drop = FALSE]
  pref <- sub("(.*)\\.\\d+$", "\\1", colnames(expr))
  tbl <- as.data.frame(expr) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(prefix = sub("(.*)\\.\\d+$", "\\1", sample)) %>%
    group_by(gene, prefix) %>%
    summarise(expression = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = prefix, values_from = expression)
  mat <- as.matrix(tbl[, -1, drop = FALSE])
  rownames(mat) <- tbl$gene
  miss <- setdiff(TARGET_PREFIX, colnames(mat))
  if (length(miss)) stop("Missing ImmGen subsets: ", paste(miss, collapse = ", "))
  mat[, TARGET_PREFIX, drop = FALSE]
}

IMM_MEAN <- merge_immgen(imm_raw)
colnames(IMM_MEAN) <- SUBSET_ORDER
say("ImmGen merged matrix: %d genes x %d subsets", nrow(IMM_MEAN), ncol(IMM_MEAN))

remap_symbol <- function(g) {
  alt <- c(Cd1d = "Cd1d1", CR2 = "Cr2", FCER2A = "Fcer2a")
  ifelse(g %in% names(alt), unname(alt[g]), g)
}

genes_in_immgen <- function(genes) {
  genes <- unique(remap_symbol(as.character(genes)))
  genes[genes %in% rownames(IMM_MEAN)]
}

zscore_genes <- function(genes) {
  genes <- genes_in_immgen(genes)
  if (length(genes) < 1) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(IMM_MEAN),
                  dimnames = list(NULL, colnames(IMM_MEAN))))
  }
  mat <- IMM_MEAN[genes, , drop = FALSE]
  keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  t(scale(t(mat)))
}

peak_subset <- function(mat) {
  if (!nrow(mat)) return(character())
  colnames(mat)[max.col(mat, ties.method = "first")]
}

union_blocks <- function(genes_a, genes_b, lab_a, lab_b) {
  a <- unique(genes_in_immgen(genes_a))
  b <- unique(genes_in_immgen(genes_b))
  shared <- intersect(a, b)
  a_only <- setdiff(a, shared)
  b_only <- setdiff(b, shared)
  za <- zscore_genes(a_only)
  zs <- zscore_genes(shared)
  zb <- zscore_genes(b_only)
  # After variance filter, blocks can shrink.
  a_only <- rownames(za); shared <- rownames(zs); b_only <- rownames(zb)
  mat <- rbind(za, zs, zb)
  blk_lab <- c(
    sprintf("%s only (%d)", lab_a, length(a_only)),
    sprintf("Shared (%d)", length(shared)),
    sprintf("%s only (%d)", lab_b, length(b_only))
  )
  blk_n <- c(length(a_only), length(shared), length(b_only))
  keep_blk <- blk_n > 0
  memb <- factor(rep(blk_lab[keep_blk], blk_n[keep_blk]),
                 levels = blk_lab[keep_blk])
  delta <- p_raw <- p_holm <- setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  if (length(a_only) >= 2 && length(b_only) >= 2) {
    delta <- colMeans(za) - colMeans(zb)
    p_raw <- sapply(colnames(mat), function(s) {
      suppressWarnings(wilcox.test(za[, s], zb[, s], exact = FALSE)$p.value)
    })
    p_holm <- p.adjust(p_raw, method = "holm")
  }
  cor_full <- if (length(a) >= 3 && length(b) >= 3) {
    za_all <- colMeans(zscore_genes(a))
    zb_all <- colMeans(zscore_genes(b))
    suppressWarnings(cor(za_all, zb_all, use = "complete.obs"))
  } else NA_real_
  list(mat = mat, memb = memb, delta = delta, p = p_raw, p_holm = p_holm,
       a_only = a_only, b_only = b_only, shared = shared,
       n_a = length(a), n_b = length(b), cor_full = cor_full,
       lab_a = lab_a, lab_b = lab_b)
}

draw_split_ht <- function(u, title, row_fs = 5.2, col_fs = 6.2,
                          show_names = NULL) {
  mat <- u$mat
  if (!nrow(mat)) stop("empty heatmap matrix: ", title)
  if (is.null(show_names)) show_names <- nrow(mat) <= 45
  zlim <- max(2.5, ceiling(max(abs(mat), na.rm = TRUE) * 2) / 2)
  col_fun <- colorRamp2(c(-zlim, 0, zlim), c(DIV_LO, DIV_MID, DIV_HI))
  has_delta <- !all(is.na(u$delta))
  ta <- NULL
  if (has_delta) {
    dmax <- max(0.4, ceiling(max(abs(u$delta), na.rm = TRUE) * 10) / 10)
    ta <- HeatmapAnnotation(
      delta = anno_barplot(
        as.numeric(u$delta), baseline = 0, border = FALSE, bar_width = 0.72,
        gp = gpar(fill = ifelse(u$delta > 0, CL_COL[["1"]], CL_COL[["0"]]),
                  col = NA),
        ylim = c(-dmax, dmax),
        axis_param = list(gp = gpar(fontsize = 5, fontfamily = FAM))
      ),
      annotation_label = "\u0394 mean z",
      annotation_name_gp = gpar(fontsize = 5.5, fontfamily = FAM),
      annotation_name_side = "left",
      height = unit(14, "mm")
    )
  }
  Heatmap(
    mat, name = "z", col = col_fun,
    cluster_rows = TRUE, cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    row_split = u$memb, row_gap = unit(1.2, "mm"),
    row_title_gp = gpar(fontsize = 5.8, fontfamily = FAM, fontface = "bold"),
    row_title_rot = 0,
    show_row_names = show_names,
    row_names_gp = gpar(fontsize = row_fs, fontfamily = FAM, fontface = "italic"),
    column_names_gp = gpar(fontsize = col_fs, fontfamily = FAM),
    column_names_rot = 45,
    show_row_dend = FALSE,
    heatmap_legend_param = list(
      title = "z", title_gp = gpar(fontsize = 5.8, fontfamily = FAM),
      labels_gp = gpar(fontsize = 5.2, fontfamily = FAM),
      legend_height = unit(18, "mm"), grid_width = unit(2.4, "mm")
    ),
    top_annotation = ta,
    border = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 7, fontfamily = FAM, fontface = "bold"),
    width = unit(88, "mm"),
    height = unit(max(38, min(110, nrow(mat) * 1.55 + 22)), "mm")
  )
}

grab_ht <- function(ht, width_in, height_in) {
  tmp <- tempfile(fileext = ".svg")
  svglite::svglite(tmp, width = width_in, height = height_in)
  on.exit({ if (dev.cur() > 1L) dev.off(); unlink(tmp) }, add = TRUE)
  draw(ht, merge_legend = TRUE, padding = unit(c(1, 1, 1, 1), "mm"))
  grid.grabExpr({
    draw(ht, merge_legend = TRUE, padding = unit(c(1, 1, 1, 1), "mm"))
  }, width = width_in, height = height_in)
}

write_union_csv <- function(u, path) {
  out <- data.frame(gene = rownames(u$mat), gene_set = as.character(u$memb),
                    u$mat, check.names = FALSE)
  write.csv(out, path, row.names = FALSE)
}

write_delta_csv <- function(u, path) {
  write.csv(data.frame(
    subset = names(u$delta),
    delta_mean_z = as.numeric(u$delta),
    wilcoxon_p = as.numeric(u$p),
    wilcoxon_p_holm = as.numeric(u$p_holm)
  ), path, row.names = FALSE)
}

export_gg <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = FAM)
  print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = dpi, compression = "lzw")
  print(plot); dev.off()
}

# =============================================================================
# P0-1 identity: cluster 1 vs 0 markers (no KEGG)
# =============================================================================
say("=== P0-1 cluster markers ===")
b <- readRDS(file.path(FIGDIR, "Bcell_reclustered.rds"))
DefaultAssay(b) <- "RNA"
md <- slot(b, "meta.data")
md$cluster <- factor(as.character(md$seurat_clusters), levels = c("0", "1", "2"))
md$IV_group <- factor(md$IV_group, levels = c("IntraNeg", "IntraPos"))
Idents(b) <- md$cluster
say("cells %d; clusters %s", ncol(b),
    paste(sprintf("%s=%d", names(table(md$cluster)), table(md$cluster)), collapse = ", "))

mk_path <- file.path(SRC, "SourceData_P0-1_FindMarkers_cluster1_vs_0_unfiltered.csv")
if (file.exists(mk_path)) {
  mk10 <- read.csv(mk_path, check.names = FALSE)
  if (!"gene" %in% names(mk10) && "X" %in% names(mk10)) mk10$gene <- mk10$X
  say("cached FindMarkers: %d genes", nrow(mk10))
} else {
  mk10 <- FindMarkers(b, ident.1 = "1", ident.2 = "0",
                      logfc.threshold = 0.25, min.pct = 0.1,
                      test.use = "wilcox", verbose = FALSE)
  mk10$gene <- rownames(mk10)
  write.csv(mk10, mk_path, row.names = FALSE)
}

PADJ <- 0.001
LFC  <- 0.5
mk_sig <- mk10 %>% filter(p_val_adj < PADJ, abs(avg_log2FC) >= LFC)
cl1_up <- mk_sig$gene[mk_sig$avg_log2FC > 0]
cl0_up <- mk_sig$gene[mk_sig$avg_log2FC < 0]
say("cluster markers padj<%g |log2FC|>=%g: cl1-up %d, cl0-up %d (of %d tested)",
    PADJ, LFC, length(cl1_up), length(cl0_up), nrow(mk10))
STATS$n_mk_tested <- nrow(mk10)
STATS$n_cl1_up <- length(cl1_up)
STATS$n_cl0_up <- length(cl0_up)

write.csv(mk_sig, file.path(SRC, "SourceData_P0-1_cluster_markers_sig.csv"),
          row.names = FALSE)

u_id <- union_blocks(cl1_up, cl0_up, "Subcluster 1", "Subcluster 0")
say("identity ImmGen: cl1-only %d, shared %d, cl0-only %d, cor=%.3f, Holm sig %d/8",
    length(u_id$a_only), length(u_id$shared), length(u_id$b_only),
    u_id$cor_full, sum(u_id$p_holm < 0.05, na.rm = TRUE))
if (length(u_id$delta)) {
  i <- which.max(abs(u_id$delta))
  say("  largest |delta|: %s = %+.3f (P=%s, Holm=%s)",
      names(u_id$delta)[i], u_id$delta[i],
      fmt_p(u_id$p[i]), fmt_p(u_id$p_holm[i]))
}
write_union_csv(u_id, file.path(SRC, "SourceData_P0-1_identity_zscore.csv"))
write_delta_csv(u_id, file.path(SRC, "SourceData_P0-1_identity_delta.csv"))
# Display heatmap: top 40 per side by |avg_log2FC| so gene names stay >= 5 pt.
# Full-set delta (545 vs 228) is the inferential result and is plotted separately.
top_n_named <- function(genes, n = 40, positive = TRUE) {
  sub <- mk_sig %>% filter(gene %in% genes)
  sub <- if (positive) arrange(sub, desc(avg_log2FC)) else arrange(sub, avg_log2FC)
  head(sub$gene, n)
}
u_id_plot <- union_blocks(top_n_named(cl1_up, 40, TRUE),
                          top_n_named(cl0_up, 40, FALSE),
                          "Subcluster 1 top40", "Subcluster 0 top40")
say("identity heatmap display: cl1 %d, cl0 %d (full-set stats use %d vs %d)",
    length(u_id_plot$a_only), length(u_id_plot$b_only),
    length(u_id$a_only), length(u_id$b_only))
write_union_csv(u_id_plot, file.path(SRC, "SourceData_P0-1_identity_heatmap_top40_zscore.csv"))

# Intra DEG without KEGG (B-cluster IntraPos vs IntraNeg full tables)
say("=== P0-1 Intra DEG (no KEGG) ===")
intra1 <- read.csv(file.path(FIGDIR, "Fig2C_Bcluster1_IntraPos_vs_IntraNeg_DEG.csv"),
                   row.names = 1, check.names = FALSE)
intra0 <- read.csv(file.path(FIGDIR, "Fig2D_Bcluster0_IntraPos_vs_IntraNeg_DEG.csv"),
                   row.names = 1, check.names = FALSE)
# Direction: positive log2FC = IntraPos up (Fig2C/D convention, matching S1).
i1 <- intra1 %>% mutate(gene = rownames(intra1)) %>%
  filter(p_val_adj < PADJ, abs(avg_log2FC) >= LFC)
i0 <- intra0 %>% mutate(gene = rownames(intra0)) %>%
  filter(p_val_adj < PADJ, abs(avg_log2FC) >= LFC)
say("Intra DEG same threshold: cluster1 %d / %d tested; cluster0 %d / %d tested",
    nrow(i1), nrow(intra1), nrow(i0), nrow(intra0))
write.csv(i1, file.path(SRC, "SourceData_P0-1_intra_cluster1_sig.csv"), row.names = FALSE)
write.csv(i0, file.path(SRC, "SourceData_P0-1_intra_cluster0_sig.csv"), row.names = FALSE)
u_in <- union_blocks(i1$gene, i0$gene, "Cluster1 Intra DEG", "Cluster0 Intra DEG")
say("Intra ImmGen: c1-only %d, shared %d, c0-only %d, cor=%.3f, Holm sig %d/8",
    length(u_in$a_only), length(u_in$shared), length(u_in$b_only),
    u_in$cor_full, sum(u_in$p_holm < 0.05, na.rm = TRUE))
write_union_csv(u_in, file.path(SRC, "SourceData_P0-1_intra_zscore.csv"))
write_delta_csv(u_in, file.path(SRC, "SourceData_P0-1_intra_delta.csv"))

# =============================================================================
# P0-2 positive controls
# =============================================================================
say("=== P0-2 S1a 21 genes and textbook panel ===")
s1a <- read.csv(file.path(S1DIR, "S1a_MZvsFo_x_BcellDEG.csv"),
                check.names = FALSE)
s1a_genes <- s1a$gene
s1a_mz <- s1a$gene[grepl("^MZ", s1a$MZ_Fo_direction)]
s1a_fo <- s1a$gene[grepl("^Fo", s1a$MZ_Fo_direction)]
say("S1a genes %d (MZ-up %d, Fo-up %d); in ImmGen %d",
    length(s1a_genes), length(s1a_mz), length(s1a_fo),
    sum(genes_in_immgen(s1a_genes) %in% rownames(IMM_MEAN)))
u_s1 <- union_blocks(s1a_mz, s1a_fo, "MZ up", "Fo up")
say("S1a ImmGen: MZ-only %d, shared %d, Fo-only %d, Holm sig %d/8",
    length(u_s1$a_only), length(u_s1$shared), length(u_s1$b_only),
    sum(u_s1$p_holm < 0.05, na.rm = TRUE))
if (length(u_s1$delta)) {
  i <- which.max(abs(u_s1$delta))
  say("  largest |delta|: %s = %+.3f (Holm=%s)",
      names(u_s1$delta)[i], u_s1$delta[i], fmt_p(u_s1$p_holm[i]))
}
write_union_csv(u_s1, file.path(SRC, "SourceData_P0-2_S1a_zscore.csv"))
write_delta_csv(u_s1, file.path(SRC, "SourceData_P0-2_S1a_delta.csv"))

textbook <- c("Cr2", "Fcer2a", "Cd9", "Cd1d", "Cd36",
              "Notch2", "Dtx1", "Hes1", "Cd80", "Cd86")
tb_map <- genes_in_immgen(textbook)
say("textbook requested %d, in ImmGen after remap %d: %s",
    length(textbook), length(tb_map), paste(tb_map, collapse = ", "))
miss_tb <- setdiff(remap_symbol(textbook), tb_map)
if (length(miss_tb)) say("  missing from ImmGen: %s", paste(miss_tb, collapse = ", "))
z_tb <- zscore_genes(tb_map)
# Split textbook into MZ-canonical vs Fo-canonical for a delta, plus the rest.
mz_canon <- intersect(c("Cr2", "Cd9", "Cd1d1", "Cd1d", "Cd36", "Notch2",
                        "Dtx1", "Hes1", "Cd80", "Cd86"), rownames(z_tb))
fo_canon <- intersect(c("Fcer2a"), rownames(z_tb))
u_tb <- list(
  mat = z_tb,
  memb = factor(ifelse(rownames(z_tb) %in% fo_canon, "Fo canonical",
                       "MZ canonical"),
                levels = c("MZ canonical", "Fo canonical")),
  delta = setNames(rep(NA_real_, ncol(z_tb)), colnames(z_tb)),
  p = NA, p_holm = NA,
  a_only = mz_canon, b_only = fo_canon, shared = character(),
  lab_a = "MZ canonical", lab_b = "Fo canonical"
)
write.csv(data.frame(gene = rownames(z_tb), gene_set = as.character(u_tb$memb),
                     z_tb, check.names = FALSE),
          file.path(SRC, "SourceData_P0-2_textbook_zscore.csv"), row.names = FALSE)

# =============================================================================
# P0-3 drop T genes from original KEGG union
# =============================================================================
say("=== P0-3 drop T-lineage genes ===")
orig <- read.csv(file.path(FIG2V2, "SourceData_Fig2e_zscore_union.csv"),
                 check.names = FALSE)
# Recreate exclusive sets from the plotted union (z already computed).
sc1_all <- orig$gene[grepl("Subcluster 1 only|Shared", orig$gene_set)]
sc0_all <- orig$gene[grepl("Subcluster 0 only|Shared", orig$gene_set)]
sc1_ex  <- orig$gene[grepl("Subcluster 1 only", orig$gene_set)]
sc0_ex  <- orig$gene[grepl("Subcluster 0 only", orig$gene_set)]
say("original exclusive: sc1 %d, sc0 %d; T genes in sc1-only %s; in full sc1 set %s",
    length(sc1_ex), length(sc0_ex),
    paste(intersect(sc1_ex, TCELL_P03), collapse = ","),
    paste(intersect(sc1_all, TCELL_P03), collapse = ","))

# Sensitivity A: drop the 12-gene list from the exclusive 35 (genes present there).
sc1_drop35 <- setdiff(sc1_ex, TCELL_P03)
u_drop35 <- union_blocks(sc1_drop35, sc0_ex, "SC1 no-T exclusive", "SC0 exclusive")
# Sensitivity B: drop the 12-gene list from the full 51-gene cluster-1 set.
u_drop51 <- union_blocks(setdiff(sc1_all, TCELL_P03),
                         setdiff(sc0_all, TCELL_P03),
                         "SC1 set no-T", "SC0 set no-T")

# Original delta from the v2 source table (do not recompute from a new z).
orig_delta <- read.csv(file.path(FIG2V2, "SourceData_Fig2e_delta.csv"),
                       check.names = FALSE)
cmp <- data.frame(
  subset = factor(orig_delta$subset, levels = SUBSET_ORDER),
  original = orig_delta$delta_mean_z,
  original_p = orig_delta$wilcoxon_p,
  original_holm = orig_delta$wilcoxon_p_holm,
  drop_from_35 = as.numeric(u_drop35$delta[as.character(orig_delta$subset)]),
  drop_from_35_holm = as.numeric(u_drop35$p_holm[as.character(orig_delta$subset)]),
  drop_from_51 = as.numeric(u_drop51$delta[as.character(orig_delta$subset)]),
  drop_from_51_holm = as.numeric(u_drop51$p_holm[as.character(orig_delta$subset)])
)
write.csv(cmp, file.path(SRC, "SourceData_P0-3_delta_sensitivity.csv"),
          row.names = FALSE)
write_union_csv(u_drop35, file.path(SRC, "SourceData_P0-3_dropT_from35_zscore.csv"))
write_delta_csv(u_drop35, file.path(SRC, "SourceData_P0-3_dropT_from35_delta.csv"))
write_union_csv(u_drop51, file.path(SRC, "SourceData_P0-3_dropT_from51_zscore.csv"))
write_delta_csv(u_drop51, file.path(SRC, "SourceData_P0-3_dropT_from51_delta.csv"))
say("P0-3 MZ delta original %.3f; drop-from-35 %.3f (Holm %s); drop-from-51 %.3f (Holm %s)",
    cmp$original[cmp$subset == "MZ"],
    cmp$drop_from_35[cmp$subset == "MZ"], fmt_p(cmp$drop_from_35_holm[cmp$subset == "MZ"]),
    cmp$drop_from_51[cmp$subset == "MZ"], fmt_p(cmp$drop_from_51_holm[cmp$subset == "MZ"]))
say("P0-3 Mem delta original %.3f; drop-from-35 %.3f (Holm %s)",
    cmp$original[cmp$subset == "Mem"],
    cmp$drop_from_35[cmp$subset == "Mem"],
    fmt_p(cmp$drop_from_35_holm[cmp$subset == "Mem"]))

# =============================================================================
# P0-4 ImmGen MZ/Fo/Mem/GC signatures -> module scores on B cells
# =============================================================================
say("=== P0-4 ImmGen signatures and module scores ===")
# Expression floor: keep genes with mean count across the 8 subsets >= 20
# so rare/noisy genes do not dominate fold-changes. Ig / ribosome / mt / Gm
# are excluded because they dominate B-cell DE and are not subset-identity
# programmes. Before/after counts are recorded.
expr_mean <- rowMeans(IMM_MEAN)
ok_expr <- expr_mean >= 20
n_before <- nrow(IMM_MEAN)
drop_pat <- "^(Igh[vdj]|Igk[vdj]|Igl[vdj]|Ighm$|Ighd$|Ighe$|Rp[sl]|mt-|Gm[0-9])"
ok_name <- !grepl(drop_pat, rownames(IMM_MEAN))
keep_sig <- ok_expr & ok_name
say("signature universe: %d -> %d after mean>=20 and Ig/Rpl/mt/Gm filter",
    n_before, sum(keep_sig))
M <- IMM_MEAN[keep_sig, , drop = FALSE]
lfc <- function(a, b) log2((a + 1) / (b + 1))
fc_mz  <- lfc(M[, "MZ"],  M[, "Fo"])
fc_mem <- lfc(M[, "Mem"], M[, "Fo"])
fc_gc  <- lfc(rowMeans(M[, c("GC-CB", "GC-CC"), drop = FALSE]), M[, "Fo"])
fc_fo  <- -fc_mz

top_n <- function(x, n = 30, up = TRUE) {
  x <- x[is.finite(x)]
  x <- if (up) sort(x, decreasing = TRUE) else sort(x, decreasing = FALSE)
  # require |log2FC| >= 1 so the tail is not near-zero noise
  x <- x[abs(x) >= 1]
  head(names(x), n)
}
sig_mz  <- top_n(fc_mz,  30, TRUE)
sig_fo  <- top_n(fc_mz,  30, FALSE)
sig_mem <- top_n(fc_mem, 30, TRUE)
sig_gc  <- top_n(fc_gc,  30, TRUE)
# Keep only genes present in the Seurat RNA assay.
rna_genes <- rownames(b[["RNA"]])
to_seurat <- function(g) intersect(g, rna_genes)
sig_mz  <- to_seurat(sig_mz)
sig_fo  <- to_seurat(sig_fo)
sig_mem <- to_seurat(sig_mem)
sig_gc  <- to_seurat(sig_gc)
say("modules in Seurat: MZ %d, Fo %d, Mem %d, GC %d",
    length(sig_mz), length(sig_fo), length(sig_mem), length(sig_gc))

sig_tbl <- rbind(
  data.frame(module = "MZ",  gene = sig_mz,  log2FC_vs_Fo = unname(fc_mz[sig_mz])),
  data.frame(module = "Fo",  gene = sig_fo,  log2FC_vs_Fo = unname(fc_mz[sig_fo])),
  data.frame(module = "Mem", gene = sig_mem, log2FC_vs_Fo = unname(fc_mem[sig_mem])),
  data.frame(module = "GC",  gene = sig_gc,  log2FC_vs_Fo = unname(fc_gc[sig_gc]))
)
write.csv(sig_tbl, file.path(SRC, "SourceData_P0-4_ImmGen_module_genes.csv"),
          row.names = FALSE)

b <- AddModuleScore(b, features = list(sig_mz, sig_fo, sig_mem, sig_gc),
                    name = "P0mod", assay = "RNA", search = FALSE)
# AddModuleScore names columns P0mod1..4
md <- slot(b, "meta.data")
md$cluster  <- factor(as.character(md$seurat_clusters), levels = c("0", "1", "2"))
md$IV_group <- factor(md$IV_group, levels = c("IntraNeg", "IntraPos"))
md$MZ_score  <- md$P0mod1
md$Fo_score  <- md$P0mod2
md$Mem_score <- md$P0mod3
md$GC_score  <- md$P0mod4
md$MZ_minus_Fo <- md$MZ_score - md$Fo_score

score_long <- md %>%
  transmute(cell = rownames(md), cluster, IV_group,
            MZ = MZ_score, Fo = Fo_score, Mem = Mem_score, GC = GC_score,
            MZ_minus_Fo) %>%
  pivot_longer(c(MZ, Fo, Mem, GC, MZ_minus_Fo),
               names_to = "module", values_to = "score")
write.csv(score_long, file.path(SRC, "SourceData_P0-4_module_scores_per_cell.csv"),
          row.names = FALSE)

# Descriptive two-sided Wilcoxon; cells are the observation unit (n = 1 pooled
# sample per compartment). Holm across 4 modules x 2 primary contrasts.
mod4 <- c("MZ", "Fo", "Mem", "GC")
stat_mod <- lapply(mod4, function(m) {
  x <- score_long %>% filter(module == m)
  cl1 <- x$score[x$cluster == "1"]
  cl0 <- x$score[x$cluster == "0"]
  ivp <- x$score[x$IV_group == "IntraPos"]
  ivn <- x$score[x$IV_group == "IntraNeg"]
  data.frame(
    module = m,
    mean_cl0 = mean(cl0), mean_cl1 = mean(cl1),
    delta_cl1_minus_cl0 = mean(cl1) - mean(cl0),
    p_cluster = wilcox.test(cl1, cl0)$p.value,
    mean_IV_neg = mean(ivn), mean_IV_pos = mean(ivp),
    delta_IVpos_minus_IVneg = mean(ivp) - mean(ivn),
    p_IV = wilcox.test(ivp, ivn)$p.value
  )
}) %>% bind_rows()
stat_mod$p_cluster_holm <- p.adjust(stat_mod$p_cluster, method = "holm")
stat_mod$p_IV_holm <- p.adjust(stat_mod$p_IV, method = "holm")
write.csv(stat_mod, file.path(SRC, "SourceData_P0-4_module_score_tests.csv"),
          row.names = FALSE)
print(stat_mod)

mean_se <- score_long %>%
  filter(module %in% mod4) %>%
  group_by(module, cluster, IV_group) %>%
  summarise(n = n(), mean = mean(score), se = sd(score) / sqrt(n()),
            .groups = "drop")
write.csv(mean_se, file.path(SRC, "SourceData_P0-4_module_score_mean_se.csv"),
          row.names = FALSE)

# =============================================================================
# Figures
# =============================================================================
say("=== drawing figures ===")
ht_id <- draw_split_ht(u_id_plot, "P0-1  Cluster 1 vs 0 markers (top 40 per side)")
delta_id_df <- data.frame(
  subset = factor(names(u_id$delta), levels = SUBSET_ORDER),
  delta = as.numeric(u_id$delta),
  p_holm = as.numeric(u_id$p_holm)
)
p_id_delta <- ggplot(delta_id_df, aes(subset, delta, fill = delta > 0)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#888888") +
  geom_col(width = 0.72, colour = NA) +
  geom_text(aes(label = ifelse(p_holm < 0.05, sprintf("q=%s", fmt_p(p_holm)), "n.s.")),
          vjust = ifelse(delta_id_df$delta >= 0, -0.3, 1.2), size = 5 / .pt, family = FAM) +
  scale_fill_manual(values = c("TRUE" = CL_COL[["1"]], "FALSE" = CL_COL[["0"]]),
                    guide = "none") +
  labs(x = NULL, y = "\u0394 mean z (all sig. markers)",
       title = "Full-set \u0394 (545 vs 228 genes)",
       subtitle = "Holm-adjusted Wilcoxon across 8 subsets") +
  theme_p0()
ht_s1 <- draw_split_ht(u_s1, "P0-2  S1a MZ \u2229 Intra (21 genes)", row_fs = 5.5)
ht_tb <- Heatmap(
  u_tb$mat, name = "z",
  col = colorRamp2(c(-2.5, 0, 2.5), c(DIV_LO, DIV_MID, DIV_HI)),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = u_tb$memb, row_gap = unit(1.2, "mm"),
  row_title_gp = gpar(fontsize = 5.8, fontfamily = FAM, fontface = "bold"),
  row_title_rot = 0,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6.2, fontfamily = FAM, fontface = "italic"),
  column_names_gp = gpar(fontsize = 6.2, fontfamily = FAM),
  column_names_rot = 45,
  show_row_dend = FALSE,
  heatmap_legend_param = list(
    title = "z", title_gp = gpar(fontsize = 5.8, fontfamily = FAM),
    labels_gp = gpar(fontsize = 5.2, fontfamily = FAM),
    legend_height = unit(16, "mm"), grid_width = unit(2.4, "mm")
  ),
  border = FALSE,
  column_title = "P0-2  Textbook MZ / Fo panel",
  column_title_gp = gpar(fontsize = 7, fontfamily = FAM, fontface = "bold"),
  width = unit(78, "mm"),
  height = unit(42, "mm")
)
ht_in <- draw_split_ht(u_in, "P0-1  Intra DEG, no KEGG filter",
                       show_names = FALSE)

# P0-3 grouped bars
cmp_long <- cmp %>%
  select(subset, original, drop_from_35, drop_from_51) %>%
  pivot_longer(-subset, names_to = "method", values_to = "delta") %>%
  mutate(method = recode(method,
                         original = "Original KEGG set",
                         drop_from_35 = "Drop 12 T genes from exclusive 35",
                         drop_from_51 = "Drop 12 T genes from full SC1 set"),
         method = factor(method, levels = c(
           "Original KEGG set",
           "Drop 12 T genes from exclusive 35",
           "Drop 12 T genes from full SC1 set")))
p_d <- ggplot(cmp_long, aes(subset, delta, fill = method)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#888888") +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  scale_fill_manual(values = c("#7A7A7A", "#4C6E8A", "#D1495B"), name = NULL) +
  labs(x = NULL, y = "\u0394 mean z (SC1-only \u2212 SC0-only)",
       title = "P0-3  Dropping T-lineage genes does not create an MZ peak",
       subtitle = "Positive \u0394: subcluster-1 exclusive genes higher in that ImmGen subset. Holm-adjusted tests remain non-significant at MZ.") +
  theme_p0() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.x = element_text(size = 6.2))

# Module violins
score_plot <- score_long %>%
  filter(module %in% mod4) %>%
  mutate(module = factor(module, levels = mod4,
                         labels = c("MZ signature", "Fo signature",
                                    "Mem signature", "GC signature")),
         IV = recode(IV_group, !!!IV_LAB))
p_mod <- ggplot(score_plot, aes(cluster, score, fill = IV)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, width = 0.86,
              position = position_dodge(width = 0.85)) +
  geom_boxplot(width = 0.14, outlier.size = 0.15, outlier.alpha = 0.25,
               position = position_dodge(width = 0.85),
               linewidth = 0.25) +
  facet_wrap(~ module, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("IV\u2212 (extravascular)" = IV_COL[["IntraNeg"]],
                               "IV+ (intravascular)" = IV_COL[["IntraPos"]]),
                    name = NULL) +
  labs(x = "B-cell subcluster", y = "Module score",
       title = "P0-4  ImmGen MZ / Fo / Mem / GC scores on reclustered B cells",
       subtitle = "Signatures: top 30 ImmGen genes with |log2FC| \u2265 1 versus Fo. Cells are the observation unit (one pooled sample per compartment).") +
  theme_p0() +
  theme(legend.position = "bottom")

p_delta_mod <- ggplot(mean_se %>%
                        filter(module == "MZ") %>%
                        mutate(IV = recode(IV_group, !!!IV_LAB)),
                      aes(cluster, mean, fill = IV)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.75), width = 0.18, linewidth = 0.3) +
  scale_fill_manual(values = c("IV\u2212 (extravascular)" = IV_COL[["IntraNeg"]],
                               "IV+ (intravascular)" = IV_COL[["IntraPos"]]),
                    name = NULL) +
  labs(x = "B-cell subcluster", y = "MZ module score (mean \u00b1 s.e.)",
       title = "MZ signature mean by subcluster and compartment") +
  theme_p0() +
  theme(legend.position = "bottom")

# ---- export heatmaps (ComplexHeatmap devices) --------------------------------
export_ht <- function(ht, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = FAM)
  draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600)
  draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
}

export_ht(ht_id, file.path(OUTDIR, "Fig_P0-1_cluster_markers_ImmGen"), 180, 170)
export_gg(p_id_delta, file.path(OUTDIR, "Fig_P0-1_identity_delta_fullset"), 120, 62)
export_ht(ht_s1, file.path(OUTDIR, "Fig_P0-2_S1a21_ImmGen"), 180, 95)
export_ht(ht_tb, file.path(OUTDIR, "Fig_P0-2_textbook_ImmGen"), 180, 70)
export_ht(ht_in, file.path(OUTDIR, "Fig_P0-1_IntraDEG_ImmGen"), 180, 150)
export_gg(p_d, file.path(OUTDIR, "Fig_P0-3_delta_sensitivity"), 180, 85)
export_gg(p_mod, file.path(OUTDIR, "Fig_P0-4_module_scores"), 180, 145)
export_gg(p_delta_mod, file.path(OUTDIR, "Fig_P0-4_MZ_score_means"), 90, 70)

# Combined identity page: S1a + textbook + P0-3 (ggplot page)
p_page2 <- (plot_spacer() / p_d) + plot_layout(heights = c(0.08, 1))
# Combined modules page already faceted.

# Combined PDF with all heatmaps as pages
cairo_pdf(file.path(OUTDIR, "Figure_P0_all_panels.pdf"),
          width = 180 / 25.4, height = 220 / 25.4, family = FAM, onefile = TRUE)
draw(ht_id, merge_legend = TRUE,
     column_title = "a  Cluster 1 vs 0 markers (top 40/side) on ImmGen",
     column_title_gp = gpar(fontsize = 8, fontfamily = FAM, fontface = "bold"))
print(p_id_delta + labs(title = "a'  Full marker set delta (545 vs 228 genes)"))
draw(ht_s1, merge_legend = TRUE,
     column_title = "b  S1a MZ \u2229 Intra genes (positive control)",
     column_title_gp = gpar(fontsize = 8, fontfamily = FAM, fontface = "bold"))
draw(ht_tb, merge_legend = TRUE,
     column_title = "c  Textbook MZ / Fo panel",
     column_title_gp = gpar(fontsize = 8, fontfamily = FAM, fontface = "bold"))
print(p_d + labs(title = "d  P0-3  T-gene drop does not rescue an MZ identity"))
draw(ht_in, merge_legend = TRUE,
     column_title = "e  Within-cluster Intra DEGs without KEGG filter",
     column_title_gp = gpar(fontsize = 8, fontfamily = FAM, fontface = "bold"))
print(p_mod)
print(p_delta_mod)
dev.off()

# =============================================================================
# Stats + legend
# =============================================================================
wl <- function(...) say(...)
sink_file <- file.path(OUTDIR, "Figure_P0_stats.txt")
con <- file(sink_file, open = "wt")
writeLines(c(
  "Fig. 2e P0 reanalysis — statistics",
  "================================================================",
  "",
  "Output folder does not overwrite Fig2_Nature_v2.",
  "Threshold for cluster markers and Intra DEG: Wilcoxon; |avg_log2FC| >= 0.5; adjusted P < 0.001.",
  "Cells treated as observations; one pooled sample per IV compartment.",
  "",
  sprintf("[P0-1 identity] FindMarkers cluster 1 vs 0"),
  sprintf("  genes tested            : %d", nrow(mk10)),
  sprintf("  significant             : %d (cl1-up %d, cl0-up %d)",
          nrow(mk_sig), length(cl1_up), length(cl0_up)),
  sprintf("  in ImmGen after z       : cl1-only %d, shared %d, cl0-only %d",
          length(u_id$a_only), length(u_id$shared), length(u_id$b_only)),
  sprintf("  profile correlation     : r = %s", fmt_p(u_id$cor_full)),
  sprintf("  Holm-significant subsets: %d / 8", sum(u_id$p_holm < 0.05, na.rm = TRUE)),
  paste(sprintf("    %s  delta=%+.3f  P=%s  Holm=%s",
                names(u_id$delta), u_id$delta, fmt_p(u_id$p), fmt_p(u_id$p_holm)),
        collapse = "\n"),
  "",
  sprintf("[P0-1 Intra DEG] no KEGG filter, same threshold"),
  sprintf("  cluster 1 Intra sig     : %d / %d", nrow(i1), nrow(intra1)),
  sprintf("  cluster 0 Intra sig     : %d / %d", nrow(i0), nrow(intra0)),
  sprintf("  in ImmGen after z       : c1-only %d, shared %d, c0-only %d",
          length(u_in$a_only), length(u_in$shared), length(u_in$b_only)),
  sprintf("  profile correlation     : r = %s", fmt_p(u_in$cor_full)),
  sprintf("  Holm-significant subsets: %d / 8", sum(u_in$p_holm < 0.05, na.rm = TRUE)),
  paste(sprintf("    %s  delta=%+.3f  P=%s  Holm=%s",
                names(u_in$delta), u_in$delta, fmt_p(u_in$p), fmt_p(u_in$p_holm)),
        collapse = "\n"),
  "",
  sprintf("[P0-2 S1a] n=%d genes; MZ-up %d; Fo-up %d",
          length(s1a_genes), length(s1a_mz), length(s1a_fo)),
  sprintf("  Holm-significant subsets: %d / 8", sum(u_s1$p_holm < 0.05, na.rm = TRUE)),
  paste(sprintf("    %s  delta=%+.3f  P=%s  Holm=%s",
                names(u_s1$delta), u_s1$delta, fmt_p(u_s1$p), fmt_p(u_s1$p_holm)),
        collapse = "\n"),
  "",
  sprintf("[P0-2 textbook] requested %s", paste(textbook, collapse = ", ")),
  sprintf("  plotted                 : %s", paste(rownames(z_tb), collapse = ", ")),
  sprintf("  missing                 : %s",
          if (length(miss_tb)) paste(miss_tb, collapse = ", ") else "none"),
  "",
  "[P0-3] original KEGG exclusive sets vs drop-T",
  sprintf("  T genes in exclusive 35 : %s",
          paste(intersect(sc1_ex, TCELL_P03), collapse = ", ")),
  sprintf("  T genes in full SC1 51  : %s",
          paste(intersect(sc1_all, TCELL_P03), collapse = ", ")),
  paste(apply(cmp, 1, function(r) {
    sprintf("    %s  orig=%s  drop35=%s (Holm %s)  drop51=%s (Holm %s)",
            r[["subset"]],
            sprintf("%+.3f", as.numeric(r[["original"]])),
            sprintf("%+.3f", as.numeric(r[["drop_from_35"]])),
            fmt_p(as.numeric(r[["drop_from_35_holm"]])),
            sprintf("%+.3f", as.numeric(r[["drop_from_51"]])),
            fmt_p(as.numeric(r[["drop_from_51_holm"]])))
  }), collapse = "\n"),
  "",
  "[P0-4] ImmGen module scores",
  sprintf("  signature universe      : %d -> %d (mean>=20, drop Ig/Rpl/mt/Gm)",
          n_before, sum(keep_sig)),
  sprintf("  module sizes            : MZ %d, Fo %d, Mem %d, GC %d",
          length(sig_mz), length(sig_fo), length(sig_mem), length(sig_gc)),
  sprintf("  Mem reference n         : 1 ImmGen sample (B.mem.Sp.1); FC only"),
  capture.output(print(stat_mod, row.names = FALSE, digits = 4)),
  "",
  "INTERPRETATION LOCK",
  "  Identity markers (P0-1) and S1a/textbook (P0-2) are the tests of whether",
  "  ImmGen can separate MZ from Fo. The original Fig. 2e KEGG set is not.",
  "  P0-3 shows the Mem peak is not rescued by dropping T genes.",
  "  P0-4 tests the same claim in cell space rather than gene-set space.",
  "",
  "Log:",
  log_lines
), con)
close(con)

writeLines(c(
  "Fig. P0 | Identity gene sets, not KEGG core, mapped onto ImmGen splenic B subsets.",
  "",
  "a, ImmGen z-scores for genes marking B-cell subcluster 1 versus 0",
  sprintf("(Wilcoxon, |avg_log2FC| \u2265 0.5, adjusted P < 0.001; %d versus %d genes).",
          length(cl1_up), length(cl0_up)),
  "Columns keep developmental order (T1\u2013GC-CC) and are not clustered; genes are",
  "clustered within membership blocks. Right-hand bars: \u0394 = mean z(subcluster-1-only)",
  "\u2212 mean z(subcluster-0-only); Wilcoxon per subset, Holm-adjusted across eight subsets.",
  "",
  "b, Same layout for the 21 genes shared by B-cell IntraPos versus IntraNeg DEG",
  "and CITE MZ versus Fo DEG (Fig. S1a). Split by the CITE identity direction.",
  "",
  "c, Textbook MZ/Fo panel (Cr2, Fcer2a, Cd9, Cd1d1, Cd36, Notch2, Dtx1, Hes1, Cd80, Cd86)",
  "after mapping onto the same ImmGen matrix. Fcer2a is the Fo canonical row.",
  "",
  "d, Sensitivity of the original Fig. 2e \u0394 statistic after removing 12 T-lineage",
  "pathway genes from the exclusive 35-gene block or from the full subcluster-1 set.",
  "",
  "e, Within-cluster IntraPos versus IntraNeg DEG (no KEGG filter) mapped as in a.",
  "",
  "f, ImmGen-derived MZ, Fo, Mem and GC module scores (top 30 genes with |log2FC| \u2265 1",
  "versus Fo on the merged ULI matrix; Ig/ribosomal/mitochondrial/Gm genes excluded)",
  "on 5,975 reclustered B cells. Violin: kernel density; box: median and IQR.",
  "n = 1 pooled sample per compartment; cells are the observation unit.",
  "",
  "Source data are provided as Source Data files in this folder."
), file.path(OUTDIR, "Figure_P0_legend.txt"))

writeLines(c(
  "Fig2e_P0_ImmGen_reanalysis",
  "==========================",
  "New folder. Does not overwrite Fig2_Nature_v2.",
  "Rebuild:",
  "  \"C:\\Program Files\\R\\R-4.4.1\\bin\\x64\\Rscript.exe\" scripts\\run_P0_reanalysis.R"
), file.path(OUTDIR, "README.txt"))

say("DONE. Wrote figures and source data to %s", OUTDIR)
