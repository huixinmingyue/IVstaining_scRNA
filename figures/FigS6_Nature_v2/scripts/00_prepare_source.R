#!/usr/bin/env Rscript
# CITE T cells (same experiment as Fig. 5 B): drop QC clusters, recluster
# remaining cells, then subset for UMAP panels. New folder only.
options(warn = 1)
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUT <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS6_Nature_v2")
}
SRC <- file.path(OUT, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

HI <- 6
LO <- 1
FC_CUT <- 0.25
PCA_DIMS <- 1:15
RES <- 0.2

OBJ <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/CITE_Tcell_reclustered.rds")
stopifnot(file.exists(OBJ))
cat("Loading CITE T cells ...\n")
obj0 <- readRDS(OBJ)
DefaultAssay(obj0) <- "RNA"
stopifnot(ncol(obj0) == 4314)
stopifnot("Proteinfeatures" %in% colnames(slot(obj0, "meta.data")))

# First-pass labels on the parent T object (used only to drop QC clusters).
clus_map0 <- c(
  "0" = "Naive CD4 T",
  "1" = "Naive CD8 T",
  "2" = "Unassigned",
  "3" = "Memory/activated CD4 T",
  "4" = "B-T doublet",
  "5" = "Myeloid-contaminated",
  "6" = "Memory-like CD8 T"
)
qc_drop <- c("B-T doublet", "Myeloid-contaminated")
cl0 <- as.character(obj0$seurat_clusters)
stopifnot(all(cl0 %in% names(clus_map0)))
lab0 <- unname(clus_map0[cl0])
n_before <- ncol(obj0)
keep <- !lab0 %in% qc_drop
cat("QC drop (first-pass clusters 4 and 5): keep", sum(keep), "of", n_before, "\n")

obj <- subset(obj0, cells = colnames(obj0)[keep])
rm(obj0)
gc(verbose = FALSE)

recluster_rna <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = PCA_DIMS, verbose = FALSE)
  seu <- FindClusters(seu, resolution = RES, verbose = FALSE)
  seu <- RunUMAP(seu, dims = PCA_DIMS, verbose = FALSE)
  seu
}

cat("Reclustering QC-filtered T cells (PCA dims 1:15, res", RES, ") ...\n")
obj <- recluster_rna(obj)
cat("new clusters\n")
print(table(obj$seurat_clusters))

rna <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                error = function(e) GetAssayData(obj, assay = "RNA"))
marker_genes <- c(
  "Cd4", "Cd8a", "Cd8b1", "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r",
  "S100a4", "Cxcr3", "Foxp3", "Il2ra",
  "Eomes", "Il2rb", "Nkg7", "Ccl5", "Gzma",
  "Klrb1c", "Ms4a1", "Cd79a", "Adgre1", "Cd74", "Trdc"
)
marker_genes <- marker_genes[marker_genes %in% rownames(rna)]

pct_by_cluster <- function(mat, cl, genes) {
  cl <- as.character(cl)
  out <- lapply(sort(unique(cl)), function(k) {
    cells <- which(cl == k)
    v <- vapply(genes, function(g) mean(mat[g, cells] > 0), numeric(1))
    data.frame(cluster = k, n = length(cells), t(v), stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

ann_tbl <- pct_by_cluster(rna, obj$seurat_clusters, marker_genes)
write.csv(ann_tbl, file.path(SRC, "recluster_marker_pct.csv"), row.names = FALSE)
cat("marker pct\n")
print(ann_tbl)

annotate_cluster <- function(row) {
  p <- function(g) {
    if (!g %in% names(row)) return(NA_real_)
    as.numeric(row[[g]])
  }
  if (!is.na(p("Ms4a1")) && p("Ms4a1") >= 0.50 &&
      !is.na(p("Cd79a")) && p("Cd79a") >= 0.50) {
    return("B-T doublet")
  }
  if (!is.na(p("Adgre1")) && p("Adgre1") >= 0.25) {
    return("Myeloid-contaminated")
  }
  if (!is.na(p("Klrb1c")) && p("Klrb1c") >= 0.50) {
    return("NKT-like T")
  }
  cd4 <- p("Cd4"); cd8 <- p("Cd8b1"); cd8a <- p("Cd8a")
  # Memory/activated CD4 before naive CD4 (S100a4 / Foxp3 / Cxcr3).
  if (!is.na(cd4) && cd4 >= 0.50 && (is.na(cd8) || cd8 < 0.35)) {
    if ((!is.na(p("S100a4")) && p("S100a4") >= 0.15) ||
        (!is.na(p("Foxp3")) && p("Foxp3") >= 0.10) ||
        (!is.na(p("Cxcr3")) && p("Cxcr3") >= 0.12)) {
      return("Memory/activated CD4 T")
    }
    return("Naive CD4 T")
  }
  # Memory-like CD8: Eomes/Il2rb programme (Fig. 3 CD8-3a), even if Cd8b1 < 50%.
  if (!is.na(cd4) && cd4 < 0.35 &&
      !is.na(p("Eomes")) && p("Eomes") >= 0.40 &&
      !is.na(p("Il2rb")) && p("Il2rb") >= 0.40) {
    return("Memory-like CD8 T")
  }
  is_cd8 <- !is.na(cd4) && cd4 < 0.35 &&
    ((!is.na(cd8) && cd8 >= 0.50) || (!is.na(cd8a) && cd8a >= 0.50))
  if (is_cd8) {
    if (!is.na(p("Ccl5")) && p("Ccl5") >= 0.45 &&
        (is.na(p("Ccr7")) || p("Ccr7") < 0.40)) {
      return("Effector/cytotoxic CD8 T")
    }
    return("Naive CD8 T")
  }
  "Unassigned"
}

ann_tbl$Tcelltype <- vapply(seq_len(nrow(ann_tbl)), function(i) {
  annotate_cluster(ann_tbl[i, , drop = FALSE])
}, character(1))
write.csv(ann_tbl, file.path(SRC, "recluster_annotation.csv"), row.names = FALSE)
cat("cluster annotation\n")
print(ann_tbl[, c("cluster", "n", "Tcelltype")])

map_new <- setNames(ann_tbl$Tcelltype, as.character(ann_tbl$cluster))
lab_new <- unname(map_new[as.character(obj$seurat_clusters)])
if (any(lab_new %in% qc_drop)) {
  cat("Residual QC clusters after recluster; subset and recluster again ...\n")
  obj <- subset(obj, cells = colnames(obj)[!lab_new %in% qc_drop])
  obj <- recluster_rna(obj)
  rna <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                  error = function(e) GetAssayData(obj, assay = "RNA"))
  marker_genes <- marker_genes[marker_genes %in% rownames(rna)]
  ann_tbl <- pct_by_cluster(rna, obj$seurat_clusters, marker_genes)
  ann_tbl$Tcelltype <- vapply(seq_len(nrow(ann_tbl)), function(i) {
    annotate_cluster(ann_tbl[i, , drop = FALSE])
  }, character(1))
  write.csv(ann_tbl, file.path(SRC, "recluster_annotation.csv"), row.names = FALSE)
  map_new <- setNames(ann_tbl$Tcelltype, as.character(ann_tbl$cluster))
  lab_new <- unname(map_new[as.character(obj$seurat_clusters)])
  stopifnot(!any(lab_new %in% qc_drop))
}

fun_pref <- c(
  "Naive CD4 T", "Naive CD8 T", "Memory/activated CD4 T",
  "Memory-like CD8 T", "Effector/cytotoxic CD8 T", "NKT-like T", "Unassigned"
)
fun_lv <- fun_pref[fun_pref %in% unique(lab_new)]
obj$Tcelltype <- factor(lab_new, levels = fun_lv)
stopifnot(!anyNA(obj$Tcelltype))

saveRDS(obj, file.path(SRC, "CITE_Tcell_QCfiltered_reclustered.rds"))

um <- as.data.frame(Embeddings(obj, "umap"))
colnames(um) <- c("umap_1", "umap_2")
md <- slot(obj, "meta.data")
df <- data.frame(
  cell = colnames(obj),
  umap_1 = um$umap_1,
  umap_2 = um$umap_2,
  cluster = as.character(md$seurat_clusters),
  Tcelltype = as.character(obj$Tcelltype),
  Proteinfeatures = as.numeric(md$Proteinfeatures),
  nCount_RNA = as.numeric(md$nCount_RNA),
  percent.mt = as.numeric(md$percent.mt),
  stringsAsFactors = FALSE
)
df$Tcelltype <- factor(df$Tcelltype, levels = fun_lv)
stopifnot(nrow(df) == ncol(obj), all(is.finite(df$Proteinfeatures)))

n_all <- as.data.frame(table(df$Tcelltype), stringsAsFactors = FALSE)
names(n_all) <- c("Tcelltype", "n")
write.csv(n_all, file.path(SRC, "n_by_Tcelltype.csv"), row.names = FALSE)
cat("kept after recluster / subset:\n")
print(n_all)

saveRDS(df, file.path(SRC, "panel_a_umap.rds"))
write.csv(df, file.path(SRC, "panel_a_umap.csv"), row.names = FALSE)

df$bin_mid <- pmin(pmax(ceiling(df$Proteinfeatures), 1), 20)
df$bin <- factor(df$bin_mid)
prop_d <- df %>%
  count(bin, bin_mid, Tcelltype, name = "freq") %>%
  group_by(bin, bin_mid) %>%
  mutate(n = sum(freq), prop = freq / n, pct = 100 * prop) %>%
  ungroup() %>%
  rename(group = Tcelltype)
keep_bins <- unique(prop_d$bin[prop_d$n >= 20])
prop_d <- prop_d[prop_d$bin %in% keep_bins, , drop = FALSE]
write.csv(prop_d, file.path(SRC, "panel_c_bin_prop.csv"), row.names = FALSE)

cor_d <- prop_d %>%
  group_by(group) %>%
  summarise(
    spearman_rho = suppressWarnings(cor(bin_mid, pct, method = "spearman")),
    spearman_p = tryCatch(cor.test(bin_mid, pct, method = "spearman",
                                   exact = FALSE)$p.value, error = function(e) NA_real_),
    n_bins = n(),
    .groups = "drop"
  )
write.csv(cor_d, file.path(SRC, "panel_c_spearman.csv"), row.names = FALSE)
cat("Spearman\n"); print(as.data.frame(cor_d))

df$PF_class <- ifelse(df$Proteinfeatures >= HI, "High (PF>=6)",
               ifelse(df$Proteinfeatures <= LO, "Low (PF<=1)", "Mid"))
df$PF_class <- factor(df$PF_class,
                      levels = c("Low (PF<=1)", "Mid", "High (PF>=6)"))
# Panel d: subset oligo-CD45 tails on the reclustered embedding.
df_split <- df[df$PF_class %in% c("Low (PF<=1)", "High (PF>=6)"), ]
df_split$PF_class <- droplevels(df_split$PF_class)
saveRDS(df_split, file.path(SRC, "panel_d_umap_split.rds"))
write.csv(df_split, file.path(SRC, "panel_d_umap_split.csv"), row.names = FALSE)
cat("split n\n"); print(table(df_split$PF_class))

prop_e <- df_split %>%
  count(PF_class, Tcelltype, name = "n") %>%
  group_by(PF_class) %>%
  mutate(pct = 100 * n / sum(n), total = sum(n)) %>%
  ungroup()
write.csv(prop_e, file.path(SRC, "panel_e_proportion.csv"), row.names = FALSE)
saveRDS(prop_e, file.path(SRC, "panel_e_proportion.rds"))

n_tail <- as.data.frame(table(df_split$Tcelltype, df_split$PF_class),
                        stringsAsFactors = FALSE)
names(n_tail) <- c("group", "PF_class", "n")
n_low <- n_tail$n[n_tail$PF_class == "Low (PF<=1)"]
names(n_low) <- as.character(n_tail$group[n_tail$PF_class == "Low (PF<=1)"])
n_high <- n_tail$n[n_tail$PF_class == "High (PF>=6)"]
names(n_high) <- as.character(n_tail$group[n_tail$PF_class == "High (PF>=6)"])
cor_fun <- cor_d %>%
  filter(as.character(group) != "Unassigned") %>%
  mutate(n_low = unname(n_low[as.character(group)]),
         n_high = unname(n_high[as.character(group)])) %>%
  filter(n_low >= 50, n_high >= 50)
if (!nrow(cor_fun)) stop("no functional T class has >= 50 cells in both tails")
deg_class <- as.character(cor_fun$group[which.max(cor_fun$spearman_rho)])
cat("DEG class:", deg_class, "\n")
print(as.data.frame(cor_fun))

cells_sub <- df$cell[as.character(df$Tcelltype) == deg_class]
obj_sub <- obj[, cells_sub]
obj_sub$PF_tail <- ifelse(obj_sub$Proteinfeatures >= HI, "High",
                   ifelse(obj_sub$Proteinfeatures <= LO, "Low", "Mid"))
obj_sub <- obj_sub[, obj_sub$PF_tail %in% c("High", "Low")]
Idents(obj_sub) <- factor(obj_sub$PF_tail, levels = c("Low", "High"))
cat("DEG n Low/High:", paste(table(Idents(obj_sub)), collapse = "/"), "\n")
stopifnot(min(table(Idents(obj_sub))) >= 20)

DefaultAssay(obj_sub) <- "RNA"
deg <- FindMarkers(
  obj_sub, ident.1 = "High", ident.2 = "Low",
  test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25,
  verbose = FALSE
)
deg$gene <- rownames(deg)
deg <- deg[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
write.csv(deg, file.path(SRC, "panel_f_DEG_unfiltered.csv"), row.names = FALSE)
deg_f <- deg[is.finite(deg$avg_log2FC) & abs(deg$avg_log2FC) >= FC_CUT, ]
write.csv(deg_f, file.path(SRC, "panel_f_DEG.csv"), row.names = FALSE)
cat("DEG unfiltered", nrow(deg), " |log2FC|>=", FC_CUT, " ", nrow(deg_f),
    " padj<0.05 ", sum(deg_f$p_val_adj < 0.05, na.rm = TRUE), "\n", sep = "")

write.csv(data.frame(deg_class = deg_class, n_low = sum(obj_sub$PF_tail == "Low"),
                     n_high = sum(obj_sub$PF_tail == "High"),
                     stringsAsFactors = FALSE),
          file.path(SRC, "panel_f_deg_class.csv"), row.names = FALSE)

ann_lines <- c(
  "CITE T cells: QC drop then recluster then subset (Fig. S6 v2)",
  paste("Parent T object n =", n_before),
  paste("After dropping first-pass B-T doublet and myeloid-contaminated:", sum(keep)),
  paste("After recluster / residual QC subset:", ncol(obj)),
  "Recluster: NormalizeData, 2000 HVG, ScaleData, PCA 30, neighbors/UMAP dims 1:15, RNA_snn res 0.2",
  "UMAP panels a and d use this reclustered embedding; d is the high/low oligo-CD45 subset of it.",
  "",
  "New cluster -> Fig. 3 v4 taxonomy (cluster-level majority markers):"
)
for (i in seq_len(nrow(ann_tbl))) {
  ann_lines <- c(ann_lines, paste0(
    ann_tbl$cluster[i], "  n=", ann_tbl$n[i], "  ", ann_tbl$Tcelltype[i]
  ))
}
ann_lines <- c(ann_lines, paste("DEG class:", deg_class))
writeLines(ann_lines, file.path(SRC, "cluster_annotation.txt"))

writeLines(c(
  "CITE-seq T cells from the same Cd45-TotalSeqB experiment as Fig. 5 B cells.",
  "Pipeline: drop first-pass QC clusters, recluster remaining T cells, then subset for UMAP.",
  paste("n_parent_T", n_before),
  paste("n_after_first_QC", sum(keep)),
  paste("n_plotted", nrow(df)),
  paste("deg_class", deg_class),
  paste("HI", HI, "LO", LO, "FC_CUT", FC_CUT, "RES", RES)
), file.path(SRC, "README.txt"))

write.csv(data.frame(fun_lv = fun_lv, stringsAsFactors = FALSE),
          file.path(SRC, "display_levels.csv"), row.names = FALSE)

cat("DONE ->", SRC, "\n")
