#!/usr/bin/env Rscript
# Extract lineage markers and ImmGen T/NK module scores for Fig. 3 v2.
# New folder only. Does not overwrite Fig3_Nature_v1.
options(warn = 1)
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

OUT <- file.path(IV_FIG_ROOT, "Fig3_Nature_v2")
SRC <- file.path(OUT, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
md <- read.csv(file.path(SRC, "cell_metadata_reannotated.csv"),
               stringsAsFactors = FALSE)
stopifnot(nrow(md) == 6098)

cat("Loading CD4T_reclustered.rds ...\n")
t4 <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/CD4T_reclustered.rds"))
DefaultAssay(t4) <- "RNA"
genes <- c("Cd4", "Cd8a", "Cd8b1", "Klrb1c", "Ms4a1", "Cd3e")
stopifnot(all(genes %in% rownames(t4)))
stopifnot(all(md$barcode %in% colnames(t4)))

get_layer <- function(layer) {
  tryCatch(GetAssayData(t4, assay = "RNA", layer = layer),
           error = function(e) GetAssayData(t4, assay = "RNA"))
}
dat <- get_layer("data")
expr <- as.data.frame(t(as.matrix(dat[genes, md$barcode])))
colnames(expr) <- genes

# ImmGen identity modules: same rule as Fig. 2d / P0-4
# top 30 genes with |log2FC| >= 1 on ImmGen pairwise axes vs naive CD4.
pw_path <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS1b_reannotate_ImmGen_v1/source_data/ImmGen_pairwise_allgenes.csv")
pw <- read.csv(pw_path, stringsAsFactors = FALSE)
drop_pat <- "^(Rp[sl]|Mrpl|Mrps|mt-|Mt-|Gm[0-9]|Igh[vdj]|Igk[vdj]|Igl[vdj])"
rna <- rownames(t4[["RNA"]])
pick_mod <- function(contrast, direction, n = 30) {
  d <- pw[pw$contrast == contrast & pw$direction == direction, ]
  d <- d[is.finite(d$log2FC) & abs(d$log2FC) >= 1, ]
  d <- d[!grepl(drop_pat, d$gene), ]
  d <- d[order(-abs(d$log2FC)), ]
  intersect(d$gene, rna)[seq_len(min(n, length(intersect(d$gene, rna))))]
}
sig <- list(
  Naive_CD4 = pick_mod("CD8naive_vs_CD4naive", "CD4_naive up"),
  Naive_CD8 = pick_mod("CD8naive_vs_CD4naive", "CD8_naive up"),
  NKT       = pick_mod("NKT_vs_CD4naive", "NKT up"),
  NK        = pick_mod("NKmature_vs_CD4naive", "NK_mature up")
)
sig_tbl <- bind_rows(lapply(names(sig), function(m) {
  data.frame(module = m, gene = sig[[m]], stringsAsFactors = FALSE)
}))
print(table(sig_tbl$module))
write.csv(sig_tbl, file.path(SRC, "SourceData_Fig3_ImmGen_module_genes.csv"),
          row.names = FALSE)

# Score on the same barcode order as metadata
cells_use <- md$barcode
obj <- subset(t4, cells = cells_use)
obj <- AddModuleScore(obj, features = unname(sig), name = "IGmod", seed = 42)
mod_cols <- paste0("IGmod", seq_along(sig))
stopifnot(all(mod_cols %in% colnames(obj[[]])))
mod <- as.data.frame(obj[[]][cells_use, mod_cols, drop = FALSE])
colnames(mod) <- names(sig)

df <- cbind(
  data.frame(
    cell = md$barcode,
    UMAP1 = md$UMAP_1,
    UMAP2 = md$UMAP_2,
    t_class = factor(md$t_class, levels = CLASS_LV),
    IV_group = factor(md$IV_group, levels = c("IntraNeg", "IntraPos")),
    cluster_ref = md$cluster,
    stringsAsFactors = FALSE
  ),
  expr,
  mod
)
stopifnot(nrow(df) == 6098)
print(table(df$t_class, df$IV_group))

write.csv(df, file.path(SRC, "SourceData_Fig3_T_cells.csv"), row.names = FALSE)
saveRDS(df, file.path(SRC, "Fig3_T_cells.rds"))

sc <- df %>%
  tidyr::pivot_longer(all_of(names(sig)), names_to = "module", values_to = "score")
write.csv(sc[, c("cell", "t_class", "IV_group", "module", "score")],
          file.path(SRC, "SourceData_Fig3_module_scores_per_cell.csv"),
          row.names = FALSE)
cat("Wrote", nrow(df), "cells\nDONE\n")
