#!/usr/bin/env Rscript
# Extract subtype labels and nine ImmGen T/NK module scores for Fig. 3 v3.
# New folder only. Does not overwrite Fig3_Nature_v2.
options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
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
  file.path(IV_FIG_ROOT, "Fig3_Nature_v4")
}
SRC <- file.path(OUT, "source_data")
S3  <- file.path(IV_FIG_ROOT, "FigS3_Nature_v1")
S3S <- file.path(S3, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)
dir.create(S3S, showWarnings = FALSE, recursive = TRUE)

CLASS_LV <- c("true_CD4", "true_CD8", "NKT", "DP", "unassigned", "doublet_c2")
SUB_LV <- c(
  "CD8-0", "CD4-1a", "CD4-1b", "Doublet", "CD8-3a", "CD8-3b",
  "CD4-other", "CD8-other", "NKT", "DP", "Unassigned"
)

md <- read.csv(file.path(SRC, "cell_metadata_reannotated.csv"),
               stringsAsFactors = FALSE)
stopifnot(nrow(md) == 6098)

assign_subtype <- function(t_class, cluster) {
  n <- length(t_class)
  out <- character(n)
  for (i in seq_len(n)) {
    tc <- t_class[i]; cl <- as.character(cluster[i])
    if (tc == "true_CD4") {
      out[i] <- if (cl %in% c("1a", "1b")) paste0("CD4-", cl) else "CD4-other"
    } else if (tc == "true_CD8") {
      out[i] <- if (cl %in% c("0", "3a", "3b")) paste0("CD8-", cl) else "CD8-other"
    } else if (tc == "doublet_c2") {
      out[i] <- "Doublet"
    } else if (tc == "unassigned") {
      out[i] <- "Unassigned"
    } else {
      out[i] <- tc
    }
  }
  out
}
md$subtype <- assign_subtype(md$t_class, md$cluster)
md$subtype <- factor(md$subtype, levels = SUB_LV)
stopifnot(!anyNA(md$subtype))
cat("subtype x class\n")
print(table(md$subtype, md$t_class))
cat("subtype x cluster\n")
print(table(md$subtype, md$cluster))

write.csv(
  as.data.frame.matrix(table(md$subtype, md$t_class)),
  file.path(SRC, "subtype_by_class.csv")
)
write.csv(
  as.data.frame.matrix(table(md$subtype, md$cluster)),
  file.path(SRC, "subtype_by_cluster.csv")
)

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

# ImmGen identity modules for every heatmap column (GSE109125 spleen T/NK).
# Naive CD4 vs naive CD8; every other subset vs naive CD4.
# Same filters as Fig. 2d / Fig. 3 v2: |log2FC| >= 1, up-side mean >= 10,
# drop ribo/mt/Gm/Ig, top 30 genes present in the RNA assay.
GZ <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS1b_reannotate_ImmGen_v1/source_data/GSE109125_Genes_count_table.tsv.gz")
stopifnot(file.exists(GZ))
imm <- read.delim(gzfile(GZ), check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- colnames(imm)[1]
rownames(imm) <- make.unique(as.character(imm[[gene_col]]))
imm[[gene_col]] <- NULL
mat <- as.matrix(imm)
storage.mode(mat) <- "numeric"
keep_pat <- c(
  "^T\\.4\\.Nve\\.Sp#",
  "^T\\.4\\.19-8-TCRb\\+CD4\\+\\.Sp#",
  "^Treg\\.4\\.25hi\\.Sp#",
  "^T\\.8\\.Nve\\.Sp#",
  "^NKT\\.19-8-TCRb\\+CD1daGalCerTet\\+\\.Sp#",
  "^Tgd\\.Sp#",
  "^NK\\.27\\+11b-\\.Sp#",
  "^NK\\.27\\+11b\\+\\.Sp#",
  "^NK\\.27-11b\\+\\.Sp#"
)
keep_cols <- unique(unlist(lapply(keep_pat, function(p) grep(p, colnames(mat), value = TRUE))))
cat("kept ImmGen columns:", length(keep_cols), "\n")
mat <- mat[, keep_cols, drop = FALSE]
prefix_of <- function(x) sub("#\\d+$", "", x)
pref <- prefix_of(colnames(mat))
merged <- sapply(unique(pref), function(p) rowMeans(mat[, pref == p, drop = FALSE]))
name_map <- c(
  "T.4.Nve.Sp" = "CD4_naive",
  "T.4.19-8-TCRb+CD4+.Sp" = "CD4_total",
  "Treg.4.25hi.Sp" = "Treg",
  "T.8.Nve.Sp" = "CD8_naive",
  "NKT.19-8-TCRb+CD1daGalCerTet+.Sp" = "NKT",
  "Tgd.Sp" = "Tgd",
  "NK.27+11b-.Sp" = "NK_immature",
  "NK.27+11b+.Sp" = "NK_int",
  "NK.27-11b+.Sp" = "NK_mature"
)
colnames(merged) <- ifelse(colnames(merged) %in% names(name_map),
                           name_map[colnames(merged)], colnames(merged))
col_order <- c("CD4_naive", "CD4_total", "Treg", "CD8_naive",
               "NKT", "Tgd", "NK_immature", "NK_int", "NK_mature")
stopifnot(all(col_order %in% colnames(merged)))
merged <- merged[, col_order, drop = FALSE]
write.csv(data.frame(immgen_col = keep_cols, prefix = prefix_of(keep_cols),
                     mapped = unname(name_map[prefix_of(keep_cols)]),
                     stringsAsFactors = FALSE),
          file.path(SRC, "ImmGen_column_map.csv"), row.names = FALSE)

FC_ID <- 1.0
MIN_MEAN <- 10
drop_pat <- "^(Rp[sl]|Mrpl|Mrps|mt-|Mt-|Gm[0-9]|Igh[vdj]|Igk[vdj]|Igl[vdj])"
rna <- rownames(t4[["RNA"]])
mod_spec <- list(
  Naive_CD4    = c("CD4_naive", "CD8_naive"),
  Total_CD4    = c("CD4_total", "CD4_naive"),
  Treg         = c("Treg", "CD4_naive"),
  Naive_CD8    = c("CD8_naive", "CD4_naive"),
  NKT          = c("NKT", "CD4_naive"),
  gdT          = c("Tgd", "CD4_naive"),
  NK_immature  = c("NK_immature", "CD4_naive"),
  NK_int       = c("NK_int", "CD4_naive"),
  NK_mature    = c("NK_mature", "CD4_naive")
)
pick_mod <- function(a, b, n = 30) {
  m1 <- merged[, a]; m2 <- merged[, b]
  lfc <- log2((m1 + 1) / (m2 + 1))
  keep <- is.finite(lfc) & lfc >= FC_ID & m1 >= MIN_MEAN
  g <- names(m1)[keep]
  g <- g[!grepl(drop_pat, g)]
  g <- g[order(-lfc[g])]
  hit <- intersect(g, rna)
  hit[seq_len(min(n, length(hit)))]
}
sig <- lapply(mod_spec, function(ab) pick_mod(ab[1], ab[2]))
names(sig) <- names(mod_spec)
sig_tbl <- bind_rows(lapply(names(sig), function(m) {
  ab <- mod_spec[[m]]
  m1 <- merged[, ab[1]]; m2 <- merged[, ab[2]]
  lfc <- log2((m1 + 1) / (m2 + 1))
  data.frame(module = m, contrast = paste0(ab[1], "_vs_", ab[2]),
             gene = sig[[m]], log2FC = unname(lfc[sig[[m]]]),
             stringsAsFactors = FALSE)
}))
print(table(sig_tbl$module))
stopifnot(all(lengths(sig) >= 10))
write.csv(sig_tbl, file.path(SRC, "SourceData_Fig3_ImmGen_module_genes.csv"),
          row.names = FALSE)

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
    subtype = factor(md$subtype, levels = SUB_LV),
    IV_group = factor(md$IV_group, levels = c("IntraNeg", "IntraPos")),
    cluster_ref = md$cluster,
    stringsAsFactors = FALSE
  ),
  expr,
  mod
)
stopifnot(nrow(df) == 6098)
print(table(df$subtype, df$IV_group))

write.csv(df, file.path(SRC, "SourceData_Fig3_T_cells.csv"), row.names = FALSE)
saveRDS(df, file.path(SRC, "Fig3_T_cells.rds"))

comp <- df %>%
  count(subtype, IV_group, name = "n") %>%
  group_by(IV_group) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup()
write.csv(comp, file.path(SRC, "n_by_subtype_IV.csv"), row.names = FALSE)

sc <- df %>%
  pivot_longer(all_of(names(sig)), names_to = "module", values_to = "score")
write.csv(sc[, c("cell", "t_class", "subtype", "IV_group", "module", "score")],
          file.path(SRC, "SourceData_Fig3_module_scores_per_cell.csv"),
          row.names = FALSE)

file.copy(file.path(SRC, "Fig3_T_cells.rds"),
          file.path(S3S, "Fig3_T_cells.rds"), overwrite = TRUE)
file.copy(file.path(SRC, "cell_metadata_reannotated.csv"),
          file.path(S3S, "cell_metadata_reannotated.csv"), overwrite = TRUE)
file.copy(file.path(SRC, "SourceData_Fig3_T_cells.csv"),
          file.path(S3S, "SourceData_Fig3_T_cells.csv"), overwrite = TRUE)
cat("Wrote", nrow(df), "cells; copied T-cell table to FigS3\nDONE\n")
