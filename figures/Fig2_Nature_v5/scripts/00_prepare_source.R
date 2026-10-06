#!/usr/bin/env Rscript
# Extract B-cell UMAP coordinates and Cr2 / Fcer2a expression for Fig. 2 a-b.
# New folder only. Does not overwrite Fig2_Nature_v4.
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

OUT <- file.path(IV_FIG_ROOT, "Fig2_Nature_v5")
SRC <- file.path(OUT, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

cat("Loading Bcell_reclustered.rds ...\n")
b <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/Bcell_reclustered.rds"))
DefaultAssay(b) <- "RNA"
stopifnot(all(c("Cr2", "Fcer2a") %in% rownames(b)))

emb <- as.data.frame(Embeddings(b, "umap"))
colnames(emb) <- c("UMAP1", "UMAP2")
md <- b@meta.data
get_layer <- function(layer) {
  tryCatch(GetAssayData(b, assay = "RNA", layer = layer),
           error = function(e) GetAssayData(b, assay = "RNA"))
}
dat <- get_layer("data")

df <- data.frame(
  cell = colnames(b),
  UMAP1 = emb$UMAP1,
  UMAP2 = emb$UMAP2,
  cluster = factor(as.character(md$seurat_clusters), levels = c("0", "1", "2")),
  IV_group = factor(md$IV_group, levels = c("IntraNeg", "IntraPos")),
  Cr2 = as.numeric(dat["Cr2", colnames(b)]),
  Fcer2a = as.numeric(dat["Fcer2a", colnames(b)]),
  stringsAsFactors = FALSE
)
stopifnot(nrow(df) == 5975)
print(table(df$cluster, df$IV_group))

write.csv(df, file.path(SRC, "SourceData_Fig2ab_cells.csv"), row.names = FALSE)
saveRDS(df, file.path(SRC, "Fig2ab_cells.rds"))
cat("Wrote", nrow(df), "cells -> source_data/SourceData_Fig2ab_cells.csv\n")
cat("DONE\n")
