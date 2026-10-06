#!/usr/bin/env Rscript
# One-time extract of log-normalised gene expression for Fig. S2a.
# Writes source_data/FigS2a_violin_expr.csv so make_figureS2.R does not
# need Seurat at plot time.

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
ROOT <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  normalizePath(getwd())
}

rds <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/Bcell_reclustered.rds")
stopifnot(file.exists(rds))

p3_genes <- c("Cr2", "Fcer2a", "Cd93", "Cd24a", "Mki67", "Aicda", "Bcl6",
              "Sdc1", "Jchain", "Xbp1", "Cd1d1", "Ighd", "Ighm", "Sell", "Adgre1")
p6_genes <- c("Adgre5", "S1pr1", "Cxcr5", "Itgb1", "Itgb2", "Notch2",
              "Cr2", "Fcer2a", "Cd1d1", "Dtx1")
genes <- unique(c(p3_genes, p6_genes))

b <- readRDS(rds)
stopifnot(all(genes %in% rownames(b)))

fd <- FetchData(b, vars = c(genes, "seurat_clusters"))
fd$cell <- rownames(fd)
n_cell <- nrow(fd)
cat(sprintf("cells: %d; genes: %d; clusters:\n", n_cell, length(genes)))
print(table(fd$seurat_clusters))

wide <- fd[, c("cell", "seurat_clusters", genes)]
out_wide <- file.path(ROOT, "source_data", "FigS2a_violin_expr_wide.csv")
write.csv(wide, out_wide, row.names = FALSE)

long <- wide %>%
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
  mutate(
    seurat_clusters = as.character(seurat_clusters),
    block = ifelse(gene %in% p3_genes, "Query markers", "Positioning genes")
  )
# Shared genes Cr2, Fcer2a, Cd1d1 stay in Query markers (P3a), not duplicated.
stopifnot(n_distinct(long$cell) == n_cell)
stopifnot(nrow(long) == n_cell * length(genes))

out_long <- file.path(ROOT, "source_data", "FigS2a_violin_expr.csv")
write.csv(long, out_long, row.names = FALSE)
cat("wrote", out_long, "\n")
cat("wrote", out_wide, "\n")
