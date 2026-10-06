#!/usr/bin/env Rscript
# Extract table-marker RNA for Fig. S3 v4 functional identities.
# New folder only. Does not overwrite FigS3_Nature_v3.
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUT <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS3_Nature_v4")
}
SRC <- file.path(OUT, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

cells <- readRDS(file.path(SRC, "Fig3_T_cells.rds"))
stopifnot(nrow(cells) == 6098)

genes <- c(
  "Cd4", "Cd3d", "Cd3e", "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r",
  "S100a4", "Cxcr3", "Cxcr6",
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2",
  "Cd8b1", "Cd44", "Eomes", "Il2rb",
  "Nkg7", "Gzma", "Ccl5", "Prf1", "Gzmb",
  "Trac", "Klrb1c",
  "Trdc", "Trgc1", "Trgc2",
  "Ms4a1", "Cd79a", "Lyz2", "S100a8"
)

cat("Loading CD4T_reclustered.rds ...\n")
t4 <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/CD4T_reclustered.rds"))
DefaultAssay(t4) <- "RNA"
present <- genes[genes %in% rownames(t4)]
missing <- setdiff(genes, present)
if (length(missing)) cat("missing genes:", paste(missing, collapse = ", "), "\n")
stopifnot(length(present) >= 20)
stopifnot(all(cells$cell %in% colnames(t4)))

dat <- tryCatch(GetAssayData(t4, assay = "RNA", layer = "data"),
                error = function(e) GetAssayData(t4, assay = "RNA"))
expr <- as.data.frame(t(as.matrix(dat[present, cells$cell])))
# Keep UMAP / identity columns; replace overlapping gene columns.
id_cols <- c("cell", "UMAP1", "UMAP2", "t_class", "subtype", "IV_group", "cluster_ref")
out <- cbind(cells[, id_cols, drop = FALSE], expr)
stopifnot(nrow(out) == 6098)
saveRDS(out, file.path(SRC, "FigS3_annot_expr.rds"))
write.csv(out, file.path(SRC, "SourceData_FigS3_annot_expr.csv"), row.names = FALSE)
writeLines(c(
  paste("requested", paste(genes, collapse = ",")),
  paste("present", paste(present, collapse = ",")),
  paste("missing", paste(missing, collapse = ","))
), file.path(SRC, "annot_marker_genes.txt"))
cat("Wrote", nrow(out), "cells x", length(present), "genes\nDONE\n")
