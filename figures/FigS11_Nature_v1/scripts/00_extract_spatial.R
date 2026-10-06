#!/usr/bin/env Rscript
# One-time extract of the plotted Visium matrices.
# Does not draw a figure. Reads the Seurat object used for the source plate
# and writes coordinates, prediction scores and four-gene expression.

suppressPackageStartupMessages(library(Seurat))
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
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else {
  normalizePath(file.path(getwd(), ".."))
}
SRC <- file.path(ROOT, "source_data")
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)

sce_path <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig7_S7_redraw/data/sce_immune5.rds")
sce <- readRDS(sce_path)
coords <- GetTissueCoordinates(sce)
stopifnot(all(c("imagecol", "imagerow") %in% colnames(coords)))

# Match make_figureS8.R: x = imagecol; y = -imagerow, then y = -y, so y = imagerow.
bc <- rownames(coords)
coord_df <- data.frame(
  barcode = bc,
  x = coords$imagecol,
  y = coords$imagerow,
  stringsAsFactors = FALSE
)
write.csv(coord_df, file.path(SRC, "FigS11_coords.csv"), row.names = FALSE)

pred <- as.matrix(sce[["predictions"]]@data)
stopifnot(all(bc %in% colnames(pred)))
pred <- pred[, bc, drop = FALSE]
pred_df <- as.data.frame(t(pred), stringsAsFactors = FALSE)
pred_df <- data.frame(barcode = rownames(pred_df), pred_df, check.names = FALSE)
write.csv(pred_df, file.path(SRC, "FigS11_predictions.csv"), row.names = FALSE)

genes <- c("Lgals3", "Prdx1", "Mki67", "Cxcl10")
DefaultAssay(sce) <- "Spatial"
expr <- tryCatch(
  GetAssayData(sce, assay = "Spatial", layer = "data"),
  error = function(e) GetAssayData(sce, assay = "Spatial", slot = "data")
)
missing <- setdiff(genes, rownames(expr))
if (length(missing)) stop("Missing genes: ", paste(missing, collapse = ", "))
expr <- expr[genes, bc, drop = FALSE]
expr_df <- as.data.frame(t(as.matrix(expr)), stringsAsFactors = FALSE)
expr_df <- data.frame(barcode = rownames(expr_df), expr_df, check.names = FALSE)
write.csv(expr_df, file.path(SRC, "FigS11_expression.csv"), row.names = FALSE)

cat("spots", nrow(coord_df), "\n")
cat("pred features", paste(setdiff(colnames(pred_df), "barcode"), collapse = ","), "\n")
cat("expr genes", paste(genes, collapse = ","), "\n")
cat("x", paste(range(coord_df$x), collapse = " "), " y", paste(range(coord_df$y), collapse = " "), "\n")
cat("wrote", SRC, "\n")
