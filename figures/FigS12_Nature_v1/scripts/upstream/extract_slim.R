#!/usr/bin/env Rscript
# Extract slim plotting payloads from Seurat RDS (coords + pred + genes)
options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: extract_slim.R <in.rds> <out.rds>")

fix_images <- function(sce) {
  for (nm in names(sce@images)) {
    img <- sce@images[[nm]]
    ok <- tryCatch({ slot(img, "misc"); TRUE }, error = function(e) FALSE)
    if (!ok) {
      slot(img, "misc", check = FALSE) <- list()
      sce@images[[nm]] <- img
    }
  }
  sce
}

in_path  <- args[[1]]
out_path <- args[[2]]
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

cat("Loading", in_path, "...\n")
sce <- fix_images(readRDS(in_path))
cat("dim:", paste(dim(sce), collapse = "x"), "\n")

coords <- GetTissueCoordinates(sce)
cn <- colnames(coords)
if (all(c("imagerow", "imagecol") %in% cn)) {
  coords <- data.frame(
    x = coords$imagecol, y = -coords$imagerow,
    row.names = rownames(coords)
  )
} else if (all(c("x", "y") %in% cn)) {
  coords <- coords[, c("x", "y"), drop = FALSE]
} else {
  num <- vapply(coords, is.numeric, logical(1))
  coords <- data.frame(
    x = coords[[which(num)[1]]], y = coords[[which(num)[2]]],
    row.names = rownames(coords)
  )
}
coords$y <- -coords$y

pred <- as.matrix(sce[["predictions"]]@data)
DefaultAssay(sce) <- "Spatial"
genes <- c("Lgals3", "Prdx1", "Mki67")
miss <- setdiff(genes, rownames(sce))
if (length(miss)) stop("Missing genes: ", paste(miss, collapse = ", "))
expr <- as.matrix(GetAssayData(sce, assay = "Spatial", layer = "data")[genes, , drop = FALSE])

need_pred <- c("n-5", "p-5", "n-3", "p-3", "n-0", "p-0")
miss_p <- setdiff(need_pred, rownames(pred))
if (length(miss_p)) stop("Missing predictions: ", paste(miss_p, collapse = ", "))

slim <- list(
  coords = coords,
  pred = pred[need_pred, , drop = FALSE],
  expr = expr,
  n_spots = ncol(sce),
  source = basename(in_path)
)
saveRDS(slim, out_path)
cat("Saved", out_path, "spots=", slim$n_spots, "\n")
