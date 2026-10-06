#!/usr/bin/env Rscript
# Extract Visium H&E raster + aligned tissue coordinates for ggplot backgrounds
options(warn = 1)
suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: extract_he.R <in.rds> <out.rds>")

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

in_path <- args[[1]]
out_path <- args[[2]]
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

cat("Loading", in_path, "...\n")
sce <- fix_images(readRDS(in_path))
img_name <- names(sce@images)[[1]]
img_obj <- sce@images[[img_name]]
cat("image:", img_name, "class:", class(img_obj)[1], "\n")

# Raw histology array (rows x cols x RGB), values usually 0-1 or 0-255
he <- img_obj@image
if (max(he, na.rm = TRUE) > 1.5) he <- he / 255
he <- pmin(pmax(he, 0), 1)
cat("HE dim:", paste(dim(he), collapse = "x"), "\n")

# Scale factors for Visium
sf <- tryCatch(img_obj@scale.factors, error = function(e) NULL)
if (!is.null(sf)) {
  cat("scale.factors:", paste(names(sf), unlist(sf), sep = "=", collapse = ", "), "\n")
}

# Coordinates in the same system used by SpatialFeaturePlot / GetTissueCoordinates
coords_raw <- GetTissueCoordinates(sce)
cat("raw coord cols:", paste(colnames(coords_raw), collapse = ", "), "\n")

# For Visium, Seurat image is indexed [row, col, channel]
# Tissue coordinates typically imagecol (x) / imagerow (y) in full-res pixels,
# but plotted coordinates after GetTissueCoordinates are often already scaled.
# We store BOTH raw GetTissueCoordinates and a ggplot-aligned version.

cn <- colnames(coords_raw)
if (all(c("imagerow", "imagecol") %in% cn)) {
  # Full-res pixel coords
  x_full <- coords_raw$imagecol
  y_full <- coords_raw$imagerow
} else if (all(c("x", "y") %in% cn)) {
  x_full <- coords_raw$x
  y_full <- coords_raw$y
} else {
  num <- which(vapply(coords_raw, is.numeric, logical(1)))
  x_full <- coords_raw[[num[1]]]
  y_full <- coords_raw[[num[2]]]
}

# Determine whether coords are already in lowres image pixel space
nr <- dim(he)[1]
nc <- dim(he)[2]
# If coords exceed image dims substantially, scale by lowres factor
max_x <- max(x_full, na.rm = TRUE)
max_y <- max(y_full, na.rm = TRUE)
cat("coord max x/y:", max_x, max_y, "; HE ncol/nrow:", nc, nr, "\n")

scale_x <- 1
scale_y <- 1
if (!is.null(sf) && !is.null(sf$lowres)) {
  # If coords look like full-res, map to lowres image pixels
  if (max_x > nc * 1.2 || max_y > nr * 1.2) {
    scale_x <- sf$lowres
    scale_y <- sf$lowres
    cat("Applying lowres scale factor:", sf$lowres, "\n")
  }
}

x_img <- x_full * scale_x
y_img <- y_full * scale_y

# ggplot: x right, y up. Image row 1 is top -> flip y for annotation_raster
coords <- data.frame(
  x = x_img,
  y = -y_img,
  row.names = rownames(coords_raw)
)

# Image extent in the same flipped coordinate system:
# x from 0.5..nc+0.5 or 1..nc; use 0.5 to nc+0.5 centered pixels
# Seurat typically uses image pixels with origin top-left.
xmin <- 0.5
xmax <- nc + 0.5
ymin <- -(nr + 0.5)
ymax <- -0.5

he_out <- list(
  image = he,
  xmin = xmin,
  xmax = xmax,
  ymin = ymin,
  ymax = ymax,
  coords = coords,
  image_name = img_name,
  source = basename(in_path),
  dim = dim(he)
)
saveRDS(he_out, out_path)
cat("Saved", out_path, "\n")
