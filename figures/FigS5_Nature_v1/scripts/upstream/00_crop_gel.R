suppressPackageStartupMessages(library(magick))
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
img <- image_read(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Figures_Nature_R_20260824/FigureS3/FigureS3.png"))
print(image_info(img))
info <- image_info(img)
w <- info$width
h <- info$height
# top-left gel (panel a) with lane labels
crop <- image_crop(img, geometry = sprintf("%dx%d+%d+%d",
                                           as.integer(w * 0.40),
                                           as.integer(h * 0.36),
                                           0, 0))
outdir <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-1_NC_PFsplit")
dir.create(file.path(outdir, "assets"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "scripts"), recursive = TRUE, showWarnings = FALSE)
image_write(crop, file.path(outdir, "assets", "panel_a_from_FigS3_crop1.png"))
# also try slightly different crops for QA
crop2 <- image_crop(img, geometry = sprintf("%dx%d+%d+%d",
                                            as.integer(w * 0.38),
                                            as.integer(h * 0.34),
                                            as.integer(w * 0.01),
                                            as.integer(h * 0.01)))
image_write(crop2, file.path(outdir, "assets", "panel_a_from_FigS3_crop2.png"))
cat("done crops\n")
