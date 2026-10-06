#!/usr/bin/env Rscript
# Per-panel crops of Figure10.png for visual QA. R / magick only.
suppressPackageStartupMessages({ library(magick) })
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

OUT <- file.path(IV_FIG_ROOT, "Fig10_Nature_v1")
QA  <- file.path(OUT, "qa")

FIG_W_MM <- 183
R1 <- 42; R2 <- 49; R3 <- 54; R4 <- 66; PAD_B <- 3
FIG_H_MM <- R1 + R2 + R3 + R4 + PAD_B

img <- image_read(file.path(OUT, "Figure10.png"))
info <- image_info(img)
PX_MM_X <- info$width / FIG_W_MM
PX_MM_Y <- info$height / FIG_H_MM
cat("page", info$width, "x", info$height, "px for", FIG_W_MM, "x", FIG_H_MM, "mm\n")

crop_mm <- function(name, x, y, w, h, out_w = 1600) {
  g <- sprintf("%dx%d+%d+%d", round(w * PX_MM_X), round(h * PX_MM_Y),
               round(x * PX_MM_X), round(y * PX_MM_Y))
  im <- image_crop(img, g)
  image_write(image_resize(im, paste0(out_w, "x")), file.path(QA, paste0(name, ".png")))
  cat(sprintf("  %-18s %6.1f x %6.1f mm\n", name, w, h))
}

Y1 <- 0; Y2 <- R1; Y3 <- R1 + R2; Y4 <- R1 + R2 + R3

crop_mm("panel_a", 0,   Y1, 52,  R1)
crop_mm("panel_b", 51,  Y1, 132, R1)
crop_mm("panel_ab", 0,  Y1, 183, R1)
crop_mm("panel_c", 0,   Y2, 52,  R2)
crop_mm("panel_d", 51,  Y2, 63,  R2)
crop_mm("panel_e", 113, Y2, 70,  R2)
crop_mm("panel_cde", 0, Y2, 183, R2)
crop_mm("panel_f", 0,   Y3, 183, R3)
crop_mm("panel_g", 0,   Y4, 65,  R4 + PAD_B)
crop_mm("panel_h", 64,  Y4, 53,  R4 + PAD_B)
crop_mm("panel_i", 116, Y4, 67,  R4 + PAD_B)
crop_mm("panel_ghi", 0, Y4, 183, R4 + PAD_B)

image_write(image_resize(img, "1400x"), file.path(QA, "_full_preview.png"))
cat("done\n")
