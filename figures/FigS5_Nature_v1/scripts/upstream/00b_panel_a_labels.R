suppressPackageStartupMessages({
  library(magick)
  library(ggplot2)
  library(cowplot)
  library(png)
  library(grid)
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

outdir <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-1_NC_PFsplit")
dir.create(file.path(outdir, "assets"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "source_data"), showWarnings = FALSE, recursive = TRUE)

# --- 1) Tight crop of labeled gel from Nature FigureS3 ---
img <- image_read(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Figures_Nature_R_20260824/FigureS3/FigureS3.png"))
info <- image_info(img)
w <- info$width
h <- info$height
# iterative tight crop: gel only (exclude violin panels)
crop <- image_crop(img, geometry = sprintf("%dx%d+%d+%d",
                                           as.integer(w * 0.335),
                                           as.integer(h * 0.295),
                                           as.integer(w * 0.008),
                                           as.integer(h * 0.008)))
image_write(crop, file.path(outdir, "assets", "panel_a_labeled_crop.png"))
print(image_info(crop))

# --- 2) Also build annotated 图片2 (dual gel) with lane labels ---
fs <- list.files(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_Bcelltype"),
                 pattern = "\\.tif$", full.names = TRUE)
idx <- which(file.info(fs)$size > 600000)[1]
gel2 <- image_read(fs[idx])
image_write(gel2, file.path(outdir, "assets", "panel_a_raw_图片2.png"), format = "png")

# Prefer the Nature-labeled crop if it looks gel-only; copy to source_data
file.copy(file.path(outdir, "assets", "panel_a_labeled_crop.png"),
          file.path(outdir, "source_data", "panel_a_schematic.png"),
          overwrite = TRUE)

# Build a clean ggplot annotation plate on 图片2 as fallback / improved a
raw_png <- file.path(outdir, "assets", "panel_a_raw_图片2.png")
img_r <- readPNG(raw_png)
# lane labels matching FigureS3 panel a (ladder + 6 samples)
labs <- c("M1/42-Oligo", "30-F11-Oligo", "M1/42", "30-F11", "M1/42-Oligo", "30-F11-Oligo")
# x positions for 6 sample lanes (ladder is leftmost); fractions of panel width
# MW labels already in image; place angled labels above sample lanes
xs <- seq(0.22, 0.92, length.out = 6)
lab_df <- data.frame(x = xs, y = 0.985, lab = labs, stringsAsFactors = FALSE)
grp <- data.frame(
  x = c(0.22, 0.92, 0.37, 0.77),
  y = c(1.03, 1.03, 1.06, 1.06),
  lab = c("", "", "1", "2"),
  stringsAsFactors = FALSE
)
# bracket segments
brk <- data.frame(
  x = c(0.18, 0.48, 0.55, 0.95),
  xend = c(0.48, 0.48, 0.95, 0.95),
  y = c(1.02, 1.02, 1.02, 1.02),
  yend = c(1.02, 1.045, 1.02, 1.045)
)

p_ann <- ggdraw() +
  draw_image(raw_png, x = 0, y = 0, width = 1, height = 0.92) +
  draw_line(x = c(0.18, 0.48), y = c(0.955, 0.955), color = "black", size = 0.4) +
  draw_line(x = c(0.55, 0.95), y = c(0.955, 0.955), color = "black", size = 0.4) +
  draw_label("1", x = 0.33, y = 0.975, size = 8, fontface = "bold", fontfamily = "Arial") +
  draw_label("2", x = 0.75, y = 0.975, size = 8, fontface = "bold", fontfamily = "Arial") +
  draw_label(labs[1], x = xs[1], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[2], x = xs[2], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[3], x = xs[3], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[4], x = xs[4], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[5], x = xs[5], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[6], x = xs[6], y = 0.935, size = 4.5, angle = 40, fontfamily = "Arial", hjust = 0)

ggsave(file.path(outdir, "assets", "panel_a_annotated_图片2.png"), p_ann,
       width = 70, height = 110, units = "mm", dpi = 600, bg = "white")
# Use annotated 图片2 as primary panel a (full dual gel + labels)
file.copy(file.path(outdir, "assets", "panel_a_annotated_图片2.png"),
          file.path(outdir, "source_data", "panel_a_schematic.png"),
          overwrite = TRUE)
cat("panel a written\n")
