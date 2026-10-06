suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
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
raw_png <- file.path(outdir, "assets", "panel_a_raw_图片2.png")

# Lane order matches FigureS3.pdf gel labels (ladder implicit at left)
labs <- c("M1/42-Oligo", "30-F11-Oligo", "M1/42", "30-F11", "M1/42-Oligo", "30-F11-Oligo")
# approximate lane centers over sample lanes (skip MW ladder ~0.12)
xs <- c(0.28, 0.40, 0.52, 0.64, 0.76, 0.88)

p_ann <- ggdraw() +
  draw_image(raw_png, x = 0.02, y = 0.00, width = 0.96, height = 0.88) +
  # group brackets: 1 = lanes 1-4, 2 = lanes 5-6 (FigureS3 style)
  draw_line(x = c(0.24, 0.68), y = c(0.94, 0.94), color = "black", linewidth = 0.45) +
  draw_line(x = c(0.72, 0.92), y = c(0.94, 0.94), color = "black", linewidth = 0.45) +
  draw_label("1", x = 0.46, y = 0.965, size = 9, fontface = "bold", fontfamily = "Arial") +
  draw_label("2", x = 0.82, y = 0.965, size = 9, fontface = "bold", fontfamily = "Arial") +
  draw_label(labs[1], x = xs[1], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[2], x = xs[2], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[3], x = xs[3], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[4], x = xs[4], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[5], x = xs[5], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0) +
  draw_label(labs[6], x = xs[6], y = 0.915, size = 5.2, angle = 40, fontfamily = "Arial", hjust = 0)

out1 <- file.path(outdir, "assets", "panel_a_annotated_图片2.png")
out2 <- file.path(outdir, "source_data", "panel_a_schematic.png")
ggsave(out1, p_ann, width = 72, height = 115, units = "mm", dpi = 600, bg = "white")
file.copy(out1, out2, overwrite = TRUE)
cat("rewrote panel a labels\n")
