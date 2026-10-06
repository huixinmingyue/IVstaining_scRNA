#!/usr/bin/env Rscript
# Export plotted-ready tables for Fig. S5 v1. New folder only.
options(warn = 1)
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
  file.path(IV_FIG_ROOT, "FigS5_Nature_v1")
}
SRC <- file.path(OUT, "source_data")
SDP <- file.path(OUT, "source_data_plotted")
dir.create(SDP, showWarnings = FALSE)

need <- c("panel_a_raw.png", "panel_c_dotplot.csv", "panel_d_logistic.csv",
          "panel_b_umap.rds")
for (f in need) {
  p <- file.path(SRC, f)
  if (!file.exists(p)) stop("Missing ", p)
  cat("ok", f, "\n")
}

um <- readRDS(file.path(SRC, "panel_b_umap.rds"))
stopifnot(all(c("umap_1", "umap_2", "Bcelltype") %in% names(um)))
write.csv(um, file.path(SRC, "panel_b_umap.csv"), row.names = FALSE)
n_tab <- as.data.frame(table(um$Bcelltype), stringsAsFactors = FALSE)
names(n_tab) <- c("Bcelltype", "n")
n_tab$n <- as.integer(n_tab$n)
write.csv(n_tab, file.path(SDP, "SourceData_FigS5_n_Bcelltype.csv"), row.names = FALSE)
cat("UMAP cells:", nrow(um), "\n")
print(n_tab)
cat("DONE\n")
