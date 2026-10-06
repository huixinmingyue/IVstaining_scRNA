#!/usr/bin/env Rscript
# Verify Fig. S4 v2 source tables. New folder only.
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
  file.path(IV_FIG_ROOT, "FigS4_Nature_v3")
}
SRC <- file.path(OUT, "source_data")
need <- c(
  "SourceData_section_markers_radial.csv",
  "SourceData_radial_rnorm.csv",
  "SourceData_transfer_scores.csv"
)
for (f in need) {
  p <- file.path(SRC, f)
  if (!file.exists(p)) stop("Missing ", p)
  cat("ok", f, "\n")
}
map <- data.frame(
  source_map_title = c("true CD4 1a (Maf-high)", "true CD4 1b (Tcf7-high)",
                       "TNK CD8 cluster 0", "TNK CD8 cluster 3a",
                       "TNK CD8 cluster 3b"),
  transfer_col = c("T.CD4_1a", "T.CD4_1b", "T.CD8_0", "T.CD8_3a", "T.CD8_3b"),
  fig3_v3_label = c("CD4 1a", "CD4 1b", "CD8 0", "CD8 3a", "CD8 3b"),
  fig3_v4_label = c("Memory/activated CD4 T", "Naive CD4 T",
                    "Naive CD8 T", "Memory-like CD8 T",
                    "Effector/cytotoxic CD8 T"),
  s4_role = c("c", "c", "c", "c", "c"),
  stringsAsFactors = FALSE
)
write.csv(map, file.path(SRC, "cluster_name_map_Fig3v4.csv"), row.names = FALSE)
cat("DONE\n")
