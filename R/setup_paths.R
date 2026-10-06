# Shared path helpers for IVstaining_scRNA figure scripts.
# Each plotting script also embeds a copy of this logic so that
# `Rscript figures/<Fig>/scripts/make_figure*.R` works from any cwd.

iv_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)[1])
  if (!length(f) || is.na(f) || !nzchar(f)) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  dirname(normalizePath(f, winslash = "/", mustWork = TRUE))
}

iv_figure_dir <- function(script_dir = iv_script_dir()) {
  d <- script_dir
  if (basename(d) %in% c("upstream", "qa", "original_scripts")) d <- dirname(d)
  if (identical(basename(d), "scripts")) dirname(d) else d
}

iv_repo_root <- function(start = iv_figure_dir()) {
  d <- start
  repeat {
    if (dir.exists(file.path(d, "figures")) &&
        file.exists(file.path(d, "README.md"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Could not find the IVstaining_scRNA repository root from ", start)
    }
    d <- parent
  }
}

iv_fig_root <- function() file.path(iv_repo_root(), "figures")
