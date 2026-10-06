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
b <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_noC4_PF_PT_corr/data/CITE_Bcell_noC4_annotated.rds"))
md <- slot(b, "meta.data")
keep <- as.character(md[["cluster_raw"]]) == "2"
cells <- colnames(b)[keep]
pf <- as.numeric(md[["Proteinfeatures"]][keep])
names(pf) <- cells
hi <- names(pf)[pf >= 6]
lo <- names(pf)[pf <= 1]
cat("nH", length(hi), "nL", length(lo), "\n")
mat <- tryCatch(
  GetAssayData(b, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(b, assay = "RNA", slot = "data")
)
g <- "Cd274"
x1 <- as.numeric(mat[g, hi])
x2 <- as.numeric(mat[g, lo])
wt <- wilcox.test(x1, x2)
fc <- log2(mean(expm1(x1)) + 1) - log2(mean(expm1(x2)) + 1)
row <- data.frame(
  gene = g, avg_log2FC = fc, p_val = wt$p.value,
  pct.1 = mean(x1 > 0), pct.2 = mean(x2 > 0),
  p_val_adj = wt$p.value, stringsAsFactors = FALSE
)
print(row)
write.csv(row, file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-1_NC_PFsplit/source_data/panel_g_Cd274_added.csv"),
          row.names = FALSE)
cat("wrote Cd274\n")
