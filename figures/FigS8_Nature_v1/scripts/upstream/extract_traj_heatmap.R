#!/usr/bin/env Rscript
# Extract slim payloads for native-ggplot Fig.S6 (trajectory + heatmap + spatial)
options(warn = 1)
suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
})

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS6_Nature_vector/source_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 1) Monocle trajectory slim ----
cat("=== monocle trajectory ===\n")
suppressPackageStartupMessages(library(monocle))
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
e <- new.env()
load("D:/多发性硬化项目/EAE-Ab-inject/monocle/monocle.RData", envir = e)
cds <- get("cds", envir = e)
pd <- as.data.frame(pData(cds))
cat("cells:", ncol(cds), " States:", paste(names(table(pd$State)), collapse=","), "\n")
cat("group:", paste(names(table(pd$group)), collapse=","), "\n")

# DDRTree cell coords (reducedDimS) and graph nodes (reducedDimK)
cell_xy <- t(as.matrix(monocle::reducedDimS(cds)))
colnames(cell_xy) <- c("comp1", "comp2")
node_xy <- t(as.matrix(monocle::reducedDimK(cds)))
colnames(node_xy) <- c("comp1", "comp2")
mst <- monocle::minSpanningTree(cds)
ed <- igraph::as_data_frame(mst, what = "edges")
edge_df <- data.frame(
  x = node_xy[ed$from, "comp1"], y = node_xy[ed$from, "comp2"],
  xend = node_xy[ed$to, "comp1"], yend = node_xy[ed$to, "comp2"],
  stringsAsFactors = FALSE
)
# Branch points: nodes with degree > 2 (monocle marks these)
deg <- igraph::degree(mst)
bp_names <- names(deg)[deg > 2]
bp_df <- data.frame(
  label = as.character(seq_along(bp_names)),
  comp1 = node_xy[bp_names, "comp1"],
  comp2 = node_xy[bp_names, "comp2"],
  stringsAsFactors = FALSE
)
traj <- list(
  cells = data.frame(
    cell = rownames(cell_xy),
    comp1 = cell_xy[, "comp1"],
    comp2 = cell_xy[, "comp2"],
    State = as.character(pd$State),
    group = as.character(pd$group),
    Pseudotime = as.numeric(pd$Pseudotime),
    celltype = as.character(if ("celltype" %in% names(pd)) pd$celltype else pd$seurat_clusters),
    stringsAsFactors = FALSE
  ),
  edges = edge_df,
  branch_points = bp_df
)
saveRDS(traj, file.path(OUT, "traj_slim.rds"))
cat("traj cells:", nrow(traj$cells), " edges:", nrow(traj$edges),
    " branch_points:", nrow(traj$branch_points), "\n")
rm(list = ls(envir = e), envir = e); gc(verbose = FALSE)

# ---- 2) Heatmap matrix from s6c export + gene metadata order ----
cat("=== S6C heatmap matrix ===\n")
meta <- read.csv(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig6_S6_redraw/data/s6c_export/meta.csv"),
                 stringsAsFactors = FALSE)
genes_all <- readLines(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig6_S6_redraw/data/s6c_export/genes.txt"))
cells_all <- readLines(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig6_S6_redraw/data/s6c_export/cells.txt"))
X <- Matrix::readMM(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig6_S6_redraw/data/s6c_export/exprs.mtx"))
# export wrote genes x cells; Python used .T so cells x genes — check dims
cat("MM dim:", paste(dim(X), collapse="x"),
    " n_genes:", length(genes_all), " n_cells:", length(cells_all), "\n")
if (nrow(X) == length(genes_all) && ncol(X) == length(cells_all)) {
  # genes x cells
} else if (nrow(X) == length(cells_all) && ncol(X) == length(genes_all)) {
  X <- t(X)
} else stop("exprs.mtx dim mismatch")
rownames(X) <- genes_all
colnames(X) <- cells_all
meta <- meta[match(cells_all, meta$cell_id), ]
stopifnot(all(meta$cell_id == cells_all))

# size-factor normalize + log1p
sf <- as.numeric(meta$Size_Factor)
sf[!is.finite(sf) | sf <= 0] <- 1
X <- Matrix::t(Matrix::t(X) / sf)  # divide each column by sf
X <- log1p(X)

gmeta <- read.csv(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "Fig6_S6_redraw/figS6C_results/S6C_all_gene_metadata.csv"),
                  stringsAsFactors = FALSE)
genes_use <- gmeta$gene
miss <- setdiff(genes_use, rownames(X))
if (length(miss)) cat("WARN missing genes:", paste(miss, collapse=","), "\n")
genes_use <- intersect(genes_use, rownames(X))
gmeta <- gmeta[match(genes_use, gmeta$gene), ]

# orient pseudotime: homeostatic early
pt <- as.numeric(meta$Pseudotime)
cl <- as.character(meta$cluster)
med0 <- median(pt[cl == "0"], na.rm = TRUE)
med1 <- median(pt[cl == "1"], na.rm = TRUE)
if (is.finite(med0) && is.finite(med1) && med0 > med1) {
  pt <- max(pt, na.rm = TRUE) - pt
  cat("reversed Pseudotime for homeostatic-early\n")
}
pt <- (pt - min(pt, na.rm = TRUE)) / (diff(range(pt, na.rm = TRUE)) + 1e-12)
n_bins <- 8
edges <- seq(0, 1, length.out = n_bins + 1)
bin_id <- as.integer(cut(pt, breaks = edges, include.lowest = TRUE, labels = FALSE)) - 1L

mat <- matrix(NA_real_, nrow = length(genes_use), ncol = n_bins,
              dimnames = list(genes_use, sprintf("%.2f", edges[seq_len(n_bins)])))
Xu <- as.matrix(X[genes_use, , drop = FALSE])
for (i in seq_along(genes_use)) {
  mb <- tapply(Xu[i, ], bin_id, mean, na.rm = TRUE)
  mb <- as.numeric(mb[as.character(0:(n_bins - 1))])
  mb[is.na(mb)] <- 0
  s <- sd(mb)
  mat[i, ] <- if (is.finite(s) && s > 0) (mb - mean(mb)) / s else mb - mean(mb)
}
hm <- list(
  mat = mat,
  gene_meta = gmeta,
  bin_edges = edges,
  n_bins = n_bins
)
saveRDS(hm, file.path(OUT, "heatmap_slim.rds"))
cat("heatmap genes:", nrow(mat), " bins:", ncol(mat), "\n")

cat("DONE extract monocle+heatmap\n")