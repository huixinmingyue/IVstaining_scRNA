#!/usr/bin/env Rscript
# =============================================================================
# 02  Label transfer: Fig. 2 Fo vs MZ B; hard-gated CD4 1a/1b and TNK CD8 0/3a/3b
# One FindTransferAnchors per reference; two TransferData calls each.
# Does not overwrite Fig. 4 or the volcano folder.
# =============================================================================
options(warn = 1)
options(future.globals.maxSize = 16 * 1024^3)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_spatial_anatomy_reprojection")
SRC <- file.path(OUT, "source_data")
SD  <- file.path(OUT, "source_data_plotted")
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)
dir.create(SD,  recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT, "run_log_02.txt")
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("start 02 %s", format(Sys.time()))

B_OBJ  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "figures/Bcell_reclustered.rds")
T_OBJ  <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS2_Nature_split_c1c3/source_data/CD4T_reclustered_split_c1c3.rds")
IDX    <- file.path(SRC, "spatial_index.rds")
if (!file.exists(IDX)) stop("Run 01 first (missing spatial_index.rds)")
idx <- readRDS(IDX)
OBJ <- idx$spatial_path
if (!file.exists(OBJ)) stop("spatial object missing")

assay_use <- idx$assay_use
sym <- idx$symbol_map

# ---- helpers ----
layer_or_slot <- function(obj, assay, what) {
  tryCatch(GetAssayData(obj, assay = assay, layer = what),
           error = function(e) GetAssayData(obj, assay = assay, slot = what))
}

ensure_pca <- function(obj, assay = "RNA", npcs = 30) {
  DefaultAssay(obj) <- assay
  if (!"pca" %in% Reductions(obj)) {
    say("  computing PCA on %d cells", ncol(obj))
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  } else if (!length(VariableFeatures(obj))) {
    obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  }
  obj
}

collapse_dup_symbols <- function(mat, symbols) {
  symbols <- as.character(symbols)
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  key <- tolower(symbols)
  if (anyDuplicated(key)) {
    say("  collapsing %d duplicate symbol rows", sum(duplicated(key)))
    mat <- Matrix::rowsum(mat, group = key)
    # restore a display-case symbol per group (first)
    first <- !duplicated(key)
    map_case <- setNames(symbols[first], key[first])
    rownames(mat) <- unname(map_case[rownames(mat)])
  } else {
    rownames(mat) <- symbols
  }
  mat
}

make_query <- function(sp, features_symbol, assay_use, sym) {
  counts <- layer_or_slot(sp, assay_use, "counts")
  ids <- names(sym)[tolower(sym) %in% tolower(features_symbol)]
  ids <- intersect(ids, rownames(counts))
  say("  query features mapped: %d / %d requested symbols", length(ids), length(unique(features_symbol)))
  sub <- counts[ids, , drop = FALSE]
  sub <- collapse_dup_symbols(sub, unname(sym[ids]))
  q <- CreateSeuratObject(counts = sub, min.cells = 0, min.features = 0)
  q <- NormalizeData(q, verbose = FALSE)
  q
}

run_transfer <- function(ref, query, refdata_list, dims = 1:30, tag = "x") {
  common <- intersect(rownames(ref), rownames(query))
  feats <- intersect(VariableFeatures(ref), common)
  if (length(feats) < 200) feats <- common[seq_len(min(2000L, length(common)))]
  say("%s common genes %d  features used %d", tag, length(common), length(feats))
  if (length(common) < 200) stop(tag, ": too few common genes")
  npcs <- min(30L, ncol(Embeddings(ref, "pca")))
  dims_use <- seq_len(npcs)
  say("%s FindTransferAnchors dims=1:%d ...", tag, npcs)
  t0 <- Sys.time()
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    features = feats,
    dims = dims_use,
    reference.reduction = "pca",
    normalization.method = "LogNormalize",
    verbose = TRUE
  )
  n_anc <- tryCatch(nrow(slot(anchors, "anchors")), error = function(e) NA)
  say("%s anchors in %.1f min; n_anchors=%s", tag,
      as.numeric(difftime(Sys.time(), t0, units = "mins")), n_anc)
  out <- list()
  for (nm in names(refdata_list)) {
    rd <- as.character(refdata_list[[nm]])
    names(rd) <- colnames(ref)
    rd <- unname(rd[colnames(ref)])
    if (anyNA(rd) || !length(rd)) stop(tag, " ", nm, ": NA refdata")
    say("%s TransferData %s labels: %s", tag, nm, paste(sort(unique(rd)), collapse = ", "))
    pred <- tryCatch(
      TransferData(
        anchorset = anchors,
        refdata = rd,
        dims = dims_use,
        verbose = TRUE
      ),
      error = function(e) {
        say("%s TransferData %s FAILED: %s", tag, nm, conditionMessage(e))
        NULL
      }
    )
    out[[nm]] <- pred
  }
  out
}

pred_to_df <- function(pred) {
  if (is.null(pred)) return(NULL)
  if (is.data.frame(pred)) {
    df <- pred
    if (!"barcode" %in% names(df) && !is.null(rownames(df))) {
      # keep as is; barcodes added later from rownames
    }
    say("  prediction (data.frame) cols: %s", paste(utils::head(colnames(df), 12), collapse = ", "))
    return(df)
  }
  mat <- tryCatch(
    as.matrix(GetAssayData(pred, layer = "data")),
    error = function(e) as.matrix(GetAssayData(pred, slot = "data"))
  )
  df <- as.data.frame(t(mat), stringsAsFactors = FALSE)
  say("  prediction (assay) cols: %s", paste(utils::head(colnames(df), 12), collapse = ", "))
  df
}

# ---- load spatial once ----
say("loading spatial")
sp <- readRDS(OBJ)
DefaultAssay(sp) <- assay_use

# =============================================================================
# B: Fo (cluster 0) vs MZ (cluster 1)
# =============================================================================
say("loading B reference")
bref <- readRDS(B_OBJ)
DefaultAssay(bref) <- "RNA"
stopifnot(all(c("seurat_clusters", "IV_group") %in% colnames(bref[[]])))
bref$b_subset <- dplyr::recode(
  as.character(bref$seurat_clusters),
  "0" = "Fo", "1" = "MZ", "2" = "Other",
  .default = "Other"
)
iv <- ifelse(as.character(bref$IV_group) == "IntraPos", "IV+", "IV-")
bref$b_iv <- paste(bref$b_subset, iv, sep = "_")
say("B n=%d", ncol(bref))
print(table(bref$b_subset, bref$IV_group))
print(table(bref$b_iv))
bref <- ensure_pca(bref)

# =============================================================================
# T: hard-gated true CD4 1a/1b and TNK-gated CD8 0/3a/3b
# =============================================================================
say("loading T reference and applying gates")
tobj <- readRDS(T_OBJ)
DefaultAssay(tobj) <- "RNA"
stopifnot(all(c("IV_group", "seurat_clusters") %in% colnames(tobj[[]])))
tobj$cluster_ref <- as.character(tobj$seurat_clusters)
cnt <- layer_or_slot(tobj, "RNA", "counts")
need <- c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Cd8b1", "Klrb1c", "Ms4a1", "Cd79a")
stopifnot(all(need %in% rownames(cnt)))
det <- function(g) as.integer(cnt[g, ] > 0)
T_pos <- (det("Cd3e") + det("Cd3d") + det("Cd3g")) > 0
CD4_pos <- det("Cd4") == 1
CD8_pos <- det("Cd8b1") == 1
B_strong <- det("Ms4a1") == 1 & det("Cd79a") == 1
NKT_pos <- T_pos & det("Klrb1c") == 1
tnk <- rep("unresolved", ncol(tobj))
tnk[B_strong & T_pos] <- "B-T doublet"
tnk[B_strong & !T_pos] <- "B contaminant"
rest <- tnk == "unresolved"
tnk[rest & T_pos & NKT_pos & !CD4_pos] <- "NKT-like"
rest <- tnk == "unresolved"
tnk[rest & T_pos & CD4_pos & CD8_pos] <- "CD4/CD8 DP"
rest <- tnk == "unresolved"
tnk[rest & T_pos & CD8_pos & !CD4_pos] <- "CD8 T"
rest <- tnk == "unresolved"
tnk[rest & T_pos & CD4_pos & !CD8_pos] <- "CD4 T"
tobj$tnk_class <- tnk
cd8_any <- det("Cd8a") == 1 | det("Cd8b1") == 1
klr <- det("Klrb1c") == 1
c2 <- tobj$cluster_ref == "2"
hgate <- rep("unassigned", ncol(tobj))
hgate[c2] <- "doublet_c2"
hgate[!c2 & klr] <- "NKT"
hgate[!c2 & !klr & (det("Cd4") == 1) & cd8_any] <- "DP"
hgate[!c2 & !klr & (det("Cd4") == 1) & !cd8_any] <- "true_CD4"
hgate[!c2 & !klr & (det("Cd4") == 0) & cd8_any] <- "true_CD8"
tobj$hgate <- hgate

cd4_ok <- tobj$hgate == "true_CD4" & tobj$cluster_ref %in% c("1a", "1b")
cd8_ok <- tobj$tnk_class == "CD8 T" & tobj$cluster_ref %in% c("0", "3a", "3b")
keep <- cd4_ok | cd8_ok
say("T keep %d / %d  (CD4 %d, CD8 %d)", sum(keep), ncol(tobj), sum(cd4_ok), sum(cd8_ok))

tref <- subset(tobj, cells = colnames(tobj)[keep])
tref$t_subset <- ifelse(
  tref$hgate == "true_CD4",
  paste0("CD4_", tref$cluster_ref),
  paste0("CD8_", tref$cluster_ref)
)
ivt <- ifelse(as.character(tref$IV_group) == "IntraPos", "IV+", "IV-")
tref$t_lineage <- ifelse(grepl("^CD4_", tref$t_subset), "trueCD4", "TNKcd8")
tref$t_iv <- paste(tref$t_lineage, ivt, sep = "_")
print(table(tref$t_subset, tref$IV_group))
print(table(tref$t_iv))
n_by <- as.data.frame(table(tref$t_subset, tref$IV_group), stringsAsFactors = FALSE)
write.csv(n_by, file.path(SRC, "transfer_T_reference_n.csv"), row.names = FALSE)
n_b <- as.data.frame(table(bref$b_subset, bref$IV_group), stringsAsFactors = FALSE)
write.csv(n_b, file.path(SRC, "transfer_B_reference_n.csv"), row.names = FALSE)

# drop PCA from parent object; recompute on gated cells
if ("pca" %in% Reductions(tref)) tref@reductions$pca <- NULL
VariableFeatures(tref) <- character(0)
tref <- ensure_pca(tref)

# ---- query genes = union of both references' variable features ----
feat <- unique(c(VariableFeatures(bref), VariableFeatures(tref),
                 c("Cr2", "Fcer2a", "Cd4", "Cd8b1", "Ms4a1", "Cd79a", "Jchain", "Marco", "Adgre1")))
say("building query with up to %d symbols", length(feat))
q <- make_query(sp, feat, assay_use, sym)
# stash coords
q$x <- idx$x[match(colnames(q), idx$barcodes)]
q$y <- idx$y[match(colnames(q), idx$barcodes)]
q$louvain <- idx$louvain[match(colnames(q), idx$barcodes)]
rm(sp); gc(verbose = FALSE)

# ---- module scores (always; TransferData can fail on Stereo-seq) ----
say("B Fo vs MZ markers")
Idents(bref) <- "b_subset"
mk_b <- FindMarkers(
  bref, ident.1 = "MZ", ident.2 = "Fo", assay = "RNA",
  min.pct = 0.2, logfc.threshold = 0.25, verbose = FALSE
)
mk_b$gene <- rownames(mk_b)
write.csv(mk_b, file.path(SRC, "markers_MZ_vs_Fo.csv"), row.names = FALSE)
top_mz <- intersect(mk_b$gene[mk_b$avg_log2FC > 0][seq_len(min(40, sum(mk_b$avg_log2FC > 0)))], rownames(q))
top_fo <- intersect(mk_b$gene[mk_b$avg_log2FC < 0][seq_len(min(40, sum(mk_b$avg_log2FC < 0)))], rownames(q))
say("  module genes MZ %d  Fo %d", length(top_mz), length(top_fo))
write.csv(data.frame(module = c(rep("MZ", length(top_mz)), rep("Fo", length(top_fo))),
                     gene = c(top_mz, top_fo)),
          file.path(SRC, "module_genes_MZ_Fo.csv"), row.names = FALSE)

say("T subtype markers")
Idents(tref) <- "t_subset"
mk_t <- FindAllMarkers(
  tref, assay = "RNA", only.pos = TRUE, min.pct = 0.2,
  logfc.threshold = 0.25, verbose = FALSE
)
write.csv(mk_t, file.path(SRC, "markers_T_subsets.csv"), row.names = FALSE)
t_lev <- sort(unique(as.character(tref$t_subset)))
t_mods <- lapply(t_lev, function(lv) {
  g <- mk_t$gene[mk_t$cluster == lv]
  intersect(g[seq_len(min(30, length(g)))], rownames(q))
})
names(t_mods) <- t_lev
for (lv in names(t_mods)) say("  %s module n=%d", lv, length(t_mods[[lv]]))

mod_list <- c(list(Fo = top_fo, MZ = top_mz), t_mods)
mod_list <- mod_list[vapply(mod_list, length, 1L) >= 5]
mod_names <- paste0("mod", gsub("[^A-Za-z0-9]", "", names(mod_list)))
q <- AddModuleScore(q, features = unname(mod_list), name = "tmpMod", assay = "RNA")
# AddModuleScore names tmpMod1, tmpMod2, ...
mod_cols <- grep("^tmpMod", colnames(utils::head(q[[]])), value = TRUE)
# more reliable:
mod_cols <- grep("^tmpMod", colnames(q[[]]), value = TRUE)
say("module score columns: %s", paste(mod_cols, collapse = ", "))
if (length(mod_cols) != length(mod_list)) {
  say("WARNING module column count %d != %d", length(mod_cols), length(mod_list))
}
mod_df <- q[[]][, mod_cols, drop = FALSE]
colnames(mod_df) <- paste0("mod.", names(mod_list)[seq_along(mod_cols)])
mod_df$barcode <- colnames(q)
if (all(c("mod.MZ", "mod.Fo") %in% names(mod_df))) {
  mod_df$mod_MZ_minus_Fo <- mod_df$mod.MZ - mod_df$mod.Fo
}

# ---- transfers ----
bref2 <- subset(bref, b_subset %in% c("Fo", "MZ"))
b_pred <- run_transfer(
  bref2, q,
  refdata_list = list(
    subset = as.character(bref2$b_subset),
    subset_iv = as.character(bref2$b_iv)
  ),
  tag = "B"
)
t_pred <- run_transfer(
  tref, q,
  refdata_list = list(
    subset = as.character(tref$t_subset),
    lineage_iv = as.character(tref$t_iv)
  ),
  tag = "T"
)

safe_pred_df <- function(pred, prefix) {
  if (is.null(pred)) return(NULL)
  d <- pred_to_df(pred)
  cn <- colnames(d)
  cn <- sub("^prediction\\.score\\.", paste0(prefix, "."), cn)
  cn <- sub("^prediction.score.", paste0(prefix, "."), cn, fixed = TRUE)
  cn[cn == "predicted.id"] <- paste0(prefix, ".id")
  cn[cn == "predicted.id.score"] <- paste0(prefix, ".max")
  colnames(d) <- make.names(cn, unique = TRUE)
  if (is.null(d$barcode) || !("barcode" %in% names(d))) {
    d$barcode <- rownames(d)
  }
  d
}
b_sub <- safe_pred_df(b_pred$subset, "B")
b_iv  <- safe_pred_df(b_pred$subset_iv, "Biv")
t_sub <- safe_pred_df(t_pred$subset, "T")
t_iv  <- safe_pred_df(t_pred$lineage_iv, "Tiv")

out <- data.frame(
  barcode = idx$barcodes,
  x = idx$x, y = idx$y, louvain = idx$louvain,
  stringsAsFactors = FALSE
)
join_by_bc <- function(a, b) {
  if (is.null(b)) return(a)
  m <- match(a$barcode, b$barcode)
  cbind(a, b[m, setdiff(names(b), "barcode"), drop = FALSE])
}
out <- join_by_bc(out, mod_df)
out <- join_by_bc(out, b_sub)
out <- join_by_bc(out, b_iv)
out <- join_by_bc(out, t_sub)
out <- join_by_bc(out, t_iv)

# derived deltas (make.names converts IV+ to IV.)
if (all(c("B.MZ", "B.Fo") %in% names(out))) {
  out$B_MZ_minus_Fo <- as.numeric(out$B.MZ) - as.numeric(out$B.Fo)
}
mz_hi <- grep("^Biv\\.MZ_IV", names(out), value = TRUE)
# prefer plus then minus via original: MZ_IV. is plus (make.names of MZ_IV+), MZ_IV- stays or MZ_IV.
say("Biv columns: %s", paste(grep("^Biv\\.", names(out), value = TRUE), collapse = ", "))
say("Tiv columns: %s", paste(grep("^Tiv\\.", names(out), value = TRUE), collapse = ", "))
say("T columns: %s", paste(grep("^T\\.", names(out), value = TRUE), collapse = ", "))
say("mod columns: %s", paste(grep("^mod", names(out), value = TRUE), collapse = ", "))

write.csv(out, file.path(SD, "SourceData_transfer_scores.csv"), row.names = FALSE)

# paired tests in niches
paired_stat <- function(high, low, label, n_note) {
  ok <- is.finite(high) & is.finite(low)
  high <- high[ok]; low <- low[ok]
  wt <- wilcox.test(high, low, paired = TRUE, conf.int = TRUE, exact = FALSE)
  data.frame(
    contrast = label, n = length(high), note = n_note,
    median_high = median(high), median_low = median(low),
    median_paired_delta = unname(wt$estimate),
    ci_low = wt$conf.int[1], ci_high = wt$conf.int[2],
    p_value = wt$p.value,
    frac_high_gt_low = mean(high > low),
    test = "two-sided Wilcoxon signed-rank (paired by spatial bin)",
    stringsAsFactors = FALSE
  )
}

stats <- list()
bnic <- out$louvain %in% c("1", "5")
if (all(c("B.MZ", "B.Fo") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(out$B.MZ[bnic], out$B.Fo[bnic],
                            "MZ vs Fo TransferData in Louvain 1+5", "B niche bins")
}
if (all(c("mod.MZ", "mod.Fo") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(out$mod.MZ[bnic], out$mod.Fo[bnic],
                            "MZ vs Fo module score in Louvain 1+5", "B niche bins")
}
if ("B_MZ_minus_Fo" %in% names(out) && "B_IV_delta" %in% names(out)) {
  # join B IV delta from marker table if present
}
if (all(c("Biv.MZ_IV+", "Biv.MZ_IV-") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out[["Biv.MZ_IV+"]][bnic], out[["Biv.MZ_IV-"]][bnic],
    "MZ IV+ vs MZ IV- in Louvain 1+5", "B niche bins")
}
if (all(c("Biv.Fo_IV+", "Biv.Fo_IV-") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out[["Biv.Fo_IV+"]][bnic], out[["Biv.Fo_IV-"]][bnic],
    "Fo IV+ vs Fo IV- in Louvain 1+5", "B niche bins")
}
tnic <- out$louvain %in% c("1", "4", "5", "8")
if (all(c("Tiv.trueCD4_IV+", "Tiv.trueCD4_IV-") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out[["Tiv.trueCD4_IV+"]][tnic], out[["Tiv.trueCD4_IV-"]][tnic],
    "trueCD4 IV+ vs IV- in Louvain 1+4+5+8", "B/T niche bins")
}
if (all(c("Tiv.TNKcd8_IV+", "Tiv.TNKcd8_IV-") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out[["Tiv.TNKcd8_IV+"]][tnic], out[["Tiv.TNKcd8_IV-"]][tnic],
    "TNKcd8 IV+ vs IV- in Louvain 1+4+5+8", "B/T niche bins")
}

# CD4_1a vs CD4_1b paired in T-ish bins
if (all(c("T.CD4_1a", "T.CD4_1b") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out$T.CD4_1a[tnic], out$T.CD4_1b[tnic],
    "CD4_1a vs CD4_1b in Louvain 1+4+5+8", "B/T niche bins")
}
if (all(c("T.CD8_0", "T.CD4_1b") %in% names(out))) {
  stats[[length(stats) + 1]] <- paired_stat(
    out$T.CD8_0[tnic], out$T.CD4_1b[tnic],
    "CD8_0 vs CD4_1b in Louvain 1+4+5+8", "B/T niche bins")
}

stat_df <- if (length(stats)) do.call(rbind, stats) else data.frame()
if (nrow(stat_df)) {
  write.csv(stat_df, file.path(SRC, "transfer_paired_stats.csv"), row.names = FALSE)
  print(stat_df)
} else {
  say("no paired stats (transfer columns missing)")
}

# merge radial
rad_path <- file.path(SD, "SourceData_section_markers_radial.csv")
delta_col <- if ("B_MZ_minus_Fo" %in% names(out)) {
  "B_MZ_minus_Fo"
} else if ("mod_MZ_minus_Fo" %in% names(out)) {
  "mod_MZ_minus_Fo"
} else {
  NA_character_
}
if (file.exists(rad_path) && !is.na(delta_col)) {
  rad <- read.csv(rad_path, stringsAsFactors = FALSE)
  m <- match(rad$barcode, out$barcode)
  rad$B_MZ_minus_Fo <- out[[delta_col]][m]
  rad$B.MZ <- if ("B.MZ" %in% names(out)) out$B.MZ[m] else if ("mod.MZ" %in% names(out)) out$mod.MZ[m] else NA
  rad$B.Fo <- if ("B.Fo" %in% names(out)) out$B.Fo[m] else if ("mod.Fo" %in% names(out)) out$mod.Fo[m] else NA
  write.csv(rad, rad_path, row.names = FALSE)
  rb2 <- rad %>%
    filter(is.finite(r_norm), r_norm <= 3) %>%
    mutate(r_bin = floor(r_norm * 10) / 10) %>%
    group_by(r_bin) %>%
    summarise(
      n = n(),
      B_MZ_minus_Fo_mean = mean(B_MZ_minus_Fo, na.rm = TRUE),
      B_MZ_minus_Fo_sem  = sd(B_MZ_minus_Fo, na.rm = TRUE) / sqrt(sum(is.finite(B_MZ_minus_Fo))),
      B_Fo_mean = mean(B.Fo, na.rm = TRUE),
      B_MZ_mean = mean(B.MZ, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rb2, file.path(SD, "SourceData_radial_transfer.csv"), row.names = FALSE)
}

say("done 02 %s", format(Sys.time()))
