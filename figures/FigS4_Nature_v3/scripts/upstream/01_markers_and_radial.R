#!/usr/bin/env Rscript
# =============================================================================
# 01  Marker overlay tables + follicle-centred radial profiles
# Writes only into Fig4_spatial_anatomy_reprojection/. Does not touch Fig. 4.
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
dir.create(file.path(OUT, "scripts"), showWarnings = FALSE)

LOG <- file.path(OUT, "run_log_01.txt")
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("start 01 %s", format(Sys.time()))

UM_PER_UNIT <- 0.5
STEP <- 50
MIN_FOLL_BINS <- 40L

cands <- unique(c(
  Sys.glob("D:/*/spleen_spatial/SP_stereoSeq/sp_seurat.rds"),
  "D:/多发性硬化项目/spleen_spatial/SP_stereoSeq/sp_seurat.rds"
))
cands <- cands[file.exists(cands)]
if (!length(cands)) stop("No Stereo-seq object found")
OBJ <- cands[which.max(file.info(cands)$size)]
say("spatial object %.1f MB", file.info(OBJ)$size / 1e6)

sp <- readRDS(OBJ)
say("n bins %d  n genes %d  assays %s", ncol(sp), nrow(sp),
    paste(Assays(sp), collapse = ", "))
assay_use <- if ("SCT" %in% Assays(sp)) "SCT" else if ("RNA" %in% Assays(sp)) "RNA" else Assays(sp)[1]
DefaultAssay(sp) <- assay_use
say("assay %s", assay_use)

mf <- tryCatch(sp[[assay_use]]@meta.features, error = function(e) NULL)
if (!is.null(mf) && "real_gene_name" %in% colnames(mf)) {
  sym <- as.character(mf$real_gene_name)
  names(sym) <- rownames(sp)
  empty <- is.na(sym) | !nzchar(sym)
  sym[empty] <- names(sym)[empty]
} else {
  sym <- setNames(rownames(sp), rownames(sp))
}

map_ids <- function(wanted) {
  hit <- names(sym)[tolower(sym) %in% tolower(wanted)]
  unmatched <- wanted[!tolower(wanted) %in% tolower(sym)]
  if (length(unmatched)) say("  unmatched among %s: %s",
                             paste(wanted, collapse = "/"),
                             paste(unmatched, collapse = ", "))
  unique(hit)
}

# Display name -> candidate symbols (mouse)
gene_alias <- list(
  Cr2      = c("Cr2"),
  Fcer2a   = c("Fcer2a"),
  Cd4      = c("Cd4"),
  Jchain   = c("Jchain", "Igj", "JCHAIN"),
  Marco    = c("Marco"),
  Adgre1   = c("Adgre1", "Emr1"),
  Ms4a1    = c("Ms4a1"),
  Cd79a    = c("Cd79a"),
  Cd3e     = c("Cd3e"),
  Cd8b1    = c("Cd8b1")
)

mapped <- lapply(gene_alias, map_ids)
for (nm in names(mapped)) {
  say("map %s -> %s", nm, if (length(mapped[[nm]])) paste(unname(sym[mapped[[nm]]]), collapse = ",") else "NONE")
}
needed_disp <- c("Cr2", "Fcer2a", "Cd4", "Jchain", "Marco", "Adgre1")
miss <- needed_disp[vapply(needed_disp, function(g) !length(mapped[[g]]), logical(1))]
if (length(miss)) stop("Required genes not mapped: ", paste(miss, collapse = ", "))

get_mat <- function(layer) {
  tryCatch(
    GetAssayData(sp, assay = assay_use, layer = layer),
    error = function(e) GetAssayData(sp, assay = assay_use, slot = layer)
  )
}
dat <- get_mat("data")

coords <- NULL
if (length(sp@images)) {
  coords <- tryCatch(
    as.data.frame(GetTissueCoordinates(sp)),
    error = function(e) {
      img <- names(sp@images)[1]
      as.data.frame(GetTissueCoordinates(sp, image = img))
    }
  )
}
md0 <- sp[[]]
if (is.null(coords) && all(c("x", "y") %in% colnames(md0))) {
  say("coords from metadata x/y")
  coords <- data.frame(x = as.numeric(md0$x), y = as.numeric(md0$y),
                       barcode = colnames(sp), stringsAsFactors = FALSE)
  rownames(coords) <- colnames(sp)
} else if (is.null(coords) && "spatial" %in% Reductions(sp)) {
  say("coords from reduction spatial")
  em <- as.data.frame(Embeddings(sp, "spatial"))
  coords <- data.frame(x = em[[1]], y = em[[2]], barcode = rownames(em),
                       stringsAsFactors = FALSE)
} else if (is.null(coords)) {
  stop("No spatial coordinates (no images, no x/y metadata, no spatial reduction)")
}
cn <- tolower(colnames(coords))
colnames(coords) <- cn
if (!all(c("x", "y") %in% colnames(coords))) {
  coords$x <- coords[[1]]
  coords$y <- coords[[2]]
}
coords$barcode <- rownames(coords)
# align to object cells
coords <- coords[match(colnames(sp), coords$barcode), , drop = FALSE]
if (anyNA(coords$x)) {
  coords$barcode <- colnames(sp)
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)
}

md <- sp[[]]
louv <- if ("louvain" %in% colnames(md)) as.character(md$louvain) else {
  clcol <- grep("louvain|seurat_clusters|ident", colnames(md), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(clcol)) stop("No louvain column")
  say("using cluster column %s", clcol)
  as.character(md[[clcol]])
}

pull_gene <- function(ids) {
  if (!length(ids)) return(rep(NA_real_, ncol(sp)))
  ids <- intersect(ids, rownames(dat))
  if (!length(ids)) return(rep(NA_real_, ncol(sp)))
  if (length(ids) == 1L) as.numeric(dat[ids, ]) else as.numeric(Matrix::colMeans(dat[ids, , drop = FALSE]))
}

df <- data.frame(
  barcode = colnames(sp),
  x = as.numeric(coords$x),
  y = as.numeric(coords$y),
  louvain = louv,
  stringsAsFactors = FALSE
)
for (nm in names(gene_alias)) df[[nm]] <- pull_gene(mapped[[nm]])

# Join Fig. 4 B IV scores (paired delta) without modifying that folder
fig4_b <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_Nature_v2/source_data_plotted/SourceData_Fig4bc_Bcell_bins.csv")
fig4_t <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_Nature_v2/source_data_plotted/SourceData_Fig4hi_CD4T_bins.csv")
fig4_p <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_Nature_v2/source_data_plotted/SourceData_Fig4ef_Plasma_bins.csv")
if (file.exists(fig4_b)) {
  bsc <- read.csv(fig4_b, stringsAsFactors = FALSE)
  m <- match(df$barcode, as.character(bsc$barcode))
  df$B_IV_high <- bsc$Bcell.spcd452.high[m]
  df$B_IV_low  <- bsc$Bcell.spcd452.low[m]
  df$B_IV_delta <- df$B_IV_high - df$B_IV_low
  say("joined Fig4 B scores: %d / %d bins", sum(is.finite(df$B_IV_delta)), nrow(df))
}
if (file.exists(fig4_t)) {
  tsc <- read.csv(fig4_t, stringsAsFactors = FALSE)
  m <- match(df$barcode, as.character(tsc$barcode))
  df$T_IV_high <- tsc$CD4T.spcd45.high[m]
  df$T_IV_low  <- tsc$CD4T.spcd45.low[m]
  df$T_IV_delta <- df$T_IV_high - df$T_IV_low
}
if (file.exists(fig4_p)) {
  psc <- read.csv(fig4_p, stringsAsFactors = FALSE)
  m <- match(df$barcode, as.character(psc$barcode))
  df$P_IV_high <- psc$Plasma.spcd452.high[m]
  df$P_IV_low  <- psc$Plasma.spcd452.low[m]
  df$P_IV_delta <- df$P_IV_high - df$P_IV_low
}

# ---- follicles: connected components of Louvain 5 ----
find_components <- function(x, y, step = STEP) {
  n <- length(x)
  key <- paste(x, y, sep = ":")
  imap <- setNames(seq_len(n), key)
  seen <- logical(n)
  comps <- integer(n)
  cid <- 0L
  for (i in seq_len(n)) {
    if (seen[i]) next
    cid <- cid + 1L
    stack <- i
    while (length(stack)) {
      j <- stack[length(stack)]
      stack <- stack[-length(stack)]
      if (seen[j]) next
      seen[j] <- TRUE
      comps[j] <- cid
      xx <- x[j]; yy <- y[j]
      cand <- paste(c(xx - step, xx + step, xx, xx),
                    c(yy, yy, yy - step, yy + step), sep = ":")
      nb <- unname(imap[cand])
      nb <- nb[!is.na(nb) & !seen[nb]]
      if (length(nb)) stack <- c(stack, nb)
    }
  }
  comps
}

f_idx <- which(df$louvain == "5" & is.finite(df$x) & is.finite(df$y))
say("cluster 5 bins: %d", length(f_idx))
comp <- find_components(df$x[f_idx], df$y[f_idx])
tab <- as.data.frame(table(comp), stringsAsFactors = FALSE)
colnames(tab) <- c("comp", "n")
tab$comp <- as.integer(as.character(tab$comp))
keep_c <- tab$comp[tab$n >= MIN_FOLL_BINS]
say("follicles with n>=%d: %d / %d components", MIN_FOLL_BINS, length(keep_c), nrow(tab))
print(tab[order(-tab$n), , drop = FALSE])

df$follicle_id <- NA_integer_
df$follicle_id[f_idx] <- ifelse(comp %in% keep_c, comp, NA_integer_)

foll <- df %>%
  filter(!is.na(follicle_id)) %>%
  group_by(follicle_id) %>%
  summarise(
    n = n(),
    cx = mean(x), cy = mean(y),
    # equivalent radius of n squares of side STEP
    r_eq_units = sqrt(n * STEP * STEP / pi),
    r_p90_units = as.numeric(quantile(sqrt((x - mean(x))^2 + (y - mean(y))^2), 0.90)),
    .groups = "drop"
  )
foll$radius_units <- pmax(foll$r_eq_units, 1)
foll$radius_um <- foll$radius_units * UM_PER_UNIT
write.csv(foll, file.path(SRC, "follicle_centroids.csv"), row.names = FALSE)
say("using %d follicles; radius um min=%.0f median=%.0f max=%.0f",
    nrow(foll), min(foll$radius_um), median(foll$radius_um), max(foll$radius_um))

# nearest follicle for every bin
cx <- foll$cx; cy <- foll$cy; rid <- foll$follicle_id; r0 <- foll$radius_units
nn <- nrow(df)
nearest <- integer(nn)
dmin <- numeric(nn)
chunk <- 4000L
for (s in seq(1L, nn, by = chunk)) {
  e <- min(s + chunk - 1L, nn)
  dx <- outer(df$x[s:e], cx, "-")
  dy <- outer(df$y[s:e], cy, "-")
  d2 <- dx * dx + dy * dy
  j <- max.col(-d2, ties.method = "first")
  nearest[s:e] <- rid[j]
  dmin[s:e] <- sqrt(d2[cbind(seq_len(e - s + 1L), j)])
}
df$nearest_follicle <- nearest
df$dist_centroid_units <- dmin
df$dist_centroid_um <- dmin * UM_PER_UNIT
df$radius_ref_units <- r0[match(df$nearest_follicle, rid)]
df$r_norm <- df$dist_centroid_units / df$radius_ref_units

# signed distance to cluster-5 edge: 0 on the B-niche boundary
# approx: dist to nearest cluster-5 bin; negative if inside cluster 5
c5 <- df[df$louvain == "5", c("x", "y"), drop = FALSE]
# subsample grid already 25 um; compute to all c5 via chunked nearest
c5x <- c5$x; c5y <- c5$y
d5 <- numeric(nn)
for (s in seq(1L, nn, by = chunk)) {
  e <- min(s + chunk - 1L, nn)
  dx <- outer(df$x[s:e], c5x, "-")
  dy <- outer(df$y[s:e], c5y, "-")
  d2 <- dx * dx + dy * dy
  j <- max.col(-d2, ties.method = "first")
  d5[s:e] <- sqrt(d2[cbind(seq_len(e - s + 1L), j)])
}
df$dist_c5_units <- d5
df$dist_c5_um <- d5 * UM_PER_UNIT
df$signed_c5_um <- ifelse(df$louvain == "5", -df$dist_centroid_um * 0, df$dist_c5_um)
# inside cluster 5, signed distance is negative distance to the component rim:
# approximate rim as radius_ref - dist_centroid for members
inside <- df$louvain == "5"
df$signed_c5_um[inside] <- (df$dist_centroid_units[inside] - df$radius_ref_units[inside]) * UM_PER_UNIT

write.csv(
  data.frame(
    display = names(gene_alias),
    aliases = vapply(gene_alias, paste, "", collapse = "/"),
    n_features = vapply(mapped, length, 1L),
    ensembl = vapply(mapped, function(x) paste(x, collapse = ";"), ""),
    symbol = vapply(mapped, function(x) paste(unique(unname(sym[x])), collapse = ";"), "")
  ),
  file.path(SRC, "gene_id_map.csv"), row.names = FALSE
)

# slim plotted table (all bins)
keep_cols <- c("barcode", "x", "y", "louvain",
               names(gene_alias),
               "B_IV_high", "B_IV_low", "B_IV_delta",
               "T_IV_high", "T_IV_low", "T_IV_delta",
               "P_IV_high", "P_IV_low", "P_IV_delta",
               "follicle_id", "nearest_follicle",
               "dist_centroid_um", "r_norm", "dist_c5_um", "signed_c5_um")
keep_cols <- intersect(keep_cols, names(df))
write.csv(df[, keep_cols, drop = FALSE],
          file.path(SD, "SourceData_section_markers_radial.csv"),
          row.names = FALSE)

# radial bins in r_norm
rb <- df %>%
  filter(is.finite(r_norm), r_norm <= 3) %>%
  mutate(r_bin = floor(r_norm * 10) / 10) %>%  # 0.1 radius units
  group_by(r_bin) %>%
  summarise(
    n = n(),
    across(any_of(c(names(gene_alias), "B_IV_delta", "T_IV_delta", "P_IV_delta")),
           list(mean = ~mean(.x, na.rm = TRUE),
                sem  = ~{ v <- .x[is.finite(.x)]; if (!length(v)) NA_real_ else stats::sd(v) / sqrt(length(v)) })),
    .groups = "drop"
  )
write.csv(rb, file.path(SD, "SourceData_radial_rnorm.csv"), row.names = FALSE)

# Spearman: does B IV delta rise toward r_norm ~ 1?
if ("B_IV_delta" %in% names(df)) {
  ok <- is.finite(df$r_norm) & is.finite(df$B_IV_delta) & df$r_norm <= 2.5
  if (sum(ok) > 50) {
    ct <- cor.test(df$r_norm[ok], df$B_IV_delta[ok], method = "spearman", exact = FALSE)
    # peak: which r_bin has max mean B_IV_delta among bins with n>=50
    pk <- rb %>% filter(n >= 50, is.finite(B_IV_delta_mean))
    peak_row <- if (nrow(pk)) pk[which.max(pk$B_IV_delta_mean), , drop = FALSE] else NULL
    peak_cr2 <- rb %>% filter(n >= 50, is.finite(Cr2_mean))
    peak_cr2 <- if (nrow(peak_cr2)) peak_cr2[which.max(peak_cr2$Cr2_mean), , drop = FALSE] else NULL
    stats <- data.frame(
      test = "Spearman B_IV_delta vs r_norm (r_norm<=2.5)",
      n = sum(ok),
      rho = unname(ct$estimate),
      p_value = ct$p.value,
      peak_B_IV_rbin = if (!is.null(peak_row)) peak_row$r_bin else NA,
      peak_B_IV_mean = if (!is.null(peak_row)) peak_row$B_IV_delta_mean else NA,
      peak_Cr2_rbin = if (!is.null(peak_cr2)) peak_cr2$r_bin else NA,
      peak_Cr2_mean = if (!is.null(peak_cr2)) peak_cr2$Cr2_mean else NA,
      n_follicles = nrow(foll),
      min_follicle_bins = MIN_FOLL_BINS,
      um_per_unit = UM_PER_UNIT,
      stringsAsFactors = FALSE
    )
    write.csv(stats, file.path(SRC, "radial_BIV_stats.csv"), row.names = FALSE)
    print(stats)
  }
}

# save a lightweight query-ready feature map for script 02
saveRDS(list(
  barcodes = df$barcode,
  x = df$x, y = df$y, louvain = df$louvain,
  symbol_map = sym,
  assay_use = assay_use,
  spatial_path = OBJ,
  um_per_unit = UM_PER_UNIT,
  step = STEP
), file.path(SRC, "spatial_index.rds"))

say("wrote marker/radial tables")
say("done 01 %s", format(Sys.time()))
