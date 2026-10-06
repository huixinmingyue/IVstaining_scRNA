#!/usr/bin/env Rscript
# =============================================================================
# 02  ImmGen spleen T/NK overlay for hard-gated Intra DEGs
# Identity axes include NKT and CD8, not only CD4 subsets.
# Gene sets = Intra DEG (tech/ribo removed), NOT KEGG leading edge.
# =============================================================================
options(warn = 1)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newSfigures/FigS1b_reannotate_ImmGen_v1")
SRC <- file.path(OUT, "source_data")
GZ  <- file.path(SRC, "GSE109125_Genes_count_table.tsv.gz")
FC_ID <- 1.0          # bulk identity-axis fold-change
MIN_MEAN <- 10        # min mean counts on the up side

say <- function(...) cat(sprintf(...), "\n")

imm <- read.delim(gzfile(GZ), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(imm) > 2)
gene_col <- colnames(imm)[1]
rownames(imm) <- make.unique(as.character(imm[[gene_col]]))
imm[[gene_col]] <- NULL
mat <- as.matrix(imm)
storage.mode(mat) <- "numeric"

# Keep steady-state spleen T/NK; drop LPS / aCD3 activation.
keep_pat <- c(
  "^T\\.4\\.Nve\\.Sp#",
  "^T\\.4\\.19-8-TCRb\\+CD4\\+\\.Sp#",
  "^Treg\\.4\\.25hi\\.Sp#",
  "^T\\.8\\.Nve\\.Sp#",
  "^NKT\\.19-8-TCRb\\+CD1daGalCerTet\\+\\.Sp#",
  "^Tgd\\.Sp#",
  "^NK\\.27\\+11b-\\.Sp#",
  "^NK\\.27\\+11b\\+\\.Sp#",
  "^NK\\.27-11b\\+\\.Sp#"
)
keep_cols <- unique(unlist(lapply(keep_pat, function(p) grep(p, colnames(mat), value = TRUE))))
say("kept ImmGen columns: %d", length(keep_cols))
print(keep_cols)
mat <- mat[, keep_cols, drop = FALSE]

prefix_of <- function(x) sub("#\\d+$", "", x)
pref <- prefix_of(colnames(mat))
merged <- sapply(unique(pref), function(p) {
  rowMeans(mat[, pref == p, drop = FALSE])
})
# readable names
name_map <- c(
  "T.4.Nve.Sp" = "CD4_naive",
  "T.4.19-8-TCRb+CD4+.Sp" = "CD4_total",
  "Treg.4.25hi.Sp" = "Treg",
  "T.8.Nve.Sp" = "CD8_naive",
  "NKT.19-8-TCRb+CD1daGalCerTet+.Sp" = "NKT",
  "Tgd.Sp" = "Tgd",
  "NK.27+11b-.Sp" = "NK_immature",
  "NK.27+11b+.Sp" = "NK_int",
  "NK.27-11b+.Sp" = "NK_mature"
)
colnames(merged) <- ifelse(colnames(merged) %in% names(name_map),
                           name_map[colnames(merged)], colnames(merged))
col_order <- c("CD4_naive", "CD4_total", "Treg", "CD8_naive",
               "NKT", "Tgd", "NK_immature", "NK_int", "NK_mature")
col_order <- intersect(col_order, colnames(merged))
merged <- merged[, col_order, drop = FALSE]
write.csv(data.frame(immgen_col = keep_cols, prefix = prefix_of(keep_cols),
                     mapped = name_map[prefix_of(keep_cols)]),
          file.path(SRC, "ImmGen_column_map.csv"), row.names = FALSE)

# Pairwise identity axes (S1a analogue)
contrasts <- list(
  NKT_vs_CD4naive = c("NKT", "CD4_naive"),
  CD8naive_vs_CD4naive = c("CD8_naive", "CD4_naive"),
  Treg_vs_CD4naive = c("Treg", "CD4_naive"),
  Tgd_vs_CD4naive = c("Tgd", "CD4_naive"),
  NKmature_vs_CD4naive = c("NK_mature", "CD4_naive")
)
pairwise <- bind_rows(lapply(names(contrasts), function(cn) {
  a <- contrasts[[cn]][1]; b <- contrasts[[cn]][2]
  if (!all(c(a, b) %in% colnames(merged))) return(NULL)
  m1 <- merged[, a]; m2 <- merged[, b]
  lfc <- log2((m1 + 1) / (m2 + 1))
  up_mean <- ifelse(lfc >= 0, m1, m2)
  data.frame(
    contrast = cn, gene = rownames(merged),
    log2FC = lfc, mean_A = m1, mean_B = m2,
    direction = ifelse(lfc >= FC_ID & up_mean >= MIN_MEAN, paste0(a, " up"),
                       ifelse(lfc <= -FC_ID & up_mean >= MIN_MEAN, paste0(b, " up"), "NS")),
    stringsAsFactors = FALSE
  )
}))
write.csv(pairwise, file.path(SRC, "ImmGen_pairwise_allgenes.csv"), row.names = FALSE)

read_lab <- function(f) {
  p <- file.path(SRC, f)
  if (!file.exists(p)) return(NULL)
  read.csv(p, stringsAsFactors = FALSE)
}
labs <- list(
  true_CD4 = read_lab("S1b_true_CD4_sig_labelled.csv"),
  true_CD8 = read_lab("S1b_true_CD8_sig_labelled.csv"),
  NKT      = read_lab("S1b_NKT_sig_labelled.csv")
)

# Intersection tables (S1a-style)
inter_list <- list()
for (nm in names(labs)) {
  dg <- labs[[nm]]
  if (is.null(dg) || !nrow(dg)) next
  hits <- pairwise %>%
    filter(gene %in% dg$gene, direction != "NS") %>%
    left_join(dg %>% select(gene, intra_direction = direction, avg_log2FC, p_val_adj),
              by = "gene")
  hits$query_set <- nm
  inter_list[[nm]] <- hits
}
inter <- bind_rows(inter_list)
write.csv(inter, file.path(SRC, "ImmGen_pairwise_x_IntraDEG.csv"), row.names = FALSE)

# Compact: one winning ImmGen identity per gene (max |log2FC| among non-NS)
if (nrow(inter)) {
  win <- inter %>%
    group_by(query_set, gene, intra_direction) %>%
    slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
    ungroup()
  write.csv(win, file.path(SRC, "ImmGen_winning_identity.csv"), row.names = FALSE)
}

# Heatmap matrix: union of labelled Intra DEGs, log1p z-score
genes_use <- unique(unlist(lapply(labs, function(d) if (is.null(d)) NULL else d$gene)))
genes_use <- intersect(genes_use, rownames(merged))
say("heatmap genes: %d / labelled union", length(genes_use))
expr <- log1p(merged[genes_use, , drop = FALSE])
z <- t(scale(t(expr)))
z[!is.finite(z)] <- 0
# membership annotation
mem <- bind_rows(lapply(names(labs), function(nm) {
  d <- labs[[nm]]
  if (is.null(d)) return(NULL)
  d %>% transmute(gene, query_set = nm, intra_direction = direction)
}))
z_df <- as.data.frame(z)
z_df$gene <- rownames(z_df)
z_long <- z_df %>%
  left_join(mem, by = "gene") %>%
  distinct()
write.csv(z, file.path(SRC, "ImmGen_zscore_matrix.csv"))
write.csv(z_long, file.path(SRC, "ImmGen_zscore_long.csv"), row.names = FALSE)

# Peak population per gene
peak <- data.frame(
  gene = rownames(z),
  peak_subset = colnames(z)[max.col(z, ties.method = "first")],
  peak_z = apply(z, 1, max),
  stringsAsFactors = FALSE
)
peak <- peak %>% left_join(mem %>% distinct(gene, .keep_all = TRUE), by = "gene")
write.csv(peak, file.path(SRC, "ImmGen_peak_subset.csv"), row.names = FALSE)

# Peak enrichment: IntraPos vs IntraNeg genes vs NKT/CD8 vs CD4/Treg
peak_tab <- with(peak[!is.na(peak$intra_direction), ],
                 table(intra_direction, peak_subset))
write.csv(as.data.frame.matrix(peak_tab), file.path(SRC, "ImmGen_peak_crosstab.csv"))

say("DONE 02 ImmGen")
print(table(inter$contrast, inter$direction))
print(peak_tab)
