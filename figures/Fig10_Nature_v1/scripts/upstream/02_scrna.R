#!/usr/bin/env Rscript
# MCAO microglia identity: scDblFinder, neutrophil/BAM modules, within-cluster
# DEG, DAM/ICAM/SAMC/IFN/G2M scores, Hallmark glycolysis/OXPHOS/IFNa/G2M.
# NEW folder only.
options(warn = 1)
options(future.globals.maxSize = 8 * 1024^3)
set.seed(42)

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig9_S8S9_depth_decontam_reanalysis")
SRC <- file.path(OUT, "source_data")
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT, "run_log_02_scrna.txt")
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("start scrna %s", format(Sys.time()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(matrixTests)
  library(msigdbr)
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

get_layer <- function(obj, assay, layer) {
  x <- tryCatch(GetAssayData(obj, assay = assay, layer = layer), error = function(e) NULL)
  if (!is.null(x) && nrow(x) > 0) return(x)
  GetAssayData(obj, assay = assay, slot = ifelse(layer == "counts", "counts", "data"))
}

zmean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) < 3) {
    say("  skip set n=%d (%s)", length(genes), paste(genes, collapse = ","))
    return(rep(NA_real_, ncol(mat)))
  }
  s <- numeric(ncol(mat))
  n_ok <- 0L
  for (g in genes) {
    x <- as.numeric(mat[g, ])
    sdx <- stats::sd(x)
    if (!is.finite(sdx) || sdx == 0) next
    s <- s + (x - mean(x)) / sdx
    n_ok <- n_ok + 1L
  }
  if (n_ok < 3L) return(rep(NA_real_, ncol(mat)))
  s / n_ok
}

fast_markers <- function(mat, pos, neg, min.pct = 0.05, max_cells = 2000L) {
  pos <- intersect(pos, colnames(mat))
  neg <- intersect(neg, colnames(mat))
  if (length(pos) < 10L || length(neg) < 10L) {
    say("  skip DEG nPOS=%d nNEG=%d", length(pos), length(neg))
    return(data.frame())
  }
  X1 <- mat[, pos, drop = FALSE]
  X2 <- mat[, neg, drop = FALSE]
  pct1 <- Matrix::rowMeans(X1 > 0)
  pct2 <- Matrix::rowMeans(X2 > 0)
  keep <- which(pct1 >= min.pct | pct2 >= min.pct)
  if (!length(keep)) return(data.frame())
  genes <- rownames(mat)[keep]
  mu1 <- Matrix::rowMeans(expm1(X1[keep, , drop = FALSE]))
  mu2 <- Matrix::rowMeans(expm1(X2[keep, , drop = FALSE]))
  lfc <- log2((mu1 + 1) / (mu2 + 1))
  rm(X1, X2)
  set.seed(42)
  pos_s <- if (length(pos) > max_cells) sample(pos, max_cells) else pos
  neg_s <- if (length(neg) > max_cells) sample(neg, max_cells) else neg
  say("    Wilcoxon on %d genes, %d vs %d cells (full n %d vs %d)",
      length(genes), length(pos_s), length(neg_s), length(pos), length(neg))
  M1 <- as.matrix(mat[genes, pos_s, drop = FALSE])
  M2 <- as.matrix(mat[genes, neg_s, drop = FALSE])
  wt <- matrixTests::row_wilcoxon_twosample(M1, M2, alternative = "two.sided", exact = FALSE)
  rm(M1, M2); gc(verbose = FALSE)
  p <- wt$pvalue
  p[!is.finite(p)] <- 1
  data.frame(gene = genes, p_val = p, avg_log2FC = as.numeric(lfc),
             pct.1 = as.numeric(pct1[keep]), pct.2 = as.numeric(pct2[keep]),
             p_val_adj = p.adjust(p, method = "BH"),
             n_pos = length(pos), n_neg = length(neg),
             n_pos_wilcox = length(pos_s), n_neg_wilcox = length(neg_s),
             stringsAsFactors = FALSE)
}

NEU <- c("S100a8", "S100a9", "Csf3r", "Cxcr2")
BAM <- c("Ms4a7", "Mrc1", "Lyve1", "Ccr2", "Ly6c2")
MG_CL <- c("0", "3", "5")
LOCKED <- c("Lgals3", "Spp1", "Ifit1", "Ifit2", "Ifit3", "Isg15", "Irf7",
            "Cxcl10", "Cxcl2", "Ccl5", "P2ry12", "Tmem119", "Hexb", "Siglech",
            "Mki67", "Trem2", "Cst7", "Apoe", "Clec7a", "Itgax", "Fth1",
            "Prdx1", "S100a8", "S100a9", "Ms4a7", "Ccr2")

SETS <- list(
  DAM = c("Clec7a", "Itgax", "Cst7", "Trem2", "Apoe", "Lpl", "Tyrobp", "Axl"),
  ICAM = c("Lgals3", "Plau", "Srxn1", "Ankrd33b", "Gas2l3", "Nes", "Edn1",
           "Fam129b", "Vat1", "Fmn1"),
  SAMC = c("Spp1", "Fabp5", "Gpnmb", "Lgals3", "Lpl", "Fth1", "Cd63", "Ctsb", "Mmp12"),
  IFN = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Irf7", "Rsad2", "Oasl1", "Oasl2", "Stat1"),
  G2M = c("Mki67", "Top2a", "Birc5", "Cdk1", "Pcna", "Tyms", "Pclaf", "Ube2c")
)

say("== load MCAO.Project2.rds ==")
obj <- readRDS("D:/MCAO/MCAO.Project2.rds")
obj$contact <- ifelse(obj$orig.ident %in% c("p1", "p2"), "POS", "NEG")
obj$cluster <- as.character(obj$seurat_clusters)
obj$is_mg <- obj$cluster %in% MG_CL
say("n=%d Mg=%d", ncol(obj), sum(obj$is_mg))
print(table(obj$cluster, obj$contact))

dat <- get_layer(obj, "RNA", "data")
cts <- get_layer(obj, "RNA", "counts")

# cluster identity check
id_genes <- unique(c(NEU, BAM, "Hexb", "P2ry12", "Siglech", "Mki67", "Lgals3"))
id_genes <- intersect(id_genes, rownames(dat))
cl_mean <- sapply(sort(unique(obj$cluster)), function(cl) {
  Matrix::rowMeans(dat[id_genes, obj$cluster == cl, drop = FALSE])
})
write.csv(cl_mean, file.path(SRC, "cluster_marker_means.csv"))
say("cluster marker means written")
print(round(cl_mean, 3))

meta_path <- file.path(SRC, "mcao_cell_meta.csv")
resume <- file.exists(meta_path)
if (resume) {
  meta_try <- read.csv(meta_path, stringsAsFactors = FALSE)
  resume <- nrow(meta_try) == ncol(obj) &&
    all(meta_try$cell %in% colnames(obj)) &&
    all(c("keep", "drop", "dbl_class", "neu_high", "bam_high") %in% names(meta_try))
}
if (resume) {
  say("== resume from existing cell meta; skip scDblFinder ==")
  meta <- meta_try[match(colnames(obj), meta_try$cell), ]
  rownames(meta) <- NULL
  for (nm in c("is_mg", "neu_high", "bam_high", "drop", "keep")) {
    meta[[nm]] <- as.logical(meta[[nm]])
  }
  obj$dbl_class <- meta$dbl_class
  obj$dbl_score <- meta$dbl_score
  obj$neu_score <- meta$neu_score
  obj$bam_score <- meta$bam_score
  obj$neu_high <- meta$neu_high
  obj$bam_high <- meta$bam_high
  obj$drop <- meta$drop
  obj$keep <- meta$keep
  say("Mg drop any=%d keep=%d", sum(meta$is_mg & meta$drop), sum(meta$keep))
  rm(meta_try)
} else {
  # scDblFinder on counts, by sample
  say("== scDblFinder ==")
  sce <- SingleCellExperiment(list(counts = cts))
  colData(sce)$orig.ident <- obj$orig.ident
  sce <- scDblFinder(sce, samples = "orig.ident")
  obj$dbl_class <- as.character(sce$scDblFinder.class)
  obj$dbl_score <- as.numeric(sce$scDblFinder.score)
  say("doublets %d / %d (%.2f%%)", sum(obj$dbl_class == "doublet"), ncol(obj),
      100 * mean(obj$dbl_class == "doublet"))
  rm(sce); gc(verbose = FALSE)

  neu_use <- intersect(NEU, rownames(dat))
  bam_use <- intersect(BAM, rownames(dat))
  say("neu genes: %s", paste(neu_use, collapse = ","))
  say("bam genes: %s", paste(bam_use, collapse = ","))
  obj$neu_score <- zmean(dat, neu_use)
  obj$bam_score <- zmean(dat, bam_use)

  s100 <- intersect(c("S100a8", "S100a9"), rownames(dat))
  s100_cl <- sapply(sort(unique(obj$cluster)), function(cl) {
    mean(Matrix::rowMeans(dat[s100, obj$cluster == cl, drop = FALSE]))
  })
  neu_cl <- names(s100_cl)[which.max(s100_cl)]
  say("neutrophil-like cluster=%s (S100 mean %.3f)", neu_cl, max(s100_cl))
  non_neu <- obj$cluster != neu_cl
  p90_neu <- as.numeric(quantile(obj$neu_score[non_neu], 0.90, na.rm = TRUE))
  p90_bam_c0 <- as.numeric(quantile(obj$bam_score[obj$cluster == "0"], 0.90, na.rm = TRUE))
  obj$neu_high <- obj$is_mg & (obj$neu_score > p90_neu)
  obj$bam_high <- obj$is_mg & (obj$bam_score > p90_bam_c0)
  obj$drop <- obj$is_mg & (obj$dbl_class == "doublet" | obj$neu_high | obj$bam_high)
  obj$keep <- obj$is_mg & !obj$drop
  say("P90 neu (non-Neu cells)=%.3f; P90 bam (c0)=%.3f", p90_neu, p90_bam_c0)
  say("Mg drop doublet=%d neu_high=%d bam_high=%d any=%d keep=%d",
      sum(obj$is_mg & obj$dbl_class == "doublet"),
      sum(obj$neu_high), sum(obj$bam_high), sum(obj$drop), sum(obj$keep))

  meta <- data.frame(
    cell = colnames(obj),
    cluster = obj$cluster,
    contact = obj$contact,
    orig.ident = as.character(obj$orig.ident),
    is_mg = obj$is_mg,
    nCount_RNA = obj$nCount_RNA,
    dbl_class = obj$dbl_class,
    dbl_score = obj$dbl_score,
    neu_score = obj$neu_score,
    bam_score = obj$bam_score,
    neu_high = obj$neu_high,
    bam_high = obj$bam_high,
    drop = obj$drop,
    keep = obj$keep,
    stringsAsFactors = FALSE
  )
  write.csv(meta, file.path(SRC, "mcao_cell_meta.csv"), row.names = FALSE)
}

comp <- meta %>%
  filter(is_mg) %>%
  group_by(contact, cluster, orig.ident) %>%
  summarise(n_all = n(), n_keep = sum(keep), .groups = "drop")
write.csv(comp, file.path(SRC, "composition_per_mouse_before_after.csv"), row.names = FALSE)

arm_comp <- function(keep_only) {
  d <- meta[meta$is_mg, ]
  if (keep_only) d <- d[d$keep, ]
  d %>%
    dplyr::count(contact, cluster) %>%
    dplyr::group_by(contact) %>%
    dplyr::mutate(pct = 100 * n / sum(n), total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cleaned = keep_only)
}
arm <- bind_rows(arm_comp(FALSE), arm_comp(TRUE))
write.csv(arm, file.path(SRC, "composition_arm_before_after.csv"), row.names = FALSE)
print(as.data.frame(arm))

enr <- meta %>%
  filter(is_mg) %>%
  group_by(cluster, contact) %>%
  summarise(
    n = n(),
    pct_dbl = 100 * mean(dbl_class == "doublet"),
    pct_neu = 100 * mean(neu_high),
    pct_bam = 100 * mean(bam_high),
    pct_drop = 100 * mean(drop),
    n_keep = sum(keep),
    .groups = "drop"
  )
write.csv(enr, file.path(SRC, "contamination_by_cluster_contact.csv"), row.names = FALSE)
print(as.data.frame(enr))

# ---- DEG: cluster 5 vs 0, and within-cluster POS vs NEG ----
say("== DEG ==")
run_pair <- function(pos, neg, label, cleaned) {
  say("  %s cleaned=%s nPOS=%d nNEG=%d", label, cleaned, length(pos), length(neg))
  fm <- fast_markers(dat, pos, neg, min.pct = 0.01)
  if (!nrow(fm)) return(fm)
  fm$contrast <- label
  fm$cleaned <- cleaned
  fm$sig <- fm$p_val_adj < 0.05 & abs(fm$avg_log2FC) >= 0.5
  fm
}
cells_mg <- colnames(obj)[obj$is_mg]
deg <- list()
for (cleaned in c(FALSE, TRUE)) {
  ok <- if (cleaned) meta$keep else meta$is_mg
  id <- meta$cell[ok]
  c5 <- id[meta$cluster[match(id, meta$cell)] == "5"]
  c0 <- id[meta$cluster[match(id, meta$cell)] == "0"]
  deg[[length(deg) + 1]] <- run_pair(c5, c0, "cluster5_vs_cluster0", cleaned)
  for (cl in MG_CL) {
    in_cl <- id[meta$cluster[match(id, meta$cell)] == cl]
    pos <- in_cl[meta$contact[match(in_cl, meta$cell)] == "POS"]
    neg <- in_cl[meta$contact[match(in_cl, meta$cell)] == "NEG"]
    deg[[length(deg) + 1]] <- run_pair(pos, neg, paste0("within_c", cl, "_POS_vs_NEG"), cleaned)
  }
}
deg_df <- bind_rows(deg)
write.csv(deg_df, file.path(SRC, "DEG_all_contrasts.csv"), row.names = FALSE)

locked_df <- deg_df %>% filter(gene %in% LOCKED)
write.csv(locked_df, file.path(SRC, "DEG_locked_genes.csv"), row.names = FALSE)
print(as.data.frame(
  locked_df %>%
    filter(gene %in% c("Lgals3", "Ifit1", "Cxcl10", "P2ry12", "S100a8", "Ms4a7")) %>%
    select(contrast, cleaned, gene, avg_log2FC, p_val_adj, n_pos, n_neg, sig)
))

# ---- custom + Hallmark module scores on Mg ----
say("== gene-set scores ==")
md <- NULL
gmt_path <- file.path(SRC, "mh.all.v2024.1.Mm.symbols.gmt")
pick_h <- function(name) {
  if (file.exists(gmt_path)) {
    ln <- grep(paste0("^", name, "\\b"), readLines(gmt_path, warn = FALSE), value = TRUE)
    if (length(ln)) {
      p <- strsplit(ln[1], "\\s+")[[1]]
      return(setdiff(p[-(1:2)], c("", "NA")))
    }
  }
  if (is.null(md) || !"gs_name" %in% names(md)) return(character())
  unique(md$gene_symbol[md$gs_name == name])
}
SETS$H_glyc <- pick_h("HALLMARK_GLYCOLYSIS")
SETS$H_oxph <- pick_h("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
SETS$H_ifna <- pick_h("HALLMARK_INTERFERON_ALPHA_RESPONSE")
SETS$H_g2m  <- pick_h("HALLMARK_G2M_CHECKPOINT")
set_used <- lapply(SETS, function(g) intersect(g, rownames(dat)))
saveRDS(set_used, file.path(SRC, "genesets_used.rds"))
for (nm in names(set_used)) say("  %s %d/%d", nm, length(set_used[[nm]]), length(SETS[[nm]]))

mg_cells <- colnames(obj)[obj$is_mg]
mdat <- dat[, mg_cells, drop = FALSE]
score_tab <- meta[meta$is_mg, ]
for (nm in names(set_used)) {
  sc <- zmean(mdat, set_used[[nm]])
  score_tab[[nm]] <- sc
}
score_tab$H_glyc_minus_oxph <- score_tab$H_glyc - score_tab$H_oxph
write.csv(score_tab, file.path(SRC, "mg_module_scores.csv"), row.names = FALSE)

score_long <- tidyr::pivot_longer(
  score_tab,
  cols = c(names(SETS), "H_glyc_minus_oxph"),
  names_to = "set", values_to = "score"
)
tests <- score_long %>%
  filter(is.finite(score)) %>%
  group_by(set) %>%
  summarise(
    n_pos = sum(contact == "POS"),
    n_neg = sum(contact == "NEG"),
    mean_pos = mean(score[contact == "POS"]),
    mean_neg = mean(score[contact == "NEG"]),
    delta = mean_pos - mean_neg,
    p_all = wilcox.test(score[contact == "POS"], score[contact == "NEG"], exact = FALSE)$p.value,
    .groups = "drop"
  )
# cleaned POS vs NEG
tests2 <- score_long %>%
  filter(keep, is.finite(score)) %>%
  group_by(set) %>%
  summarise(
    mean_pos_clean = mean(score[contact == "POS"]),
    mean_neg_clean = mean(score[contact == "NEG"]),
    delta_clean = mean_pos_clean - mean_neg_clean,
    p_clean = wilcox.test(score[contact == "POS"], score[contact == "NEG"], exact = FALSE)$p.value,
    n_pos_clean = sum(contact == "POS"),
    n_neg_clean = sum(contact == "NEG"),
    .groups = "drop"
  )
tests <- left_join(tests, tests2, by = "set") %>%
  mutate(q_all = p.adjust(p_all, "BH"), q_clean = p.adjust(p_clean, "BH"))
write.csv(tests, file.path(SRC, "module_POS_vs_NEG_tests.csv"), row.names = FALSE)
print(as.data.frame(tests))

clx <- score_long %>%
  filter(is.finite(score), cluster %in% MG_CL) %>%
  group_by(set, cluster, contact) %>%
  summarise(n = n(), mean = mean(score), median = median(score), .groups = "drop")
write.csv(clx, file.path(SRC, "module_cluster_by_contact.csv"), row.names = FALSE)

# within-cluster contact tests for hero sets
hero <- c("DAM", "ICAM", "SAMC", "IFN", "G2M", "H_glyc", "H_oxph", "H_ifna", "H_g2m",
          "H_glyc_minus_oxph")
wcl <- score_long %>%
  filter(set %in% hero, cluster %in% MG_CL, is.finite(score)) %>%
  group_by(set, cluster) %>%
  summarise(
    n_pos = sum(contact == "POS"), n_neg = sum(contact == "NEG"),
    mean_pos = mean(score[contact == "POS"]),
    mean_neg = mean(score[contact == "NEG"]),
    delta = mean_pos - mean_neg,
    p = wilcox.test(score[contact == "POS"], score[contact == "NEG"], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(set) %>%
  mutate(q = p.adjust(p, "BH")) %>%
  ungroup()
write.csv(wcl, file.path(SRC, "module_within_cluster_contact.csv"), row.names = FALSE)

say("DONE scrna %s", format(Sys.time()))
