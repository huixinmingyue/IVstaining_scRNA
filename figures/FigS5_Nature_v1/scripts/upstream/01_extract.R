#!/usr/bin/env Rscript
# Extract slim data for revised Fig.S3-1 (PF>=6 vs <=1 + volcano G)
options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS3-1_NC_PFsplit/source_data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

HI <- 6
LO <- 1

# ---------- res0.2 for b/c/d ----------
rds0 <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_res0.2_PF_corr/data/CITE_Bcell_res0.2_annotated.rds")
cat("Loading", rds0, "\n")
b <- readRDS(rds0)
umap <- as.data.frame(Embeddings(b, "umap"))
colnames(umap) <- c("umap_1", "umap_2")
md <- slot(b, "meta.data")
df_b <- data.frame(
  cell = colnames(b),
  umap_1 = umap$umap_1,
  umap_2 = umap$umap_2,
  cluster = as.character(md$cluster_raw),
  Bcelltype = as.character(md$Bcelltype),
  Proteinfeatures = as.numeric(md$Proteinfeatures),
  stringsAsFactors = FALSE
)
saveRDS(df_b, file.path(OUT, "panel_b_umap.rds"))

# DotPlot genes: Bcelltype marker DotPlot + myeloid C4 markers
genes_bt <- c("Cr2", "Fcer2a", "Cd19", "Ms4a1", "Cd79a", "Vpreb3", "Igll1",
              "Cd93", "Sox4", "Cd24a", "Vpreb1", "Rag1")
genes_my <- c("Slpi", "Axl", "Mafb", "Hpgd", "Fcna")
top_genes <- c(genes_bt, genes_my)
rn <- rownames(b[["RNA"]])
top_genes <- top_genes[top_genes %in% rn]
cat("DotPlot genes kept:", paste(top_genes, collapse = ","), "\n")

cl_ord <- sort(unique(df_b$cluster))
mat <- tryCatch(
  GetAssayData(b, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(b, assay = "RNA", slot = "data")
)
dot_rows <- list()
for (g in top_genes) {
  for (cl in cl_ord) {
    cells <- colnames(b)[as.character(b$cluster_raw) == cl]
    x <- as.numeric(mat[g, cells])
    dot_rows[[length(dot_rows) + 1]] <- data.frame(
      gene = g, cluster = cl,
      avg.exp = mean(x), pct.exp = 100 * mean(x > 0),
      stringsAsFactors = FALSE
    )
  }
}
dot <- bind_rows(dot_rows) %>%
  group_by(gene) %>%
  mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) %>%
  ungroup()
dot$avg.exp.scaled <- pmax(pmin(dot$avg.exp.scaled, 2.5), -2.5)
dot$gene <- factor(dot$gene, levels = top_genes)
dot$cluster <- factor(dot$cluster, levels = cl_ord)
ann <- read.csv(
  file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_res0.2_PF_corr/data/Bcell_cluster_annotation.csv"),
  stringsAsFactors = FALSE
)
ann_map <- setNames(ann$Bcelltype, as.character(ann$cluster))
dot$Bcelltype <- unname(ann_map[as.character(dot$cluster)])
saveRDS(dot, file.path(OUT, "panel_c_dotplot.rds"))
write.csv(dot, file.path(OUT, "panel_c_dotplot.csv"), row.names = FALSE)

# logistic curves (same as before)
df_ok <- df_b[is.finite(df_b$Proteinfeatures), , drop = FALSE]
btype_levels <- ann$Bcelltype[ann$Bcelltype %in% unique(df_ok$Bcelltype)]
grid <- data.frame(
  Proteinfeatures = seq(0, as.numeric(quantile(df_ok$Proteinfeatures, 0.99)), length.out = 200)
)
pred_bt <- bind_rows(lapply(btype_levels, function(g) {
  y <- as.integer(df_ok$Bcelltype == g)
  fit <- glm(y ~ log1p(Proteinfeatures), data = df_ok, family = binomial())
  data.frame(
    Proteinfeatures = grid$Proteinfeatures,
    Bcelltype = g,
    prob = as.numeric(predict(fit, newdata = grid, type = "response")),
    stringsAsFactors = FALSE
  )
}))
saveRDS(pred_bt, file.path(OUT, "panel_d_logistic.rds"))
write.csv(pred_bt, file.path(OUT, "panel_d_logistic.csv"), row.names = FALSE)

rm(b, mat); gc(verbose = FALSE)

# ---------- noC4 for e/f/g ----------
rds1 <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_noC4_PF_PT_corr/data/CITE_Bcell_noC4_annotated.rds")
cat("Loading", rds1, "\n")
b2 <- readRDS(rds1)
umap2 <- as.data.frame(Embeddings(b2, "umap"))
colnames(umap2) <- c("umap_1", "umap_2")
md2 <- slot(b2, "meta.data")
df_e <- data.frame(
  cell = colnames(b2),
  umap_1 = umap2$umap_1,
  umap_2 = umap2$umap_2,
  cluster = as.character(md2$cluster_raw),
  Bcelltype = as.character(md2$Bcelltype),
  Proteinfeatures = as.numeric(md2$Proteinfeatures),
  stringsAsFactors = FALSE
)
df_e <- df_e[is.finite(df_e$Proteinfeatures), , drop = FALSE]
df_e$PF_class <- ifelse(df_e$Proteinfeatures >= HI, "High (PF>=6)",
                 ifelse(df_e$Proteinfeatures <= LO, "Low (PF<=1)", "Mid"))
df_e$PF_class <- factor(df_e$PF_class,
                        levels = c("Low (PF<=1)", "Mid", "High (PF>=6)"))
# split UMAP: keep High and Low only (classic split)
df_split <- df_e[df_e$PF_class %in% c("Low (PF<=1)", "High (PF>=6)"), , drop = FALSE]
df_split$PF_class <- droplevels(df_split$PF_class)
saveRDS(df_e, file.path(OUT, "panel_ef_noC4_pf.rds"))
saveRDS(df_split, file.path(OUT, "panel_e_umap_split.rds"))
write.csv(df_split, file.path(OUT, "panel_e_umap_split.csv"), row.names = FALSE)

# proportion by PF class (High vs Low) within noC4
prop <- df_split %>%
  count(PF_class, Bcelltype, name = "n") %>%
  group_by(PF_class) %>%
  mutate(pct = 100 * n / sum(n), total = sum(n)) %>%
  ungroup()
saveRDS(prop, file.path(OUT, "panel_f_proportion.rds"))
write.csv(prop, file.path(OUT, "panel_f_proportion.csv"), row.names = FALSE)
cat("Split counts:\n"); print(table(df_split$PF_class))

# Volcano G: C2 ge6 le1
deg_path <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig3_CITE_Bcell_noC4_PF_threshold_sweep/data/DEG_C2_ge6_le1.csv")
deg <- read.csv(deg_path, stringsAsFactors = FALSE)
write.csv(deg, file.path(OUT, "panel_g_DEG_C2_ge6_le1.csv"), row.names = FALSE)
cat("DEG C2 ge6 le1 genes:", nrow(deg), "\n")
lab <- c("Cd274", "Ezr", "Marcks", "Cr2", "Dtx1", "Klf2")
cat("label present:", paste(lab[lab %in% deg$gene], collapse = ","), "\n")
cat("label missing:", paste(lab[!lab %in% deg$gene], collapse = ","), "\n")

writeLines(c(
  "Fig.S3-1 revised source_data",
  paste0("Generated: ", Sys.time()),
  "E/F thresholds: High PF>=6, Low PF<=1, no C4",
  "C genes: Bcelltype DotPlot + myeloid Slpi/Axl/Mafb/Hpgd/Fcna",
  "G: DEG_C2_ge6_le1 volcano"
), file.path(OUT, "README.txt"))

cat("DONE ->", OUT, "\n")
