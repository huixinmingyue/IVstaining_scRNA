#!/usr/bin/env Rscript
# =============================================================================
# Build plotted source tables for Fig. 1 v3 and Fig. S1a
# New folders only — does not overwrite Fig1_Nature_v2 or FigS1_Nature_b_allDEG
#
# CD4 T / CD8 T identities come from the RNA-rule classes in
# Figure_TNK_reannotation panel b (cell_audit_table.csv). Marker statistics
# are recomputed on those display identities from the extra+intra object.
# =============================================================================

options(warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
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

FIG1 <- file.path(IV_FIG_ROOT, "Fig1_Nature_v3")
FIGS1 <- file.path(IV_FIG_ROOT, "FigS1_Nature_v2")
SRC1 <- file.path(FIG1, "source_data")
SRC1S <- file.path(FIGS1, "source_data")
dir.create(SRC1S, showWarnings = FALSE, recursive = TRUE)

CT_DISP <- c(
  "CD4 T", "CD8 T", "NKT-like", "CD4/CD8 DP", "T other",
  "B-T doublet", "B cell", "cDC2", "NK", "Plasma cell", "Neutrophil"
)
IV_LAB <- c(IntraNeg = "IV\u2212 (extravascular)",
            IntraPos = "IV+ (intravascular)")

genes <- c(
  "Il7r", "Cd4", "Cd8b1", "Cd8a", "Cd27", "Klrb1c", "Cd3e",
  "Cd19", "Ms4a1", "Cd79a", "Klf4", "Clec9a", "Cd68", "Itgax",
  "Ncr1", "Cd160", "Ccl3", "Gzmb", "Nkg7", "Jchain", "Mzb1",
  "S100a8", "S100a9"
)

map_display <- function(celltype, t_class) {
  out <- as.character(celltype)
  out[celltype == "Bcell"] <- "B cell"
  out[celltype == "cDC2"] <- "cDC2"
  out[celltype == "NK"] <- "NK"
  out[celltype == "Neutrophils"] <- "Neutrophil"
  out[celltype == "Plasma"] <- "Plasma cell"
  is_t <- celltype %in% c("CD4_T", "CD8_T")
  tcl <- as.character(t_class)
  out[is_t] <- tcl[is_t]
  out[out == "B contaminant"] <- "B cell"
  out[out == "NK-like"] <- "NK"
  out[out == "unresolved"] <- "T other"
  out
}

cat("Reading TNK panel-b audit table ...\n")
audit <- read.csv(file.path(SRC1, "cell_audit_table.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
n_in <- nrow(audit)
stopifnot(n_in == 15271)

audit$lineage <- map_display(audit$celltype, audit$t_class)
bad <- setdiff(unique(audit$lineage), CT_DISP)
if (length(bad)) stop("Unexpected display class: ", paste(bad, collapse = ", "))
audit$lineage <- factor(audit$lineage, levels = CT_DISP)
audit$IV_group <- factor(audit$IV_group, levels = c("IntraNeg", "IntraPos"))
audit$IV_lab <- factor(unname(IV_LAB[as.character(audit$IV_group)]),
                       levels = unname(IV_LAB))

# Coordinate provenance: prefer the published Fig. 1 UMAP object
umap <- readRDS(file.path(SRC1, "Fig1D_umap.rds"))
if ("cell" %in% names(umap)) {
  umap$barcode <- as.character(umap$cell)
} else {
  umap$barcode <- rownames(umap)
}
m <- match(audit$barcode, umap$barcode)
n_umap <- sum(!is.na(m))
cat("UMAP matched:", n_umap, "/", n_in, "\n")
stopifnot(n_umap == n_in)
audit$umap_1 <- umap$umap_1[m]
audit$umap_2 <- umap$umap_2[m]

# Lumping ledger
lump <- audit %>%
  mutate(
    t_class_raw = ifelse(celltype %in% c("CD4_T", "CD8_T"),
                         as.character(t_class), NA_character_),
    rule = case_when(
      t_class_raw == "B contaminant" ~ "B contaminant -> B cell",
      t_class_raw == "NK-like" ~ "NK-like -> NK",
      t_class_raw == "unresolved" ~ "unresolved -> T other",
      TRUE ~ "identity kept"
    )
  ) %>%
  count(celltype, t_class_raw, lineage, rule, name = "n") %>%
  arrange(desc(n))
write.csv(lump, file.path(SRC1, "lineage_crosstab.csv"), row.names = FALSE)
print(lump)

cat("\nDisplay class x IV:\n")
tab_iv <- with(audit, table(lineage, IV_group))
print(tab_iv)
write.csv(as.data.frame.matrix(tab_iv),
          file.path(SRC1, "lineage_by_IV_counts.csv"))

ids <- audit %>%
  transmute(
    barcode,
    original_celltype = celltype,
    t_class = ifelse(celltype %in% c("CD4_T", "CD8_T"),
                     as.character(t_class), NA_character_),
    lineage, IV_group, IV_lab, umap_1, umap_2
  )
write.csv(ids, file.path(SRC1, "cell_identities.csv"), row.names = FALSE)
saveRDS(ids, file.path(SRC1, "Fig1D_umap_reannotated.rds"))
file.copy(file.path(SRC1, "cell_identities.csv"),
          file.path(SRC1S, "cell_identities.csv"), overwrite = TRUE)

prop <- ids %>%
  count(IV_group, IV_lab, lineage, name = "Freq") %>%
  group_by(IV_group) %>%
  mutate(total = sum(Freq), Percentage = 100 * Freq / total) %>%
  ungroup()
stopifnot(all(prop$total %in% c(7519, 7752)))
write.csv(prop, file.path(SRC1, "Fig1E_proportions.csv"), row.names = FALSE)

ct_order <- prop %>%
  group_by(lineage) %>%
  summarise(mean_pct = mean(Percentage), .groups = "drop") %>%
  arrange(desc(mean_pct))
write.csv(ct_order, file.path(SRC1, "lineage_order.csv"), row.names = FALSE)
file.copy(file.path(SRC1, "lineage_order.csv"),
          file.path(SRC1S, "lineage_order.csv"), overwrite = TRUE)
cat("\nMean abundance order:\n")
print(ct_order)

# -----------------------------------------------------------------------------
# Marker matrix on the display identities (Seurat DotPlot arithmetic)
# -----------------------------------------------------------------------------
cat("\nLoading extra+intra object for marker genes ...\n")
obj <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "sp.group.rds"))
DefaultAssay(obj) <- "RNA"
missing_bc <- setdiff(ids$barcode, colnames(obj))
if (length(missing_bc)) stop("Barcodes missing from object: ", length(missing_bc))
sp <- obj[, ids$barcode]
rm(obj); gc()
if (inherits(sp[["RNA"]], "Assay5")) {
  try(sp[["RNA"]] <- JoinLayers(sp[["RNA"]]), silent = TRUE)
}

get_mat <- function(layer) {
  tryCatch(
    GetAssayData(sp, assay = "RNA", layer = layer),
    error = function(e) GetAssayData(sp, assay = "RNA")
  )
}
genes <- intersect(genes, rownames(sp))
stopifnot(length(genes) >= 20)
cat("Genes used:", paste(genes, collapse = ", "), "\n")

dat <- get_mat("data")
cnt <- get_mat("counts")
dat <- dat[genes, ids$barcode, drop = FALSE]
cnt <- cnt[genes, ids$barcode, drop = FALSE]
lineage <- as.character(ids$lineage)
stopifnot(identical(colnames(dat), ids$barcode))

avg_mat <- sapply(genes, function(g) {
  tapply(expm1(as.numeric(dat[g, ])), lineage, mean)
})
pct_mat <- sapply(genes, function(g) {
  tapply(as.numeric(dat[g, ]) > 0, lineage, function(z) 100 * mean(z))
})
# scale() is per-gene across display classes
scale_safe <- function(x) {
  if (length(x) < 2 || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}
scaled_mat <- apply(avg_mat, 2, scale_safe)

dot <- as.data.frame(as.table(avg_mat), stringsAsFactors = FALSE)
names(dot) <- c("celltype", "gene", "avg")
dot$pct <- as.numeric(as.table(pct_mat))
dot$avg_scaled <- as.numeric(as.table(scaled_mat))
dot$celltype <- as.character(dot$celltype)
dot$gene <- as.character(dot$gene)
dot <- dot[dot$celltype %in% CT_DISP, ]
write.csv(dot, file.path(SRC1S, "FigS1a_markers_dotplot.csv"), row.names = FALSE)
file.copy(file.path(SRC1S, "FigS1a_markers_dotplot.csv"),
          file.path(SRC1, "FigS1a_markers_dotplot.csv"), overwrite = TRUE)

cat("\nSanity: Cd4 / Cd8b1 mean expm1 by class\n")
print(dot %>% filter(gene %in% c("Cd4", "Cd8b1", "Cd8a", "Klrb1c", "Ms4a1")) %>%
        select(celltype, gene, avg, pct) %>%
        arrange(gene, desc(avg)))

cat("\nCells retained:", n_in, "\n")
cat("DONE prepare ->", FIG1, "and", FIGS1, "\n")
