FigS6_Nature_v2 - Nature-style Figure S6 package
================================================

Figure S6 | Oligo-CD45 labelling of splenic T cells from the same
CITE-seq experiment as Fig. 5.

This folder is new. It does not overwrite FigS6_Nature_v1.

UMAP rule
---------
Panels a and d are not leftover cells plotted on the parent T-cell
UMAP. First-pass B-T doublet and myeloid-contaminated clusters are
removed, remaining T cells are reclustered (PCA 1:15, res 0.2), and
panel d is the high/low oligo-CD45 subset of that embedding.

Panel map
---------
a    Reclustered T-cell UMAP (Fig. 3 v4 taxonomy)
b    Oligo-CD45 copies by subtype
c    Subtype share along the labelling gradient
d-e  High (>= 6) vs low (<= 1) tails: UMAP subset and composition
f    Memory/activated CD4 T differential expression at the same thresholds

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figureS6.R

Requires R >= 4.4 with Seurat (prepare only), ggplot2, dplyr, patchwork,
ggrepel, scales, svglite, ragg and systemfonts.
