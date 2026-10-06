Fig5_Nature_v1 - Nature-style Figure 5 package
==============================================

Figure 5 | In vivo oligo-CD45 labelling stratifies splenic B cells.

This folder is new. It does not overwrite any earlier figure package.

Panel map
---------
a    Delivery and CITE-seq schematic
b-c  B-cell UMAP and oligo-CD45 copies by subtype (cluster 4 removed)
d    Subtype share along the labelling gradient
e-f  High (>= 6) vs low (<= 1) tails: UMAP and composition
g    Marginal-zone B differential expression at the same thresholds

Rebuild
-------
Rscript scripts/make_figure5.R

Requires R >= 4.4 with ggplot2, dplyr, patchwork, ggrepel, scales,
svglite, ragg, png, grid and systemfonts.
