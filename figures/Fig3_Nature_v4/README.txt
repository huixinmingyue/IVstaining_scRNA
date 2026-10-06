Fig3_Nature_v4 - Nature-style Figure 3 package
==============================================

Figure 3 | Hard-gated CD4 and CD8 T cells occupy naive, memory and
effector islands and shift with intravascular labelling.

This folder is new. It does not overwrite Fig3_Nature_v3.

Panel map
---------
a  IV-split UMAP with functional labels on CD4 1a/1b and CD8 0/3a/3b
b  composition (Fig. 1e style)
c  ImmGen heatmap of true-CD4 Intra DEGs
d  nine ImmGen modules matching panel c columns

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figure3.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales, svglite,
ragg, systemfonts and Seurat (prepare script only).
