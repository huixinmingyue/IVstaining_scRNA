Fig4_Nature_v1 - Nature-style Figure 4 package
==============================================

Figure 4 | Spatial CD45-IV reference signatures in splenic B-cell,
plasma-cell and CD4 T-cell compartments.

This folder is new. It does not overwrite any earlier figure package.

Panel map
---------
a  B-cell niche (clusters 1, 5) and Bcell.spcd452 scores
b  Cluster 5 and Bcell.spcd45 scores
c  Peri-follicular niche (4, 5, 8, 9) and Plasma.spcd452 scores
d  B/T follicular niche (1, 4, 5, 8) and CD4T.spcd45 scores (cluster 8)
e  Clusters 4, 8 and CD4T.spcd45 scores (cluster 4)

Rebuild
-------
Rscript scripts/make_figure4.R

Requires R >= 4.4 with ggplot2, patchwork, scales, svglite, ragg,
systemfonts and grid.
