FigS3_Nature_v4 - Nature-style Figure S3 package
================================================

Figure S3 | Gate detection, lineage markers, and table-marker support
for the Fig. 3 functional CD4 / CD8 / NKT / doublet labels.

This folder is new. It does not overwrite FigS3_Nature_v3.

Panel map
---------
a    RNA detection after the hard gate
b-e  Cd4 / Cd8b1 / Klrb1c / Ms4a1
f    Dot plot of annotation-table genes on the seven new identities
g-j  Ccr7 / S100a4 / Eomes / Ccl5 (same format as b-e)

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figureS3.R

Requires R >= 4.4 with Seurat (prepare only), ggplot2, dplyr, tidyr,
patchwork, scales, svglite, ragg and systemfonts.
