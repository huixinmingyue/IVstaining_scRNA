FigS2_Nature_v1 - assembled Fig. S2
===================================

New folder. Does not overwrite any previous FigS2 package or Fig2_Nature_v4.

Panel map
---------
a  Query-marker violins (P3 panel a) plus unique positioning genes (P6 panel a)
b  ImmGen identity-gene heatmap (P2 panel c)
c  Within-cluster IV GSEA after nCount regression (P4 entire figure)
d  Signaling subset of the same within-cluster IV GSEA (identifiers as in c)

Contents
--------
FigureS2.pdf / .svg / .png / .tiff
FigureS2_legend.txt
FigureS2_QA.txt
figure_contract.txt
scripts/make_figureS2.R
scripts/00_extract_violin_expr.R
source_data/                 copied inputs plus cell-level violin extract
source_data_plotted/         exact plotted extracts
source_figures/              original source figures and original scripts

Rebuild
-------
Rscript scripts/00_extract_violin_expr.R   # only if FigS2a_violin_expr.csv is missing
Rscript scripts/make_figureS2.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, stringr, ggrepel,
svglite, ragg, systemfonts. The extract script also needs Seurat.
