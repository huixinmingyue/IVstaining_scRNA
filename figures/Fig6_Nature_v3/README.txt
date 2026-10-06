Fig6_Nature_v3 - Nature-style Figure 6 package
==============================================

Fig. 6 | Vascular-contact CNS immune cells rewire original outgoing
signalling and myeloid programmes.

This folder is new. It does not overwrite Fig6_Nature_v1 or
Fig6_Nature_v2.

Changed relative to v2
----------------------
c and d share one row as matched horizontal POS/NEG bar charts.
g and h show the All-cells locked-gene series only.

Panel map
---------
a  IV CD45 FACS gates (NEG / POS; Extra CD45 vs Intra CD45)
b  UMAP of CNS immune cells split by POS and NEG
c  Cell-type composition
d  CellChat outgoing POS/NEG ratio, original matrices only
   (GALECTIN to CSF)
e  Macrophage KEGG ORA, signalling terms (display filter)
f  Microglia KEGG ORA, signalling terms (display filter)
g  Locked macrophage genes, all cells
h  Locked microglia genes, all cells

Rebuild
-------
Rscript scripts/make_figure6.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, tibble, patchwork,
scales, svglite, ragg, systemfonts, png, ggrepel, scattermore and magick.

Outputs
-------
Figure6.pdf  Figure6.svg  Figure6.png  Figure6.tiff
Figure6_legend.txt  Figure6_stats.txt  Figure6_QA.txt
source_data/          tables and RDS used for plotting
source_figures/       FACS TIFF and source figure rasters
source_data_plotted/  tables actually drawn
scripts/make_figure6.R
