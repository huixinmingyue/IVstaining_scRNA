Fig7_Nature_v3 - Nature-style Figure 7 package
==============================================

Fig. 7 | Vessel-contacting microglia acquire an IFN-γ / IL-17-chemokine
programme in EAE.

This folder is new. It does not overwrite Fig7_Nature_v1 or v2.

Changed relative to v2
----------------------
Previous f is split into f (trajectory) and g (Cxcl10 / Id2).
Those two share a row with the immunofluorescence plate, now h.

Panel map
---------
a  UMAP of EAE microglia split by POS / NEG
b  Cluster composition
c  KEGG GSEA (signalling terms)
d  IFN-γ / IL-17-axis module scores, All cells only
e  Leading-edge heatmap (IL-17, NOD-like, Chemokine)
f  Pseudotime trajectory (pseudotime and cluster)
g  Cxcl10 and Id2 along pseudotime
h  Immunofluorescence of the vessel-contact niche

Rebuild
-------
Rscript scripts/make_figure7.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales,
svglite, ragg, systemfonts, png.

Canvas: 183 x 210 mm.

Outputs
-------
Figure7.pdf  Figure7.svg  Figure7.png  Figure7.tiff
Figure7_legend.txt  Figure7_stats.txt  Figure7_QA.txt
source_data/          tables and RDS used for plotting
source_figures/       source rasters
source_data_plotted/  tables actually drawn
scripts/make_figure7.R
