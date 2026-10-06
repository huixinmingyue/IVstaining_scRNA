Fig10_Nature_v1 - Nature-style Figure 10 package
================================================

Fig. 10 | Vascular-contact microglia after MCAO: state shift, day-3
territory, and cluster-by-contact modules.

This folder is new. It does not overwrite any earlier figure package.

Panel map
---------
a    UMAP of all sorted immune cells; microglial clusters 0/3/5 in colour
b    Hexb / Siglech / P2ry12 / Mki67 on the same embedding
c    Microglial state composition, MCAO_NEG versus MCAO_POS, both mice
d    Nine state-defining markers across clusters 0, 3 and 5
e    Volcano, cluster 5 versus cluster 0
f    MCAO_POS fraction plus Lgals3 / Mki67 / P2ry12 / Prdx1 maps
g    Detection rate by POS dominance (Wilson 95% CI, Fisher)
h    Spearman rho against the continuous MCAO_POS fraction
i    Mean z-scored module scores by cluster and contact

Rebuild
-------
Rscript scripts/make_figure10.R
Rscript scripts/qa_crops.R

Canvas: 183 x 214 mm. Backend: R only.

Requires R with ggplot2, dplyr, tidyr, patchwork, ggrepel, scales,
svglite, ragg, systemfonts, scattermore, cowplot and grid, plus magick
for the QA crops.

Deliverables
------------
  Figure10.pdf / .svg / .png / .tiff
  Figure10_legend.txt, Figure10_stats.txt, Figure10_QA.txt,
  Figure10_provenance.txt, figure_contract.txt
  source_data/            input tables used to draw the figure
  source_data_plotted/    per-panel CSVs written by the plotting script
  source_figures/         parent figure and module-heatmap source displays
  scripts/                plotting and QA scripts
