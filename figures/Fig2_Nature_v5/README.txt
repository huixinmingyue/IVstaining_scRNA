Fig2_Nature_v5 - Nature-style Figure 2 package
==============================================

Figure 2 | Intravascular-enriched B-cell subcluster 1 is MZ-like.
Panels: a UMAP by IV plus subcluster composition (Fig. 1e style);
b Cr2 / Fcer2a; c S1a MZ intersect Intra vs ImmGen; d ImmGen module scores.

This folder is new. It does not overwrite Fig2_Nature_v4.

Contents
--------
Figure2.pdf / .svg / .png / .tiff
Figure2_legend.txt
Figure2_QA.txt
figure_contract.txt
scripts/00_prepare_source.R
scripts/make_figure2.R
scripts/upstream/
assets/
source_data/
source_data_plotted/
qa/

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figure2.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales, svglite,
ragg, systemfonts and Seurat (prepare script only). Paths inside the scripts
are absolute; edit OUTDIR if the package is moved.
