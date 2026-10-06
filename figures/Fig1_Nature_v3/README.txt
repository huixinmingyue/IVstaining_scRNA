Fig1_Nature_v3 - Nature-style Figure 1 package
==============================================

Figure 1 | Intravascular antibody labelling separates intravascular from
extravascular splenic immune cells for single-cell transcriptomics.
Panels: a, b immunofluorescence; c flow gating and sorting-to-scRNA-seq
workflow; d UMAP with RNA-rule CD4 T / CD8 T identities; e composition.
The marker dot plot that was panel f in v2 is now Fig. S1a.

This folder is new. It does not overwrite Fig1_Nature_v2.

Contents
--------
Figure1.pdf / .svg / .png / .tiff
Figure1_legend.txt
Figure1_QA.txt
figure_contract.txt
scripts/00_prepare_source.R   rebuild identities, composition and S1a markers
scripts/make_figure1.R        rebuild this figure
assets/                      micrographs, flow plots and workflow schematic
source_data/                 input tables, UMAP, mapping rules, identities
source_data_plotted/         exact matrices behind panels d and e
qa/                          panel crops and embedded-image dpi table

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figure1.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales, png,
svglite, ragg, systemfonts and Seurat (prepare script only). Paths inside
the scripts are absolute; edit FIG1 / OUTDIR if the package is moved.
