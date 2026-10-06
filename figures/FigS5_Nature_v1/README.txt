FigS5_Nature_v1 - Nature-style Figure S5 package
================================================

Figure S5 | Oligo-CD45 conjugate QC and Harmony B-cell annotation.

This folder is new. It does not overwrite any prior figure package.

Panel map
---------
a    Oligo-CD45 conjugate gels (source panel a)
b    Harmony B-cell UMAP (source panel b)
c    B-lineage + myeloid marker dot plot (source panel c)
d    Logistic P(Bcelltype) versus Proteinfeatures (source panel d)

Source panels e-g of the original seven-panel figure are not in Fig. S5.

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figureS5.R

Requires R >= 4.4 with ggplot2, dplyr, patchwork, ggrepel, scales,
png, svglite, ragg and systemfonts.

Folder contents
---------------
FigureS5.pdf / .svg / .png / .tiff
scripts/make_figureS5.R
scripts/00_prepare_source.R
scripts/upstream/
source_data/
source_data_plotted/
source_figures/
assets/
qa/
