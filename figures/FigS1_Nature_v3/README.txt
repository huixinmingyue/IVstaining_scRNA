FigS1_Nature_v3 - assembled Fig. S1
===================================

New folder. Does not overwrite FigS1_Nature_v2, the QC audit package,
FigS1_Nature_b_allDEG, or the reannotation package.

Panel map
---------
a  QC violins by IV (from the QC/DEG/plasma figure, panel a)
b  lineage marker dot plot (from the previous Fig. S1 marker package)
c  B-cell volcano, MZ vs Fo intersection labels (from the all-DEG figure, panel a)
d  true CD4 volcano (from the hard-gate reannotation figure, panel d)
e  true CD8 volcano (from the hard-gate reannotation figure, panel e)
f  plasma genes without Aicda (from the QC/DEG/plasma figure, panel f)

Contents
--------
FigureS1.pdf / .svg / .png / .tiff
FigureS1_legend.txt
FigureS1_QA.txt
figure_contract.txt
scripts/make_figureS1.R
source_data/                 copied inputs used to draw each panel
source_data_plotted/         exact plotted extracts
source_figures/              original source figures and original scripts

Rebuild
-------
"C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe" scripts\make_figureS1.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales, cowplot,
svglite, ragg, systemfonts.
