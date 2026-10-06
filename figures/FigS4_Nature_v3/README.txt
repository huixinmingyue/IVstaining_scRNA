FigS4_Nature_v3 - Nature-style Figure S4 package
================================================

Figure S4 | Anatomical landmarks and Fig. 3 T-subset programmes
on the Stereo-seq spleen section.

This folder is new. It does not overwrite FigS4_Nature_v1,
FigS4_Nature_v2 or Fig4_Nature_v1.

Panel map
---------
a    Six anatomical gene maps (one lettered panel)
b    Follicle-centred radial gene profiles
c    Five T-subset transfer maps in one row
d    Paired B/T-niche scores

What changed versus v2
----------------------
Compact 183 x 270 mm layout (was 400 mm). Vertical colour bars on
spatial maps. Panel c is five maps in one row (no empty cell).
Radial legend sits above the traces. Display names still match
Fig. 3 v4 a/b.

Cluster-name alignment (Fig. 3 v4 a/b)
--------------------------------------
CD4-1a  ->  Memory/activated CD4 T
CD4-1b  ->  Naive CD4 T
CD8-0   ->  Naive CD8 T
CD8-3a  ->  Memory-like CD8 T
CD8-3b  ->  Effector/cytotoxic CD8 T

Rebuild
-------
Rscript scripts/00_prepare_source.R
Rscript scripts/make_figureS4.R

Requires R >= 4.4 with ggplot2, dplyr, tidyr, patchwork, scales,
svglite, ragg and systemfonts.

Folder contents
---------------
FigureS4.pdf / .svg / .png / .tiff
scripts/make_figureS4.R
scripts/00_prepare_source.R
scripts/upstream/
source_data/
source_data_plotted/
source_figures/
qa/
