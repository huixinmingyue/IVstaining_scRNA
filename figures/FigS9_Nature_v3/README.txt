FigS9_Nature_v3 - Nature-style Figure S9 package
===============================================

Fig. S9 | Spatial cluster markers and POS versus NEG immune-lineage
mapping in the EAE Visium section.

This folder is new. It does not overwrite FigS9_Nature_v1 or
FigS9_Nature_v2.

Changed relative to v2
----------------------
Panel b replaces the Ccl2 map with Cxcl10.
Cxcl10 uses the source spatial-map colour scale, fixed at 0-2.5.
The other panels are unchanged.

Panel map
---------
a  Top marker heatmap, Visium clusters 0-8
b  Spatial maps of Mki67, Lgals3, Cxcl10, Hexb
c  DC POS | NEG
d  B-cell POS | NEG
e  Neu1 POS | NEG
f  Plasma-cell POS | NEG
g  Macrophage POS | NEG
h  All-spot and window (L1–L3 / L4–L6) score boxplots

Rebuild
-------
Rscript scripts/make_figureS9.R
Rscript qa/preview.R

Canvas: 183 x 248 mm.
