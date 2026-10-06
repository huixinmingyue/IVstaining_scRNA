Fig. 9 | Nature-format multipanel figure
==========================================

New folder. Fig9_Nature_v1 and every earlier folder were left unchanged.

Backend: R only. Final size 183.0 x 246.6 mm.

Panels
------
a  T-cell maps plus the within-spot EAE_POS minus EAE_NEG difference
b  Microglia maps plus the same paired difference
c  H&E plus the four scores inside L1-L3
d  The same five-tile layout inside L4-L6
e  Within-spot difference for both cell types: Outside, L1-L3, L4-L6
f  Spatial Cxcl10 in EAE_POS and EAE_NEG microglia
g  Cxcl10-positive fraction, all spots, L1-L3 and L4-L6
h  Cxcl10 on the vessel mask, enrichment and distance gradient
i  Cell-to-vessel distance, all spots, L1-L3 and L4-L6

L1-L3 and L4-L6 are tested separately and are not pooled.
L1-L3 contains 3 EAE_NEG microglial spots. That bar and that violin
are drawn, and the sample size is stated on the figure and in the legend.

Rebuild
-------
Rscript scripts/make_figure9.R
