Fig8_Nature_v4 - Nature-style Figure 8 package
==============================================

Fig. 8 | Perivascular windows and vessel-contact scores in the EAE brain.

This folder is new. It does not overwrite Fig8_Nature_v1, v2, or v3.

Changed relative to v3
----------------------
Panel a is a vector schematic of EAE sectioning, slice collection,
and Visium spatial transcriptomics (replaces the previous raster).
Panel e keeps H&E and the cluster map; L1–L6 zooms are removed.
Panel h adds microglia POS/NEG window crops for L1–L3 and L4–L6
(Fig. S6 d layout) plus All immune spots / Windows boxplots
(Fig. 8 f–g statistics).

Panel map
---------
a    Vector sampling schematic (section, slice, Visium)
b    Sagittal Laminin/DAPI
c    Coronal IF with L1–L6 labels
d    L1–L3 IF zooms plus L4–L6 crops
e    H&E and cluster map with six windows (no zooms)
f    T-cell whole-section POS/NEG maps; all-spot and window boxplots
g    Microglia whole-section POS/NEG maps; all-spot and window boxplots
h    Microglia POS/NEG window crops L1–L6; all-spot and window boxplots

Rebuild
-------
Rscript scripts/make_figure8.R

Canvas: 183 x 470 mm.
