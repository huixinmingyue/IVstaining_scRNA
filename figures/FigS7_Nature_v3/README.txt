FigS7_Nature_v3 - Nature-style Figure S7 package
================================================

Figure S7 | Marker-defined CNS immune clusters, CellChat counts,
T/B KEGG signalling, and locked genes in all cells.

This folder is new. It does not overwrite FigS7_Nature_v1 or v2.

Panel map
---------
a    Top marker heatmap
b    Canonical lineage marker dot plot
c    CellChat interaction-number networks (was v2 g)
d    T-cell KEGG signalling GSEA
e    B-cell KEGG signalling GSEA
f    T-cell locked genes, all cells only (was v2 e)
g    B-cell locked genes, all cells only (was v2 f)

Rebuild
-------
Rscript scripts/make_figureS7.R
Rscript qa/preview.R
