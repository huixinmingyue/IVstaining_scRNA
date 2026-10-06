Source data for Figure S8 (copied; not modified)

traj_slim.rds
  Monocle DDRTree cells (comp1, comp2, State, group, Pseudotime),
  MST edges, and branch points. n = 1,807 microglia.

heatmap_slim.rds
  50 genes x 8 pseudotime bins (z-scored bin-mean expression)
  plus gene_meta (group, peak bin, Spearman rho vs pseudotime).

Fig6f_violin.csv
  Long-format normalized expression for P2ry12, Hexb, Il1b, Cxcl2,
  Cxcl10, Id2 by Seurat cluster and EAE_POS / EAE_NEG.

Fig6f_wilcoxon_BH.csv
  Two-sided Wilcoxon rank-sum tests within gene x cluster,
  Benjamini-Hochberg-adjusted across 24 tests (from Fig. 6 v3 panel f).
