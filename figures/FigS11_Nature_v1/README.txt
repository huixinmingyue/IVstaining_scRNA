Fig. S11 package

Redrawn from the supplied Figure S8 plate (copied into source_figures/).
Same observations, same panel order, same colour mapping. Figure number, folder, journal width and TIFF export are new. No previous folder was overwritten.

How to rebuild
1. scripts/00_extract_spatial.R writes coordinates, prediction scores and four-gene expression into source_data/. It reads the Visium Seurat object named in that script.
2. scripts/make_figureS11.R reads source_data/ and writes FigureS11.svg, FigureS11.pdf, FigureS11.png and FigureS11.tiff.

Contents
- FigureS11.svg / .pdf / .png / .tiff
- scripts/ plotting and extraction code
- source_data/ matrices used to draw the figure, including the H&E raster and the original spot tables
- source_data_plotted/ tables written by the plotting script
- source_figures/ the source plate and its original script
- figure_contract.txt, FigureS11_legend.txt, FigureS11_QA.txt
