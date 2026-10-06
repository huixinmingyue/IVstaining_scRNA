# IVstaining_scRNA

Figure reproduction scripts for a Nature-style manuscript on intravascular
antibody labelling, splenic immune-cell scRNA-seq, and related spatial
analyses.

This repository contains the **plotting scripts** and **tabular source data**
used to assemble main Figures 1–10 and Supplementary Figures S1–S12. Rendered
figures (PDF/PNG/TIFF), microscopy assets, and large Seurat/RDS objects are
not included.

## Repository layout

```
figures/
  Fig1_Nature_v3/          # Figure 1 package
    scripts/               # make_figure*.R and optional 00_prepare_source.R
      upstream/            # earlier versions kept for provenance
    source_data/           # CSV/TSV/TXT inputs
    source_data_plotted/   # exact matrices behind plotted panels (where present)
    README.txt             # panel map and rebuild notes
  Fig2_Nature_v5/
  ...
  FigS12_Nature_v1/
R/setup_paths.R            # shared path helpers
```

| Manuscript panel | Folder |
| --- | --- |
| Fig. 1 | `figures/Fig1_Nature_v3` |
| Fig. 2 | `figures/Fig2_Nature_v5` |
| Fig. 3 | `figures/Fig3_Nature_v4` |
| Fig. 4 | `figures/Fig4_Nature_v1` |
| Fig. 5 | `figures/Fig5_Nature_v1` |
| Fig. 6 | `figures/Fig6_Nature_v3` |
| Fig. 7 | `figures/Fig7_Nature_v3` |
| Fig. 8 | `figures/Fig8_Nature_v4` |
| Fig. 9 | `figures/Fig9_Nature_v2` |
| Fig. 10 | `figures/Fig10_Nature_v1` |
| Fig. S1 | `figures/FigS1_Nature_v3` |
| Fig. S2 | `figures/FigS2_Nature_v1` |
| Fig. S3 | `figures/FigS3_Nature_v4` |
| Fig. S4 | `figures/FigS4_Nature_v3` |
| Fig. S5 | `figures/FigS5_Nature_v1` |
| Fig. S6 | `figures/FigS6_Nature_v2` |
| Fig. S7 | `figures/FigS7_Nature_v3` |
| Fig. S8 | `figures/FigS8_Nature_v1` |
| Fig. S9 | `figures/FigS9_Nature_v3` |
| Fig. S10 | `figures/FigS10_Nature_v1` |
| Fig. S11 | `figures/FigS11_Nature_v1` |
| Fig. S12 | `figures/FigS12_Nature_v1` |

## Requirements

- R >= 4.4
- Common packages: `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `scales`,
  `svglite`, `ragg`, `systemfonts`, `png`, `grid`, `scattermore`, `ggrepel`,
  `cowplot`
- Prepare scripts (`00_prepare_source.R`) additionally need `Seurat` and
  `Matrix` when rebuilding identities from objects

Install the plotting stack:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "scales",
  "svglite", "ragg", "systemfonts", "png", "scattermore",
  "ggrepel", "cowplot"
))
```

## Rebuild a figure

From the repository root (or any working directory):

```bash
Rscript figures/Fig1_Nature_v3/scripts/make_figure1.R
```

If a folder contains `00_prepare_source.R`, run that first. Paths are resolved
from the script location; hardcoded `F:/IVstaining_scRNA/...` figure-package
paths have been replaced with `IV_FIG_ROOT`.

Microscopy plates and some spatial rasters live in local `assets/` directories
that are not in this repository, so image panels will be skipped or fail unless
those files are supplied separately.

## What is not in this repository

- Rendered figure files (PDF, SVG, PNG, TIFF)
- Immunofluorescence / flow / H&E assets
- Seurat and other RDS objects (often >10–50 MB)

Scripts that still need a local analysis tree (Seurat objects, ImmGen tables,
legacy figure folders) read

```r
file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "...")
```

Set `IV_RAW_DATA` to the root of that local tree if you have it. These
`scripts/upstream/` files document how source tables were generated; they are
not required to inspect the plotted CSVs.

## Licence and citation

Code is provided for manuscript reproducibility. Please cite the accompanying
paper when using these scripts.
