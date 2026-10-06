#!/usr/bin/env Rscript
# =============================================================================
# Fig. S2-2 — Per-cluster IntraPos vs IntraNeg volcanoes (Nature SI)
#
# Core conclusion
#   Within refined CD4 T-cell subclusters (Fig. S2c: 0, 1a, 1b, 3a, 3b),
#   IntraPos and IntraNeg cells show distinct DEG programmes after filtering
#   ribosomal / mitochondrial / Gm noise genes. Cluster 2 omitted (panel f removed).
#
# Panel map (a–e): one volcano per retained cluster
# Thresholds: padj < 0.001, |log2FC| >= 0.5; exclude Rpl/Rps/mt/Gm
# Labels: all remaining significant genes; no rings; no leader lines
# Archetype: quantitative grid | Backend: R native ggplot2 vector
# Size: 183 x 290 mm | Exports: PDF / SVG / PNG (600 dpi)
# New folder FigS2-2_Nature_volcano_noC2 — does not overwrite prior packages
# =============================================================================

options(warn = 1)
set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(systemfonts)
  library(cowplot)
  library(grid)
})
# --- path resolver (GitHub repo IVstaining_scRNA) ---
.iv <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)[1])
  if (!length(f) || is.na(f) || !nzchar(f)) {
    fig <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    list(script = NA_character_, fig = fig, root = NA_character_,
         fig_root = dirname(fig))
  } else {
    sp <- normalizePath(f, winslash = "/", mustWork = TRUE)
    d <- dirname(sp)
    if (basename(d) %in% c("upstream", "qa", "original_scripts")) d <- dirname(d)
    fig <- if (identical(basename(d), "scripts")) dirname(d) else d
    root <- fig
    repeat {
      if (dir.exists(file.path(root, "figures")) &&
          file.exists(file.path(root, "README.md"))) break
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    list(script = sp, fig = fig, root = root,
         fig_root = file.path(root, "figures"))
  }
})
if (!exists("OUTDIR") || is.null(OUTDIR)) OUTDIR <- .iv$fig
IV_ROOT     <- .iv$root
IV_FIG_ROOT <- .iv$fig_root
IV_FIG      <- .iv$fig
# --- end path resolver ---

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/FigS2-2_Nature_volcano_noC2")
SRC <- file.path(OUT, "source_data")
DEGDIR <- file.path(OUT, "DEG_tables")
dir.create(file.path(OUT, "scripts"), showWarnings = FALSE, recursive = TRUE)

fc_cut <- 0.5
padj_cut <- 0.001
BASE <- 6.2
FAM <- "Arial"
FIG_W_MM <- 183
FIG_H_MM <- 290

T_LEVELS <- c("0", "1a", "1b", "3a", "3b")
PANEL_ORDER <- c("0", "1a", "1b", "3a", "3b")
panel_letters <- setNames(letters[seq_along(PANEL_ORDER)], PANEL_ORDER)
panel_titles <- setNames(paste0("Cluster ", T_LEVELS), T_LEVELS)

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args <- list(reg()$Helvetica)
      names(args) <- FAM
      try(do.call(reg, args), silent = TRUE)
    }
  }
}

is_exclude <- function(g) {
  grepl("^(Rpl|Rps|Mrpl|Mrps|mt-|Mt-|Gm)", g, ignore.case = TRUE) |
    g %in% c("AY036118")
}

assign_dir <- function(df) {
  df %>%
    mutate(
      direction = case_when(
        !is.na(p_val_adj) & p_val_adj < padj_cut & avg_log2FC >= fc_cut ~ "IntraPos up",
        !is.na(p_val_adj) & p_val_adj < padj_cut & avg_log2FC <= -fc_cut ~ "IntraNeg up",
        TRUE ~ "NS"
      ),
      neglog10_padj = pmin(-log10(pmax(p_val_adj, 1e-300)), 50)
    )
}

theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = BASE),
      axis.text = element_text(size = BASE - 0.5, colour = "black"),
      plot.title = element_text(size = BASE + 0.3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE - 1.0, hjust = 0.5, colour = "grey35"),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.margin = margin(2, 4, 2, 2),
      panel.grid = element_blank()
    )
)

# ---- Load DEGs ----
deg_rds <- file.path(SRC, "deg_list.rds")
if (!file.exists(deg_rds)) stop("Missing ", deg_rds, " — run 01_compute_DEG.R first")
deg_raw <- readRDS(deg_rds)
n_tab <- read.csv(file.path(SRC, "n_cells_by_cluster_IV.csv"), stringsAsFactors = FALSE)
rownames(n_tab) <- n_tab$cluster

deg_list <- lapply(T_LEVELS, function(cl) {
  if (is.null(deg_raw[[cl]])) return(NULL)
  assign_dir(deg_raw[[cl]])
})
names(deg_list) <- T_LEVELS

# Volcano: coloured side-stack labels, no rings, no leader lines
volcano_panel <- function(deg, cl, tag) {
  d <- deg
  d$direction <- factor(d$direction, levels = c("IntraNeg up", "NS", "IntraPos up"))
  lab <- d %>% filter(direction != "NS", !is_exclude(gene)) %>% arrange(p_val_adj)
  nlab <- nrow(lab)
  n_left <- sum(as.character(lab$direction) == "IntraNeg up")
  n_right <- sum(as.character(lab$direction) == "IntraPos up")
  n_side_max <- max(n_left, n_right, 1)
  ncols_guess <- if (n_side_max > 70) 3L else if (n_side_max > 35) 2L else 1L
  n_row_guess <- ceiling(n_side_max / ncols_guess)
  lab_size <- max(1.8, min(2.4, 72 / max(n_row_guess, 8)))
  xmax_data <- max(3, quantile(abs(d$avg_log2FC), 0.999, na.rm = TRUE) * 1.05, fc_cut * 1.2)
  # Extra side room for dense panels (many labelled genes)
  xmax <- xmax_data * (1.35 + 0.65 * ncols_guess)
  ymax <- max(5, max(d$neglog10_padj, na.rm = TRUE) * 1.05)

  n_neg <- n_tab[cl, "IntraNeg"]
  n_pos <- n_tab[cl, "IntraPos"]
  subtit <- sprintf(
    "n Neg=%s; Pos=%s; labelled=%d",
    format(n_neg, big.mark = ","), format(n_pos, big.mark = ","), nlab
  )

  p <- ggplot(d, aes(avg_log2FC, neglog10_padj, colour = direction)) +
    geom_point(size = 0.4, stroke = 0, alpha = 0.7) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(padj_cut), linewidth = 0.25,
               linetype = "dashed", colour = "grey50") +
    scale_colour_manual(
      values = c("IntraNeg up" = "#2F6F8F", "NS" = "#B0B0B0", "IntraPos up" = "#B2473E"),
      name = NULL, drop = FALSE
    ) +
    coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), clip = "off") +
    labs(
      x = expression(log[2] * "FC (IntraPos / IntraNeg)"),
      y = expression(-log[10] * " adjusted " * italic(P)),
      title = panel_titles[[cl]],
      subtitle = subtit,
      tag = tag
    ) +
    theme(
      legend.position = "none",
      plot.tag.position = c(0.02, 0.98),
      plot.margin = margin(2, 5, 2, 2)
    )

  if (nlab > 0) {
    stack_side <- function(dat, side = c("left", "right")) {
      side <- match.arg(side)
      n <- nrow(dat)
      if (!n) return(dat[0, ])
      dat <- dat %>% arrange(desc(neglog10_padj), desc(abs(avg_log2FC)))
      ncols <- if (n > 70) 3L else if (n > 35) 2L else 1L
      n_row <- ceiling(n / ncols)
      dat$col <- rep(seq_len(ncols), each = n_row)[seq_len(n)]
      dat$row <- ave(seq_len(n), dat$col, FUN = seq_along)
      y_top <- ymax * 0.985
      y_bot <- ymax * 0.015
      gap <- (xmax - xmax_data * 1.05) / max(ncols, 1)
      if (side == "left") {
        dat$lab_x <- -xmax + 0.06 * xmax_data + (dat$col - 1) * gap
      } else {
        dat$lab_x <- xmax - 0.06 * xmax_data - (dat$col - 1) * gap
      }
      dat$lab_y <- if (n_row == 1) {
        (y_top + y_bot) / 2
      } else {
        y_top - (dat$row - 1) * (y_top - y_bot) / (n_row - 1)
      }
      dat
    }

    lab_L <- stack_side(lab %>% filter(as.character(direction) == "IntraNeg up"), "left")
    lab_R <- stack_side(lab %>% filter(as.character(direction) == "IntraPos up"), "right")

    if (nrow(lab_L)) {
      p <- p + geom_text(
        data = lab_L, aes(lab_x, lab_y, label = gene, colour = direction),
        size = lab_size, fontface = "italic", family = FAM,
        hjust = 0, vjust = 0.5, lineheight = 0.8,
        show.legend = FALSE, inherit.aes = FALSE
      )
    }
    if (nrow(lab_R)) {
      p <- p + geom_text(
        data = lab_R, aes(lab_x, lab_y, label = gene, colour = direction),
        size = lab_size, fontface = "italic", family = FAM,
        hjust = 1, vjust = 0.5, lineheight = 0.8,
        show.legend = FALSE, inherit.aes = FALSE
      )
    }
  }
  list(plot = p, labelled = lab)
}

plots <- list()
lab_rows <- list()
summary_rows <- list()

for (cl in PANEL_ORDER) {
  if (is.null(deg_list[[cl]])) {
    cat("WARN: no DEG for cluster", cl, "\n")
    next
  }
  res <- volcano_panel(deg_list[[cl]], cl, panel_letters[[cl]])
  plots[[cl]] <- res$plot
  lab <- res$labelled
  d <- deg_list[[cl]]
  n_excl <- sum(d$direction != "NS" & is_exclude(d$gene))
  if (nrow(lab) > 0) {
    lab$panel <- panel_letters[[cl]]
    lab$cluster <- cl
    lab_rows[[cl]] <- lab
  }
  summary_rows[[cl]] <- tibble(
    panel = panel_letters[[cl]],
    cluster = cl,
    n_IntraNeg = n_tab[cl, "IntraNeg"],
    n_IntraPos = n_tab[cl, "IntraPos"],
    n_genes_tested = nrow(d),
    n_IntraPos_up = sum(d$direction == "IntraPos up"),
    n_IntraNeg_up = sum(d$direction == "IntraNeg up"),
    n_sig_excluded_RplRps_mt_Gm = n_excl,
    n_labelled = nrow(lab),
    fc_cutoff = fc_cut,
    padj_cutoff = padj_cut
  )
  cat(sprintf(
    "Panel %s (cluster %s): labelled %d (Pos-up %d; Neg-up %d; excluded noise %d)\n",
    panel_letters[[cl]], cl, nrow(lab),
    if (nrow(lab)) sum(lab$direction == "IntraPos up") else 0L,
    if (nrow(lab)) sum(lab$direction == "IntraNeg up") else 0L,
    n_excl
  ))
}

lab_all <- bind_rows(lab_rows)
write.csv(
  lab_all %>%
    select(panel, cluster, gene, direction, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
  file.path(SRC, "FigS2-2_labelled_genes.csv"),
  row.names = FALSE
)
write.csv(
  lab_all %>%
    select(panel, cluster, gene, direction, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
  file.path(OUT, "FigS2-2_labelled_genes.csv"),
  row.names = FALSE
)

# Also write per-cluster filtered significant tables (before noise filter + after)
for (cl in T_LEVELS) {
  if (is.null(deg_list[[cl]])) next
  d <- deg_list[[cl]]
  sig <- d %>%
    filter(direction != "NS") %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
  write.csv(
    sig,
    file.path(DEGDIR, sprintf("DEG_sig_cluster_%s_padj001_FCge0.5.csv", cl)),
    row.names = FALSE
  )
  write.csv(
    sig %>% filter(!is_exclude(gene)),
    file.path(DEGDIR, sprintf("DEG_sig_labelled_cluster_%s.csv", cl)),
    row.names = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(SRC, "FigS2-2_DEG_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(OUT, "FigS2-2_DEG_summary.csv"), row.names = FALSE)

# Shared legend
leg_plot <- plots[["0"]] +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
shared_leg <- cowplot::get_legend(leg_plot)

# Reading order: a=0, b=1a, c=1b, d=3a, e=3b (cluster 2 / former panel f omitted)
fig <- (
  plots[["0"]] /
  plots[["1a"]] /
  (plots[["1b"]] | plots[["3a"]] | plots[["3b"]]) /
  wrap_elements(panel = shared_leg)
) +
  plot_layout(heights = c(1.25, 1.25, 1.05, 0.08))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUT, "FigureS2-2")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())

cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in, family = FAM)
print(fig); invisible(dev.off())

cat("Writing PNG ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())

writeLines(c(
  "Fig. S2-2 | IntraPos versus IntraNeg differential expression within refined CD4 T-cell subclusters.",
  "",
  "a–e, Volcano plots of IntraPos versus IntraNeg DEGs within each refined T-cell subcluster",
  "from Fig. S2c (a, cluster 0; b, 1a; c, 1b; d, 3a; e, 3b). Cluster 2 is omitted.",
  "Wilcoxon rank-sum tests (Seurat FindMarkers); Benjamini–Hochberg adjusted P.",
  "Genes with adjusted P < 0.001 and |log2 fold change| >= 0.5 are coloured",
  "(red, IntraPos up; blue, IntraNeg up). Ribosomal (Rpl/Rps), mitochondrial (mt-)",
  "and Gm* genes are excluded from labelling; all remaining significant genes are",
  "labelled without highlight rings or leader lines. n cells per IV arm are shown",
  "in each panel subtitle. Source data are provided as a Source Data file."
), file.path(OUT, "FigureS2-2_legend.txt"))

writeLines(c(
  "Figure S2-2 | Nature SI quantitative-grid volcano package (no cluster 2)",
  paste0("Generated: ", Sys.time()),
  paste0("Final size: ", FIG_W_MM, " x ", FIG_H_MM, " mm"),
  "Exports: FigureS2-2.pdf / .svg / .png (600 dpi)",
  "Backend: R exclusive; native ggplot2 vector",
  sprintf("Thresholds: padj < %g; |log2FC| >= %g", padj_cut, fc_cut),
  "Noise filter (labels): Rpl/Rps/Mrpl/Mrps/mt-/Gm*/AY036118",
  "No highlight rings; no leader lines; side-stack italic gene labels (>= ~5 pt)",
  "Layout: a=0, b=1a full-width; c=1b, d=3a, e=3b row; cluster 2 omitted",
  "New folder FigS2-2_Nature_volcano_noC2/ — does not overwrite prior packages."
), file.path(OUT, "FigureS2-2_QA.txt"))

cat("DONE ->", OUT, "\n")
print(as.data.frame(summary_df))
