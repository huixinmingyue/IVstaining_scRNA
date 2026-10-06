#!/usr/bin/env Rscript
# =============================================================================
# Fig. S4 v3 - Spatial anatomy markers, follicle radius, and T-subset transfer
#
# Core conclusion
#   Anatomical genes recover follicle / MZ / red-pulp landmarks on the
#   Stereo-seq section. Follicle-centred radius traces Cr2 at the edge and
#   Fcer2a inside. Label transfer of Fig. 3 naive / memory / effector
#   programmes onto the same section places memory/activated CD4 T over
#   naive CD4 T in B/T niches.
#
# Evidence chain
#   a  Six anatomical gene maps (former markers a-f as one panel)
#   b  Follicle-centred radial profiles (former radial panel a)
#   c  Prediction maps for Fig. 3 v4 functional T identities
#   d  Paired scores in B/T niches (former T maps f-g; Fig. 3 v4 names)
#
# Archetype : image plate + quant
# Backend   : R exclusive
# Size      : 183 x 270 mm
# New folder — does not overwrite any prior package
# =============================================================================

options(warn = 1)
set.seed(42)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(systemfonts)
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

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
OUTDIR <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  file.path(IV_FIG_ROOT, "FigS4_Nature_v3")
}
SRC <- file.path(OUTDIR, "source_data")
SDP <- file.path(OUTDIR, "source_data_plotted")
QA  <- file.path(OUTDIR, "qa")
dir.create(SDP, showWarnings = FALSE)
dir.create(QA, showWarnings = FALSE)

FIG_W_MM <- 183
FIG_H_MM <- 270
width_mm <- 183
height_mm <- 270
BASE <- 6.5
FAM  <- "Arial"
PT <- function(x) x / .pt

fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
  cat("font:", FAM, "->", fm$path, "\n")
  for (reg in list(grDevices::postscriptFonts, grDevices::pdfFonts)) {
    if (!FAM %in% names(reg())) {
      args <- list(reg()$Helvetica); names(args) <- FAM
      try(do.call(reg, args), silent = TRUE)
    }
  }
}

UM_PER_UNIT <- 0.5
STEP <- 50
SCALEBAR_UM <- 1000
score_ramp <- grDevices::hcl.colors(256, "YlGnBu", rev = TRUE)[52:256]
proj_ramp <- grDevices::hcl.colors(256, "Plasma")[1:235]
BG_COL <- "#EBEBEB"

# Fig. 3 v4 a/b functional names and colours for the five transferred islands.
wrap_fun <- function(s) sub(" (CD[48] T)$", "\n\\1", s)
FIG3_LAB <- c(
  CD4_1a = "Memory/activated CD4 T",
  CD4_1b = "Naive CD4 T",
  CD8_0  = "Naive CD8 T",
  CD8_3a = "Memory-like CD8 T",
  CD8_3b = "Effector/cytotoxic CD8 T"
)
FIG3_COL <- c(
  "Memory/activated CD4 T"   = "#E07A5F",
  "Naive CD4 T"              = "#F2B6A0",
  "Naive CD8 T"              = "#3D5A80",
  "Memory-like CD8 T"        = "#2A9D8F",
  "Effector/cytotoxic CD8 T" = "#7BC9BE"
)
GENE_LAB <- c(
  Cr2 = "Cr2 (CD21)", Fcer2a = "Fcer2a (CD23)", Cd4 = "Cd4",
  Jchain = "Jchain", Marco = "Marco", Adgre1 = "Adgre1 (F4/80)"
)
GENE_COL <- c(
  "Cr2 (CD21)" = "#0072B2", "Fcer2a (CD23)" = "#E69F00",
  Cd4 = "#CC79A7", Jchain = "#009E73",
  Marco = "#D55E00", "Adgre1 (F4/80)" = "#000000"
)

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title = element_text(size = BASE, colour = "black"),
    axis.text = element_text(size = BASE - 0.5, colour = "black"),
    legend.title = element_text(size = BASE - 0.2),
    legend.text = element_text(size = BASE - 0.6),
    legend.key.size = unit(2.4, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 0.7, hjust = 0.5, colour = "grey30"),
    plot.margin = margin(1.2, 1.5, 1.2, 1.5, "mm"),
    panel.grid = element_blank()
  )
theme_sp <- theme_void(base_family = FAM) +
  theme(
    plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5,
                              family = FAM, lineheight = 0.92,
                              margin = margin(b = 0.6, unit = "mm")),
    legend.title = element_text(size = 5.2, family = FAM, hjust = 0.5),
    legend.text = element_text(size = 5.2, family = FAM),
    plot.margin = margin(0.6, 0.4, 0.4, 0.6, "mm")
  )

tag_it <- function(p, t) {
  p + labs(tag = t) +
    theme(plot.tag = element_text(size = 8, face = "bold", family = FAM,
                                  hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}

cat("Load section ...\n")
df <- read.csv(file.path(SRC, "SourceData_section_markers_radial.csv"),
               stringsAsFactors = FALSE)
df$xp <- as.numeric(df$x)
df$yp <- -as.numeric(df$y)
PAD <- 150
XLIM <- range(df$xp, na.rm = TRUE) + c(-PAD, PAD)
YLIM <- range(df$yp, na.rm = TRUE) + c(-PAD, PAD)

runs <- do.call(rbind, lapply(split(df, df$yp), function(g) {
  xs <- sort(unique(g$xp))
  brk <- c(0, which(diff(xs) != STEP), length(xs))
  do.call(rbind, lapply(seq_len(length(brk) - 1), function(i) {
    seg <- xs[(brk[i] + 1):brk[i + 1]]
    data.frame(xmin = min(seg) - STEP / 2, xmax = max(seg) + STEP / 2,
               ymin = g$yp[1] - STEP / 2, ymax = g$yp[1] + STEP / 2)
  }))
}))
section_layer <- function() {
  geom_rect(data = runs,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = BG_COL, colour = NA, inherit.aes = FALSE)
}
sb_len <- SCALEBAR_UM / UM_PER_UNIT
sb_x1 <- XLIM[2] - 0.06 * diff(XLIM)
sb_x0 <- sb_x1 - sb_len
sb_y  <- YLIM[1] + 0.045 * diff(YLIM)

lim01 <- function(v) {
  v <- v[is.finite(v)]
  q <- unname(quantile(v, c(0.01, 0.99)))
  if (q[1] == q[2]) q[2] <- q[1] + 1e-6
  q
}

make_feat <- function(xp, yp, val, title, scalebar = FALSE,
                      legend_name = "Expr.",
                      ramp = score_ramp, size = 0.18) {
  lim <- lim01(val)
  d <- data.frame(xp = xp, yp = yp,
                  val = pmin(pmax(as.numeric(val), lim[1]), lim[2]))
  d <- d[is.finite(d$xp) & is.finite(d$yp), , drop = FALSE]
  d <- d[order(d$val), , drop = FALSE]
  br <- pretty(lim, n = 3)
  br <- br[br >= lim[1] - 1e-12 & br <= lim[2] + 1e-12]
  if (length(br) < 2) br <- lim
  p <- ggplot() +
    section_layer() +
    geom_point(data = d, aes(xp, yp, colour = val),
               size = size, stroke = 0, shape = 16) +
    scale_colour_gradientn(colours = ramp, limits = lim, oob = squish,
                           name = legend_name, breaks = br) +
    guides(colour = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(1.5, "mm"), barheight = unit(11, "mm"),
      ticks.colour = "white", frame.colour = "black", frame.linewidth = 0.2
    )) +
    coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = title) +
    theme_sp +
    theme(legend.position = "right",
          legend.margin = margin(0, 0, 0, 0.4, "mm"),
          legend.box.spacing = unit(0.4, "mm"),
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.25))
  if (scalebar) {
    p <- p +
      annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
               linewidth = 0.55, colour = "black") +
      annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y + 0.032 * diff(YLIM),
               label = "1 mm", size = PT(5.2), family = FAM)
  }
  p
}

# -----------------------------------------------------------------------------
# a  anatomical gene maps
# -----------------------------------------------------------------------------
cat("Panel a ...\n")
plist <- vector("list", length(GENE_LAB))
for (i in seq_along(GENE_LAB)) {
  g <- names(GENE_LAB)[i]
  stopifnot(g %in% names(df))
  plist[[i]] <- make_feat(df$xp, df$yp, df[[g]], GENE_LAB[[g]],
                          scalebar = (i == 1L), legend_name = "Expr.",
                          size = 0.18)
}
pA <- wrap_plots(plist, ncol = 3)

# -----------------------------------------------------------------------------
# b  follicle-centred radial profiles
# -----------------------------------------------------------------------------
cat("Panel b ...\n")
rb <- read.csv(file.path(SRC, "SourceData_radial_rnorm.csv"),
               stringsAsFactors = FALSE)
gene_long <- bind_rows(lapply(names(GENE_LAB), function(g) {
  mcol <- paste0(g, "_mean"); scol <- paste0(g, "_sem")
  if (!mcol %in% names(rb)) return(NULL)
  data.frame(r_bin = rb$r_bin, n = rb$n,
             mean = rb[[mcol]], sem = rb[[scol]],
             gene = unname(GENE_LAB[[g]]), stringsAsFactors = FALSE)
}))
gene_long <- gene_long[is.finite(gene_long$mean) & gene_long$n >= 30, ]
gene_long$gene <- factor(gene_long$gene, levels = unname(GENE_LAB))

pB <- ggplot(gene_long, aes(r_bin, mean, colour = gene, fill = gene)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.45) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = GENE_COL, name = NULL) +
  scale_fill_manual(values = GENE_COL, name = NULL) +
  labs(x = "Normalised radius (0 = follicle centroid, 1 = cluster-5 edge)",
       y = "Mean log-normalised expression") +
  theme_pub +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.box.spacing = unit(0.4, "mm"),
        legend.key.width = unit(4.5, "mm"),
        legend.key.height = unit(1.6, "mm"),
        plot.margin = margin(1.5, 2, 1, 2, "mm")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8, alpha = 1)),
         fill = "none")

# -----------------------------------------------------------------------------
# c-d  T-subset transfer maps and paired scores
# -----------------------------------------------------------------------------
cat("Panels c-d ...\n")
tf <- read.csv(file.path(SRC, "SourceData_transfer_scores.csv"),
               stringsAsFactors = FALSE)
tf$xp <- as.numeric(tf$x)
tf$yp <- -as.numeric(tf$y)
tf$louvain <- as.character(tf$louvain)
tcols <- c(
  CD4_1a = "T.CD4_1a", CD4_1b = "T.CD4_1b",
  CD8_0 = "T.CD8_0", CD8_3a = "T.CD8_3a", CD8_3b = "T.CD8_3b"
)
stopifnot(all(tcols %in% names(tf)))

# Source map order with Fig. 3 v4 titles; five maps in one row.
map_order <- c("CD4_1a", "CD4_1b", "CD8_0", "CD8_3a", "CD8_3b")
t_panels <- lapply(map_order, function(nm) {
  make_feat(tf$xp, tf$yp, tf[[tcols[[nm]]]], wrap_fun(unname(FIG3_LAB[[nm]])),
            legend_name = "Score", ramp = proj_ramp, size = 0.16)
})
pC <- wrap_plots(t_panels, ncol = 5)

paired_box <- function(high, low, lab_h, lab_l, accent) {
  ok <- is.finite(high) & is.finite(low)
  high <- high[ok]; low <- low[ok]
  wt <- wilcox.test(high, low, paired = TRUE, conf.int = TRUE, exact = FALSE)
  dd <- data.frame(
    score = c(low, high),
    group = factor(rep(c(lab_l, lab_h), each = length(low)),
                   levels = c(lab_l, lab_h))
  )
  fill_vals <- c(unname(FIG3_COL[lab_l]), accent)
  names(fill_vals) <- c(lab_l, lab_h)
  p_str <- if (!is.finite(wt$p.value) || wt$p.value < 1e-16) {
    "P < 1e-16"
  } else {
    sprintf("P = %.1e", wt$p.value)
  }
  ann <- sprintf("median \u0394 = %+.3f\n95%% CI %+.3f to %+.3f\n%s, n = %s",
                 unname(wt$estimate), wt$conf.int[1], wt$conf.int[2],
                 p_str, format(length(low), big.mark = ","))
  list(
    plot = ggplot(dd, aes(group, score, fill = group)) +
      geom_jitter(width = 0.17, height = 0, size = 0.12, alpha = 0.10,
                  stroke = 0, colour = "grey30", shape = 16) +
      geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.92,
                   linewidth = 0.3, colour = "black") +
      annotate("text", x = 1.5, y = max(dd$score) * 1.02, label = ann,
               size = PT(5.2), family = FAM, vjust = 0, lineheight = 1.05) +
      scale_x_discrete(labels = wrap_fun) +
      scale_fill_manual(values = fill_vals, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.30))) +
      labs(x = NULL, y = "Prediction score") +
      theme_pub +
      theme(axis.text.x = element_text(size = 5.5, lineheight = 0.95,
                                       margin = margin(t = 1.0, unit = "mm"))),
    n = length(low), estimate = unname(wt$estimate),
    ci_lo = wt$conf.int[1], ci_hi = wt$conf.int[2], p = wt$p.value
  )
}

tnic <- tf$louvain %in% c("1", "4", "5", "8")
lab_mem_cd4 <- unname(FIG3_LAB[["CD4_1a"]])
lab_naive_cd4 <- unname(FIG3_LAB[["CD4_1b"]])
lab_naive_cd8 <- unname(FIG3_LAB[["CD8_0"]])
box1 <- paired_box(tf[["T.CD4_1a"]][tnic], tf[["T.CD4_1b"]][tnic],
                   lab_mem_cd4, lab_naive_cd4, FIG3_COL[[lab_mem_cd4]])
box2 <- paired_box(tf[["T.CD8_0"]][tnic], tf[["T.CD4_1b"]][tnic],
                   lab_naive_cd8, lab_naive_cd4, FIG3_COL[[lab_naive_cd8]])
pD <- box1$plot | box2$plot
cat(sprintf("d left %s vs %s: n=%d median d=%.3f P=%.1e\n",
            lab_mem_cd4, lab_naive_cd4, box1$n, box1$estimate, box1$p))
cat(sprintf("d right %s vs %s: n=%d median d=%.3f P=%.1e\n",
            lab_naive_cd8, lab_naive_cd4, box2$n, box2$estimate, box2$p))

# -----------------------------------------------------------------------------
# Assemble
# -----------------------------------------------------------------------------
cat("Assemble ...\n")
row_a <- wrap_elements(full = pA)
row_c <- wrap_elements(full = pC)
row_d <- wrap_elements(full = pD)
fig <- tag_it(row_a, "a") / tag_it(pB, "b") / tag_it(row_c, "c") / tag_it(row_d, "d") +
  plot_layout(heights = c(1.22, 0.58, 0.72, 0.62)) &
  theme(plot.margin = margin(1.4, 1.4, 1.2, 2.2, "mm"))

w_in <- FIG_W_MM / 25.4
h_in <- FIG_H_MM / 25.4
base <- file.path(OUTDIR, "FigureS4")

cat("Writing SVG ...\n")
svglite::svglite(paste0(base, ".svg"), width = w_in, height = h_in, bg = "white")
print(fig); invisible(dev.off())
cat("Writing PDF ...\n")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in,
                     family = FAM, bg = "white")
print(fig); invisible(dev.off())
cat("Writing PNG ...\n")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
print(fig); invisible(dev.off())
cat("Writing TIFF ...\n")
ragg::agg_tiff(paste0(base, ".tiff"), width = w_in, height = h_in, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); invisible(dev.off())

write.csv(gene_long, file.path(SDP, "SourceData_FigS4b_radial.csv"), row.names = FALSE)
stats_d <- data.frame(
  panel = c("d_left", "d_right"),
  contrast = c("Memory/activated CD4 T vs Naive CD4 T",
               "Naive CD8 T vs Naive CD4 T"),
  n = c(box1$n, box2$n),
  median_paired_delta = c(box1$estimate, box2$estimate),
  ci_low = c(box1$ci_lo, box2$ci_lo),
  ci_high = c(box1$ci_hi, box2$ci_hi),
  p_value = c(box1$p, box2$p),
  stringsAsFactors = FALSE
)
write.csv(stats_d, file.path(SDP, "SourceData_FigS4d_paired.csv"), row.names = FALSE)
cat("DONE ->", OUTDIR, "\n")
