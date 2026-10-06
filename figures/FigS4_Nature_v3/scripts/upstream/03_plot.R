#!/usr/bin/env Rscript
# =============================================================================
# 03  Figures for anatomy overlay, radial profiles, and (if present) transfers
# R exclusive. Writes into Fig4_spatial_anatomy_reprojection/ only.
# =============================================================================
options(warn = 1)
set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(systemfonts)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig4_spatial_anatomy_reprojection")
SRC <- file.path(OUT, "source_data")
SD  <- file.path(OUT, "source_data_plotted")
FIG <- file.path(OUT, "figures")
QA  <- file.path(OUT, "qa")
for (d in c(FIG, QA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT, "run_log_03.txt")
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("start 03 %s", format(Sys.time()))
# Projection colour refresh: skip anatomy plates unless this flag is TRUE.
REMAKE_ANATOMY <- FALSE

FAM <- "Arial"
BASE <- 6.5
fm <- tryCatch(systemfonts::match_fonts(FAM), error = function(e) NULL)
if (!is.null(fm) && nzchar(fm$path)) {
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
# Transfer maps: high-chroma Plasma (dark purple -> magenta -> orange -> yellow).
# Truncate the last ~8% so the high end stays yellow, not paper-white.
proj_ramp <- grDevices::hcl.colors(256, "Plasma")[1:235]
div_ramp  <- grDevices::hcl.colors(256, "Purple-Orange")
BG_COL <- "#EBEBEB"

theme_pub <- theme_classic(base_size = BASE, base_family = FAM) +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.0, "pt"),
    axis.title = element_text(size = BASE - 0.3, colour = "black"),
    axis.text = element_text(size = BASE - 0.8, colour = "black"),
    legend.title = element_text(size = BASE - 0.5),
    legend.text = element_text(size = BASE - 0.8),
    legend.key.size = unit(2.4, "mm"),
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5, colour = "grey30"),
    plot.tag = element_text(size = 8, face = "bold", family = FAM, hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    plot.margin = margin(1.5, 2, 1.5, 3),
    panel.grid = element_blank()
  )
theme_sp <- theme_void(base_family = FAM) +
  theme(
    plot.title = element_text(size = BASE, face = "bold", hjust = 0.5, family = FAM),
    plot.subtitle = element_text(size = BASE - 1.2, hjust = 0.5, colour = "grey30"),
    plot.tag = element_text(size = 8, face = "bold", family = FAM, hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    legend.title = element_text(size = BASE - 0.5, family = FAM),
    legend.text = element_text(size = BASE - 0.8, family = FAM),
    plot.margin = margin(1, 1, 1, 4)
  )

save_pub <- function(plot, filename, width_mm = 183, height_mm = 140, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = FAM)
  print(plot); dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot); dev.off()
  say("wrote %s", basename(filename))
}

df <- read.csv(file.path(SD, "SourceData_section_markers_radial.csv"),
               stringsAsFactors = FALSE)
df$xp <- as.numeric(df$x)
df$yp <- -as.numeric(df$y)
PAD <- 150
XLIM <- range(df$xp, na.rm = TRUE) + c(-PAD, PAD)
YLIM <- range(df$yp, na.rm = TRUE) + c(-PAD, PAD)

# lossless section outline (same idea as Fig. 4)
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

make_feat <- function(val, title, tag, scalebar = FALSE, legend_name = "log-normalised\nexpression") {
  lim <- lim01(val)
  d <- data.frame(xp = df$xp, yp = df$yp,
                  val = pmin(pmax(as.numeric(val), lim[1]), lim[2]))
  d <- d[order(d$val), , drop = FALSE]
  p <- ggplot() +
    section_layer() +
    geom_point(data = d, aes(xp, yp, colour = val),
               size = 0.22, stroke = 0, shape = 16) +
    scale_colour_gradientn(colours = score_ramp, limits = lim, oob = scales::squish,
                           name = legend_name) +
    guides(colour = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(18, "mm"), barheight = unit(1.5, "mm"),
      ticks.colour = "white", frame.colour = "black", frame.linewidth = 0.25
    )) +
    coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = title, tag = tag) +
    theme_sp +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.25))
  if (scalebar) {
    p <- p +
      annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
               linewidth = 0.6, colour = "black") +
      annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y + 0.028 * diff(YLIM),
               label = "1 mm", size = (BASE - 1.0) / .pt, family = FAM)
  }
  p
}

# ---- Figure markers ----
if (REMAKE_ANATOMY) {
genes <- c(
  Cr2    = "Cr2 (CD21)",
  Fcer2a = "Fcer2a (CD23)",
  Cd4    = "Cd4",
  Jchain = "Jchain",
  Marco  = "Marco",
  Adgre1 = "Adgre1 (F4/80)"
)
tags <- letters[1:6]
plist <- vector("list", 6)
for (i in seq_along(genes)) {
  g <- names(genes)[i]
  plist[[i]] <- make_feat(df[[g]], genes[[i]], tags[i], scalebar = (i == 1L))
}
fig_m <- wrap_plots(plist, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
save_pub(fig_m, file.path(FIG, "Figure_markers"), 183, 145)

# ---- radial summaries ----
rb <- read.csv(file.path(SD, "SourceData_radial_rnorm.csv"), stringsAsFactors = FALSE)

gene_long <- bind_rows(lapply(names(genes), function(g) {
  mcol <- paste0(g, "_mean"); scol <- paste0(g, "_sem")
  if (!mcol %in% names(rb)) return(NULL)
  data.frame(r_bin = rb$r_bin, n = rb$n,
             mean = rb[[mcol]], sem = rb[[scol]],
             gene = genes[[g]], stringsAsFactors = FALSE)
}))
gene_long <- gene_long[is.finite(gene_long$mean) & gene_long$n >= 30, ]

p_g <- ggplot(gene_long, aes(r_bin, mean, colour = gene, fill = gene)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.45) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c(
    "Cr2 (CD21)" = "#0072B2", "Fcer2a (CD23)" = "#E69F00",
    "Cd4" = "#CC79A7", "Jchain" = "#009E73",
    "Marco" = "#D55E00", "Adgre1 (F4/80)" = "#000000"
  ), name = NULL) +
  scale_fill_manual(values = c(
    "Cr2 (CD21)" = "#0072B2", "Fcer2a (CD23)" = "#E69F00",
    "Cd4" = "#CC79A7", "Jchain" = "#009E73",
    "Marco" = "#D55E00", "Adgre1 (F4/80)" = "#000000"
  ), name = NULL) +
  labs(x = "Normalised radius (0 = follicle centroid, 1 = cluster-5 edge)",
       y = "Mean log-normalised expression", tag = "a",
       title = "Anatomical markers") +
  theme_pub +
  theme(legend.position = "right")

# B IV delta radial — only bins that carry Fig. 4 scores
if ("B_IV_delta_mean" %in% names(rb)) {
  biv <- rb[is.finite(rb$B_IV_delta_mean) & rb$n >= 20, ]
  p_b <- ggplot(biv, aes(r_bin, B_IV_delta_mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey40") +
    geom_ribbon(aes(ymin = B_IV_delta_mean - B_IV_delta_sem,
                    ymax = B_IV_delta_mean + B_IV_delta_sem),
                fill = "#0072B2", colour = NA, alpha = 0.18) +
    geom_line(colour = "#0072B2", linewidth = 0.5) +
    labs(x = "Normalised radius",
         y = "Mean B IV-high \u2212 IV-low\nprediction score",
         tag = "b", title = "Fig. 4 B reference, same radius") +
    theme_pub
} else {
  p_b <- ggplot() + theme_void() + labs(tag = "b")
}

# distance-to-cluster-5 alternative (less sensitive to small components)
df$dist_c5_um <- as.numeric(df$dist_c5_um)
alt <- df %>%
  filter(is.finite(dist_c5_um), dist_c5_um <= 400) %>%
  mutate(d_bin = floor(dist_c5_um / 25) * 25) %>%
  group_by(d_bin) %>%
  summarise(
    n = n(),
    Cr2 = mean(Cr2, na.rm = TRUE),
    Fcer2a = mean(Fcer2a, na.rm = TRUE),
    Marco = mean(Marco, na.rm = TRUE),
    Adgre1 = mean(Adgre1, na.rm = TRUE),
    Cd4 = mean(Cd4, na.rm = TRUE),
    Jchain = mean(Jchain, na.rm = TRUE),
    B_IV_delta = mean(B_IV_delta, na.rm = TRUE),
    B_IV_n = sum(is.finite(B_IV_delta)),
    .groups = "drop"
  )
write.csv(alt, file.path(SD, "SourceData_radial_dist_c5.csv"), row.names = FALSE)
alt_long <- alt %>%
  filter(n >= 30) %>%
  pivot_longer(c(Cr2, Fcer2a, Marco, Adgre1, Cd4, Jchain),
               names_to = "gene", values_to = "mean")
alt_long$gene <- recode(alt_long$gene,
  Cr2 = "Cr2 (CD21)", Fcer2a = "Fcer2a (CD23)", Marco = "Marco",
  Adgre1 = "Adgre1 (F4/80)", Cd4 = "Cd4", Jchain = "Jchain")
p_c <- ggplot(alt_long, aes(d_bin, mean, colour = gene)) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = c(
    "Cr2 (CD21)" = "#0072B2", "Fcer2a (CD23)" = "#E69F00",
    "Cd4" = "#CC79A7", "Jchain" = "#009E73",
    "Marco" = "#D55E00", "Adgre1 (F4/80)" = "#000000"
  ), name = NULL) +
  labs(x = "Distance to nearest cluster-5 bin (\u03bcm)",
       y = "Mean log-normalised expression", tag = "c",
       title = "Distance to cluster 5 (edge-based)") +
  theme_pub +
  theme(legend.position = "right")

p_d <- ggplot(alt %>% filter(B_IV_n >= 20), aes(d_bin, B_IV_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_line(colour = "#0072B2", linewidth = 0.5) +
  labs(x = "Distance to nearest cluster-5 bin (\u03bcm)",
       y = "Mean B IV-high \u2212 IV-low", tag = "d",
       title = "Fig. 4 B delta vs cluster-5 distance") +
  theme_pub

# transfer radial if present
tr_path <- file.path(SD, "SourceData_radial_transfer.csv")
has_tr <- file.exists(file.path(SD, "SourceData_transfer_scores.csv"))
if (file.exists(tr_path)) {
  tr <- read.csv(tr_path, stringsAsFactors = FALSE)
  tr <- tr[is.finite(tr$B_MZ_minus_Fo_mean) & tr$n >= 30, ]
  p_e <- ggplot(tr, aes(r_bin, B_MZ_minus_Fo_mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey40") +
    geom_ribbon(aes(ymin = B_MZ_minus_Fo_mean - B_MZ_minus_Fo_sem,
                    ymax = B_MZ_minus_Fo_mean + B_MZ_minus_Fo_sem),
                fill = "#009E73", colour = NA, alpha = 0.18) +
    geom_line(colour = "#009E73", linewidth = 0.5) +
    labs(x = "Normalised radius", y = "Mean MZ \u2212 Fo prediction score",
         tag = "e", title = "Fo vs MZ transfer") +
    theme_pub
} else {
  p_e <- ggplot() + theme_void() +
    labs(tag = "e", title = "Fo vs MZ transfer (pending script 02)") +
    theme_pub
}

fig_r <- (p_g | p_b) / (p_c | p_d) / p_e + plot_layout(heights = c(1, 1, 0.9))
save_pub(fig_r, file.path(FIG, "Figure_radial"), 183, 195)
} else {
  say("skip anatomy plates (REMAKE_ANATOMY = FALSE)")
}

# ---- transfer figures ----
tf_path <- file.path(SD, "SourceData_transfer_scores.csv")
if (file.exists(tf_path)) {
  tf <- read.csv(tf_path, stringsAsFactors = FALSE)
  tf$xp <- as.numeric(tf$x); tf$yp <- -as.numeric(tf$y)
  # align order to df for section layer already built from df
  pick_col <- function(d, patterns) {
    for (p in patterns) {
      hits <- grep(p, names(d), value = TRUE, ignore.case = TRUE)
      if (length(hits)) return(hits[1])
    }
    NA_character_
  }
  col_fo <- pick_col(tf, c("^B\\.Fo$", "^mod\\.Fo$", "B\\.Fo", "mod\\.Fo"))
  col_mz <- pick_col(tf, c("^B\\.MZ$", "^mod\\.MZ$", "B\\.MZ", "mod\\.MZ"))
  say("B Fo col=%s  MZ col=%s", col_fo, col_mz)

  make_score <- function(val, title, tag, lnm = "Prediction\nscore") {
    make_feat(val, title, tag, legend_name = lnm)
  }
  # rebuild make_feat using tf coords if barcodes match df
  make_tf <- function(val, title, tag, lnm = "Prediction\nscore",
                      ramp = proj_ramp, diverging = FALSE) {
    v <- as.numeric(val)
    if (diverging) {
      m <- max(abs(range(v[is.finite(v)])), 1e-6)
      lim <- c(-m, m)
      cols <- div_ramp
    } else {
      lim <- lim01(v)
      cols <- ramp
    }
    d <- data.frame(xp = tf$xp, yp = tf$yp,
                    val = pmin(pmax(v, lim[1]), lim[2]))
    d <- d[is.finite(d$xp) & is.finite(d$yp), , drop = FALSE]
    d <- d[order(d$val), , drop = FALSE]
    ggplot() +
      section_layer() +
      geom_point(data = d, aes(xp, yp, colour = val),
                 size = 0.28, stroke = 0, shape = 16) +
      scale_colour_gradientn(colours = cols, limits = lim, oob = scales::squish,
                             name = lnm) +
      guides(colour = guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = unit(18, "mm"), barheight = unit(1.5, "mm"),
        ticks.colour = "black", frame.colour = "black", frame.linewidth = 0.25
      )) +
      coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
      labs(title = title, tag = tag) +
      theme_sp +
      theme(legend.position = "bottom",
            panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.25))
  }

  paired_box <- function(high, low, lab_h, lab_l, title, tag, accent) {
    ok <- is.finite(high) & is.finite(low)
    high <- high[ok]; low <- low[ok]
    wt <- wilcox.test(high, low, paired = TRUE, conf.int = TRUE, exact = FALSE)
    dd <- data.frame(
      score = c(low, high),
      group = factor(rep(c(lab_l, lab_h), each = length(low)), levels = c(lab_l, lab_h))
    )
    ann <- sprintf("median \u0394 = %+.3f\n95%% CI %+.3f to %+.3f\nP = %.1e, n = %s",
                   unname(wt$estimate), wt$conf.int[1], wt$conf.int[2],
                   wt$p.value, format(length(low), big.mark = ","))
    ggplot(dd, aes(group, score, fill = group)) +
      geom_jitter(width = 0.17, height = 0, size = 0.14, alpha = 0.12,
                  stroke = 0, colour = "grey30", shape = 16) +
      geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.92,
                   linewidth = 0.3, colour = "black") +
      annotate("text", x = 1.5, y = max(dd$score) * 1.02, label = ann,
               size = (BASE - 1.3) / .pt, family = FAM, vjust = 0, lineheight = 1.05) +
      scale_fill_manual(values = c("#BFBFBF", accent), guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.28))) +
      labs(x = NULL, y = "Prediction score", title = title, tag = tag) +
      theme_pub
  }

  if (!is.na(col_fo) && !is.na(col_mz)) {
    delta <- as.numeric(tf[[col_mz]]) - as.numeric(tf[[col_fo]])
    pA <- make_tf(tf[[col_fo]], "Fo-like (cluster 0)", "a")
    pB <- make_tf(tf[[col_mz]], "MZ-like (cluster 1)", "b")
    pC <- make_tf(delta, "MZ minus Fo", "c", lnm = "MZ \u2212 Fo", diverging = TRUE)
    bnic <- tf$louvain %in% c("1", "5")
    pD <- paired_box(tf[[col_mz]][bnic], tf[[col_fo]][bnic],
                     "MZ", "Fo", "Louvain 1+5 bins", "d", "#009E73")
    fig_b <- (pA | pB | pC | pD) + plot_layout(widths = c(1, 1, 1, 0.85))
    save_pub(fig_b, file.path(FIG, "Figure_B_MZvsFo"), 183, 78)
  }

  tcols <- list(
    CD4_1a = pick_col(tf, c("T\\.CD4_1a", "mod\\.CD4_1a", "CD4_1a")),
    CD4_1b = pick_col(tf, c("T\\.CD4_1b", "mod\\.CD4_1b", "CD4_1b")),
    CD8_0  = pick_col(tf, c("T\\.CD8_0", "mod\\.CD8_0", "CD8_0")),
    CD8_3a = pick_col(tf, c("T\\.CD8_3a", "mod\\.CD8_3a", "CD8_3a")),
    CD8_3b = pick_col(tf, c("T\\.CD8_3b", "mod\\.CD8_3b", "CD8_3b"))
  )
  say("T columns: %s", paste(sprintf("%s=%s", names(tcols), unlist(tcols)), collapse = "; "))
  t_panels <- list()
  labs_t <- c("true CD4 1a (Maf-high)", "true CD4 1b (Tcf7-high)",
              "TNK CD8 cluster 0", "TNK CD8 cluster 3a", "TNK CD8 cluster 3b")
  for (i in seq_along(tcols)) {
    cn <- tcols[[i]]
    if (is.na(cn)) next
    t_panels[[length(t_panels) + 1]] <- make_tf(tf[[cn]], labs_t[i], letters[i])
  }
  tnic <- tf$louvain %in% c("1", "4", "5", "8")
  box_list <- list()
  if (!is.na(tcols$CD4_1a) && !is.na(tcols$CD4_1b)) {
    box_list[[1]] <- paired_box(tf[[tcols$CD4_1a]][tnic], tf[[tcols$CD4_1b]][tnic],
                                "CD4 1a", "CD4 1b", "B/T niches",
                                letters[length(t_panels) + 1], "#CC79A7")
  }
  if (!is.na(tcols$CD8_0) && !is.na(tcols$CD4_1b)) {
    box_list[[length(box_list) + 1]] <- paired_box(
      tf[[tcols$CD8_0]][tnic], tf[[tcols$CD4_1b]][tnic],
      "CD8 0", "CD4 1b", "B/T niches",
      letters[length(t_panels) + length(box_list) + 1], "#0072B2")
  }
  # lineage IV if present
  col_cd4p <- pick_col(tf, c("trueCD4_IV\\+", "trueCD4_IV\\."))
  col_cd4n <- pick_col(tf, c("trueCD4_IV-", "trueCD4_IV\\.$", "Tiv\\.trueCD4_IV"))
  # more robust: grep both
  iv_cd4_hi <- grep("trueCD4_IV\\+|trueCD4_IV\\.$", names(tf), value = TRUE)
  iv_cd4_lo <- grep("trueCD4_IV-|trueCD4_IV-$", names(tf), value = TRUE)
  say("CD4 IV cols hi=%s lo=%s", paste(iv_cd4_hi, collapse=","), paste(iv_cd4_lo, collapse=","))
  if (length(iv_cd4_hi) && length(iv_cd4_lo)) {
    box_list[[length(box_list) + 1]] <- paired_box(
      tf[[iv_cd4_hi[1]]][tnic], tf[[iv_cd4_lo[1]]][tnic],
      "CD4 IV+", "CD4 IV-", "true CD4, B/T niches",
      letters[length(t_panels) + length(box_list) + 1], "#CC79A7")
  }
  iv_cd8_hi <- grep("TNKcd8_IV\\+|TNKcd8_IV\\.$", names(tf), value = TRUE)
  iv_cd8_lo <- grep("TNKcd8_IV-", names(tf), value = TRUE)
  if (length(iv_cd8_hi) && length(iv_cd8_lo)) {
    box_list[[length(box_list) + 1]] <- paired_box(
      tf[[iv_cd8_hi[1]]][tnic], tf[[iv_cd8_lo[1]]][tnic],
      "CD8 IV+", "CD8 IV-", "TNK CD8, B/T niches",
      letters[length(t_panels) + length(box_list) + 1], "#0072B2")
  }

  if (length(t_panels)) {
    fig_t <- wrap_plots(c(t_panels, box_list), ncol = 3)
    save_pub(fig_t, file.path(FIG, "Figure_T_CD4CD8"), 183, 155)
  }
} else {
  say("transfer scores not yet present; skipped B/T transfer figures")
}

# legends
writeLines(c(
  "Figure_markers | Anatomical genes on the spleen Stereo-seq section used in Fig. 4.",
  "a–f, Log-normalised expression of Cr2, Fcer2a, Cd4, Jchain, Marco and Adgre1 (F4/80).",
  "All six maps share one coordinate frame. Colour limits are the 1st–99th percentile of each gene.",
  "Scale bar, 1 mm (0.5 um per chip unit, same constant as Fig. 4)."
), file.path(OUT, "Figure_markers_legend.txt"))

writeLines(c(
  "Figure_radial | Follicle-centred and cluster-5-edge profiles.",
  "a, Mean expression versus normalised radius (0 = centroid of Louvain-5 components with >=40 bins; 1 = equivalent edge). Ribbons, SEM.",
  "b, Fig. 4 B-cell CD45-IV-high minus IV-low score on the same radius (only bins that carry that score).",
  "c, Same genes versus Euclidean distance to the nearest Louvain-5 bin.",
  "d, Fig. 4 B IV delta versus that edge distance.",
  "e, MZ-like minus Fo-like TransferData scores versus normalised radius (script 02)."
), file.path(OUT, "Figure_radial_legend.txt"))

writeLines(c(
  "Figure_B_MZvsFo | Fo-like versus MZ-like label transfer.",
  "Reference: Fig. 2 B-cell clusters 0 (Fo-like, Fcer2a-high) and 1 (MZ-like, Cr2-high).",
  "a, Fo-like prediction score. b, MZ-like prediction score. Sequential Plasma scale (purple to yellow).",
  "c, MZ minus Fo on a symmetric Purple-Orange diverging scale.",
  "d, Paired scores in Louvain 1+5 bins; two-sided Wilcoxon signed-rank; median paired difference with 95% CI."
), file.path(OUT, "Figure_B_MZvsFo_legend.txt"))

writeLines(c(
  "Figure_T_CD4CD8 | Hard-gated true CD4 (1a/1b) and TNK-gated CD8 (0/3a/3b) label transfer.",
  "Gates match the hard-gated CD4 / TNK-gated CD8 volcano package (cluster 2 dropped; CD8 uses Cd8b1).",
  "Maps show prediction scores on the same Stereo-seq section using the Plasma sequential scale (purple to yellow).",
  "Boxplots are paired Wilcoxon signed-rank tests in Louvain 1+4+5+8 bins."
), file.path(OUT, "Figure_T_CD4CD8_legend.txt"))

say("done 03 %s", format(Sys.time()))
