#!/usr/bin/env Rscript
# Fig.7 / S5 / S6 figures after mapping corrected L4–L6. New folder only.
options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
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

OUT <- file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "newfigures/Fig7_S5_S6_spatial_reanalysis_lesion_corrected")
SRC <- file.path(OUT, "source_data")
FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

FAM <- "Arial"
BASE <- 6.5
theme_set(
  theme_classic(base_size = BASE, base_family = FAM) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.title = element_text(size = BASE, colour = "black"),
      axis.text = element_text(size = 5.6, colour = "black"),
      legend.title = element_text(size = 5.6),
      legend.text = element_text(size = 5.2),
      legend.key.size = unit(2.4, "mm"),
      plot.title = element_text(size = BASE, face = "bold", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold", family = FAM),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid = element_blank()
    )
)
set_cols <- c(old = "#6C3483", added = "#C0392B", pooled = "#1A7A4C",
              added_vs_strict_extra = "#E67E22")
zone_cols <- c(Meninges = "#4C78A8", Lesion = "#C0392B", Plasma = "#8E44AD",
               Parenchyma = "#7F8C8D", Outside = "#9E9E9E")
C_OLD <- "#6C3483"; C_NEW <- "#C0392B"

save_pub <- function(plot, filename, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  base <- file.path(FIG, filename)
  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
  print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = FAM)
  print(plot); grDevices::dev.off()
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in",
                res = 600, background = "white")
  print(plot); grDevices::dev.off()
  cat("wrote", filename, "\n")
}

dat <- readRDS(file.path(SRC, "spot_table_reanalysis.rds"))
rects <- read.csv(file.path(SRC, "lesion_rects_expanded.csv"), stringsAsFactors = FALSE)
item1 <- read.csv(file.path(SRC, "item1_lesion_vs_outside_T_Mg_CXCL.csv"), stringsAsFactors = FALSE)
item4 <- read.csv(file.path(SRC, "item4_immune_restricted_wilcoxon.csv"), stringsAsFactors = FALSE)
item5w <- read.csv(file.path(SRC, "item5_MgIV_vs_MgIVC_within_windows.csv"), stringsAsFactors = FALSE)
item5o <- read.csv(file.path(SRC, "item5_Mg_windows_vs_outside.csv"), stringsAsFactors = FALSE)
item8 <- read.csv(file.path(SRC, "item8_distance_spearman.csv"), stringsAsFactors = FALSE)
item8b <- read.csv(file.path(SRC, "item8_distance_bins.csv"), stringsAsFactors = FALSE)
item9z <- read.csv(file.path(SRC, "item9_Cxcl10pos_winner_by_zone.csv"), stringsAsFactors = FALSE)
item3 <- read.csv(file.path(SRC, "item3_S5h_signed_effects.csv"), stringsAsFactors = FALSE)
clust <- read.csv(file.path(SRC, "item3_cluster7_rename.csv"), stringsAsFactors = FALSE)
detp <- read.csv(file.path(SRC, "item10_insitu_detection_by_zone.csv"), stringsAsFactors = FALSE)
fullsamp <- readRDS(file.path(SRC, "fullsection_gene_plot_sample.rds"))
fullsum <- read.csv(file.path(SRC, "item2_fullsection_zone_summary.csv"), stringsAsFactors = FALSE)
winov <- read.csv(file.path(SRC, "S6d_vs_Fig7_window_overlap.csv"), stringsAsFactors = FALSE)
he <- readRDS(file.path(Sys.getenv("IV_RAW_DATA", unset = "F:/IVstaining_scRNA"), "IVstaining/Fig5_Nature_vector/source_data/HE_Brain_lowres.rds"))
if (max(he, na.rm = TRUE) > 1.5) he <- he / 255
nr <- dim(he)[1]; nc <- dim(he)[2]
rast <- grDevices::as.raster(he)

# ---- Figure A: map + identity flip ----
pMap <- ggplot(dat, aes(x_img, -y_img)) +
  annotation_raster(rast, xmin = 0.5, xmax = nc + 0.5, ymin = -(nr + 0.5), ymax = -0.5) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -ymax, ymax = -ymin,
                              colour = set),
            fill = NA, linewidth = 0.35, inherit.aes = FALSE) +
  scale_colour_manual(values = c(original = C_OLD, added = C_NEW), name = NULL) +
  coord_fixed(xlim = c(0.5, nc + 0.5), ylim = c(-(nr + 0.5), -0.5), expand = FALSE) +
  labs(title = "Original L1–L3 (purple) and corrected L4–L6 (red)", tag = "a") +
  theme_void(base_family = FAM) +
  theme(plot.title = element_text(size = BASE, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 8, face = "bold"),
        legend.position = "bottom", legend.text = element_text(size = 5.2))

key <- c("Tcell EAE_NEG", "Tcell_POS_bias", "Mg EAE_POS", "expr_Cxcl10",
         "Mg_clean_POS_cos")
b1 <- item1 %>%
  filter(score %in% key, lesion_set %in% c("old", "added", "pooled",
                                          "added_vs_strict_extra")) %>%
  mutate(
    score = recode(score,
                   `Tcell EAE_NEG` = "T NEG",
                   Tcell_POS_bias = "T POS-NEG bias",
                   `Mg EAE_POS` = "Mg POS",
                   expr_Cxcl10 = "Cxcl10",
                   Mg_clean_POS_cos = "Clean Mg POS cosine"),
    lesion_set = recode(lesion_set,
                        old = "Original L1–L3",
                        added = "Added vs complement",
                        pooled = "Pooled six",
                        added_vs_strict_extra = "Added vs strict extra"),
    lesion_set = factor(lesion_set, levels = c(
      "Original L1–L3", "Added vs complement", "Added vs strict extra", "Pooled six"))
  )
pBar <- ggplot(b1, aes(score, r_rb, fill = lesion_set)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(position = position_dodge(0.78), width = 0.7, colour = "black",
           linewidth = 0.15) +
  scale_fill_manual(values = c(
    "Original L1–L3" = "#6C3483",
    "Added vs complement" = "#C0392B",
    "Added vs strict extra" = "#E67E22",
    "Pooled six" = "#1A7A4C"
  ), name = NULL) +
  labs(x = NULL, y = "r (lesion vs outside)", tag = "b",
       title = "Fig. 7j cuff claims reverse in the corrected red boxes") +
  theme(legend.position = "top", axis.text.x = element_text(size = 5.3),
        legend.text = element_text(size = 5.0))

set.seed(1)
dat_v <- dat[sample.int(nrow(dat), min(nrow(dat), 8000L)), ]
v1 <- dat_v %>%
  mutate(scope = case_when(
    in_old ~ "Original",
    in_new ~ "Added",
    TRUE ~ "Extra"
  )) %>%
  select(scope, `Tcell EAE_NEG`, Tcell_POS_bias, `Mg EAE_POS`) %>%
  pivot_longer(-scope, names_to = "score", values_to = "value") %>%
  mutate(
    score = recode(score, `Tcell EAE_NEG` = "T NEG",
                   Tcell_POS_bias = "T POS-NEG bias", `Mg EAE_POS` = "Mg POS"),
    scope = factor(scope, levels = c("Extra", "Original", "Added"))
  )
pV <- ggplot(v1, aes(scope, value, fill = scope)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(Extra = "#9E9E9E", Original = "#6C3483", Added = "#C0392B"),
                    guide = "none") +
  facet_wrap(~ score, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Prediction score", tag = "c",
       title = "Same spots: original cuffs raise T NEG and drop Mg POS; added windows do the opposite") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 5.8))

figA <- (pMap | pBar) / pV + plot_layout(heights = c(1.15, 1.0))
save_pub(figA, "Figure_A_window_identity", 183, 170)

# ---- Figure B: per-window Mg + Hexb full + S5 leftover ----
o_mg <- item5o %>%
  filter(score == "MgIV", window %in% c("L1", "L2", "L3", "L4", "L5", "L6",
                                        "Original", "Added")) %>%
  mutate(window = factor(window, levels = c("L1", "L2", "L3", "Original",
                                            "L4", "L5", "L6", "Added")))
pF <- ggplot(o_mg, aes(window, r_rb, fill = window %in% c("L4", "L5", "L6", "Added"))) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%+.2f", r_rb),
                vjust = ifelse(r_rb < 0, 1.35, -0.35)),
            size = 5.0 / ggplot2::.pt, family = FAM) +
  scale_fill_manual(values = c("FALSE" = "#6C3483", "TRUE" = "#C0392B"), guide = "none") +
  labs(x = NULL, y = "r (window vs strict extra, MgIV)", tag = "a",
       title = "L1/L2 deplete MgIV; L4–L6 enrich MgIV versus extra-window meninges")

wuse <- item5w %>%
  filter(window %in% c("L1", "L2", "L3", "L4", "L5", "L6", "Outside")) %>%
  mutate(window = factor(window, levels = c("L1", "L2", "L3", "L4", "L5", "L6", "Outside")))
pE <- ggplot(wuse, aes(window, r_rb)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 2.0, colour = "#C0392B") +
  geom_text(aes(label = sprintf("%.2f", r_rb)), vjust = -0.7,
            size = 5.0 / ggplot2::.pt, family = FAM) +
  labs(x = NULL, y = "r (MgIV > MgIVC within window)", tag = "b",
       title = "Within-window MgIV > MgIVC is not lesion-specific (also true outside)") +
  coord_cartesian(ylim = c(0, 0.75))

fullsamp$zone_full <- factor(fullsamp$zone_full,
                             levels = c("Meninges", "Lesion", "Parenchyma", "Plasma"))
fullsum$zone_full <- factor(fullsum$zone_full,
                            levels = c("Meninges", "Lesion", "Parenchyma", "Plasma"))
fullsum$lab <- sprintf("mean=%.3f\n%.1f%%+", fullsum$Hexb_mean, fullsum$Hexb_pct)
pC <- ggplot(fullsamp, aes(zone_full, Hexb, fill = zone_full)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.8) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  geom_text(data = fullsum, aes(zone_full, Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.1, size = 5.0 / ggplot2::.pt, family = FAM) +
  scale_fill_manual(values = zone_cols, guide = "none") +
  labs(x = NULL, y = "Hexb", tag = "c",
       title = "Full section, pooled six windows: Hexb is highest in the expanded lesion set, not lowest") +
  theme(axis.text.x = element_text(size = 5.2))

item3$score <- factor(item3$score, levels = rev(item3$score[order(item3$r_rb)]))
pH <- ggplot(item3, aes(r_rb, score, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("POS > NEG" = "#E67E22", "NEG > POS" = "#8E44AD"),
                    name = NULL) +
  labs(x = "r (POS vs NEG scores)", y = NULL, tag = "d",
       title = "S5h is independent of windows: lineage-specific IV bias") +
  theme(legend.position = "top")

figB <- (pF | pE) / (pC | pH) + plot_layout(heights = c(1, 1.1))
save_pub(figB, "Figure_B_windows_Hexb_S5", 183, 160)

# ---- Figure C: Cxcl10 source + distance ----
item9z$zone <- factor(item9z$zone, levels = c("Meninges", "Lesion", "Plasma"))
item9z$lesion_set <- recode(item9z$lesion_set,
                            old = "Original L1–L3", added = "Added L4–L6",
                            pooled = "Pooled six")
item9z$lesion_set <- factor(item9z$lesion_set,
                            levels = c("Original L1–L3", "Added L4–L6", "Pooled six"))
pJ <- ggplot(item9z, aes(zone, Freq, fill = winner)) +
  geom_col(position = "fill", width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(MacIVC = "#E67E22", MgIV = "#1ABC9C", MgIVC = "#8E44AD"),
                    name = "Highest score") +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ lesion_set, nrow = 1) +
  labs(x = NULL, y = "Fraction of Cxcl10+ spots", tag = "a",
       title = "Cxcl10+ winner: MacIVC majority only in original cuffs; added windows are MgIV-led") +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 5.6))

bT <- item8b %>%
  filter(score %in% c("Tcell EAE_NEG", "Hexb", "expr_Cxcl10")) %>%
  mutate(
    score = recode(score, `Tcell EAE_NEG` = "T NEG", expr_Cxcl10 = "Cxcl10"),
    distance = recode(distance,
                      dist_lesion_old = "To original L1–L3",
                      dist_lesion_added = "To added L4–L6",
                      dist_lesion = "To pooled six"),
    distance = factor(distance, levels = c("To original L1–L3", "To added L4–L6",
                                           "To pooled six"))
  )
pL <- ggplot(bT, aes(dist_mid, mean, colour = score)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  facet_wrap(~ distance, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("T NEG" = "#C0392B", Hexb = "#4C78A8", Cxcl10 = "#E67E22"),
                      name = NULL) +
  labs(x = "Distance to window (image pixels)", y = "Mean score", tag = "b",
       title = "Near original cuffs: T NEG up, Hexb down. Near added windows the T NEG gradient reverses") +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(size = 5.5))

det_show <- detp %>%
  filter(lesion_set == "pooled",
         gene %in% c("Hexb", "P2ry12", "Tmem119", "Cd3e", "Ptprc", "Gfap")) %>%
  mutate(zone = factor(zone, levels = c("Meninges", "Lesion", "Plasma")))
pM <- ggplot(det_show, aes(zone, pct_pos, fill = zone)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = zone_cols, guide = "none") +
  facet_wrap(~ gene, nrow = 1) +
  labs(x = NULL, y = "% spots detected", tag = "c",
       title = "Pooled six windows: Hexb/P2ry12/Tmem119 no longer drop in the lesion class") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 5.6),
        axis.text.x = element_text(size = 5.0))

figC <- pJ / pL / pM + plot_layout(heights = c(1.05, 1.05, 0.95))
save_pub(figC, "Figure_C_source_distance_insitu", 183, 195)

# ---- Figure D: gates + clean Mg ----
gate <- item4 %>%
  filter(gate == "imm_p75",
         score %in% c("Tcell EAE_NEG", "Mg EAE_POS", "expr_Cxcl10"),
         lesion_set %in% c("old", "added", "pooled")) %>%
  mutate(
    lab = recode(score, `Tcell EAE_NEG` = "T NEG", `Mg EAE_POS` = "Mg POS",
                 expr_Cxcl10 = "Cxcl10"),
    lesion_set = recode(lesion_set, old = "Original", added = "Added",
                        pooled = "Pooled"),
    lesion_set = factor(lesion_set, levels = c("Original", "Added", "Pooled"))
  )
pt_path <- file.path(SRC, "item4_Ptprc_positive_wilcoxon.csv")
pO <- ggplot(gate, aes(lab, r_rb, fill = lesion_set)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(position = position_dodge(0.75), width = 0.68, colour = "black",
           linewidth = 0.2) +
  scale_fill_manual(values = c(Original = "#6C3483", Added = "#C0392B",
                               Pooled = "#1A7A4C"), name = NULL) +
  labs(x = NULL, y = "r (lesion vs outside), imm-sum P75", tag = "a",
       title = "Immune-high spots: original T NEG up / Mg POS down; added windows reverse both") +
  theme(legend.position = "top")

if (file.exists(pt_path)) {
  pt <- read.csv(pt_path, stringsAsFactors = FALSE)
  pt <- pt %>%
    filter(score %in% c("Tcell EAE_NEG", "Mg EAE_POS", "expr_Cxcl10")) %>%
    mutate(
      lab = recode(score, `Tcell EAE_NEG` = "T NEG", `Mg EAE_POS` = "Mg POS",
                   expr_Cxcl10 = "Cxcl10"),
      lesion_set = recode(lesion_set, old = "Original", added = "Added",
                          pooled = "Pooled"),
      lesion_set = factor(lesion_set, levels = c("Original", "Added", "Pooled"))
    )
  pP <- ggplot(pt, aes(lab, r_rb, fill = lesion_set)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_col(position = position_dodge(0.75), width = 0.68, colour = "black",
             linewidth = 0.2) +
    scale_fill_manual(values = c(Original = "#6C3483", Added = "#C0392B",
                                 Pooled = "#1A7A4C"), name = NULL) +
    labs(x = NULL, y = "r (lesion vs outside), Ptprc>0", tag = "b",
         title = "Ptprc>0 gate: original Cxcl10 up in cuffs; added windows keep a smaller Cxcl10 excess") +
    theme(legend.position = "top")
} else {
  pP <- ggplot() + labs(tag = "b")
}

if (all(c("Mg_clean_POS_cos", "Mg_clean_NEG_cos") %in% names(dat))) {
  vc <- dat_v %>%
    mutate(scope = case_when(in_old ~ "Original", in_new ~ "Added", TRUE ~ "Extra")) %>%
    select(scope, Mg_clean_POS_cos, Mg_clean_NEG_cos) %>%
    pivot_longer(-scope, names_to = "score", values_to = "value") %>%
    mutate(
      score = recode(score, Mg_clean_POS_cos = "Clean Mg POS",
                     Mg_clean_NEG_cos = "Clean Mg NEG"),
      scope = factor(scope, levels = c("Extra", "Original", "Added"))
    )
  pN <- ggplot(vc, aes(scope, value, fill = scope)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.85) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(Extra = "#9E9E9E", Original = "#6C3483",
                                 Added = "#C0392B"), guide = "none") +
    facet_wrap(~ score, nrow = 1) +
    labs(x = NULL, y = "Cosine to clean Mg mean", tag = "c",
         title = "Decontaminated Mg cosine: depleted in original cuffs, enriched in added windows") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 5.6))
} else {
  pN <- ggplot() + labs(tag = "c")
}

figD <- pO / pP / pN + plot_layout(heights = c(1, 1, 1.05))
save_pub(figD, "Figure_D_gates_cleanMg", 183, 190)

writeLines(c(
  "Fig. R1 | Fig. 7 cuff-window claims after adding the corrected L4–L6 red boxes.",
  "a, H&E with original L1–L3 (purple) and user-corrected L4–L6 (red).",
  "b, Rank-biserial effects for lesion versus outside under four contrasts. Added vs complement includes original cuffs in the outside group; added vs strict extra excludes all six windows from outside.",
  "c, Prediction scores split by original windows, added windows, and extra-window meninges/plasma.",
  "Source data are provided as a Source Data file."
), file.path(OUT, "Figure_A_legend.txt"), useBytes = TRUE)

writeLines(c(
  "Fig. R2 | Per-window microglia scores, full-section Hexb, and window-independent S5h.",
  "a, MgIV in each window versus tissue outside all six windows.",
  "b, Within-window MgIV versus MgIVC.",
  "c, Hexb on the full Visium section with pooled six windows as the lesion class.",
  "d, POS versus NEG prediction scores by lineage (independent of lesion windows).",
  "Source data are provided as a Source Data file."
), file.path(OUT, "Figure_B_legend.txt"), useBytes = TRUE)

writeLines(c(
  "Fig. R3 | Cxcl10+ winners, distance gradients, and in-situ detection.",
  "a, Highest of MacIVC / MgIV / MgIVC among Cxcl10-positive spots, with lesion defined as original, added, or pooled windows.",
  "b, Mean T NEG, Hexb and Cxcl10 versus distance to original, added, or pooled windows.",
  "c, Detection rates on the immune slice using pooled six windows as the lesion class.",
  "Source data are provided as a Source Data file."
), file.path(OUT, "Figure_C_legend.txt"), useBytes = TRUE)

writeLines(c(
  "Fig. R4 | Immune occupancy gates and clean-microglia cosine.",
  "a, Imm-sum P75 gate, lesion versus outside, for original, added and pooled window sets.",
  "b, Ptprc>0 gate, same contrasts.",
  "c, Cosine to decontaminated microglia means in original windows, added windows, and extra tissue.",
  "Source data are provided as a Source Data file."
), file.path(OUT, "Figure_D_legend.txt"), useBytes = TRUE)

cat("PLOTS DONE\n")
