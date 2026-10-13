###############################################################################
##  05_figures_main.R
##  Main Figures 1-4 for the manuscript
##  (Figure 5 is a conceptual schematic — create with a drawing tool)
##
##  Requires: results from 01_gateA_lin.R, 02_gateA_spatial.R,
##            03_gateB_gsc.R, 04_gateC_u87.R
##
##  *** BUG-AUDIT WARNING (fixed in this revision, do not skip) ***
##  This script had NEVER been checked against 01/02's Section renames before
##  this audit and was broken by them (would have crashed on ANY fresh
##  end-to-end run using the post-rename 01/02 scripts). Five breakages found
##  and fixed, all pure cross-script dependency drift, no new science:
##    1. Fig 1 data load: "gateA_per_patient_rho.csv" -> renamed to
##       "gateA_primary_per_patient_rho.csv" (01_gateA_lin.R Section 5 rename
##       from an earlier session was never propagated here).
##    2. Fig 2 data load: "gateA_contrast.csv" -> "qa_lin_contrast_mean.csv"
##       (broken by this session's 01_gateA_lin.R Section 8 rename; see that
##       script's header -- the OLD name is a 3rd instance of the stale
##       score_genes-vs-mean file confusion, 16/28 p=0.095 vs the correct
##       24/28 p=1.1e-5 quoted in the manuscript).
##    3. Fig 2C plotting code: `contrast$patient` (lowercase, x2 occurrences)
##       -> `contrast$Patient` (capital) -- broken by this session's column
##       rename in 01_gateA_lin.R Section 8.
##    4. Fig 3 Panel C data load (+ a commented-out duplicate reference):
##       "qa_spatial_null.csv" -> "qa_greenwald_mean_sensitivity.csv" (broken
##       by 02_gateA_spatial.R's rename from an earlier session, never
##       propagated here).
##  Confirmed NOT broken (checked, no fix needed): gateA_component_summary.csv
##  and gateA_per_gene_rho.csv reads (filenames unchanged; content is freshly
##  correct for an end-to-end R-only run -- staleness only affects comparison
##  against the OLD ARCHIVED Python files, see 01_gateA_lin.R header);
##  phase1_rho_matrix.csv; gateB_per_line_delta_11gene.csv/gateB_summary.csv/
##  qa_loco_gsc.csv. gateB_contrast.csv and gateA_spatial_contrast.csv are
##  intentionally not read by any panel here (no dedicated Fig 3/4 panel
##  exists for the GSC/spatial contrast -- consistent, not a bug).
###############################################################################

source("00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
})

## Font settings
theme_set(theme_classic(base_size = 10, base_family = "sans"))
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: save figure as PNG + SVG ---------------------------------------
save_figure <- function(plot, name, width = 180, height = 120, units = "mm") {
  png_path <- file.path(fig_dir, paste0(name, ".png"))
  svg_path <- file.path(fig_dir, paste0(name, ".svg"))
  ggsave(png_path, plot, width = width, height = height, units = units, dpi = 300)
  ggsave(svg_path, plot, width = width, height = height, units = units)
  cat(sprintf("Saved: %s, %s\n", png_path, svg_path))
}

###############################################################################
##  Figure 1: Discovery
##  A: Per-patient Spearman rho (waterfall)
##  B: 12-program landscape (heatmap)
###############################################################################

## Load data
## *** BUG-AUDIT FIX: 01_gateA_lin.R Section 5 now writes
## "gateA_primary_per_patient_rho.csv" (renamed to avoid colliding with the
## stale, differently-shaped archived gateA_per_patient_rho.csv -- see that
## script's header). Updated here to match. ***
patient_rho <- fread(file.path(results_dir, "gateA_primary_per_patient_rho.csv"))
landscape   <- fread(file.path(results_dir, "phase1_rho_matrix.csv"))

## Panel A: Waterfall plot of per-patient rho
patient_rho_sorted <- patient_rho[order(rho), ]
patient_rho_sorted$patient <- factor(patient_rho_sorted$patient,
                                      levels = patient_rho_sorted$patient)

fig1a <- ggplot(patient_rho_sorted, aes(x = patient, y = rho, fill = rho > 0)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "#0279EE", "TRUE" = "#E9ED4C"),
                    guide = "none") +
  labs(x = "Patient (ranked)", y = "Spearman rho",
       title = "Hypoxia vs Glutamine (11-gene)") +
  annotate("text", x = 2, y = max(patient_rho$rho) * 0.9,
           label = sprintf("%d/%d positive\nP < 0.0001",
                           sum(patient_rho$rho > 0), nrow(patient_rho)),
           hjust = 0, size = 3) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Panel B: 12-program landscape heatmap
landscape_mat <- as.matrix(landscape[, -1])
rownames(landscape_mat) <- landscape$patient

## Order programs by median rho
program_order <- names(sort(colMedians(landscape_mat), decreasing = TRUE))
landscape_mat <- landscape_mat[, program_order]

## Convert to long format for ggplot
landscape_long <- reshape2::melt(landscape_mat)
colnames(landscape_long) <- c("patient", "program", "rho")

fig1b <- ggplot(landscape_long, aes(x = program, y = patient, fill = rho)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                       midpoint = 0, name = "rho") +
  labs(x = NULL, y = NULL, title = "Metabolic program landscape") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6))

figure1 <- fig1a / fig1b + plot_layout(heights = c(1, 1.2))
save_figure(figure1, "figure1_discovery", width = 180, height = 160)

###############################################################################
##  Figure 2: Decomposition
##  A: 5-component summary (dot plot with significance)
##  B: Per-gene rho (heatmap)
##  C: Synthesis vs Catabolism contrast (paired dot plot)
###############################################################################

## *** BUG-AUDIT FIX: 01_gateA_lin.R Section 8's contrast output is now
## "qa_lin_contrast_mean.csv" (renamed -- the old name "gateA_contrast.csv"
## belongs to a stale, score_genes-based archived file with a qualitatively
## different result: 16/28 synthesis-favored, p=0.095, vs the correct 24/28,
## p=1.1e-5 that matches the manuscript). comp_summary/per_gene keep their
## original names since those ARE freshly overwritten by 01_gateA_lin.R on
## every run (no rename needed there; see that script's Section 6/7 comments
## for why the *archived* copies of these two filenames must not be used for
## validation instead). ***
comp_summary <- fread(file.path(results_dir, "gateA_component_summary.csv"))
per_gene     <- fread(file.path(results_dir, "gateA_per_gene_rho.csv"))
contrast     <- fread(file.path(results_dir, "qa_lin_contrast_mean.csv"))

## Panel A: Component dot plot
comp_summary$Component <- factor(comp_summary$Component,
                                  levels = names(fate_components))
comp_summary$significant <- comp_summary$q_BH < 0.05

fig2a <- ggplot(comp_summary, aes(x = Component, y = median_rho, fill = Component)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = 0, ymax = median_rho), width = 0.2) +
  geom_text(aes(label = sprintf("%d/%d", n_positive, n_total)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = component_colors, guide = "none") +
  labs(x = NULL, y = "Median Spearman rho",
       title = "Fate component associations") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## Panel B: Per-gene heatmap
per_gene$Gene <- factor(per_gene$Gene, levels = per_gene$Gene)
fig2b <- ggplot(per_gene, aes(x = 1, y = Gene, fill = median_rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%d/%d", n_positive, n_total)), size = 3) +
  scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                       midpoint = 0, name = "median rho") +
  labs(x = NULL, y = NULL, title = "Per-gene associations") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Panel C: Synthesis vs Catabolism contrast
## *** BUG-AUDIT FIX: contrast_df's patient identifier column is "Patient"
## (capital P, matching the frozen qa_lin_contrast_mean.csv schema), not
## "patient" -- updated both references below accordingly. ***
contrast_long <- rbind(
  data.frame(patient = contrast$Patient, component = "Synthesis",
             rho = contrast$rho_Synthesis),
  data.frame(patient = contrast$Patient, component = "Catabolism",
             rho = contrast$rho_Catabolism)
)
contrast_long$component <- factor(contrast_long$component,
                                   levels = c("Catabolism", "Synthesis"))

fig2c <- ggplot(contrast_long, aes(x = component, y = rho, group = patient)) +
  geom_line(color = "grey60", linewidth = 0.3) +
  geom_point(aes(color = component), size = 1.5) +
  scale_color_manual(values = c("Synthesis" = "#1F77B4",
                                 "Catabolism" = "#FF7F0E"),
                     guide = "none") +
  labs(x = NULL, y = "Spearman rho with hypoxia",
       title = "Synthesis vs Catabolism") +
  annotate("text", x = 1.5, y = max(contrast_long$rho) * 0.95,
           label = sprintf("%d/%d synthesis-favored\nP < 0.0001",
                           sum(contrast$D > 0), nrow(contrast)),
           hjust = 0.5, size = 3)

figure2 <- (fig2a + fig2b) / fig2c + plot_layout(heights = c(1, 0.8))
save_figure(figure2, "figure2_decomposition", width = 180, height = 160)

###############################################################################
##  Figure 3: Spatial
##  A: Representative section (zh916inf) hypoxia spatial map
##  B: Same section glutamine spatial map
##  C: Per-patient raw rho
##  D: Per-patient calibrated Z-scores
###############################################################################

## For panels A-B, you need to re-load the representative section
## and plot spot-level scores on spatial coordinates.
## This requires the Seurat object from 02_gateA_spatial.R.
## Below is the code structure; uncomment and adjust paths.

# section_df <- fread(file.path(results_dir, "qa_greenwald_mean_sensitivity.csv"))
#
# ## Load representative section (zh916inf)
# counts <- Read10X_h5("data/greenwald_visium/GSM7596591_zh916inf_filtered_feature_bc_matrix.h5")
# positions <- fread("data/greenwald_visium/GSM7596591_zh916inf_tissue_positions_list.csv.gz",
#                     header = FALSE)
# colnames(positions) <- c("barcode","in_tissue","array_row","array_col",
#                          "pixel_row","pixel_col")
# positions <- positions[in_tissue == 1]
# common <- intersect(colnames(counts), positions$barcode)
# counts <- counts[, common]
# positions <- positions[match(common, barcode), ]
# seu <- CreateSeuratObject(counts) |> NormalizeData()
# expr_mat <- t(GetAssayData(seu, slot = "data"))
# hyp <- score_module(expr_mat, hypoxia_genes, colnames(expr_mat))
# glq <- score_module(expr_mat, glutamine_11, colnames(expr_mat))
# spot_df <- data.frame(x = positions$pixel_row, y = positions$pixel_col,
#                       hypoxia = hyp, glutamine = glq)
#
# fig3a <- ggplot(spot_df, aes(x = x, y = y, color = hypoxia)) +
#   geom_point(size = 0.3) +
#   scale_color_gradient(low = "white", high = "#D62728") +
#   labs(title = "Hypoxia score", x = NULL, y = NULL) +
#   coord_fixed() + theme_void()
#
# fig3b <- ggplot(spot_df, aes(x = x, y = y, color = glutamine)) +
#   geom_point(size = 0.3) +
#   scale_color_gradient(low = "white", high = "#1F77B4") +
#   labs(title = "Glutamine score", x = NULL, y = NULL) +
#   coord_fixed() + theme_void()

## Panel C: Per-patient raw rho
## *** BUG-AUDIT FIX: 02_gateA_spatial.R now writes "qa_greenwald_mean_
## sensitivity.csv" (renamed from "qa_spatial_null.csv" -- the old name
## belongs to a stale, score_genes-based archived file that correlates only
## weakly, Spearman rho=0.39/p=0.09, with the correct arithmetic-mean file
## across the same 19 sections; see that script's header for details). ***
section_df <- fread(file.path(results_dir, "qa_greenwald_mean_sensitivity.csv"))
patient_raw <- section_df[, .(raw_rho = median(observed_rho), n_sec = .N), by = Patient]
patient_raw <- patient_raw[order(raw_rho)]
patient_raw$Patient <- factor(patient_raw$Patient, levels = patient_raw$Patient)

fig3c <- ggplot(patient_raw, aes(x = Patient, y = raw_rho, fill = raw_rho > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "#0279EE", "TRUE" = "#E9ED4C"),
                    guide = "none") +
  labs(x = "Patient (ranked)", y = "Median Spearman rho",
       title = "Raw spatial association") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

## Panel D: Per-patient calibrated Z-scores
patient_z <- section_df[, .(z = median(z_score), n_sec = .N), by = Patient]
patient_z <- patient_z[order(z)]
patient_z$Patient <- factor(patient_z$Patient, levels = patient_z$Patient)

fig3d <- ggplot(patient_z, aes(x = Patient, y = z, fill = z > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "#0279EE", "TRUE" = "#E9ED4C"),
                    guide = "none") +
  labs(x = "Patient (ranked)", y = "Calibrated Z-score",
       title = "Spatial-null-calibrated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

## Combine (uncomment fig3a/fig3b when spatial maps are ready)
# figure3 <- (fig3a + fig3b) / (fig3c + fig3d) + plot_layout(heights = c(1, 1))
figure3 <- fig3c / fig3d
save_figure(figure3, "figure3_spatial", width = 180, height = 120)

###############################################################################
##  Figure 4: Perturbation (GSC)
##  A: Per-line delta (11-gene glutamine, paired dot plot)
##  B: 5-component deltas
##  C: LOCO analysis
###############################################################################

delta_11  <- fread(file.path(results_dir, "gateB_per_line_delta_11gene.csv"))
comp_b    <- fread(file.path(results_dir, "gateB_summary.csv"))
loco      <- fread(file.path(results_dir, "qa_loco_gsc.csv"))

## Panel A: Per-line delta paired dot plot
delta_long <- data.frame(
  GSC = rep(delta_11$GSC_line, 2),
  condition = rep(c("Normoxia", "Hypoxia"), each = nrow(delta_11)),
  score = c(delta_11$Normoxia_11, delta_11$Hypoxia_11)
)

fig4a <- ggplot(delta_long, aes(x = condition, y = score, group = GSC)) +
  geom_line(color = "grey60", linewidth = 0.3) +
  geom_point(aes(color = condition), size = 2) +
  scale_color_manual(values = c("Normoxia" = "grey50", "Hypoxia" = "#D62728"),
                     guide = "none") +
  labs(x = NULL, y = "Glutamine score (11-gene)",
       title = "GSC acute hypoxia response") +
  annotate("text", x = 1.5, y = max(delta_long$score) * 1.1,
           label = sprintf("7/8 increased\nP = 0.016"), hjust = 0.5, size = 3)

## Panel B: Component deltas
comp_b$Component <- factor(comp_b$Component, levels = names(fate_components))
comp_b$significant <- comp_b$q_BH < 0.05

fig4b <- ggplot(comp_b, aes(x = Component, y = median_delta, fill = Component)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", n_positive, n_total)),
            vjust = ifelse(comp_b$median_delta > 0, -0.5, 1.5), size = 3) +
  scale_fill_manual(values = component_colors, guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(x = NULL, y = "Median delta (Hypoxia - Normoxia)",
       title = "Fate component responses") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## Panel C: LOCO analysis
loco$score_set <- factor(loco$score_set,
                          levels = loco$score_set)
loco$label <- ifelse(loco$survives == TRUE, "*", "")

fig4c <- ggplot(loco, aes(x = score_set, y = median_delta, fill = survives)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = 0, size = 6) +
  scale_fill_manual(values = c("TRUE" = "#E9ED4C", "FALSE" = "grey70"),
                    labels = c("TRUE" = "P < 0.05", "FALSE" = "NS"),
                    name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(x = NULL, y = "Median delta",
       title = "Leave-one-component-out") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

figure4 <- (fig4a + fig4b) / fig4c + plot_layout(heights = c(1, 0.8))
save_figure(figure4, "figure4_perturbation", width = 180, height = 160)

###############################################################################
##  Figure 5: Context Model (SCHEMATIC)
##  This is a conceptual figure. Create with a drawing tool (e.g., Illustrator,
##  Inkscape, or BioRender) — NOT with R code.
##
##  Content:
##  - Left panel: Steady-state tumor (scRNA + Visium)
##    Hypoxia -> synthesis-favored glutamine program (GLUL up, C/N null)
##  - Right panel: Acute perturbation (GSC, 24h)
##    Hypoxia -> C/N-redistribution-dominated response (C/N up, synthesis unchanged)
##  - Bottom: U87 established cell line — not reproduced
##  - Arrows: solid = significant, dashed = trend, double-headed = null
###############################################################################

cat("\nFigures 1-4 saved. Figure 5 (schematic) must be created manually.\n")
cat("All figures saved to figures/ directory.\n")
