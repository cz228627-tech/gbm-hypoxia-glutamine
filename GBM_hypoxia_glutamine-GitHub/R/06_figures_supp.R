###############################################################################
##  06_figures_supp.R
##  Supplementary Figures S2-S7
##  (S1 is a cohort workflow schematic — create with a drawing tool)
##
##  Requires: results from all gate scripts
##
##  *** BUG-AUDIT WARNING (fixed in this revision, do not skip) ***
##  This script had NEVER been checked before this audit. Four separate
##  problems were found across Figs S2/S3/S5/S6 (see each section's own
##  inline comment for full detail); summary:
##    - Fig S2: "qa_sensitivity_scoring.csv" is not produced by any script in
##      this R port (would require re-implementing Scanpy's sc.tl.score_genes
##      natively in R -- out of scope) AND the original plotting code assumed
##      the wrong table shape (per-patient long format with method/rho
##      columns; the real file is a 14-row per-component-per-method SUMMARY).
##      Added a file-exists guard + rewrote the plot around the real schema.
##      Numerically re-verified end-to-end against the real frozen file after
##      the fix: all 14 rows parse cleanly into Component x method, 0 NAs.
##    - Fig S3: was reading the per-patient-per-state DETAIL file
##      ("qa_malignant_state.csv", no median_rho column at all after this
##      session's 01_gateA_lin.R Section 9 restructure) instead of the
##      aggregated "qa_malignant_state_summary.csv"; also referenced a
##      nonexistent "n_total" column instead of the real "n_patients".
##    - Fig S5: stale filename "qa_spatial_null.csv" -> renamed to
##      "qa_greenwald_mean_sensitivity.csv" (same propagation gap as 05).
##    - Fig S6: GENUINELY PRE-EXISTING bug, independent of any renaming in
##      this whole audit -- code used capitalized `GSC_line` (matching the
##      unrelated ARCHIVED gateB_per_sample_scores.csv's different Score_*
##      schema) but 03_gateB_gsc.R's own per_sample data.frame (and the CSV
##      it writes under the SAME filename during a fresh run) uses lowercase
##      `gsc_line` throughout. This would crash `data.table`'s `.()`
##      selection ("object 'GSC_line' not found") on ANY fresh end-to-end run
##      of the ORIGINAL (pre-audit) R port, regardless of every other fix in
##      this audit. Fixed by using `gsc_line` (lowercase) to match what
##      03_gateB_gsc.R actually produces.
##  Confirmed NOT broken (checked, no fix needed): Fig S4 (phase1_rho_matrix.
##  csv) and Fig S7 (gateC_u87_summary.csv) -- column names verified against
##  both the real frozen files and 04_gateC_u87.R's own construction code.
###############################################################################

source("00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_figure <- function(plot, name, width = 180, height = 120, units = "mm") {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = units, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), plot,
         width = width, height = height, units = units)
  cat(sprintf("Saved: %s\n", name))
}

###############################################################################
##  Fig. S2: Scoring sensitivity comparison
##  Compares arithmetic mean vs sc.tl.score_genes (with control gene subtraction)
###############################################################################

## *** BUG-AUDIT WARNING/FIX (this revision) ***
## Two separate problems were found here:
## 1. NONE of 00_setup.R/01-04 currently produce "qa_sensitivity_scoring.csv"
##    at all -- reproducing it would require re-implementing Scanpy's
##    sc.tl.score_genes() (expression-bin-stratified random control-gene
##    sampling) natively in R, which was judged out of scope for this port.
##    This figure can therefore only be reproduced if you separately obtain
##    the archived qa_sensitivity_scoring.csv (e.g. alongside the frozen
##    Python outputs) and place it in results/ yourself; a missing-file guard
##    is added below so the rest of the script still runs without it.
## 2. The plotting code assumed a per-patient long-format table with columns
##    "method"/"rho" (for a per-patient boxplot+jitter). The ACTUAL frozen
##    file is a pre-aggregated per-component-per-method SUMMARY table with
##    columns name/n_positive/n_negative/n_total/median_rho/
##    wilcoxon_p_twosided/n_genes_found/n_genes_total (14 rows: 5 components +
##    2 overall endpoints, x2 methods) -- neither "method" nor "rho" exist as
##    columns, so the original code would error immediately if the file were
##    ever supplied. Rewritten below to parse the "name" column (e.g.
##    "Fate_Synthesis__mean__hyp_mean") into Component + method and plot the
##    aggregated median_rho comparison that the real file actually supports.
if (file.exists(file.path(results_dir, "qa_sensitivity_scoring.csv"))) {
  sensitivity <- fread(file.path(results_dir, "qa_sensitivity_scoring.csv"))
  sensitivity[, Component := sub("^(Fate_|Overall_)(.*?)__(mean|score_genes)__.*$", "\\2", name)]
  sensitivity[, method := sub("^.*__(mean|score_genes)__.*$", "\\1", name)]
  sensitivity$Component <- factor(sensitivity$Component,
    levels = c(names(fate_components), "11", "13"))

  fig_s2 <- ggplot(sensitivity, aes(x = Component, y = median_rho, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("mean" = "#0279EE", "score_genes" = "#FF9400"),
                       labels = c("mean" = "Arithmetic mean", "score_genes" = "score_genes"),
                       name = NULL) +
    labs(x = NULL, y = "Median Spearman rho",
         title = "Scoring method sensitivity (component-level)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(fig_s2, "figS2_scoring_sensitivity")
} else {
  message("Skipping Fig. S2: qa_sensitivity_scoring.csv not found in results/ ",
          "(not produced by any script in this R port -- see comment above).")
}

###############################################################################
##  Fig. S3: Within-state analysis
##  Shows hypoxia-glutamine rho within each malignant state (AC, MES, OPC, NPC)
###############################################################################

## *** BUG-AUDIT FIX (this revision): this panel needs the AGGREGATED
## per-state summary (State, n_patients, n_positive, n_negative, median_rho,
## wilcoxon_p_twosided), not the per-patient-per-state detail table.
## 01_gateA_lin.R Section 9 now produces both as two separate files;
## "qa_malignant_state.csv" itself is the per-patient detail table (columns
## Patient/State/rho/direction/n_cells) since this session's restructure, so
## reading it here would silently produce an empty/wrong plot. Use
## "qa_malignant_state_summary.csv" instead. Also note the summary's patient
## -count column is "n_patients", not "n_total". ***
state_results <- fread(file.path(results_dir, "qa_malignant_state_summary.csv"))
state_results$State <- factor(state_results$State, levels = malignant_states)

fig_s3 <- ggplot(state_results, aes(x = State, y = median_rho, fill = State)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", n_positive, n_patients)),
            vjust = -0.5, size = 3) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "Median Spearman rho",
       title = "Within-state hypoxia-glutamine association") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_figure(fig_s3, "figS3_within_state", width = 120, height = 80)

###############################################################################
##  Fig. S4: Full 12-program landscape (detailed)
###############################################################################

landscape <- fread(file.path(results_dir, "phase1_rho_matrix.csv"))
landscape_mat <- as.matrix(landscape[, -1])
rownames(landscape_mat) <- landscape$patient

## Order by glutamine rho (descending)
glq_order <- order(landscape_mat[, "glutamine"], decreasing = TRUE)
landscape_mat <- landscape_mat[glq_order, ]

## Order programs by median rho
prog_order <- names(sort(colMedians(landscape_mat), decreasing = TRUE))
landscape_mat <- landscape_mat[, prog_order]

landscape_long <- reshape2::melt(landscape_mat)
colnames(landscape_long) <- c("patient", "program", "rho")

fig_s4 <- ggplot(landscape_long, aes(x = program, y = patient, fill = rho)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                       midpoint = 0, name = "rho") +
  labs(x = NULL, y = NULL, title = "12-program metabolic landscape") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6))

save_figure(fig_s4, "figS4_landscape", width = 160, height = 140)

###############################################################################
##  Fig. S5: 19-section Z-scores (spatial null calibration)
###############################################################################

## *** BUG-AUDIT FIX: renamed from "qa_spatial_null.csv" -- see
## 02_gateA_spatial.R header for why that name is stale/deprecated. ***
spatial_null <- fread(file.path(results_dir, "qa_greenwald_mean_sensitivity.csv"))
spatial_null <- spatial_null[order(z_score)]
spatial_null$Section <- factor(spatial_null$Section,
                                levels = spatial_null$Section)

fig_s5 <- ggplot(spatial_null, aes(x = Section, y = z_score, fill = z_score > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "#0279EE", "TRUE" = "#E9ED4C"),
                    guide = "none") +
  labs(x = "Section (ranked by Z-score)", y = "Calibrated Z-score",
       title = "Per-section spatial null calibration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

save_figure(fig_s5, "figS5_spatial_zscores", width = 180, height = 100)

###############################################################################
##  Fig. S6: Per-GSC paired plots (all 5 components)
###############################################################################

per_sample <- fread(file.path(results_dir, "gateB_per_sample_scores.csv"))

## *** BUG-AUDIT FIX (this revision): 03_gateB_gsc.R's own per_sample
## data.frame (and the CSV it writes) uses lowercase "gsc_line", not
## "GSC_line" -- confirmed by grepping every other internal use of this
## column within 03_gateB_gsc.R itself (delta/LOCO/contrast sections all use
## $gsc_line). The capitalized name only exists in the unrelated ARCHIVED
## frozen gateB_per_sample_scores.csv (which additionally uses a
## "Score_<Component>" prefix and capitalized "Condition" -- do not mix the
## two schemas). Without this fix, data.table's .() selection below throws
## "object 'GSC_line' not found" and the script stops -- this would happen
## on a completely fresh end-to-end run, independent of any other fix in
## this audit. ***
comp_plots <- lapply(names(fate_components), function(comp) {
  df <- per_sample[, .(gsc_line, condition, get(comp))]
  colnames(df)[3] <- "score"
  ggplot(df, aes(x = condition, y = score, group = gsc_line)) +
    geom_line(color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = condition), size = 1.5) +
    scale_color_manual(values = c("Normoxia" = "grey50", "Hypoxia" = "#D62728"),
                       guide = "none") +
    labs(x = NULL, y = "Score", title = comp) +
    theme(plot.title = element_text(size = 8))
})

fig_s6 <- wrap_plots(comp_plots, ncol = 3) +
  plot_annotation(title = "Per-GSC component responses")
save_figure(fig_s6, "figS6_per_gsc_components", width = 180, height = 120)

###############################################################################
##  Fig. S7: U87 complete results
###############################################################################

u87_summary <- fread(file.path(results_dir, "gateC_u87_summary.csv"))
u87_summary$Component <- factor(u87_summary$Component,
                                 levels = u87_summary$Component)

fig_s7 <- ggplot(u87_summary, aes(x = Component, y = mean_diff, fill = mean_diff > 0)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("p=%.3f", t_p_twosided)),
            vjust = ifelse(u87_summary$mean_diff > 0, -0.5, 1.5), size = 2.5) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "#0279EE", "TRUE" = "#E9ED4C"),
                    guide = "none") +
  labs(x = NULL, y = "Mean diff (Hypoxia - Normoxia)",
       title = "U87MG complete results") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_figure(fig_s7, "figS7_u87_complete", width = 140, height = 90)

cat("\nSupplementary figures S2-S7 saved.\n")
cat("Fig. S1 (cohort workflow) must be created manually.\n")
