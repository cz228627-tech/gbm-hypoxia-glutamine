###############################################################################
##  03_gateB_gsc.R
##  Gate B: GSE243848 — patient-derived GSC lines under acute hypoxia
##
##  Dataset: GSE243848 (GEO)
##  File:    deseq2_normalized.txt.gz (DESeq2 normalized counts)
##  Design:  8 GSC lines x {DMSO, RSL3} x {Normoxia, Hypoxia} = 32 samples
##           We use ONLY the DMSO control arm: 8 lines x 2 conditions = 16 samples
##  Columns: G44_ctr (DMSO normoxia), G44_ctr_H (DMSO hypoxia), etc.
##
##  Analysis:
##    1. Load DESeq2 normalized counts, log2(x+1) transform
##    2. Score hypoxia, glutamine (11-gene + 13-gene), 5 components per sample
##    3. Per-line delta = hypoxia score - normoxia score
##    4. Wilcoxon signed-rank test on deltas
##    5. Component-level analysis with BH correction
##    6. Leave-one-component-out (LOCO) analysis
###############################################################################

source("00_setup.R")
suppressPackageStartupMessages({
  library(data.table)
})

# ---- 1. Load Data -----------------------------------------------------------
## Adjust path to your local copy
counts_file <- "data/gse243848/deseq2_normalized.txt.gz"

cat("Loading GSE243848 DESeq2 normalized counts...\n")
raw <- fread(counts_file)

## Column structure: ensembl, gene_Ensembl (gene symbol), gene_biotype, then 32 sample columns
## Gene symbols are in column 2 ("gene_Ensembl")
gene_symbols <- raw[[2]]
sample_cols <- colnames(raw)[4:ncol(raw)]  # skip first 3 metadata columns
expr_mat <- as.matrix(raw[, ..sample_cols])
rownames(expr_mat) <- gene_symbols

cat(sprintf("Loaded: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# ---- 2. Filter to DMSO Control Arm ------------------------------------------
## DMSO normoxia columns: *_ctr (no _H suffix)
## DMSO hypoxia columns:  *_ctr_H
dmso_normoxia <- sample_cols[grep("_ctr$", sample_cols)]
dmso_hypoxia  <- sample_cols[grep("_ctr_H$", sample_cols)]

cat(sprintf("DMSO normoxia samples: %s\n", paste(dmso_normoxia, collapse = ", ")))
cat(sprintf("DMSO hypoxia samples:  %s\n", paste(dmso_hypoxia, collapse = ", ")))

## Extract GSC line names
gsc_lines <- sub("_ctr$", "", dmso_normoxia)
cat(sprintf("GSC lines: %s\n", paste(gsc_lines, collapse = ", ")))

## Subset to DMSO arm
expr_dmso <- expr_mat[, c(dmso_normoxia, dmso_hypoxia)]

# ---- 3. Log2(x+1) Transform -------------------------------------------------
expr_log <- log2(expr_dmso + 1)

## *** BUG-AUDIT FIX (CONFIRMED): missing per-gene z-scoring step ***
## The locked deliverable (glutamine_fate_gate_report.md) states the Gate B
## method explicitly: "DESeq2 normalized counts, log-transformed, z-scored
## across 16 samples." This R port previously log2-transformed and went
## STRAIGHT to score_module()'s arithmetic mean on log2 values -- z-scoring
## was missing entirely, meaning genes with larger log2-scale variance would
## silently dominate every module average, and module scores were on the
## wrong (log2, not z) scale.
## Recipe implemented: per-gene mean-center, divide by sample SD (ddof=1,
## i.e. R's default `sd()`/`scale()`), across these 16 DMSO-arm samples.
## VALIDATION (CONFIRMED, upgraded from an earlier inconclusive check): this
## recipe, combined with the signed_rank_test() zero-handling fix in
## 00_setup.R, reproduces the frozen gateB_per_line_delta_11gene.csv,
## gateB_per_line_delta.csv, gateB_summary.csv, and qa_loco_gsc.csv (for all
## overlapping rows) to full float64 precision (max abs diff ~1e-15) -- see
## analysis notebook for the reproduction run. An earlier session's numerical
## check on this z-scoring recipe was inconclusive because it was compared
## against the wrong column (13-gene Delta_OverallGlutamine values misread as
## 11-gene); the actual 11-gene column matches exactly, and so does the
## 13-gene column when checked separately. GLS2 has exact-zero DESeq2-
## normalized counts in 4/16 samples but this does not affect the fit
## quality once the correct comparison target is used.
expr_z <- t(scale(t(expr_log), center = TRUE, scale = TRUE))

# ---- 4. Score Modules -------------------------------------------------------
cat("Scoring modules...\n")

## Transpose to samples x genes for score_module (now on z-scored values)
expr_t <- t(expr_z)
gene_names <- colnames(expr_t)

hyp_scores <- score_module(expr_t, hypoxia_genes, gene_names)
glq_11_scores <- score_module(expr_t, glutamine_11, gene_names)
glq_13_scores <- score_module(expr_t, glutamine_13, gene_names)

## 5 fate components
comp_scores <- sapply(names(fate_components), function(comp) {
  score_module(expr_t, fate_components[[comp]], gene_names)
})
colnames(comp_scores) <- names(fate_components)

## Assemble per-sample scores
per_sample <- data.frame(
  sample = rownames(expr_t),
  condition = rep(c("Normoxia", "Hypoxia"), each = length(gsc_lines)),
  gsc_line = rep(gsc_lines, 2),
  hypoxia_score = hyp_scores,
  glutamine_11 = glq_11_scores,
  glutamine_13 = glq_13_scores,
  comp_scores,
  stringsAsFactors = FALSE
)
fwrite(per_sample, file.path(results_dir, "gateB_per_sample_scores.csv"))

# ---- 5. Per-Line Delta (Hypoxia - Normoxia) ---------------------------------
cat("\n=== Per-Line Delta ===\n")

normoxia_rows <- per_sample[per_sample$condition == "Normoxia", ]
hypoxia_rows  <- per_sample[per_sample$condition == "Hypoxia", ]

## Match by GSC line
delta_df <- data.frame(
  GSC_line = gsc_lines,
  Delta_OverallGlutamine_11 = NA_real_,
  Delta_OverallGlutamine_13 = NA_real_,
  Normoxia_11 = NA_real_,
  Hypoxia_11  = NA_real_,
  Normoxia_13 = NA_real_,
  Hypoxia_13  = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(gsc_lines)) {
  line <- gsc_lines[i]
  n_row <- normoxia_rows[normoxia_rows$gsc_line == line, ]
  h_row <- hypoxia_rows[hypoxia_rows$gsc_line == line, ]

  delta_df$Delta_OverallGlutamine_11[i] <- h_row$glutamine_11 - n_row$glutamine_11
  delta_df$Delta_OverallGlutamine_13[i] <- h_row$glutamine_13 - n_row$glutamine_13
  delta_df$Normoxia_11[i] <- n_row$glutamine_11
  delta_df$Hypoxia_11[i]  <- h_row$glutamine_11
  delta_df$Normoxia_13[i] <- n_row$glutamine_13
  delta_df$Hypoxia_13[i]  <- h_row$glutamine_13
}

fwrite(delta_df, file.path(results_dir, "gateB_per_line_delta_11gene.csv"))

# ---- 6. Wilcoxon Signed-Rank Tests ------------------------------------------
cat("\n=== Primary Endpoint (11-gene) ===\n")
test_11 <- signed_rank_test(delta_df$Delta_OverallGlutamine_11)
cat(sprintf("  Positive: %d/%d, median delta = %.4f, p = %.4f\n",
            test_11$n_positive, test_11$n_total, test_11$median, test_11$p_value))

cat("\n=== 13-gene version ===\n")
test_13 <- signed_rank_test(delta_df$Delta_OverallGlutamine_13)
cat(sprintf("  Positive: %d/%d, median delta = %.4f, p = %.4f\n",
            test_13$n_positive, test_13$n_total, test_13$median, test_13$p_value))

# ---- 7. Component-Level Analysis --------------------------------------------
cat("\n=== Component-Level Analysis ===\n")

component_delta <- data.frame(
  GSC_line = gsc_lines,
  matrix(NA_real_, nrow = length(gsc_lines), ncol = length(fate_components)),
  stringsAsFactors = FALSE
)
colnames(component_delta) <- c("GSC_line", names(fate_components))

for (comp in names(fate_components)) {
  for (i in seq_along(gsc_lines)) {
    line <- gsc_lines[i]
    n_val <- normoxia_rows[normoxia_rows$gsc_line == line, comp]
    h_val <- hypoxia_rows[hypoxia_rows$gsc_line == line, comp]
    component_delta[i, comp] <- h_val - n_val
  }
}

fwrite(component_delta, file.path(results_dir, "gateB_per_line_delta.csv"))

## Component summary with BH correction
comp_summary <- data.frame(
  Component = names(fate_components),
  n_positive = NA_integer_,
  n_negative = NA_integer_,
  n_total = NA_integer_,
  median_delta = NA_real_,
  wilcoxon_p = NA_real_,
  q_BH = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(fate_components)) {
  comp <- names(fate_components)[i]
  wt <- signed_rank_test(component_delta[[comp]])
  comp_summary$n_positive[i] <- wt$n_positive
  comp_summary$n_negative[i] <- wt$n_negative
  comp_summary$n_total[i]   <- wt$n_total
  comp_summary$median_delta[i] <- wt$median
  comp_summary$wilcoxon_p[i] <- wt$p_value
  cat(sprintf("  %s: %d/%d positive, delta = %.4f, p = %.4f\n",
              comp, wt$n_positive, wt$n_total, wt$median, wt$p_value))
}
comp_summary$q_BH <- bh_correct(comp_summary$wilcoxon_p)
cat(sprintf("  BH q-values: %s\n",
            paste(sprintf("%s=%.4f", comp_summary$Component, comp_summary$q_BH),
                  collapse = ", ")))

fwrite(comp_summary, file.path(results_dir, "gateB_summary.csv"))

# ---- 8. Leave-One-Component-Out (LOCO) Analysis -----------------------------
cat("\n=== LOCO Analysis ===\n")

## *** BUG-AUDIT NOTE: added a "Minus_Transport" row to exactly reproduce the
## frozen qa_loco_gsc.csv, which used a COARSER LOCO grouping than the 5 fate
## components for this particular sensitivity check: it tested
## CN_redistribution, Synthesis, Catabolism, and a combined "Transport" group
## -- confirmed (exact float64 match) to be SLC1A5+SLC7A5+SLC7A11, i.e. the 3
## Direct_uptake/Exchange_redox genes that are actually part of the 11-gene
## PRIMARY panel (SLC38A1/SLC38A2 are 13-gene-only, never in scope for this
## 11-gene LOCO). The frozen file did NOT test a separate Exchange_redox-only
## LOCO row at all. This script keeps all 5 fate-component LOCO rows (more
## granular; the 4 that overlap the frozen file -- Full_11, Minus_Synthesis,
## Minus_Catabolism, Minus_CN_redistribution -- reproduce it to full float64
## precision) AND adds "Minus_Transport" so the frozen file's exact scope is
## also fully reproducible here.
transport_genes <- c("SLC1A5", "SLC7A5", "SLC7A11")

## For each component, remove its genes from the 11-gene module and re-score
loco_results <- data.frame(
  score_set = c("Full_11", paste0("Minus_", names(fate_components)), "Minus_Transport"),
  n_genes = NA_integer_,
  n_positive = NA_integer_,
  n_negative = NA_integer_,
  n_total = NA_integer_,
  median_delta = NA_real_,
  wilcoxon_p = NA_real_,
  survives = NA,
  stringsAsFactors = FALSE
)

## Full 11-gene
loco_results[1, "n_genes"] <- length(glutamine_11)
loco_results[1, "n_positive"] <- test_11$n_positive
loco_results[1, "n_negative"] <- test_11$n_negative
loco_results[1, "n_total"] <- test_11$n_total
loco_results[1, "median_delta"] <- test_11$median
loco_results[1, "wilcoxon_p"] <- test_11$p_value
loco_results[1, "survives"] <- test_11$p_value < 0.05

## Minus each component
for (i in seq_along(fate_components)) {
  comp <- names(fate_components)[i]
  reduced_genes <- setdiff(glutamine_11, fate_components[[comp]])
  reduced_scores <- score_module(expr_t, reduced_genes, gene_names)

  ## Per-line delta
  red_normoxia <- reduced_scores[per_sample$condition == "Normoxia"]
  red_hypoxia  <- reduced_scores[per_sample$condition == "Hypoxia"]
  red_delta <- red_hypoxia - red_normoxia

  wt <- signed_rank_test(red_delta)
  idx <- i + 1
  loco_results$n_genes[idx]    <- length(reduced_genes)
  loco_results$n_positive[idx] <- wt$n_positive
  loco_results$n_negative[idx] <- wt$n_negative
  loco_results$n_total[idx]   <- wt$n_total
  loco_results$median_delta[idx] <- wt$median
  loco_results$wilcoxon_p[idx] <- wt$p_value
  loco_results$survives[idx]   <- wt$p_value < 0.05

  cat(sprintf("  Minus %s (%d genes): %d/%d positive, p = %.4f, survives = %s\n",
              comp, length(reduced_genes), wt$n_positive, wt$n_total,
              wt$p_value, wt$p_value < 0.05))
}

## "Transport" LOCO row (matches frozen qa_loco_gsc.csv scope -- see note above)
reduced_genes <- setdiff(glutamine_11, transport_genes)
reduced_scores <- score_module(expr_t, reduced_genes, gene_names)
red_normoxia <- reduced_scores[per_sample$condition == "Normoxia"]
red_hypoxia  <- reduced_scores[per_sample$condition == "Hypoxia"]
red_delta <- red_hypoxia - red_normoxia
wt <- signed_rank_test(red_delta)
idx <- match("Minus_Transport", loco_results$score_set)
loco_results$n_genes[idx]     <- length(reduced_genes)
loco_results$n_positive[idx]  <- wt$n_positive
loco_results$n_negative[idx]  <- wt$n_negative
loco_results$n_total[idx]    <- wt$n_total
loco_results$median_delta[idx] <- wt$median
loco_results$wilcoxon_p[idx]  <- wt$p_value
loco_results$survives[idx]    <- wt$p_value < 0.05
cat(sprintf("  Minus Transport (%d genes): %d/%d positive, p = %.4f, survives = %s\n",
            length(reduced_genes), wt$n_positive, wt$n_total, wt$p_value, wt$p_value < 0.05))

fwrite(loco_results, file.path(results_dir, "qa_loco_gsc.csv"))

# ---- 9. Synthesis vs Catabolism Contrast ------------------------------------
## *** BUG-AUDIT FIX: this section was MISSING from the original R port ***
## entirely, even though a frozen `gateB_contrast.csv` exists on disk (the
## same analysis as Gate A Lin Section 8 / Gate A Spatial Section 5: per-line
## D = Delta_Synthesis - Delta_Catabolism). Added here for parity; verified
## to reproduce the frozen file exactly (it is built from `component_delta`,
## already confirmed to match to full float64 precision above).
cat("\n=== Synthesis vs Catabolism Contrast ===\n")

D_gsc <- component_delta$Synthesis - component_delta$Catabolism
contrast_gsc <- data.frame(
  GSC_line        = gsc_lines,
  Delta_Synthesis  = component_delta$Synthesis,
  Delta_Catabolism = component_delta$Catabolism,
  D               = D_gsc,
  D_direction     = ifelse(D_gsc > 0, "synthesis_favored", "catabolism_favored"),
  stringsAsFactors = FALSE
)
contrast_test_gsc <- signed_rank_test(D_gsc)
cat(sprintf("  Synthesis-favored: %d/%d lines\n", sum(D_gsc > 0), length(D_gsc)))
cat(sprintf("  Median D: %.4f, p = %.4e\n",
            contrast_test_gsc$median, contrast_test_gsc$p_value))

fwrite(contrast_gsc, file.path(results_dir, "gateB_contrast.csv"))

cat("\nGate B complete.\n")
cat(sprintf("  Primary (11-gene): %d/%d positive, p = %.4f\n",
            test_11$n_positive, test_11$n_total, test_11$p_value))
cat(sprintf("  C/N redistribution: %d/%d positive, q = %.4f\n",
            comp_summary$n_positive[5], comp_summary$n_total[5],
            comp_summary$q_BH[5]))
