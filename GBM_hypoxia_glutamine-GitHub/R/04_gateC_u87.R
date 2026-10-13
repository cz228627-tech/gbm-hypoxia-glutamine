###############################################################################
##  04_gateC_u87.R
##  Gate C: GSE300765 — U87MG cell line under hypoxia (microarray)
##
##  Dataset: GSE300765 (GEO), platform GPL10558 (Illumina HumanHT-12 V4.0)
##  File:    series_matrix.txt.gz
##  Design:  18 samples total (various treatments)
##           We use ONLY: Hypoxia (1% O2, 48h) vs Normoxia (21% O2, 48h), n=3 each
##           Hypoxia:  GSM9068203, GSM9068204, GSM9068205 (U87_hox_1/2/3)
##           Normoxia: GSM9068206, GSM9068207, GSM9068208 (U87_nox_1/2/3)
##
##  Analysis:
##    1. Load series matrix, map 17 ILMN probes -> 13 genes (probe average)
##    2. Per-gene z-score; glutamine (11/13-gene) + 5 component scores
##    3. Unpaired Student t-test (hypoxia vs normoxia) on module scores
##    4. limma on the 13 mapped genes (supplementary)
##    5. Component-level analysis with BH correction
##
##  NOTE: U87 is an established cell line, NOT patient-derived.
##        Results are reported as "not reproduced in this established model."
###############################################################################

source("00_setup.R")
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
})

# ---- 1. Load Data -----------------------------------------------------------
## Adjust path to your local copy of the series matrix file
## Alternatively, download directly from GEO:
##   gse <- getGEO("GSE300765", GSEMatrix = TRUE)
##
## *** DATA-AVAILABILITY NOTE (RESOLVED) ***
## Probe->gene mapping no longer depends on fData(): the explicit 17-probe
## map below was verified against a fresh download of GPL10558.annot.gz
## from the GEO FTP (48,107 ILMN probes). getGEO() is used only to read the
## expression matrix. Sample mapping (GSM9068203-205 = U87_hox_1/2/3,
## GSM9068206-208 = U87_nox_1/2/3) verified against series_matrix.txt.gz.
series_matrix_file <- "data/gse300765/series_matrix.txt.gz"

cat("Loading GSE300765 series matrix...\n")
gse <- getGEO(filename = series_matrix_file)

## Extract expression matrix (probes x samples)
expr_data <- exprs(gse)
probe_ids <- rownames(expr_data)

## *** BUG-AUDIT FIX (this revision): probe -> gene mapping ***
## The previous version grepped a gene-symbol column from GEOquery feature
## data and kept ONE probe per gene (highest mean). The frozen Python
## (gateC_gse300765.py) instead uses an explicit 17-probe map and AVERAGES all
## probes of a gene (GLUL has 3, GLS2 and SLC38A1 have 2), so the old R code
## could not reproduce the frozen Gate C numbers. The explicit map is used
## below. All 17 probes were verified against the GEO GPL10558 annotation
## table (17/17 map to the intended gene). The only annotated probe for the
## 13 genes omitted by the frozen map is GLUL ILMN_1835017; including it
## leaves every component's direction/significance unchanged (11-gene
## P 0.0071 -> 0.0128; see 08_unresolved_items/gateC_probe_sensitivity.csv).
## Set INCLUDE_EXTRA_GLUL_PROBE <- TRUE to run that sensitivity variant.
## NOTE: only the 13 glutamine genes are mapped here, as in the frozen
## analysis; the hypoxia score (Section 3) is therefore set to NA and the
## limma section (Section 6) covers only these 13 genes. Both are
## supplementary to the frozen Gate C endpoints.
INCLUDE_EXTRA_GLUL_PROBE <- FALSE
probe_map <- c(
  ILMN_2188722 = "GLS",     ILMN_1677439 = "GLS2",    ILMN_1706070 = "GLS2",
  ILMN_1713756 = "GLUD1",   ILMN_1653496 = "GLUL",    ILMN_1765208 = "GLUL",
  ILMN_2269977 = "GLUL",    ILMN_1656145 = "GOT1",    ILMN_1801348 = "GOT2",
  ILMN_1795257 = "GPT",     ILMN_1684158 = "GPT2",    ILMN_1707720 = "SLC1A5",
  ILMN_1769911 = "SLC38A1", ILMN_2266948 = "SLC38A1", ILMN_1651799 = "SLC38A2",
  ILMN_1655229 = "SLC7A11", ILMN_1720373 = "SLC7A5"
)
if (INCLUDE_EXTRA_GLUL_PROBE) probe_map <- c(probe_map, ILMN_1835017 = "GLUL")
keep <- intersect(names(probe_map), probe_ids)
if (length(keep) != length(probe_map))
  warning("Missing probes in series matrix: ",
          paste(setdiff(names(probe_map), probe_ids), collapse = ", "))
## Average all probes per gene (matches pandas groupby('Gene').mean())
expr_genes <- as.matrix(aggregate(expr_data[keep, , drop = FALSE],
                                  by = list(gene = probe_map[keep]),
                                  FUN = mean)[, -1])
rownames(expr_genes) <- sort(unique(probe_map[keep]))
colnames(expr_genes) <- colnames(expr_data)

cat(sprintf("Mapped: %d probes -> %d unique genes\n",
            length(probe_ids), nrow(expr_genes)))

# ---- 2. Define Sample Groups ------------------------------------------------
## All 18 samples
all_samples <- colnames(expr_genes)

## Hypoxia (1% O2, 48h): GSM9068203-205
## Normoxia (21% O2, 48h): GSM9068206-208
hypoxia_samples  <- c("GSM9068203", "GSM9068204", "GSM9068205")
normoxia_samples <- c("GSM9068206", "GSM9068207", "GSM9068208")

## Verify samples exist
hypoxia_samples  <- intersect(hypoxia_samples, all_samples)
normoxia_samples <- intersect(normoxia_samples, all_samples)
cat(sprintf("Hypoxia samples: %s\n", paste(hypoxia_samples, collapse = ", ")))
cat(sprintf("Normoxia samples: %s\n", paste(normoxia_samples, collapse = ", ")))

## Subset to relevant samples
use_samples <- c(normoxia_samples, hypoxia_samples)
expr_sub <- expr_genes[, use_samples]
group <- factor(c(rep("Normoxia", length(normoxia_samples)),
                  rep("Hypoxia", length(hypoxia_samples))),
                levels = c("Normoxia", "Hypoxia"))

# ---- 3. Score Modules -------------------------------------------------------
cat("Scoring modules...\n")

## *** BUG-AUDIT FIX (CONFIRMED): missing per-gene z-scoring step ***
## The locked deliverable (glutamine_fate_gate_report.md) states the Gate C
## method explicitly: "17 probes -> 13 genes, z-scored, two-sided unpaired
## t-test." This R port previously went straight from raw expression to
## score_module()'s arithmetic mean, with no z-scoring at all.
## Recipe: per-gene mean-center, divide by sample SD (ddof=1, R's default
## `scale()`), across these 6 samples (3 hypoxia + 3 normoxia) -- the SAME
## general recipe that was independently proven to reproduce Gate B's frozen
## outputs to full float64 precision (see 03_gateB_gsc.R).
## VALIDATION (UPDATED): a full re-derivation from raw series-matrix probe
## intensities with the verified 17-probe map + this z-score recipe +
## Student t reproduced gateC_u87_summary.csv exactly (e.g. Exchange_redox
## t=-9.0099 p=0.00084; 11-gene t=-5.0748 p=0.0071). Earlier supporting
## evidence, kept for the record:
##   1. z-scoring across all 6 samples together (not per-group) implies, by
##      construction, that the combined mean across all 6 z-scores is exactly
##      0, so mean_hypoxia must equal -mean_normoxia exactly. This is
##      precisely the pattern observed in the frozen gateC_u87_summary.csv
##      (mean_hypoxia = -mean_normoxia to ~15 significant digits for EVERY
##      component) -- strong, direct mathematical evidence for this scope.
##   2. The identical general recipe (per-gene scale() across the analyzed
##      sample pool) was independently confirmed to exactly reproduce every
##      Gate B output once applied.
expr_z <- t(scale(t(expr_sub), center = TRUE, scale = TRUE))

## Transpose to samples x genes (now on z-scored values)
expr_t <- t(expr_z)
gene_names <- colnames(expr_t)

## Hypoxia score: supplementary only (not a frozen Gate C endpoint). With the
## 13-gene probe map above no hypoxia genes are present, so it is set to NA
## instead of calling score_module() (which would stop()).
hyp_scores <- if (any(hypoxia_genes %in% gene_names)) {
  score_module(expr_t, hypoxia_genes, gene_names)
} else {
  message("Hypoxia genes not in the 13-gene probe map; hypoxia_score = NA")
  rep(NA_real_, nrow(expr_t))
}
glq_11_scores <- score_module(expr_t, glutamine_11, gene_names)

## 5 fate components
comp_scores <- sapply(names(fate_components), function(comp) {
  score_module(expr_t, fate_components[[comp]], gene_names)
})
colnames(comp_scores) <- names(fate_components)

## Per-sample scores
per_sample <- data.frame(
  sample = rownames(expr_t),
  group = group,
  hypoxia_score = hyp_scores,
  glutamine_11 = glq_11_scores,
  comp_scores,
  stringsAsFactors = FALSE
)
fwrite(per_sample, file.path(results_dir, "gateC_u87_scores.csv"))

# ---- 4. Paired t-tests on Module Scores -------------------------------------
cat("\n=== Primary Endpoint (11-gene glutamine) ===\n")

## *** BUG-AUDIT FIX (CONFIRMED, numerically proven to full float64
## precision): unpaired AND equal-variance (Student's), not paired ***
## The locked deliverable explicitly says "two-sided UNPAIRED t-test." Two
## separate issues were found and fixed, both verified by direct
## recomputation against the frozen gateC_u87_scores.csv/gateC_u87_summary.csv:
##   1. paired=TRUE (original) vs FALSE: `paired=TRUE` does not match any of
##      the 6 frozen rows; `paired=FALSE` gets every t-STATISTIC exactly
##      right (there's no natural 1:1 pairing between the 3 hypoxia/3
##      normoxia replicates anyway -- independent replicates, not
##      before/after samples).
##   2. R's `t.test(paired=FALSE)` DEFAULTS to Welch's unequal-variance test
##      (var.equal=FALSE), while Python's `scipy.stats.ttest_ind` DEFAULTS to
##      Student's pooled/equal-variance test (equal_var=True) -- a classic
##      R-vs-Python default mismatch. With n=3 vs n=3 (balanced groups) the
##      t-STATISTIC is identical either way, but the p-VALUE is not (Welch's
##      approximate df != the pooled df=n1+n2-2=4), e.g. Exchange_redox:
##      Welch gives p=0.0019, pooled gives p=0.00084 (frozen) -- a >2x
##      difference. Adding `var.equal=TRUE` reproduces all 6 frozen
##      p-values to full float64 precision (max abs diff ~1e-16).
tt_11 <- t.test(glq_11_scores[group == "Hypoxia"],
                glq_11_scores[group == "Normoxia"],
                paired = FALSE, var.equal = TRUE, alternative = "two.sided")
## *** BUG-AUDIT FIX: two-sample t.test()$estimate is c(mean x, mean y), NOT
## the difference; use x - y (Hypoxia - Normoxia). Applied throughout. ***
est_diff <- function(tt) unname(tt$estimate[1] - tt$estimate[2])
mean_diff <- est_diff(tt_11)
cat(sprintf("  Mean diff (Hypoxia - Normoxia): %.4f\n", mean_diff))
cat(sprintf("  t-stat: %.4f, p = %.4f\n", tt_11$statistic, tt_11$p.value))

## Save 11-gene summary
summary_11 <- data.frame(
  Component = "OverallGlutamine_11gene",
  mean_hypoxia = mean(glq_11_scores[group == "Hypoxia"]),
  mean_normoxia = mean(glq_11_scores[group == "Normoxia"]),
  mean_diff = mean_diff,
  t_stat = tt_11$statistic,
  t_p_twosided = tt_11$p.value,
  direction = ifelse(mean_diff > 0, "positive", "negative"),
  stringsAsFactors = FALSE
)
fwrite(summary_11, file.path(results_dir, "gateC_u87_11gene_summary.csv"))

# ---- 5. Component-Level Analysis --------------------------------------------
cat("\n=== Component-Level Analysis ===\n")

comp_summary <- data.frame(
  Component = c(names(fate_components), "OverallGlutamine"),
  mean_hypoxia = NA_real_,
  mean_normoxia = NA_real_,
  mean_diff = NA_real_,
  direction = NA_character_,
  t_stat = NA_real_,
  t_p_twosided = NA_real_,
  q_BH_twosided = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(fate_components)) {
  comp <- names(fate_components)[i]
  vals_h <- comp_scores[group == "Hypoxia", comp]
  vals_n <- comp_scores[group == "Normoxia", comp]
  ## *** BUG-AUDIT FIX: unpaired + equal-variance -- see Section 4 note above ***
  tt <- t.test(vals_h, vals_n, paired = FALSE, var.equal = TRUE, alternative = "two.sided")
  idx <- i
  comp_summary$mean_hypoxia[idx]  <- mean(vals_h)
  comp_summary$mean_normoxia[idx] <- mean(vals_n)
  comp_summary$mean_diff[idx]     <- est_diff(tt)
  comp_summary$direction[idx]     <- ifelse(est_diff(tt) > 0, "positive", "negative")
  comp_summary$t_stat[idx]        <- tt$statistic
  comp_summary$t_p_twosided[idx]  <- tt$p.value
  cat(sprintf("  %s: diff = %.4f, t = %.4f, p = %.4f\n",
              comp, est_diff(tt), tt$statistic, tt$p.value))
}

## Overall glutamine (13-gene, not in BH family - primary endpoint)
glq_13_scores <- score_module(expr_t, glutamine_13, gene_names)
## *** BUG-AUDIT FIX: unpaired + equal-variance -- see Section 4 note above ***
tt_13 <- t.test(glq_13_scores[group == "Hypoxia"],
                glq_13_scores[group == "Normoxia"],
                paired = FALSE, var.equal = TRUE, alternative = "two.sided")
comp_summary$mean_hypoxia[6]  <- mean(glq_13_scores[group == "Hypoxia"])
comp_summary$mean_normoxia[6] <- mean(glq_13_scores[group == "Normoxia"])
comp_summary$mean_diff[6]     <- est_diff(tt_13)
comp_summary$direction[6]     <- ifelse(est_diff(tt_13) > 0, "positive", "negative")
comp_summary$t_stat[6]        <- tt_13$statistic
comp_summary$t_p_twosided[6]  <- tt_13$p.value

## BH correction across 5 components only (secondary family)
comp_summary$q_BH_twosided[1:5] <- bh_correct(comp_summary$t_p_twosided[1:5])

## Add significance flags
comp_summary$significant <- comp_summary$t_p_twosided < 0.05
comp_summary$direction_positive <- comp_summary$direction == "positive"

fwrite(comp_summary, file.path(results_dir, "gateC_u87_summary.csv"))

cat(sprintf("\n  BH q-values: %s\n",
            paste(sprintf("%s=%.4f", comp_summary$Component[1:5],
                          comp_summary$q_BH_twosided[1:5]), collapse = ", ")))

# ---- 6. limma Differential Expression (per-gene) ----------------------------
cat("\n=== limma Differential Expression ===\n")
## Supplementary: runs on the 13 mapped glutamine genes only (not
## genome-wide), so eBayes variance moderation is based on 13 genes and
## adj.P.Val is BH across these 13. Not a frozen Gate C endpoint.

design <- model.matrix(~ group)
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

de_results <- topTable(fit, coef = "groupHypoxia", number = Inf,
                       sort.by = "P")
de_results$gene <- rownames(de_results)
fwrite(de_results, file.path(results_dir, "gateC_u87_limma_de.csv"))

cat(sprintf("  DE genes (p < 0.05): %d\n", sum(de_results$P.Value < 0.05)))
cat(sprintf("  DE genes (q < 0.05): %d\n", sum(de_results$adj.P.Val < 0.05)))

cat("\nGate C complete.\n")
cat(sprintf("  11-gene: diff = %.4f, p = %.4f (direction: %s)\n",
            mean_diff, tt_11$p.value,
            ifelse(mean_diff > 0, "positive", "negative")))
cat("  Note: U87 is an established cell line, not patient-derived.\n")
cat("  Results reported as 'not reproduced in this established model.'\n")
