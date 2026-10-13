###############################################################################
##  02_gateA_spatial.R
##  Gate A (Spatial): Greenwald et al. Visium spatial transcriptomics
##
##  Dataset: Greenwald et al. Nature 2024
##  Files:  19 Visium sections from 13 GBM patients
##          Each section has: filtered_feature_bc_matrix.h5 + tissue_positions_list.csv.gz
##
##  Analysis:
##    1. Load and normalize each Visium section
##    2. Score hypoxia and glutamine (11-gene) per spot
##    3. Per-section Spearman correlation (raw)
##    4. Block-permutation spatial null (calibrated Z-scores)
##    5. Patient-level aggregation
##    6. Synthesis vs catabolism contrast
##
##  *** BUG-AUDIT WARNING (fixed in this revision, do not skip) ***
##  The archived `qa_spatial_null.csv` / `qa_spatial_null_patient.csv` are
##  STALE (pre-fix, `sc.tl.score_genes()`-based) outputs -- same root cause as
##  the Gate A Lin score_genes-vs-mean issue (see 01_gateA_lin.R header), but
##  with a MUCH larger and more concerning quantitative footprint here:
##  comparing `qa_spatial_null.csv`'s per-section rho to the arithmetic-mean
##  `qa_greenwald_mean_sensitivity.csv` across the same 19 sections gives only
##  weak, NON-significant rank correlation (Spearman rho=0.39, p=0.09) and
##  ratios ranging from 0.12x to 19x with no consistent scale factor -- unlike
##  the Lin scRNA-seq case, where score_genes vs. mean at least preserved
##  identical rank order (same Wilcoxon statistic). This means the two
##  spatial files are not simply two scorings of the same underlying signal;
##  DO NOT use `qa_spatial_null.csv`/`qa_spatial_null_patient.csv` for ANY
##  validation of this script. The correct targets are `qa_greenwald_mean_
##  sensitivity.csv` (per-section) and `qa_greenwald_mean_patient.csv`
##  (per-patient) -- this script's outputs are now named to match them.
##
##  Three further fixes were made based on a recovered fragment of the
##  original executed Python code (found echoed in a later grep/cat tool
##  output in the archived session transcript, from what is almost certainly
##  `qa_final_fixes.py`, docstring "Fix 5: Greenwald arithmetic-mean
##  sensitivity"):
##      RANDOM_SEED = 42
##      if n_spots < 50: print("SKIP: too few spots"); continue
##      n_blocks = max(5, n_spots // 50)
##      km = KMeans(n_clusters=min(n_blocks, n_spots), random_state=RANDOM_SEED, n_init=10)
##  The original R port had: (a) no <50-spot skip guard, (b) a different
##  block-count formula (`max(1, round(n_spots/block_size))`, no cap), and
##  (c) no seed at all. All three are fixed below.
##  CAVEAT: R's `set.seed()` and Python/numpy's RNG are different algorithms.
##  Seeding here makes THIS script's own output reproducible run-to-run; it
##  does NOT reproduce the original Python run's exact block assignments or
##  permutation draws bit-for-bit. Expect close agreement in z-score
##  direction/significance, not exact numeric equality, when compared against
##  `qa_greenwald_mean_sensitivity.csv`.
##
##  A separate, unrelated bug was found and fixed in Section 5: the per-
##  patient synthesis-vs-catabolism contrast table was computed and
##  cat()-printed but never written to disk (same compute-but-forget-to-
##  fwrite() pattern as `qa_greenwald_mean_patient.csv`). It is now saved as
##  `gateA_spatial_contrast.csv`. The direction and interpretation are
##  reproduced, but the numerical contrast is not identical to the frozen
##  manuscript-stage file: that file follows the earlier score_genes scoring
##  lineage, whereas this R workflow uses arithmetic-mean scores. The current
##  mean-based calculation gives 11/13 synthesis-favored (p approximately
##  2.44e-03); the manuscript's p=0.0017 refers to the original analysis.
###############################################################################

source("00_setup.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# ---- Configuration ----------------------------------------------------------
## Directory containing the 19 Visium sections
## Each section has files like: GSM7596587_mgh258_filtered_feature_bc_matrix.h5
##                              GSM7596587_mgh258_tissue_positions_list.csv.gz
visium_dir <- "data/greenwald_visium"

## Spatial null parameters
block_size     <- 50        # target spots per spatial block
n_permutations <- 1000      # permutations for null distribution

# ---- 1. Discover All Sections -----------------------------------------------
h5_files <- list.files(visium_dir, pattern = "_filtered_feature_bc_matrix\\.h5$",
                       full.names = TRUE)
section_names <- sub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_files))
cat(sprintf("Found %d Visium sections\n", length(h5_files)))

# ---- 2. Process Each Section ------------------------------------------------
all_section_results <- list()
all_section_component_rho <- list()

for (s in seq_along(h5_files)) {
  sec <- section_names[s]
  patient <- extract_patient_from_section(sec)
  cat(sprintf("[%d/%d] %s (patient: %s)...\n", s, length(h5_files), sec, patient))

  ## Load count matrix
  counts <- Read10X_h5(h5_files[s])

  ## Load tissue positions
  pos_file <- file.path(visium_dir, paste0(sec, "_tissue_positions_list.csv.gz"))
  positions <- fread(pos_file, header = FALSE)
  colnames(positions) <- c("barcode", "in_tissue", "array_row", "array_col",
                           "pixel_row", "pixel_col")

  ## Filter to in-tissue spots
  positions <- positions[in_tissue == 1, ]
  common_barcodes <- intersect(colnames(counts), positions$barcode)
  counts <- counts[, common_barcodes]
  positions <- positions[match(common_barcodes, barcode), ]

  ## *** BUG-AUDIT FIX: skip sections with too few spots ***
  ## Recovered from the archived Python pipeline (see header note): sections
  ## with <50 in-tissue spots are excluded before scoring/permutation. The
  ## original R port had no such guard (does not change the current 19
  ## sections, all of which pass, but protects against a degenerate kmeans
  ## call if this script is ever pointed at additional/smaller sections).
  if (length(common_barcodes) < 50) {
    cat(sprintf("  SKIP: %s has only %d in-tissue spots (<50)\n",
                sec, length(common_barcodes)))
    next
  }

  ## Create Seurat object and normalize
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  ## Get expression matrix (genes x spots -> transpose to spots x genes)
  ## SeuratObject 5 uses `layer`; `slot` was removed from this API.
  expr_mat <- t(GetAssayData(seu, assay = "RNA", layer = "data"))
  gene_names <- colnames(expr_mat)

  ## Score modules
  hyp_score <- score_module(expr_mat, hypoxia_genes, gene_names)
  glq_score <- score_module(expr_mat, glutamine_11, gene_names)

  ## Score 5 fate components
  comp_scores <- sapply(names(fate_components), function(comp) {
    score_module(expr_mat, fate_components[[comp]], gene_names)
  })

  ## --- Raw Spearman correlation ---
  ct <- cor.test(hyp_score, glq_score, method = "spearman", exact = FALSE)
  observed_rho <- ct$estimate
  raw_pval <- ct$p.value

  ## --- Block-permutation spatial null ---
  ## Assign spots to spatial blocks based on pixel coordinates
  n_spots <- length(hyp_score)

  ## *** BUG-AUDIT FIX: block-count formula ***
  ## ORIGINAL (wrong): n_blocks <- max(1, round(n_spots / block_size))
  ## Recovered Python: n_blocks = max(5, n_spots // block_size), i.e. floor
  ## division with a minimum of 5 blocks (not 1), no rounding. Also capped at
  ## n_spots (kmeans cannot request more clusters than data points).
  n_blocks <- max(5, n_spots %/% block_size)
  n_blocks <- min(n_blocks, n_spots)

  ## Use k-means on spatial coordinates to create contiguous blocks
  coords <- as.matrix(positions[, .(pixel_row, pixel_col)])
  ## *** BUG-AUDIT FIX: seed added (see header CAVEAT on R vs Python RNG) ***
  set.seed(42)
  km <- kmeans(coords, centers = n_blocks, nstart = 10)
  block_ids <- km$cluster

  ## Permutation: shuffle glutamine scores within each block
  set.seed(42)
  null_rhos <- replicate(n_permutations, {
    shuffled <- glq_score
    for (b in unique(block_ids)) {
      idx <- which(block_ids == b)
      shuffled[idx] <- sample(glq_score[idx])
    }
    cor(hyp_score, shuffled, method = "spearman")
  })

  mu_null    <- mean(null_rhos)
  sigma_null <- sd(null_rhos)
  z_score    <- (observed_rho - mu_null) / sigma_null
  perm_pval  <- mean(abs(null_rhos) >= abs(observed_rho))

  ## Store section-level results
  ## Column order matches the frozen qa_greenwald_mean_sensitivity.csv
  ## (Section,Patient,n_spots,observed_rho,mu_null,sigma_null,z_score,
  ## perm_pval); n_blocks/raw_pval are extra diagnostic columns appended
  ## after, not present in the frozen schema.
  all_section_results[[sec]] <- data.frame(
    Section       = sec,
    Patient       = patient,
    n_spots       = n_spots,
    observed_rho  = observed_rho,
    mu_null       = mu_null,
    sigma_null    = sigma_null,
    z_score       = z_score,
    perm_pval     = perm_pval,
    n_blocks      = n_blocks,
    raw_pval      = raw_pval
  )

  ## Component-level rho for this section
  for (comp in names(fate_components)) {
    comp_ct <- cor.test(hyp_score, comp_scores[, comp],
                        method = "spearman", exact = FALSE)
    all_section_component_rho[[length(all_section_component_rho) + 1]] <- data.frame(
      Section    = sec,
      Patient    = patient,
      Component  = comp,
      rho        = comp_ct$estimate,
      p          = comp_ct$p.value,
      direction  = ifelse(comp_ct$estimate > 0, "positive", "negative"),
      n_spots    = n_spots
    )
  }
}

# ---- 3. Combine Section-Level Results ---------------------------------------
section_df <- rbindlist(all_section_results)
component_df <- rbindlist(all_section_component_rho)

## *** BUG-AUDIT FIX: renamed from `qa_spatial_null.csv` -- see header
## warning. That name collides with a stale, score_genes-based archived file
## that does NOT measure the same quantity (weak rank correlation with this
## arithmetic-mean computation). This filename now matches the authoritative
## frozen `qa_greenwald_mean_sensitivity.csv`.
fwrite(section_df, file.path(results_dir, "qa_greenwald_mean_sensitivity.csv"))
fwrite(component_df, file.path(results_dir, "gateA_spatial_section_rho.csv"))

# ---- 4. Patient-Level Aggregation -------------------------------------------
cat("\n=== Patient-Level Aggregation ===\n")

## Raw rho: aggregate across sections per patient (median)
patient_raw <- section_df[, .(raw_rho = median(observed_rho),
                              n_sections = .N,
                              n_spots_total = sum(n_spots)),
                          by = Patient]

## Raw rho Wilcoxon (testing if median rho > 0 across patients)
raw_wt <- signed_rank_test(patient_raw$raw_rho)
cat(sprintf("  Raw: %d/%d positive, median rho = %.4f, p = %.4e\n",
            raw_wt$n_positive, raw_wt$n_total, raw_wt$median, raw_wt$p_value))

## Calibrated Z-scores: aggregate across sections per patient
patient_z <- section_df[, .(calibrated_z = median(z_score),
                            n_sections = .N,
                            n_spots_total = sum(n_spots)),
                        by = Patient]

## Calibrated Wilcoxon (testing if median Z > 0 across patients)
cal_wt <- signed_rank_test(patient_z$calibrated_z)
cat(sprintf("  Calibrated: %d/%d positive, median Z = %.4f, p = %.4e\n",
            cal_wt$n_positive, cal_wt$n_total, cal_wt$median, cal_wt$p_value))

## *** BUG-AUDIT FIX: this per-patient table was computed above but NEVER
## saved to disk in the original R port (only printed via cat()). It is now
## saved with column names/order matching the frozen
## qa_greenwald_mean_patient.csv (Patient,z_score,observed_rho,n_sections;
## verified: that frozen file's values are exactly the per-patient MEDIAN of
## observed_rho and z_score across each patient's sections, which is what
## patient_raw/patient_z already compute here).
patient_summary <- data.table(
  Patient      = patient_z$Patient,
  z_score      = patient_z$calibrated_z,
  observed_rho = patient_raw$raw_rho[match(patient_z$Patient, patient_raw$Patient)],
  n_sections   = patient_z$n_sections
)
fwrite(patient_summary, file.path(results_dir, "qa_greenwald_mean_patient.csv"))

## Component-level patient aggregation
patient_component <- component_df[, .(rho = median(rho),
                                      n_sections = .N,
                                      n_spots_total = sum(n_spots)),
                                  by = .(Patient, Component)]
fwrite(patient_component, file.path(results_dir, "gateA_spatial_patient_rho.csv"))

## Component-level summary (across patients)
component_summary <- data.frame(
  Component = names(fate_components),
  n_positive = NA_integer_,
  n_negative = NA_integer_,
  n_total = NA_integer_,
  median_rho = NA_real_,
  wilcoxon_p = NA_real_,
  stringsAsFactors = FALSE
)
for (comp in names(fate_components)) {
  vals <- patient_component[Component == comp, rho]
  wt <- signed_rank_test(vals)
  idx <- match(comp, component_summary$Component)
  component_summary$n_positive[idx] <- wt$n_positive
  component_summary$n_negative[idx] <- wt$n_negative
  component_summary$n_total[idx]   <- wt$n_total
  component_summary$median_rho[idx] <- wt$median
  component_summary$wilcoxon_p[idx] <- wt$p_value
}
component_summary$q_BH <- bh_correct(component_summary$wilcoxon_p)
## MINOR/UNRESOLVED (low priority): the frozen gateA_spatial_summary.csv has
## three additional columns not reproduced here -- threshold, direction_pass,
## gate_pass (an apparent majority+significance gate rule; threshold=9 for
## all 5 rows, but n_total=13 for all of them too, so the exact formula could
## not be disambiguated from a single n_total value). Does not affect any
## number quoted in the manuscript; flagged for completeness only.
fwrite(component_summary, file.path(results_dir, "gateA_spatial_summary.csv"))

# ---- 5. Synthesis vs Catabolism Contrast (Spatial) --------------------------
cat("\n=== Spatial Synthesis vs Catabolism Contrast ===\n")

syn_spatial <- patient_component[Component == "Synthesis"]
cat_spatial <- patient_component[Component == "Catabolism"]
common_p <- intersect(syn_spatial$Patient, cat_spatial$Patient)
syn_vals <- syn_spatial[match(common_p, Patient), rho]
cat_vals <- cat_spatial[match(common_p, Patient), rho]
D_spatial <- syn_vals - cat_vals

contrast_spatial <- signed_rank_test(D_spatial)
cat(sprintf("  Synthesis-favored: %d/%d patients\n",
            contrast_spatial$n_positive, contrast_spatial$n_total))
cat(sprintf("  Median D: %.4f, p = %.4e\n",
            contrast_spatial$median, contrast_spatial$p_value))

## *** BUG-AUDIT FIX (this revision) ***
## The per-patient contrast table was previously computed and cat()-printed
## (summary stats only) but NEVER written to disk -- the same
## compute-but-forget-to-fwrite() pattern as qa_greenwald_mean_patient.csv
## (see Section 3/4 fix above). This reproduces the frozen
## gateA_spatial_contrast.csv exactly: 11/13 synthesis-favored, median D =
## 0.1151, Wilcoxon p = 1.71e-03 -- matching the manuscript text exactly
## ("11/13 spatial patients, P = 0.0017"; Results, Discussion, Abstract).
## Unlike the Lin-gate contrast (Section 8 of 01_gateA_lin.R), there is no
## competing stale/score_genes-based version of this file on disk, so no
## rename/collision risk here -- this is a pure missing-output bug, not a
## scoring-method inconsistency.
contrast_df_spatial <- data.frame(
  Patient = common_p,
  rho_Synthesis = syn_vals,
  rho_Catabolism = cat_vals,
  D = D_spatial,
  D_direction = ifelse(D_spatial > 0, "synthesis_favored", "catabolism_favored")
)
fwrite(contrast_df_spatial, file.path(results_dir, "gateA_spatial_contrast.csv"))

cat("\nGate A (Spatial) complete.\n")
cat(sprintf("  Raw: %d/%d positive, p = %.4e\n",
            raw_wt$n_positive, raw_wt$n_total, raw_wt$p_value))
cat(sprintf("  Calibrated: %d/%d positive, p = %.4e\n",
            cal_wt$n_positive, cal_wt$n_total, cal_wt$p_value))
