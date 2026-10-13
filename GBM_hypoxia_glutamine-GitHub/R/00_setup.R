###############################################################################
##  00_setup.R
##  Context-Dependent Glutamine Transcriptional Remodeling Under Hypoxia in GBM
##
##  Gene modules, color palette, and helper functions.
##  Source this file before running any gate script.
###############################################################################

# ---- Packages ---------------------------------------------------------------
# Install if missing (uncomment the lines below on first run)
#
# install.packages(c("BiocManager", "data.table", "ggplot2", "patchwork",
#                    "RColorBrewer", "ggrepel", "scales", "matrixStats", "reshape2"))
# BiocManager::install(c("zellkonverter", "Seurat", "SeuratObject",
#                        "GEOquery", "limma", "ComplexHeatmap", "circlize"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(matrixStats)
})

# ---- Gene Modules (frozen) --------------------------------------------------

## 16-gene hypoxia signature
hypoxia_genes <- c(
  "CA9", "VEGFA", "NDRG1", "P4HA1", "EGLN3", "BNIP3", "ADM", "ANKRD37",
  "CDKN1A", "GYS1", "DDIT3", "F3", "JMY", "LRP1", "LOX", "ANGPTL4"
)

## 11-gene primary glutamine module (primary endpoint)
glutamine_11 <- c(
  "GLS", "GLS2", "GLUL", "GLUD1", "GOT1", "GOT2", "GPT", "GPT2",
  "SLC1A5", "SLC7A5", "SLC7A11"
)

## 13-gene glutamine module (5 fate components, 13 genes)
fate_components <- list(
  Synthesis         = "GLUL",
  Catabolism        = c("GLS", "GLS2"),
  Direct_uptake     = c("SLC1A5", "SLC38A1", "SLC38A2"),
  Exchange_redox    = c("SLC7A5", "SLC7A11"),
  CN_redistribution = c("GLUD1", "GOT1", "GOT2", "GPT", "GPT2")
)

## Flatten to 13-gene vector
glutamine_13 <- unlist(fate_components, use.names = FALSE)

## 12 metabolic programs for landscape analysis (Figure 1B)
metabolic_programs <- list(
  glycolysis = c("HK2","PFKP","ALDOA","GAPDH","PGAM1","ENO1","PKM","LDHA","SLC2A3"),
  oxphos_tca = c("NDUFA1","NDUFB1","NDUFS1","NDUFV1","SDHA","SDHB","UQCRC1","UQCRC2",
                 "COX4I1","COX5A","COX5B","COX6A1","COX7A1","COX8A","ATP5F1A","ATP5F1B",
                 "CS","ACO2","IDH2","IDH3A","IDH3B","OGDH","DLST","MDH2"),
  fatty_acid_oxidation = c("CPT1A","CPT1B","CPT1C","CPT2","SLC25A20","ACOX1","ACOX2",
                           "ACOX3","ACADM","ACADL","ACADS","ACADVL","ECHS1","HADH",
                           "HADHA","HADHB"),
  fatty_acid_synthesis = c("FASN","ACACA","ACACB","OXSM"),
  fatty_acid_desaturation = c("SCD","SCD5","FADS1","FADS2","FADS3"),
  cholesterol_mevalonate = c("HMGCS1","HMGCR","MVK","FDPS","FDFT1","SQLE","LSS",
                             "DHCR7","DHCR24","LDLR","SREBF2"),
  lipid_uptake = c("CD36","FABP4","FABP5","FABP6","FABP7","LPL","APOE","VLDLR",
                   "LDLR","SCARB1"),
  nadph_generation = c("G6PD","PGD","TKT","TALDO1","ME2","ME3","IDH1","MTHFD1",
                       "MTHFD2","ALDH1L1","ALDH1L2"),
  gsh_redox = c("GCLC","GCLM","GSS","GSR","GPX1","GPX4","SOD1","SOD2","SOD3",
                "CAT","PRDX1","PRDX2","PRDX3","PRDX4","PRDX5","PRDX6","TXN",
                "TXNRD1","NQO1"),
  one_carbon = c("MTHFR","MTR","SHMT1","SHMT2","MTHFD1","MTHFD2","DHFR","TYMS",
                 "ATIC","GART","SLC19A1"),
  glutamine = glutamine_11,
  nucleotide = c("RRM1","RRM2","RRM2B","TYMS","DHFR","ATIC","CAD","DHODH",
                 "UMPS","TK1","DTYMK","DCK","DUT")
)

## Malignant cell state annotations (Neftel et al. classification)
malignant_states <- c("AC-like", "MES-like", "OPC-like", "NPC-like")

# ---- Color Palette ----------------------------------------------------------

## Fate component colors (frozen, used in all figures)
component_colors <- c(
  Synthesis         = "#1F77B4",
  Catabolism        = "#FF7F0E",
  Direct_uptake     = "#2CA02C",
  Exchange_redox    = "#D62728",
  CN_redistribution = "#9467BD"
)

## Direction colors
direction_colors <- c(
  positive = "#E9ED4C",
  negative = "#0279EE",
  null     = "grey70"
)

# ---- Helper Functions -------------------------------------------------------

## Arithmetic mean scoring (frozen primary method)
## Returns one score per cell/sample (row of expr matrix)
score_module <- function(expr, genes, gene_names) {
  ## expr:       matrix or sparse matrix (n_cells x n_genes)
  ## genes:      character vector, module gene symbols
  ## gene_names: character vector, colnames of expr (gene symbols)
  present <- genes[genes %in% gene_names]
  if (length(present) == 0)
    stop("None of the module genes found in expression data")
  if (length(present) < length(genes))
    warning(sprintf("Missing %d/%d module genes: %s",
                    sum(!genes %in% gene_names), length(genes),
                    paste(setdiff(genes, gene_names), collapse = ", ")))
  idx <- match(present, gene_names)
  sub <- as.matrix(expr[, idx, drop = FALSE])
  rowMeans(sub, na.rm = TRUE)
}

## Per-patient Spearman correlation between two score vectors
## Returns a data.frame with one row per patient
##
## NOTE on `min_cells` (BUG AUDIT FINDING, see README "Known Issues"):
## For the main 28-patient Gate A analysis, patients are pre-filtered to
## >=50 malignant cells, so this threshold never binds there. It DOES bind
## for the within-malignant-state analysis (qa_malignant_state.csv), where
## some patients have very few cells of a given state. RESOLVED: the original
## Python (qa_steps4_5.py) uses MIN_STATE_CELLS = 20, and 01_gateA_lin.R
## Section 9 passes min_cells = 20 explicitly; this reproduces the frozen
## qa_malignant_state_summary.csv exactly (AC 25, MES 25, OPC 19, NPC 18).
## The function default of 10 is kept only for backward compatibility; it
## never binds in the main analysis (patients pre-filtered to >=50 cells).
per_patient_spearman <- function(score_x, score_y, patient_ids, min_cells = 10) {
  patients <- unique(patient_ids)
  results <- data.frame(
    patient = patients,
    rho = NA_real_,
    p_value = NA_real_,
    n_cells = NA_integer_,
    direction = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(patients)) {
    p <- patients[i]
    mask <- patient_ids == p
    x <- score_x[mask]
    y <- score_y[mask]
    n <- sum(mask)
    if (n < min_cells) next  # skip patients with too few cells
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    results$rho[i]       <- ct$estimate
    results$p_value[i]   <- ct$p.value
    results$n_cells[i]   <- n
    results$direction[i] <- ifelse(ct$estimate > 0, "positive", "negative")
  }
  results
}

## Wilcoxon signed-rank test on a vector of rho/delta values
## (two-sided, testing whether median != 0)
##
## BUG-AUDIT FIX (CONFIRMED, high-confidence): match scipy.stats.wilcoxon's
## default zero-handling. scipy's default `zero_method="wilcox"` DROPS exact-
## zero differences before ranking; a plain R `wilcox.test(values)` does NOT
## drop them the same way and gives a DIFFERENT p-value whenever exact zeroes
## are present (verified with a continuity-corrected normal approximation
## discrepancy). Reproduced exactly: Gate B's Synthesis component (a single-
## gene, GLUL-only score) has an exact-zero z-scored hypoxia-normoxia delta
## for 2/8 GSC lines; `wilcox.test(values)` gives p=0.4017 vs the frozen
## p=0.4375 (matches scipy's default). Dropping zeroes first reproduces
## 0.4375 exactly. n_total/n_positive/n_negative/median are still computed
## on the FULL vector (zeroes included) -- only the test statistic/p-value
## input excludes them, matching scipy's behavior and the frozen CSVs.
signed_rank_test <- function(values) {
  values <- values[!is.na(values)]
  n <- length(values)
  n_pos <- sum(values > 0)
  n_neg <- sum(values < 0)
  nonzero <- values[values != 0]
  if (length(nonzero) == 0) nonzero <- values  # degenerate all-zero fallback
  wt <- wilcox.test(nonzero, alternative = "two.sided")
  list(
    n_total   = n,
    n_positive = n_pos,
    n_negative = n_neg,
    median    = median(values),
    statistic = wt$statistic,
    p_value   = wt$p.value
  )
}

## BH correction across a family of p-values
bh_correct <- function(p_values) {
  p.adjust(p_values, method = "BH")
}

## Extract patient ID from Greenwald Visium section name
## e.g. "GSM7596588_zh881inf" -> "zh881"
extract_patient_from_section <- function(section_name) {
  ## Remove GSM accession prefix
  name <- sub("^GSM\\d+_", "", section_name)
  ## Remove known suffixes
  name <- sub("(inf|t1|bulk|nec|a|b)$", "", name)
  name
}

# ---- Output Directory -------------------------------------------------------
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat("Setup complete. Gene modules loaded:\n")
cat(sprintf("  Hypoxia: %d genes\n", length(hypoxia_genes)))
cat(sprintf("  Glutamine 11-gene (primary): %d genes\n", length(glutamine_11)))
cat(sprintf("  Glutamine 13-gene (5 components): %d genes\n", length(glutamine_13)))
cat(sprintf("  Metabolic programs: %d programs\n", length(metabolic_programs)))
