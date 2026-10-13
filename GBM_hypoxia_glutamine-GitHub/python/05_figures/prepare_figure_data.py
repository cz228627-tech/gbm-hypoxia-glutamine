#!/usr/bin/env python3
"""
Data preparation for figure revisions:
1. Recompute 11-gene GSC per-line delta (replacing 13-gene)
2. Prepare Visium spatial data for representative section
"""

import warnings
warnings.filterwarnings('ignore')
import os

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Recompute 11-gene GSC per-line delta
# ============================================================
print("=== Recomputing 11-gene GSC per-line delta ===")

DATA_FILE = Path(os.environ.get("GBM_GSC_FILE", DATA_DIR / "gse243848" / "deseq2_normalized.txt.gz")).expanduser()
GSC_LINES = ['G44', 'G146', 'G157', 'G816', 'G33', 'G34', 'G62', 'G88']

GLUTAMINE_11 = ['GLS','GLS2','GLUL','GLUD1','GOT1','GOT2','GPT','GPT2',
                'SLC1A5','SLC7A5','SLC7A11']
GLUTAMINE_13 = GLUTAMINE_11 + ['SLC38A1', 'SLC38A2']

df = pd.read_csv(DATA_FILE, sep='\t')
gene_col = df.columns[1]  # gene_Ensembl actually contains symbols
df = df.rename(columns={gene_col: 'symbol'})
df = df[df['symbol'].notna()]
df = df.drop_duplicates(subset='symbol')
df = df.set_index('symbol')

# Extract DMSO control samples
control_samples = []
for line in GSC_LINES:
    norm_col = f"{line}_ctr"
    hyp_col = f"{line}_ctr_H"
    if norm_col in df.columns:
        control_samples.append(norm_col)
    if hyp_col in df.columns:
        control_samples.append(hyp_col)

print(f"  Control samples: {len(control_samples)}")

expr = df[control_samples].copy()
expr = expr[(expr > 0).any(axis=1)]
expr_log = np.log1p(expr)
expr_z = expr_log.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1)

# Score 11-gene overall glutamine
glut11_present = [g for g in GLUTAMINE_11 if g in expr_z.index]
glut13_present = [g for g in GLUTAMINE_13 if g in expr_z.index]
print(f"  11-gene present: {glut11_present} ({len(glut11_present)}/11)")
print(f"  13-gene present: {glut13_present} ({len(glut13_present)}/13)")

glut11_score = expr_z.loc[glut11_present].mean(axis=0)
glut13_score = expr_z.loc[glut13_present].mean(axis=0)

# Per-line delta
per_line = []
for line in GSC_LINES:
    norm_col = f"{line}_ctr"
    hyp_col = f"{line}_ctr_H"
    if norm_col in glut11_score.index and hyp_col in glut11_score.index:
        d11 = glut11_score[hyp_col] - glut11_score[norm_col]
        d13 = glut13_score[hyp_col] - glut13_score[norm_col]
        per_line.append({
            'GSC_line': line,
            'Delta_OverallGlutamine_11': d11,
            'Delta_OverallGlutamine_13': d13,
            'Normoxia_11': glut11_score[norm_col],
            'Hypoxia_11': glut11_score[hyp_col],
            'Normoxia_13': glut13_score[norm_col],
            'Hypoxia_13': glut13_score[hyp_col],
        })

df_delta = pd.DataFrame(per_line)
df_delta.to_csv(OUT / "gateB_per_line_delta_11gene.csv", index=False)

n_pos_11 = (df_delta['Delta_OverallGlutamine_11'] > 0).sum()
n_pos_13 = (df_delta['Delta_OverallGlutamine_13'] > 0).sum()
_, p_11 = stats.wilcoxon(df_delta['Delta_OverallGlutamine_11'])
_, p_13 = stats.wilcoxon(df_delta['Delta_OverallGlutamine_13'])
median_11 = df_delta['Delta_OverallGlutamine_11'].median()

print(f"  11-gene: {n_pos_11}/8+, median Δ={median_11:.3f}, p={p_11:.4f}")
print(f"  13-gene: {n_pos_13}/8+, median Δ={df_delta['Delta_OverallGlutamine_13'].median():.3f}, p={p_13:.4f}")
print(f"  Saved to gateB_per_line_delta_11gene.csv")

# ============================================================
# 2. Prepare Visium spatial data for representative section
# ============================================================
print("\n=== Preparing Visium spatial data ===")

VISIUM_DIR = Path(os.environ.get("GBM_GREENWALD_DIR", DATA_DIR / "greenwald_visium")).expanduser()
# Use GSM7596590_zh916bulk — strong synthesis ρ=0.328, negative catabolism, 2199 spots
SECTION = "GSM7596590_zh916bulk"

h5_path = VISIUM_DIR / f"{SECTION}_filtered_feature_bc_matrix.h5"
pos_path = VISIUM_DIR / f"{SECTION}_tissue_positions_list.csv.gz"

adata = sc.read_visium(VISIUM_DIR, count_file=f"{SECTION}_filtered_feature_bc_matrix.h5")
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

# Load tissue positions
positions = pd.read_csv(pos_path, header=None,
    names=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
positions = positions[positions['in_tissue'] == 1]
positions = positions.set_index('barcode')

# Filter to in-tissue spots
adata = adata[adata.obs_names.isin(positions.index)].copy()
adata.obs = adata.obs.join(positions[['pxl_row','pxl_col','array_row','array_col']])

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Score hypoxia and glutamine
HYPOXIA_GENES = ['CA9','VEGFA','NDRG1','P4HA1','EGLN3','BNIP3','ADM','ANKRD37',
                 'CDKN1A','GYS1','DDIT3','F3','JMY','LRP1','LOX','ANGPTL4']
GLUTAMINE_11 = ['GLS','GLS2','GLUL','GLUD1','GOT1','GOT2','GPT','GPT2',
                'SLC1A5','SLC7A5','SLC7A11']

def arith_mean(adata, genes):
    present = [g for g in genes if g in adata.var_names]
    expr = adata[:, present].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return expr.mean(axis=1)

adata.obs['HypoxiaScore'] = arith_mean(adata, HYPOXIA_GENES)
adata.obs['GlutamineScore'] = arith_mean(adata, GLUTAMINE_11)
adata.obs['GLUL_expr'] = adata[:, 'GLUL'].X.toarray().ravel() if 'GLUL' in adata.var_names else 0

# Save spot data
spot_data = adata.obs[['pxl_row','pxl_col','array_row','array_col',
                        'HypoxiaScore','GlutamineScore','GLUL_expr']].copy()
spot_data['Section'] = SECTION
spot_data['Patient'] = 'zh916'
spot_data.to_csv(OUT / "visium_representative_section.csv", index=True)

# Compute section-level rho
rho, p = stats.spearmanr(adata.obs['HypoxiaScore'], adata.obs['GlutamineScore'])
print(f"  Section: {SECTION}")
print(f"  Spots: {len(spot_data)}")
print(f"  Hypoxia-Glutamine ρ: {rho:.4f}, p={p:.2e}")
print(f"  GLUL range: {adata.obs['GLUL_expr'].min():.2f} - {adata.obs['GLUL_expr'].max():.2f}")
print(f"  Saved to visium_representative_section.csv")

# Also prepare a second section for optional dual-panel
# GSM7596597_zh8811a — strong synthesis ρ=0.261, 2486 spots
SECTION2 = "GSM7596597_zh8811a"
adata2 = sc.read_visium(VISIUM_DIR, count_file=f"{SECTION2}_filtered_feature_bc_matrix.h5")
adata2.var_names = adata2.var_names.astype(str)
adata2.var_names_make_unique()

pos2 = pd.read_csv(VISIUM_DIR / f"{SECTION2}_tissue_positions_list.csv.gz", header=None,
    names=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
pos2 = pos2[pos2['in_tissue'] == 1].set_index('barcode')
adata2 = adata2[adata2.obs_names.isin(pos2.index)].copy()
adata2.obs = adata2.obs.join(pos2[['pxl_row','pxl_col','array_row','array_col']])
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)
adata2.obs['HypoxiaScore'] = arith_mean(adata2, HYPOXIA_GENES)
adata2.obs['GlutamineScore'] = arith_mean(adata2, GLUTAMINE_11)
adata2.obs['GLUL_expr'] = adata2[:, 'GLUL'].X.toarray().ravel() if 'GLUL' in adata2.var_names else 0

spot_data2 = adata2.obs[['pxl_row','pxl_col','array_row','array_col',
                          'HypoxiaScore','GlutamineScore','GLUL_expr']].copy()
spot_data2['Section'] = SECTION2
spot_data2['Patient'] = 'zh8811'
spot_data2.to_csv(OUT / "visium_representative_section2.csv", index=True)

rho2, p2 = stats.spearmanr(adata2.obs['HypoxiaScore'], adata2.obs['GlutamineScore'])
print(f"\n  Section 2: {SECTION2}")
print(f"  Spots: {len(spot_data2)}")
print(f"  Hypoxia-Glutamine ρ: {rho2:.4f}, p={p2:.2e}")
print(f"  Saved to visium_representative_section2.csv")

print("\n=== Data preparation complete ===")
