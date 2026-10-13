#!/usr/bin/env python3
"""Recompute U87 11-gene overall glutamine (two-sided t-test)."""
import warnings; warnings.filterwarnings('ignore')
import os
import numpy as np, pandas as pd, gzip
from scipy import stats
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
DATA_FILE = Path(os.environ.get("GBM_U87_FILE", DATA_DIR / "gse300765" / "series_matrix.txt.gz")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

GLUTAMINE_11 = ['GLS','GLS2','GLUL','GLUD1','GOT1','GOT2','GPT','GPT2',
                'SLC1A5','SLC7A5','SLC7A11']

PROBE_MAP = {
    'ILMN_2188722':'GLS','ILMN_1677439':'GLS2','ILMN_1706070':'GLS2',
    'ILMN_1713756':'GLUD1','ILMN_1653496':'GLUL','ILMN_1765208':'GLUL',
    'ILMN_2269977':'GLUL','ILMN_1656145':'GOT1','ILMN_1801348':'GOT2',
    'ILMN_1795257':'GPT','ILMN_1684158':'GPT2','ILMN_1707720':'SLC1A5',
    'ILMN_1769911':'SLC38A1','ILMN_2266948':'SLC38A1','ILMN_1651799':'SLC38A2',
    'ILMN_1655229':'SLC7A11','ILMN_1720373':'SLC7A5',
}

HYPOXIA_GSM = ['GSM9068203','GSM9068204','GSM9068205']
NORMOXIA_GSM = ['GSM9068206','GSM9068207','GSM9068208']

# Parse series matrix
with gzip.open(DATA_FILE, 'rt') as f:
    in_table = False; data_lines = []; col_names = None
    for line in f:
        line = line.strip()
        if line == '!series_matrix_table_begin': in_table = True; continue
        if line == '!series_matrix_table_end': break
        if in_table:
            if col_names is None:
                col_names = [c.strip('"') for c in line.split('\t')]
            else:
                data_lines.append(line)

rows = []; indices = []
for line in data_lines:
    fields = [c.strip('"') for c in line.split('\t')]
    indices.append(fields[0])
    rows.append([float(v) if v else np.nan for v in fields[1:]])
expr = pd.DataFrame(rows, index=indices, columns=col_names[1:])

# Filter to 11-gene probes only
probes_11 = [p for p in PROBE_MAP if PROBE_MAP[p] in GLUTAMINE_11 and p in expr.index]
expr_fate = expr.loc[probes_11].copy()
expr_fate['Gene'] = expr_fate.index.map(PROBE_MAP)

all_samples = HYPOXIA_GSM + NORMOXIA_GSM
expr_gene = expr_fate.groupby('Gene')[all_samples].mean()
print(f"11-gene genes found: {sorted(expr_gene.index)} ({len(expr_gene)}/11)")

# Z-score each gene across 6 samples
expr_z = expr_gene.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1)

# Score 11-gene overall
glut_score = expr_z.mean(axis=0)

hyp_vals = glut_score[HYPOXIA_GSM].values
norm_vals = glut_score[NORMOXIA_GSM].values
mean_diff = np.mean(hyp_vals) - np.mean(norm_vals)
t_stat, t_p = stats.ttest_ind(hyp_vals, norm_vals)  # two-sided

print(f"\n11-gene Overall Glutamine:")
print(f"  mean_diff (hypoxia - normoxia) = {mean_diff:.4f}")
print(f"  t_stat = {t_stat:.4f}")
print(f"  p (two-sided) = {t_p:.6f}")
print(f"  direction = {'negative' if mean_diff < 0 else 'positive'}")

# Save
df_out = pd.DataFrame([{
    'Component': 'OverallGlutamine_11gene',
    'mean_hypoxia': np.mean(hyp_vals),
    'mean_normoxia': np.mean(norm_vals),
    'mean_diff': mean_diff,
    't_stat': t_stat,
    't_p_twosided': t_p,
    'direction': 'negative' if mean_diff < 0 else 'positive',
}])
df_out.to_csv(OUT / "gateC_u87_11gene_summary.csv", index=False)
print(f"\nSaved to gateC_u87_11gene_summary.csv")
