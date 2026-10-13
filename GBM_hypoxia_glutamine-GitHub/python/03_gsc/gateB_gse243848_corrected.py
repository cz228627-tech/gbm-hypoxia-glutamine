#!/usr/bin/env python3
"""
Gate B CORRECTED: GSE243848 perturbation — 8 patient-derived GSC lines.
DMSO control arm: Normoxia (21% O2) vs Hypoxia (1% O2), 24h.

Two separate tests per component:
  1. Replication test: >=6/8 in Lin-predicted direction + two-sided Wilcoxon p < 0.05
  2. Acute-response discovery: >=6/8 in either direction + two-sided Wilcoxon p < 0.05

All tests are TWO-SIDED (unified statistical definition).
"""

import warnings
warnings.filterwarnings('ignore')
import os

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import json

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
DATA_FILE = Path(os.environ.get("GBM_GSC_FILE", DATA_DIR / "gse243848" / "deseq2_normalized.txt.gz")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

FATE_COMPONENTS = {
    'Synthesis':              ['GLUL'],
    'Catabolism':             ['GLS', 'GLS2'],
    'Direct_uptake':          ['SLC1A5', 'SLC38A1', 'SLC38A2'],
    'Exchange_redox':         ['SLC7A5', 'SLC7A11'],
    'CN_redistribution':      ['GLUD1', 'GOT1', 'GOT2', 'GPT', 'GPT2'],
}

ALL_FATE_GENES = sorted(set(g for gs in FATE_COMPONENTS.values() for g in gs))
GLUTAMINE_ALL = ALL_FATE_GENES

GSC_LINES = ['G44', 'G146', 'G157', 'G816', 'G33', 'G34', 'G62', 'G88']

# Lin Gate A predicted directions (from gateA_component_summary.csv)
# Positive = majority of patients showed positive ρ
LIN_PREDICTED = {
    'Synthesis':         'positive',   # 27/28+
    'Catabolism':        'positive',   # 24/28+
    'Direct_uptake':     'positive',   # 27/28+
    'Exchange_redox':    'positive',   # 26/28+
    'CN_redistribution': 'negative',   # 10/28+ (i.e., 18/28 negative)
    'OverallGlutamine':  'positive',   # 28/28+ (from original Phase 1)
}

def bh_correction(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    bh = ranked * n / np.arange(1, n + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    result = np.empty(n)
    result[order] = bh
    return result

def main():
    print("Loading GSE243848 DESeq2 normalized counts...")
    df = pd.read_csv(DATA_FILE, sep='\t')
    gene_col = df.columns[1]
    df = df.rename(columns={gene_col: 'symbol'})
    df = df[df['symbol'].notna()]
    df = df.drop_duplicates(subset='symbol')
    df = df.set_index('symbol')

    # Extract DMSO control arm
    ctr_normoxia = {}
    ctr_hypoxia = {}
    for line in GSC_LINES:
        norm_col = f"{line}_ctr"
        hyp_col = f"{line}_ctr_H"
        if norm_col in df.columns:
            ctr_normoxia[line] = norm_col
        if hyp_col in df.columns:
            ctr_hypoxia[line] = hyp_col

    control_samples = [ctr_normoxia[l] for l in GSC_LINES if l in ctr_normoxia] + \
                      [ctr_hypoxia[l] for l in GSC_LINES if l in ctr_hypoxia]

    expr = df[control_samples].copy()
    expr = expr[(expr > 0).any(axis=1)]
    expr_log = np.log1p(expr)
    expr_z = expr_log.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1)

    # Score components
    comp_scores = {}
    for comp, genes in FATE_COMPONENTS.items():
        present = [g for g in genes if g in expr_z.index]
        if len(present) == 0:
            comp_scores[comp] = pd.Series(0, index=control_samples)
        else:
            comp_scores[comp] = expr_z.loc[present].mean(axis=0)

    glut_present = [g for g in GLUTAMINE_ALL if g in expr_z.index]
    glut_score = expr_z.loc[glut_present].mean(axis=0)

    # Per-GSC-line Δ
    delta_records = []
    per_sample_records = []

    for sample in control_samples:
        row = {'Sample': sample}
        for comp in FATE_COMPONENTS:
            row[f'Score_{comp}'] = comp_scores[comp].get(sample, np.nan)
        row['Score_OverallGlutamine'] = glut_score.get(sample, np.nan)
        is_hypoxia = '_ctr_H' in sample
        line = sample.replace('_ctr_H', '').replace('_ctr', '')
        row['GSC_line'] = line
        row['Condition'] = 'hypoxia' if is_hypoxia else 'normoxia'
        per_sample_records.append(row)

    for line in GSC_LINES:
        norm_col = ctr_normoxia.get(line)
        hyp_col = ctr_hypoxia.get(line)
        if not norm_col or not hyp_col:
            continue
        row = {'GSC_line': line}
        for comp in FATE_COMPONENTS:
            delta = comp_scores[comp][hyp_col] - comp_scores[comp][norm_col]
            row[f'Delta_{comp}'] = delta
        delta_glut = glut_score[hyp_col] - glut_score[norm_col]
        row['Delta_OverallGlutamine'] = delta_glut
        delta_records.append(row)

    df_delta = pd.DataFrame(delta_records)
    df_samples = pd.DataFrame(per_sample_records)

    # ── Per-component: TWO tests ──
    test_components = list(FATE_COMPONENTS.keys()) + ['OverallGlutamine']
    summary_rows = []

    for comp in test_components:
        col = f'Delta_{comp}'
        deltas = df_delta[col].values
        n_pos = (deltas > 0).sum()
        n_neg = (deltas < 0).sum()
        n_total = len(deltas)
        median_delta = np.median(deltas)

        # Two-sided Wilcoxon signed-rank test
        try:
            w_stat, w_p_twosided = stats.wilcoxon(deltas)
        except ValueError:
            w_stat, w_p_twosided = np.nan, 1.0

        observed_direction = 'positive' if n_pos > n_neg else 'negative'
        predicted_direction = LIN_PREDICTED.get(comp, 'unknown')

        # ── Replication test ──
        if predicted_direction == 'positive':
            n_in_predicted = n_pos
        elif predicted_direction == 'negative':
            n_in_predicted = n_neg
        else:
            n_in_predicted = 0

        replication_direction_pass = n_in_predicted >= 6
        replication_pass = replication_direction_pass and (w_p_twosided < 0.05)

        # ── Acute-response discovery ──
        n_in_observed = max(n_pos, n_neg)
        acute_direction_pass = n_in_observed >= 6
        acute_pass = acute_direction_pass and (w_p_twosided < 0.05)

        # Directional concordance with Lin
        concordant = (observed_direction == predicted_direction) if predicted_direction != 'unknown' else None

        summary_rows.append({
            'Component': comp,
            'n_positive': n_pos,
            'n_negative': n_neg,
            'n_total': n_total,
            'median_delta': median_delta,
            'wilcoxon_p_twosided': w_p_twosided,
            'observed_direction': observed_direction,
            'lin_predicted_direction': predicted_direction,
            'concordant_with_lin': concordant,
            'n_in_predicted_direction': n_in_predicted,
            'replication_direction_pass': replication_direction_pass,
            'replication_pass': replication_pass,
            'n_in_observed_direction': n_in_observed,
            'acute_direction_pass': acute_direction_pass,
            'acute_response_pass': acute_pass,
        })

    df_summary = pd.DataFrame(summary_rows)

    # BH correction across 5 fate components (two-sided p)
    fate_mask = df_summary['Component'].isin(FATE_COMPONENTS.keys())
    df_summary['q_BH_twosided'] = np.nan
    df_summary.loc[fate_mask, 'q_BH_twosided'] = bh_correction(
        df_summary.loc[fate_mask, 'wilcoxon_p_twosided'].values
    )

    # Update pass flags with BH
    df_summary['replication_pass'] = (df_summary['replication_direction_pass']) & \
        (df_summary['wilcoxon_p_twosided'] < 0.05) & \
        ((df_summary['q_BH_twosided'] < 0.05) | (df_summary['q_BH_twosided'].isna()))
    df_summary['acute_response_pass'] = (df_summary['acute_direction_pass']) & \
        (df_summary['wilcoxon_p_twosided'] < 0.05) & \
        ((df_summary['q_BH_twosided'] < 0.05) | (df_summary['q_BH_twosided'].isna()))

    # ── Contrast: Synthesis vs Catabolism (two-sided) ──
    synth_d = df_delta['Delta_Synthesis'].values
    catab_d = df_delta['Delta_Catabolism'].values
    D = synth_d - catab_d
    n_D_pos = (D > 0).sum()
    median_D = np.median(D)
    try:
        d_stat, d_p_twosided = stats.wilcoxon(D)
    except ValueError:
        d_stat, d_p_twosided = np.nan, 1.0
    contrast_pass = (median_D > 0) and (n_D_pos >= 6) and (d_p_twosided < 0.05)

    df_contrast = pd.DataFrame({
        'GSC_line': df_delta['GSC_line'],
        'Delta_Synthesis': synth_d,
        'Delta_Catabolism': catab_d,
        'D': D,
        'D_direction': ['synthesis_favored' if d > 0 else 'catabolism_favored' for d in D],
    })

    # ── Verdict ──
    replication_passes = df_summary[fate_mask & df_summary['replication_pass']]['Component'].tolist()
    acute_passes = df_summary[fate_mask & df_summary['acute_response_pass']]['Component'].tolist()
    overall_replication = df_summary[df_summary['Component'] == 'OverallGlutamine']['replication_pass'].iloc[0]
    overall_acute = df_summary[df_summary['Component'] == 'OverallGlutamine']['acute_response_pass'].iloc[0]

    # Directional failures (opposite to predicted)
    directional_failures = df_summary[fate_mask & ~df_summary['concordant_with_lin']]['Component'].tolist()

    verdict = {
        'gate': 'B',
        'dataset': 'GSE243848',
        'n_gsc_lines': len(GSC_LINES),
        'design': 'paired normoxia vs hypoxia, DMSO control, 24h',
        'statistical_test': 'two-sided Wilcoxon signed-rank',
        'replication_test': {
            'description': 'Does GSE243848 replicate the Lin Gate A direction?',
            'criterion': '>=6/8 in Lin-predicted direction + two-sided Wilcoxon p < 0.05 + q_BH < 0.05',
            'components_passing': replication_passes,
            'n_passing': len(replication_passes),
            'overall_glutamine_pass': bool(overall_replication),
        },
        'acute_response_discovery': {
            'description': 'Is there a stable transcriptional change under 24h hypoxia?',
            'criterion': '>=6/8 in either direction + two-sided Wilcoxon p < 0.05 + q_BH < 0.05',
            'components_passing': acute_passes,
            'n_passing': len(acute_passes),
            'overall_glutamine_pass': bool(overall_acute),
        },
        'directional_failures': directional_failures,
        'directional_failures_note': 'These components showed opposite direction to Lin steady-state pattern',
        'contrast_pass': bool(contrast_pass),
        'contrast_n_D_positive': int(n_D_pos),
        'contrast_median_D': float(median_D),
        'contrast_wilcoxon_p_twosided': float(d_p_twosided),
    }

    # ── Save ──
    df_samples.to_csv(OUT / 'gateB_per_sample_scores.csv', index=False)
    df_delta.to_csv(OUT / 'gateB_per_line_delta.csv', index=False)
    df_summary.to_csv(OUT / 'gateB_summary.csv', index=False)
    df_contrast.to_csv(OUT / 'gateB_contrast.csv', index=False)
    with open(OUT / 'gateB_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    # ── Print ──
    print("\n" + "="*80)
    print("GATE B CORRECTED: GSE243848 PERTURBATION (8 GSC lines, DMSO control)")
    print("Statistical test: TWO-SIDED Wilcoxon signed-rank")
    print("="*80)

    print(f"\n{'Component':<22} {'n+':>3} {'n-':>3} {'median Δ':>10} {'p(2-sided)':>12} "
          f"{'Lin dir':>8} {'Obs dir':>8} {'Concord':>8} {'Repl':>6} {'Acute':>6}")
    print("-"*100)
    for _, row in df_summary.iterrows():
        conc = str(row['concordant_with_lin']) if row['concordant_with_lin'] is not None else 'N/A'
        print(f"{row['Component']:<22} {row['n_positive']:>3} {row['n_negative']:>3} "
              f"{row['median_delta']:>10.4f} {row['wilcoxon_p_twosided']:>12.4f} "
              f"{row['lin_predicted_direction']:>8} {row['observed_direction']:>8} "
              f"{conc:>8} "
              f"{'PASS' if row['replication_pass'] else 'fail':>6} "
              f"{'PASS' if row['acute_response_pass'] else 'fail':>6}")

    print(f"\nReplication test (vs Lin direction):")
    print(f"  Components passing: {replication_passes}")
    print(f"  Overall glutamine: {'PASS' if overall_replication else 'FAIL'}")

    print(f"\nAcute-response discovery (any direction):")
    print(f"  Components passing: {acute_passes}")
    print(f"  Overall glutamine: {'PASS' if overall_acute else 'FAIL'}")

    print(f"\nDirectional failures (opposite to Lin): {directional_failures}")

    print(f"\nContrast (Synthesis - Catabolism, two-sided):")
    print(f"  n_D_positive: {n_D_pos}/{len(D)}")
    print(f"  median_D: {median_D:.4f}")
    print(f"  Wilcoxon p (two-sided): {d_p_twosided:.4f}")
    print(f"  Contrast PASS: {contrast_pass}")
    print("="*80)

if __name__ == '__main__':
    main()
