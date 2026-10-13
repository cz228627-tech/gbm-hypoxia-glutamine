#!/usr/bin/env python3
"""
Recompute Gate A, A-spatial, and C p-values with TWO-SIDED tests.
Reads saved per-patient/per-section ρ values, recomputes Wilcoxon/t-test
with alternative='two-sided', updates q_BH and verdicts.
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
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()

FATE_COMPONENTS = [
    'Synthesis', 'Catabolism', 'Direct_uptake', 'Exchange_redox', 'CN_redistribution'
]

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

def recompute_gate_a():
    """Recompute Gate A (Lin scRNA) with two-sided Wilcoxon."""
    print("="*70)
    print("RECOMPUTING GATE A (Lin scRNA) — TWO-SIDED")
    print("="*70)

    df = pd.read_csv(OUT / "gateA_per_patient_rho.csv")
    threshold = 19  # 2/3 × 28 = 18.67 → 19
    n_patients = df['Patient'].nunique()

    summary_rows = []
    for comp in FATE_COMPONENTS:
        sub = df[df['Component'] == comp]
        rhos = sub['rho'].values
        n_pos = (rhos > 0).sum()
        n_neg = (rhos < 0).sum()
        median_rho = np.median(rhos)

        # Two-sided Wilcoxon signed-rank test
        try:
            w_stat, w_p = stats.wilcoxon(rhos)
        except ValueError:
            w_stat, w_p = np.nan, 1.0

        direction_pass = n_pos >= threshold
        summary_rows.append({
            'Component': comp,
            'n_positive': n_pos,
            'n_negative': n_neg,
            'n_total': len(rhos),
            'threshold': threshold,
            'direction_pass': direction_pass,
            'median_rho': median_rho,
            'wilcoxon_stat': w_stat,
            'wilcoxon_p_twosided': w_p,
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary['q_BH_twosided'] = bh_correction(df_summary['wilcoxon_p_twosided'].values)
    df_summary['gate_pass'] = (df_summary['q_BH_twosided'] < 0.05) & (df_summary['direction_pass'])

    # Contrast (two-sided)
    synth = df[df['Component'] == 'Synthesis'].set_index('Patient')['rho']
    catab = df[df['Component'] == 'Catabolism'].set_index('Patient')['rho']
    D = synth - catab
    n_D_pos = (D > 0).sum()
    median_D = D.median()
    try:
        d_stat, d_p = stats.wilcoxon(D.values)
    except ValueError:
        d_stat, d_p = np.nan, 1.0
    contrast_pass = (median_D > 0) and (n_D_pos >= threshold) and (d_p < 0.05)

    n_pass = df_summary['gate_pass'].sum()
    if contrast_pass:
        label = "synthesis_favored_transcriptional_shift"
    elif n_pass > 0:
        label = "glutamine_associated_state"
    else:
        label = "no_stable_component"

    verdict = {
        'gate': 'A',
        'statistical_test': 'two-sided Wilcoxon signed-rank',
        'n_patients': n_patients,
        'n_components_pass': int(n_pass),
        'components_pass': df_summary[df_summary['gate_pass']]['Component'].tolist(),
        'contrast_pass': bool(contrast_pass),
        'contrast_n_D_positive': int(n_D_pos),
        'contrast_median_D': float(median_D),
        'contrast_wilcoxon_p_twosided': float(d_p),
        'gate_a_label': label,
        'proceed_to_gateB': True,
    }

    df_summary.to_csv(OUT / 'gateA_component_summary.csv', index=False)
    df_contrast = pd.DataFrame({
        'Patient': D.index,
        'rho_Synthesis': synth.values,
        'rho_Catabolism': catab.values,
        'D': D.values,
        'D_direction': ['synthesis_favored' if d > 0 else 'catabolism_favored' for d in D.values],
    })
    df_contrast.to_csv(OUT / 'gateA_contrast.csv', index=False)
    with open(OUT / 'gateA_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    print(f"\n{'Component':<22} {'n+':>3} {'n-':>3} {'median ρ':>10} {'p(2-sided)':>12} {'q_BH':>10} {'PASS':>6}")
    print("-"*70)
    for _, row in df_summary.iterrows():
        print(f"{row['Component']:<22} {row['n_positive']:>3} {row['n_negative']:>3} "
              f"{row['median_rho']:>10.4f} {row['wilcoxon_p_twosided']:>12.6f} "
              f"{row['q_BH_twosided']:>10.6f} {'PASS' if row['gate_pass'] else 'fail':>6}")
    print(f"\nContrast: n_D+={n_D_pos}/{len(D)}, median D={median_D:.4f}, p(2-sided)={d_p:.4f}")
    print(f"Contrast PASS: {contrast_pass}")
    print(f"Label: {label}")

def recompute_gate_a_spatial():
    """Recompute Gate A-spatial (Visium) with two-sided Wilcoxon."""
    print("\n" + "="*70)
    print("RECOMPUTING GATE A-SPATIAL (Visium) — TWO-SIDED")
    print("="*70)

    df = pd.read_csv(OUT / "gateA_spatial_patient_rho.csv")
    threshold = 9  # 2/3 × 13 = 8.67 → 9
    n_patients = df['Patient'].nunique()

    summary_rows = []
    for comp in FATE_COMPONENTS:
        sub = df[df['Component'] == comp]
        rhos = sub['rho'].values
        n_pos = (rhos > 0).sum()
        n_neg = (rhos < 0).sum()
        median_rho = np.median(rhos)

        try:
            w_stat, w_p = stats.wilcoxon(rhos)
        except ValueError:
            w_stat, w_p = np.nan, 1.0

        direction_pass = n_pos >= threshold
        summary_rows.append({
            'Component': comp,
            'n_positive': n_pos,
            'n_negative': n_neg,
            'n_total': len(rhos),
            'threshold': threshold,
            'direction_pass': direction_pass,
            'median_rho': median_rho,
            'wilcoxon_stat': w_stat,
            'wilcoxon_p_twosided': w_p,
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary['q_BH_twosided'] = bh_correction(df_summary['wilcoxon_p_twosided'].values)
    df_summary['gate_pass'] = (df_summary['q_BH_twosided'] < 0.05) & (df_summary['direction_pass'])

    # Contrast (two-sided)
    synth = df[df['Component'] == 'Synthesis'].set_index('Patient')['rho']
    catab = df[df['Component'] == 'Catabolism'].set_index('Patient')['rho']
    common = synth.index.intersection(catab.index)
    synth = synth.loc[common]
    catab = catab.loc[common]
    D = synth - catab
    n_D_pos = (D > 0).sum()
    median_D = D.median()
    try:
        d_stat, d_p = stats.wilcoxon(D.values)
    except ValueError:
        d_stat, d_p = np.nan, 1.0
    contrast_pass = (median_D > 0) and (n_D_pos >= threshold) and (d_p < 0.05)

    n_pass = df_summary['gate_pass'].sum()
    if contrast_pass:
        label = "synthesis_favored_transcriptional_shift"
    elif n_pass > 0:
        label = "glutamine_associated_state"
    else:
        label = "no_stable_component"

    verdict = {
        'gate': 'A_spatial',
        'statistical_test': 'two-sided Wilcoxon signed-rank',
        'n_patients': n_patients,
        'threshold': threshold,
        'n_components_pass': int(n_pass),
        'components_pass': df_summary[df_summary['gate_pass']]['Component'].tolist(),
        'contrast_pass': bool(contrast_pass),
        'contrast_n_D_positive': int(n_D_pos),
        'contrast_median_D': float(median_D),
        'contrast_wilcoxon_p_twosided': float(d_p),
        'gate_a_spatial_label': label,
    }

    df_summary.to_csv(OUT / 'gateA_spatial_summary.csv', index=False)
    df_contrast = pd.DataFrame({
        'Patient': D.index,
        'rho_Synthesis': synth.values,
        'rho_Catabolism': catab.values,
        'D': D.values,
        'D_direction': ['synthesis_favored' if d > 0 else 'catabolism_favored' for d in D.values],
    })
    df_contrast.to_csv(OUT / 'gateA_spatial_contrast.csv', index=False)
    with open(OUT / 'gateA_spatial_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    print(f"\n{'Component':<22} {'n+':>3} {'n-':>3} {'median ρ':>10} {'p(2-sided)':>12} {'q_BH':>10} {'PASS':>6}")
    print("-"*70)
    for _, row in df_summary.iterrows():
        print(f"{row['Component']:<22} {row['n_positive']:>3} {row['n_negative']:>3} "
              f"{row['median_rho']:>10.4f} {row['wilcoxon_p_twosided']:>12.6f} "
              f"{row['q_BH_twosided']:>10.6f} {'PASS' if row['gate_pass'] else 'fail':>6}")
    print(f"\nContrast: n_D+={n_D_pos}/{len(D)}, median D={median_D:.4f}, p(2-sided)={d_p:.4f}")
    print(f"Contrast PASS: {contrast_pass}")
    print(f"Label: {label}")

def recompute_gate_c():
    """Recompute Gate C (U87) with two-sided t-test."""
    print("\n" + "="*70)
    print("RECOMPUTING GATE C (U87) — TWO-SIDED")
    print("="*70)

    df = pd.read_csv(OUT / "gateC_u87_scores.csv")

    HYPOXIA = ['GSM9068203', 'GSM9068204', 'GSM9068205']
    NORMOXIA = ['GSM9068206', 'GSM9068207', 'GSM9068208']

    test_components = FATE_COMPONENTS + ['OverallGlutamine']
    summary_rows = []

    for comp in test_components:
        col = f'Score_{comp}'
        hyp_vals = df[df['Condition'] == 'hypoxia'][col].values
        norm_vals = df[df['Condition'] == 'normoxia'][col].values

        mean_diff = np.mean(hyp_vals) - np.mean(norm_vals)
        # Two-sided t-test
        t_stat, t_p = stats.ttest_ind(hyp_vals, norm_vals)

        direction = 'positive' if mean_diff > 0 else 'negative'

        summary_rows.append({
            'Component': comp,
            'mean_hypoxia': np.mean(hyp_vals),
            'mean_normoxia': np.mean(norm_vals),
            'mean_diff': mean_diff,
            'direction': direction,
            't_stat': t_stat,
            't_p_twosided': t_p,
        })

    df_summary = pd.DataFrame(summary_rows)
    fate_mask = df_summary['Component'].isin(FATE_COMPONENTS)
    df_summary['q_BH_twosided'] = np.nan
    df_summary.loc[fate_mask, 'q_BH_twosided'] = bh_correction(
        df_summary.loc[fate_mask, 't_p_twosided'].values
    )
    df_summary['direction_positive'] = df_summary['mean_diff'] > 0
    df_summary['significant'] = (df_summary['t_p_twosided'] < 0.05) | (df_summary['q_BH_twosided'] < 0.05)

    # Contrast
    synth_hyp = df[df['Condition'] == 'hypoxia']['Score_Synthesis'].values
    synth_norm = df[df['Condition'] == 'normoxia']['Score_Synthesis'].values
    catab_hyp = df[df['Condition'] == 'hypoxia']['Score_Catabolism'].values
    catab_norm = df[df['Condition'] == 'normoxia']['Score_Catabolism'].values
    synth_diff = np.mean(synth_hyp) - np.mean(synth_norm)
    catab_diff = np.mean(catab_hyp) - np.mean(catab_norm)
    D = synth_diff - catab_diff

    fate_positive = df_summary[fate_mask & df_summary['direction_positive']]['Component'].tolist()
    fate_significant = df_summary[fate_mask & df_summary['significant']]['Component'].tolist()
    overall_positive = df_summary[df_summary['Component'] == 'OverallGlutamine']['direction_positive'].iloc[0]
    overall_significant = df_summary[df_summary['Component'] == 'OverallGlutamine']['significant'].iloc[0]

    verdict = {
        'gate': 'C',
        'statistical_test': 'two-sided unpaired t-test',
        'dataset': 'GSE300765',
        'cell_line': 'U87MG',
        'design': 'hypoxia (1% O2, 48h) vs normoxia (21% O2, 48h), 3 replicates each',
        'fate_components_positive': fate_positive,
        'fate_components_significant': fate_significant,
        'overall_glutamine_positive': bool(overall_positive),
        'overall_glutamine_significant': bool(overall_significant),
        'contrast_D': float(D),
        'contrast_synthesis_favored': bool(D > 0),
        'note': 'Supportive analysis only. n=3 per group, limited power. '
                'Not reproduced in this established cell-line model.'
    }

    df_summary.to_csv(OUT / 'gateC_u87_summary.csv', index=False)
    with open(OUT / 'gateC_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    print(f"\n{'Component':<22} {'mean diff':>10} {'direction':>10} {'t_p(2-sided)':>12} {'q_BH':>10} {'sig':>5}")
    print("-"*75)
    for _, row in df_summary.iterrows():
        q = f"{row['q_BH_twosided']:.4f}" if not np.isnan(row['q_BH_twosided']) else "N/A"
        print(f"{row['Component']:<22} {row['mean_diff']:>10.4f} {row['direction']:>10} "
              f"{row['t_p_twosided']:>12.4f} {q:>10} {'yes' if row['significant'] else 'no':>5}")
    print(f"\nContrast D = {D:.4f} (synthesis-favored: {D > 0})")
    print(f"Overall glutamine: positive={overall_positive}, significant={overall_significant}")

if __name__ == '__main__':
    recompute_gate_a()
    recompute_gate_a_spatial()
    recompute_gate_c()
