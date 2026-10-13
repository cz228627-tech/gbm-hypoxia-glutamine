#!/usr/bin/env python3
"""
Gate C: GSE300765 U87MG supportive analysis.
Microarray (Illumina HumanHT-12 V4.0), hypoxia (1% O2, 48h) vs normoxia (21% O2, 48h).
3 replicates each, unpaired t-test.

5 pre-frozen fate components (same as Gate A/B).
"""

import warnings
warnings.filterwarnings('ignore')
import os

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import json
import gzip

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
DATA_FILE = Path(os.environ.get("GBM_U87_FILE", DATA_DIR / "gse300765" / "series_matrix.txt.gz")).expanduser()
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

# Probe-to-gene mapping (from GPL10558 annotation)
PROBE_MAP = {
    'ILMN_2188722': 'GLS',
    'ILMN_1677439': 'GLS2',
    'ILMN_1706070': 'GLS2',
    'ILMN_1713756': 'GLUD1',
    'ILMN_1653496': 'GLUL',
    'ILMN_1765208': 'GLUL',
    'ILMN_2269977': 'GLUL',
    'ILMN_1656145': 'GOT1',
    'ILMN_1801348': 'GOT2',
    'ILMN_1795257': 'GPT',
    'ILMN_1684158': 'GPT2',
    'ILMN_1707720': 'SLC1A5',
    'ILMN_1769911': 'SLC38A1',
    'ILMN_2266948': 'SLC38A1',
    'ILMN_1651799': 'SLC38A2',
    'ILMN_1655229': 'SLC7A11',
    'ILMN_1720373': 'SLC7A5',
}

# Sample mapping (from series matrix column order)
# Columns: GSM9068191-GSM9068208
# 13-15 (0-based 12-14): U87_hox_1/2/3 (hypoxia)
# 16-18 (0-based 15-17): U87_nox_1/2/3 (normoxia)
HYPOXIA_GSM = ['GSM9068203', 'GSM9068204', 'GSM9068205']
NORMOXIA_GSM = ['GSM9068206', 'GSM9068207', 'GSM9068208']

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

def parse_series_matrix(filepath):
    """Parse GEO series matrix file, return expression DataFrame."""
    with gzip.open(filepath, 'rt') as f:
        in_table = False
        data_lines = []
        col_names = None
        for line in f:
            line = line.strip()
            if line == '!series_matrix_table_begin':
                in_table = True
                continue
            if line == '!series_matrix_table_end':
                break
            if in_table:
                if col_names is None:
                    # Header row
                    col_names = [c.strip('"') for c in line.split('\t')]
                else:
                    data_lines.append(line)

    # Build DataFrame
    rows = []
    indices = []
    for line in data_lines:
        fields = [c.strip('"') for c in line.split('\t')]
        indices.append(fields[0])
        rows.append([float(v) if v else np.nan for v in fields[1:]])

    df = pd.DataFrame(rows, index=indices, columns=col_names[1:])
    return df

def main():
    print("Parsing GSE300765 series matrix...")
    expr = parse_series_matrix(DATA_FILE)
    print(f"  {expr.shape[0]} probes x {expr.shape[1]} samples")
    print(f"  Samples: {list(expr.columns)}")

    # Filter to probes we need
    fate_probes = [p for p in PROBE_MAP.keys() if p in expr.index]
    print(f"  Fate probes found: {len(fate_probes)}/{len(PROBE_MAP)}")

    expr_fate = expr.loc[fate_probes].copy()
    # Map probe to gene
    expr_fate['Gene'] = expr_fate.index.map(PROBE_MAP)

    # Average across probes for the same gene
    hypoxia_samples = HYPOXIA_GSM
    normoxia_samples = NORMOXIA_GSM
    all_samples = hypoxia_samples + normoxia_samples

    # Average probes per gene
    expr_gene = expr_fate.groupby('Gene')[all_samples].mean()
    print(f"  Gene-level matrix: {expr_gene.shape[0]} genes x {len(all_samples)} samples")
    print(f"  Genes: {list(expr_gene.index)}")

    # Z-score each gene across the 6 samples
    expr_z = expr_gene.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1)

    # Score each fate component: mean of z-scored expression
    comp_scores = {}
    for comp, genes in FATE_COMPONENTS.items():
        present = [g for g in genes if g in expr_z.index]
        missing = set(genes) - set(present)
        if missing:
            print(f"  WARNING: {comp} missing: {missing}")
        if len(present) == 0:
            comp_scores[comp] = pd.Series(0, index=all_samples)
        else:
            comp_scores[comp] = expr_z.loc[present].mean(axis=0)
        print(f"  {comp}: {len(present)}/{len(genes)} genes")

    # Overall glutamine
    glut_present = [g for g in ALL_FATE_GENES if g in expr_z.index]
    glut_score = expr_z.loc[glut_present].mean(axis=0)

    # ── Per-sample scores ──
    per_sample = []
    for sample in all_samples:
        row = {'Sample': sample}
        row['Condition'] = 'hypoxia' if sample in hypoxia_samples else 'normoxia'
        for comp in FATE_COMPONENTS:
            row[f'Score_{comp}'] = comp_scores[comp].get(sample, np.nan)
        row['Score_OverallGlutamine'] = glut_score.get(sample, np.nan)
        per_sample.append(row)
    df_samples = pd.DataFrame(per_sample)

    # ── Per-component t-test (3 vs 3) ──
    summary_rows = []
    test_components = list(FATE_COMPONENTS.keys()) + ['OverallGlutamine']

    for comp in test_components:
        if comp == 'OverallGlutamine':
            scores = glut_score
        else:
            scores = comp_scores[comp]

        hyp_vals = scores[hypoxia_samples].values
        norm_vals = scores[normoxia_samples].values

        mean_diff = np.mean(hyp_vals) - np.mean(norm_vals)
        t_stat, t_p = stats.ttest_ind(hyp_vals, norm_vals, alternative='greater')

        # Direction: positive if hypoxia > normoxia
        direction = 'positive' if mean_diff > 0 else 'negative'

        summary_rows.append({
            'Component': comp,
            'mean_hypoxia': np.mean(hyp_vals),
            'mean_normoxia': np.mean(norm_vals),
            'mean_diff': mean_diff,
            'direction': direction,
            't_stat': t_stat,
            't_p': t_p,
        })

    df_summary = pd.DataFrame(summary_rows)

    # BH correction across 5 fate components
    fate_mask = df_summary['Component'].isin(FATE_COMPONENTS.keys())
    df_summary['q_BH'] = np.nan
    df_summary.loc[fate_mask, 'q_BH'] = bh_correction(
        df_summary.loc[fate_mask, 't_p'].values
    )

    # Gate: direction positive + p < 0.05 (for supportive, we report but don't gate)
    df_summary['direction_positive'] = df_summary['mean_diff'] > 0
    df_summary['significant'] = (df_summary['t_p'] < 0.05) | (df_summary['q_BH'] < 0.05)

    # ── Contrast: Synthesis vs Catabolism ──
    synth_hyp = comp_scores['Synthesis'][hypoxia_samples].values
    synth_norm = comp_scores['Synthesis'][normoxia_samples].values
    catab_hyp = comp_scores['Catabolism'][hypoxia_samples].values
    catab_norm = comp_scores['Catabolism'][normoxia_samples].values

    synth_diff = np.mean(synth_hyp) - np.mean(synth_norm)
    catab_diff = np.mean(catab_hyp) - np.mean(catab_norm)
    D = synth_diff - catab_diff

    # For contrast with n=3, we can't do a paired test. Report the difference.
    contrast = {
        'synthesis_diff': synth_diff,
        'catabolism_diff': catab_diff,
        'D': D,
        'synthesis_favored': D > 0,
    }

    # ── Verdict ──
    fate_positive = df_summary[fate_mask & df_summary['direction_positive']]['Component'].tolist()
    fate_significant = df_summary[fate_mask & df_summary['significant']]['Component'].tolist()
    overall_positive = df_summary[df_summary['Component'] == 'OverallGlutamine']['direction_positive'].iloc[0]
    overall_significant = df_summary[df_summary['Component'] == 'OverallGlutamine']['significant'].iloc[0]

    # Direction consistency with Gate B
    gate_b_directions = {
        'Synthesis': 'negative',  # 4/8 positive → negative
        'Catabolism': 'negative',  # 3/8 positive → negative
        'Direct_uptake': 'positive',  # 7/8 positive
        'Exchange_redox': 'positive',  # 6/8 positive
        'CN_redistribution': 'positive',  # 8/8 positive
        'OverallGlutamine': 'positive',  # 7/8 positive
    }

    consistency = {}
    for comp in test_components:
        u87_dir = df_summary[df_summary['Component'] == comp]['direction'].iloc[0]
        gsc_dir = gate_b_directions.get(comp, 'unknown')
        consistency[comp] = {
            'u87_direction': u87_dir,
            'gsc_direction': gsc_dir,
            'concordant': u87_dir == gsc_dir,
        }

    verdict = {
        'gate': 'C',
        'dataset': 'GSE300765',
        'cell_line': 'U87MG',
        'design': 'hypoxia (1% O2, 48h) vs normoxia (21% O2, 48h), 3 replicates each',
        'platform': 'Illumina HumanHT-12 V4.0 (microarray)',
        'fate_components_positive': fate_positive,
        'fate_components_significant': fate_significant,
        'overall_glutamine_positive': bool(overall_positive),
        'overall_glutamine_significant': bool(overall_significant),
        'contrast_D': float(D),
        'contrast_synthesis_favored': bool(D > 0),
        'direction_consistency_with_gateB': consistency,
        'note': 'Supportive analysis only. n=3 per group, limited power. '
                'Reported for direction consistency, not formal gating.'
    }

    # ── Save ──
    df_samples.to_csv(OUT / 'gateC_u87_scores.csv', index=False)
    df_summary.to_csv(OUT / 'gateC_u87_summary.csv', index=False)
    with open(OUT / 'gateC_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    # ── Print ──
    print("\n" + "="*70)
    print("GATE C: GSE300765 U87MG SUPPORTIVE ANALYSIS")
    print("="*70)
    print(f"\nPer-component results (unpaired t-test, 3v3):")
    print(df_summary[['Component','mean_diff','direction','t_stat','t_p','q_BH',
                       'direction_positive','significant']].to_string(index=False))
    print(f"\nContrast (Synthesis - Catabolism):")
    print(f"  Synthesis diff: {synth_diff:.4f}")
    print(f"  Catabolism diff: {catab_diff:.4f}")
    print(f"  D = {D:.4f} (synthesis-favored: {D > 0})")
    print(f"\nDirection consistency with Gate B (GSE243848):")
    for comp, info in consistency.items():
        print(f"  {comp}: U87={info['u87_direction']}, GSC={info['gsc_direction']}, "
              f"concordant={info['concordant']}")
    print("="*70)

if __name__ == '__main__':
    main()
