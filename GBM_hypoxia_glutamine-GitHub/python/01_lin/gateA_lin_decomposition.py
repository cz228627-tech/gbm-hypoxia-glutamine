#!/usr/bin/env python3
"""
Gate A: Lin scRNA glutamine fate decomposition (28 patients).

5 pre-frozen fate components:
  1. Synthesis:              GLUL
  2. Catabolism:             GLS, GLS2
  3. Direct uptake:          SLC1A5, SLC38A1, SLC38A2
  4. Exchange/redox handling: SLC7A5, SLC7A11
  5. C/N redistribution:    GLUD1, GOT1, GOT2, GPT, GPT2

Per-component PASS: q_BH < 0.05 AND >=19/28 patients same direction.
Contrast (synthesis vs catabolism): D_i = rho(Synth) - rho(Catab),
  claim "synthesis-favored" only if median(D)>0, >=19/28 D>0, paired Wilcoxon p<0.05.
Single-gene transparency: per-patient rho for all 13 genes.
"""

import warnings
warnings.filterwarnings('ignore')
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from pathlib import Path
import json

# ── Paths ──
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
H5AD = Path(os.environ.get("GBM_LIN_H5AD", DATA_DIR / "lin" / "gbm_integrated.h5ad")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

# ── Hypoxia signature (16 genes, frozen) ──
HYPOXIA_GENES = ['CA9','VEGFA','NDRG1','P4HA1','EGLN3','BNIP3','ADM','ANKRD37',
                 'CDKN1A','GYS1','DDIT3','F3','JMY','LRP1','LOX','ANGPTL4']

# ── 5 fate components (frozen) ──
FATE_COMPONENTS = {
    'Synthesis':              ['GLUL'],
    'Catabolism':             ['GLS', 'GLS2'],
    'Direct_uptake':          ['SLC1A5', 'SLC38A1', 'SLC38A2'],
    'Exchange_redox':         ['SLC7A5', 'SLC7A11'],
    'CN_redistribution':      ['GLUD1', 'GOT1', 'GOT2', 'GPT', 'GPT2'],
}

# All 13 genes for single-gene transparency
ALL_FATE_GENES = sorted(set(g for gs in FATE_COMPONENTS.values() for g in gs))

# ── Neoplastic annotations ──
NEOPLASTIC = ['AC-like', 'MES-like', 'OPC-like', 'NPC-like']

# ── Eligibility: >=50 neoplastic cells, IDH-WT ──
MIN_CELLS = 50

def bh_correction(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    bh = ranked * n / np.arange(1, n + 1)
    # enforce monotonicity from the end
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    result = np.empty(n)
    result[order] = bh
    return result

def main():
    print("Loading Lin h5ad...")
    adata = sc.read_h5ad(H5AD)
    adata.var_names = adata.var['features'].astype(str).values
    adata.var_names_make_unique()
    print(f"  {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Filter to neoplastic cells
    mask_neo = adata.obs['annotation'].isin(NEOPLASTIC)
    adata = adata[mask_neo].copy()
    print(f"  Neoplastic: {adata.shape[0]} cells")

    # Eligible patients: >=50 neoplastic cells
    counts = adata.obs['Patient'].value_counts()
    eligible = counts[counts >= MIN_CELLS].index.tolist()
    # IDH-WT filter: Lin dataset is IDH-WT GBM, but check if there's an IDH column
    if 'IDH_status' in adata.obs.columns:
        idh_wt = adata.obs[adata.obs['IDH_status'] == 'WT']['Patient'].unique()
        eligible = [p for p in eligible if p in idh_wt]
    eligible = sorted(eligible)
    print(f"  Eligible patients: {len(eligible)}")

    adata = adata[adata.obs['Patient'].isin(eligible)].copy()

    # Score hypoxia
    hypoxia_present = [g for g in HYPOXIA_GENES if g in adata.var_names]
    print(f"  Hypoxia genes present: {len(hypoxia_present)}/16")
    sc.tl.score_genes(adata, gene_list=hypoxia_present, score_name='HypoxiaScore',
                      use_raw=False, random_state=42)

    # Score each fate component
    for comp, genes in FATE_COMPONENTS.items():
        present = [g for g in genes if g in adata.var_names]
        missing = set(genes) - set(present)
        if missing:
            print(f"  WARNING: {comp} missing genes: {missing}")
        if len(present) == 1:
            # Single gene: use expression directly
            g = present[0]
            adata.obs[f'FateScore_{comp}'] = adata[:, g].X.toarray().ravel()
            print(f"  {comp}: single gene {g} (direct expression)")
        else:
            sc.tl.score_genes(adata, gene_list=present, score_name=f'FateScore_{comp}',
                              use_raw=False, random_state=42)
            print(f"  {comp}: {len(present)} genes (score_genes)")

    # ── Per-patient Spearman ρ ──
    records_comp = []
    records_gene = []

    for pat in eligible:
        pat_mask = adata.obs['Patient'] == pat
        pat_data = adata[pat_mask]

        hyp = pat_data.obs['HypoxiaScore'].values

        # Component-level
        for comp in FATE_COMPONENTS:
            fate = pat_data.obs[f'FateScore_{comp}'].values
            rho, p = stats.spearmanr(hyp, fate)
            records_comp.append({
                'Patient': pat,
                'Component': comp,
                'rho': rho,
                'p': p,
                'direction': 'positive' if rho > 0 else 'negative',
                'n_cells': pat_mask.sum()
            })

        # Single-gene level
        for gene in ALL_FATE_GENES:
            if gene in adata.var_names:
                expr = pat_data[:, gene].X.toarray().ravel()
                rho, p = stats.spearmanr(hyp, expr)
                records_gene.append({
                    'Patient': pat,
                    'Gene': gene,
                    'rho': rho,
                    'p': p,
                    'direction': 'positive' if rho > 0 else 'negative',
                })

    df_comp = pd.DataFrame(records_comp)
    df_gene = pd.DataFrame(records_gene)

    # ── Per-component summary ──
    summary_rows = []
    for comp in FATE_COMPONENTS:
        sub = df_comp[df_comp['Component'] == comp]
        n_pos = (sub['rho'] > 0).sum()
        n_neg = (sub['rho'] < 0).sum()
        n_total = len(sub)
        median_rho = sub['rho'].median()
        # One-sample Wilcoxon signed-rank test on rho values vs 0
        try:
            w_stat, w_p = stats.wilcoxon(sub['rho'].values, alternative='greater')
        except ValueError:
            # All zeros or insufficient data
            w_stat, w_p = np.nan, 1.0
        summary_rows.append({
            'Component': comp,
            'n_positive': n_pos,
            'n_negative': n_neg,
            'n_total': n_total,
            'threshold': int(np.ceil(2/3 * n_total)),
            'direction_pass': n_pos >= int(np.ceil(2/3 * n_total)),
            'median_rho': median_rho,
            'wilcoxon_stat': w_stat,
            'wilcoxon_p': w_p,
        })

    df_summary = pd.DataFrame(summary_rows)
    # BH correction across 5 components
    df_summary['q_BH'] = bh_correction(df_summary['wilcoxon_p'].values)
    df_summary['gate_pass'] = (df_summary['q_BH'] < 0.05) & (df_summary['direction_pass'])

    # ── Contrast: Synthesis vs Catabolism ──
    synth = df_comp[df_comp['Component'] == 'Synthesis'].set_index('Patient')['rho']
    catab = df_comp[df_comp['Component'] == 'Catabolism'].set_index('Patient')['rho']
    D = synth - catab
    n_D_pos = (D > 0).sum()
    median_D = D.median()
    try:
        d_stat, d_p = stats.wilcoxon(D.values, alternative='greater')
    except ValueError:
        d_stat, d_p = np.nan, 1.0

    contrast_pass = (median_D > 0) and (n_D_pos >= 19) and (d_p < 0.05)

    df_contrast = pd.DataFrame({
        'Patient': D.index,
        'rho_Synthesis': synth.values,
        'rho_Catabolism': catab.values,
        'D': D.values,
        'D_direction': ['synthesis_favored' if d > 0 else 'catabolism_favored' for d in D.values],
    })

    # ── Verdict ──
    n_components_pass = df_summary['gate_pass'].sum()
    any_component_pass = n_components_pass > 0

    if contrast_pass:
        gate_a_label = "synthesis_favored_transcriptional_shift"
    elif any_component_pass:
        gate_a_label = "glutamine_associated_state"
    else:
        gate_a_label = "no_stable_component"

    verdict = {
        'gate': 'A',
        'n_patients': len(eligible),
        'n_components_pass': int(n_components_pass),
        'components_pass': df_summary[df_summary['gate_pass']]['Component'].tolist(),
        'contrast_pass': bool(contrast_pass),
        'contrast_n_D_positive': int(n_D_pos),
        'contrast_median_D': float(median_D),
        'contrast_wilcoxon_p': float(d_p),
        'gate_a_label': gate_a_label,
        'proceed_to_gateB': True,  # always proceed
        'note': 'Gate A is an explanation gate. Always proceed to Gate B. '
                'Label determines title only.'
    }

    # ── Save outputs ──
    df_comp.to_csv(OUT / 'gateA_per_patient_rho.csv', index=False)
    df_gene.to_csv(OUT / 'gateA_per_gene_rho.csv', index=False)
    df_summary.to_csv(OUT / 'gateA_component_summary.csv', index=False)
    df_contrast.to_csv(OUT / 'gateA_contrast.csv', index=False)
    with open(OUT / 'gateA_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    # ── Print summary ──
    print("\n" + "="*70)
    print("GATE A: LIN scRNA GLUTAMINE FATE DECOMPOSITION")
    print("="*70)
    print(f"Patients: {len(eligible)}")
    print(f"\nPer-component results:")
    print(df_summary[['Component','n_positive','n_negative','threshold',
                       'direction_pass','median_rho','wilcoxon_p','q_BH','gate_pass']].to_string(index=False))
    print(f"\nContrast (Synthesis - Catabolism):")
    print(f"  n_D_positive: {n_D_pos}/{len(D)}")
    print(f"  median_D: {median_D:.4f}")
    print(f"  Wilcoxon p: {d_p:.4f}")
    print(f"  Contrast PASS: {contrast_pass}")
    print(f"\nGate A label: {gate_a_label}")
    print(f"Components passing: {df_summary[df_summary['gate_pass']]['Component'].tolist()}")
    print(f"\nAlways proceeding to Gate B.")
    print("="*70)

    # Also print single-gene summary
    print("\nSingle-gene per-patient direction counts:")
    gene_summary = df_gene.groupby('Gene').apply(
        lambda x: pd.Series({
            'n_positive': (x['rho'] > 0).sum(),
            'n_negative': (x['rho'] < 0).sum(),
            'median_rho': x['rho'].median(),
        })
    ).sort_values('median_rho', ascending=False)
    print(gene_summary.to_string())

if __name__ == '__main__':
    main()
