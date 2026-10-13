#!/usr/bin/env python3
"""
Gate A-spatial: Greenwald Visium glutamine fate decomposition.
19 sections from 13 IDH-WT GBM patients.

5 pre-frozen fate components (same as Gate A):
  1. Synthesis:              GLUL
  2. Catabolism:             GLS, GLS2
  3. Direct uptake:          SLC1A5, SLC38A1, SLC38A2
  4. Exchange/redox handling: SLC7A5, SLC7A11
  5. C/N redistribution:    GLUD1, GOT1, GOT2, GPT, GPT2

Per-component PASS: q_BH < 0.05 AND >=9/13 patients same direction.
Contrast: D_i = rho(Synth) - rho(Catab), median>0, >=9/13, paired Wilcoxon p<0.05.
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

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
VISIUM_DIR = Path(os.environ.get("GBM_GREENWALD_DIR", DATA_DIR / "greenwald_visium")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
OUT.mkdir(parents=True, exist_ok=True)

HYPOXIA_GENES = ['CA9','VEGFA','NDRG1','P4HA1','EGLN3','BNIP3','ADM','ANKRD37',
                 'CDKN1A','GYS1','DDIT3','F3','JMY','LRP1','LOX','ANGPTL4']

FATE_COMPONENTS = {
    'Synthesis':              ['GLUL'],
    'Catabolism':             ['GLS', 'GLS2'],
    'Direct_uptake':          ['SLC1A5', 'SLC38A1', 'SLC38A2'],
    'Exchange_redox':         ['SLC7A5', 'SLC7A11'],
    'CN_redistribution':      ['GLUD1', 'GOT1', 'GOT2', 'GPT', 'GPT2'],
}

ALL_FATE_GENES = sorted(set(g for gs in FATE_COMPONENTS.values() for g in gs))

SAMPLE_TO_PATIENT = {
    'GSM7596587_mgh258': 'mgh258',
    'GSM7596588_zh881inf': 'zh881',
    'GSM7596589_zh881t1': 'zh881',
    'GSM7596590_zh916bulk': 'zh916',
    'GSM7596591_zh916inf': 'zh916',
    'GSM7596592_zh916t1': 'zh916',
    'GSM7596593_zh1007inf': 'zh1007',
    'GSM7596594_zh1007nec': 'zh1007',
    'GSM7596595_zh1019inf': 'zh1019',
    'GSM7596596_zh1019t1': 'zh1019',
    'GSM7596597_zh8811a': 'zh8811',
    'GSM7596598_zh8811b': 'zh8811',
    'GSM7596599_zh8812': 'zh8812',
    'GSM7596600_bwh23': 'bwh23',
    'GSM7596601_bwh24': 'bwh24',
    'GSM7596602_bwh25': 'bwh25',
    'GSM7596603_bwh28': 'bwh28',
    'GSM7596604_bwh35': 'bwh35',
    'GSM7596605_mgh259': 'mgh259',
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

def load_visium(sample_key):
    """Load a Visium section from h5 file + tissue positions."""
    h5_path = VISIUM_DIR / f"{sample_key}_filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Load tissue positions
    tissue_file = VISIUM_DIR / f"{sample_key}_tissue_positions_list.csv.gz"
    tissue_pos = pd.read_csv(tissue_file, header=None,
                             names=['barcode','in_tissue','array_row','array_col',
                                    'pxl_row','pxl_col'], compression='gzip')
    tissue_pos = tissue_pos[tissue_pos['in_tissue'] == 1]
    tissue_pos.set_index('barcode', inplace=True)

    # Filter to in-tissue spots
    adata = adata[adata.obs_names.isin(tissue_pos.index)].copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def main():
    section_records = []
    gene_records = []

    for sample_key, patient in SAMPLE_TO_PATIENT.items():
        h5_path = VISIUM_DIR / f"{sample_key}_filtered_feature_bc_matrix.h5"
        if not h5_path.exists():
            print(f"  SKIP {sample_key}: file not found")
            continue

        print(f"  Loading {sample_key} (patient {patient})...")
        try:
            adata = load_visium(sample_key)
            print(f"    {adata.shape[0]} spots x {adata.shape[1]} genes")
        except Exception as e:
            print(f"    ERROR: {e}")
            continue

        # Score hypoxia
        hyp_present = [g for g in HYPOXIA_GENES if g in adata.var_names]
        if len(hyp_present) < 10:
            print(f"    WARNING: only {len(hyp_present)}/16 hypoxia genes present")
        sc.tl.score_genes(adata, gene_list=hyp_present, score_name='HypoxiaScore',
                          use_raw=False, random_state=42)

        # Score fate components
        for comp, genes in FATE_COMPONENTS.items():
            present = [g for g in genes if g in adata.var_names]
            missing = set(genes) - set(present)
            if missing:
                print(f"    {comp} missing: {missing}")
            if len(present) == 0:
                print(f"    {comp}: NO genes present, skipping")
                continue
            if len(present) == 1:
                g = present[0]
                adata.obs[f'FateScore_{comp}'] = adata[:, g].X.toarray().ravel()
            else:
                sc.tl.score_genes(adata, gene_list=present, score_name=f'FateScore_{comp}',
                                  use_raw=False, random_state=42)

        hyp = adata.obs['HypoxiaScore'].values

        # Component-level ρ
        for comp in FATE_COMPONENTS:
            col = f'FateScore_{comp}'
            if col not in adata.obs.columns:
                section_records.append({
                    'Section': sample_key, 'Patient': patient,
                    'Component': comp, 'rho': np.nan, 'n_spots': adata.shape[0]
                })
                continue
            fate = adata.obs[col].values
            rho, p = stats.spearmanr(hyp, fate)
            section_records.append({
                'Section': sample_key, 'Patient': patient,
                'Component': comp, 'rho': rho, 'p': p,
                'direction': 'positive' if rho > 0 else 'negative',
                'n_spots': adata.shape[0]
            })

        # Single-gene ρ
        for gene in ALL_FATE_GENES:
            if gene in adata.var_names:
                expr = adata[:, gene].X.toarray().ravel()
                rho, p = stats.spearmanr(hyp, expr)
                gene_records.append({
                    'Section': sample_key, 'Patient': patient,
                    'Gene': gene, 'rho': rho,
                    'direction': 'positive' if rho > 0 else 'negative',
                })

    df_section = pd.DataFrame(section_records)
    df_gene = pd.DataFrame(gene_records)

    # ── Aggregate to patient level (median ρ across sections) ──
    df_patient = df_section.groupby(['Patient', 'Component']).agg(
        rho=('rho', 'median'),
        n_sections=('Section', 'count'),
        n_spots_total=('n_spots', 'sum'),
    ).reset_index()

    patients = sorted(df_patient['Patient'].unique())
    n_patients = len(patients)
    threshold = int(np.ceil(2/3 * n_patients))
    print(f"\n  {n_patients} patients, direction threshold: {threshold}")

    # ── Per-component summary ──
    summary_rows = []
    for comp in FATE_COMPONENTS:
        sub = df_patient[df_patient['Component'] == comp]
        n_pos = (sub['rho'] > 0).sum()
        n_neg = (sub['rho'] < 0).sum()
        median_rho = sub['rho'].median()
        try:
            w_stat, w_p = stats.wilcoxon(sub['rho'].values, alternative='greater')
        except ValueError:
            w_stat, w_p = np.nan, 1.0
        summary_rows.append({
            'Component': comp,
            'n_positive': n_pos,
            'n_negative': n_neg,
            'n_total': len(sub),
            'threshold': threshold,
            'direction_pass': n_pos >= threshold,
            'median_rho': median_rho,
            'wilcoxon_stat': w_stat,
            'wilcoxon_p': w_p,
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary['q_BH'] = bh_correction(df_summary['wilcoxon_p'].values)
    df_summary['gate_pass'] = (df_summary['q_BH'] < 0.05) & (df_summary['direction_pass'])

    # ── Contrast ──
    synth = df_patient[df_patient['Component'] == 'Synthesis'].set_index('Patient')['rho']
    catab = df_patient[df_patient['Component'] == 'Catabolism'].set_index('Patient')['rho']
    common = synth.index.intersection(catab.index)
    synth = synth.loc[common]
    catab = catab.loc[common]
    D = synth - catab
    n_D_pos = (D > 0).sum()
    median_D = D.median()
    try:
        d_stat, d_p = stats.wilcoxon(D.values, alternative='greater')
    except ValueError:
        d_stat, d_p = np.nan, 1.0
    contrast_pass = (median_D > 0) and (n_D_pos >= threshold) and (d_p < 0.05)

    df_contrast = pd.DataFrame({
        'Patient': D.index,
        'rho_Synthesis': synth.values,
        'rho_Catabolism': catab.values,
        'D': D.values,
        'D_direction': ['synthesis_favored' if d > 0 else 'catabolism_favored' for d in D.values],
    })

    n_components_pass = df_summary['gate_pass'].sum()
    if contrast_pass:
        label = "synthesis_favored_transcriptional_shift"
    elif n_components_pass > 0:
        label = "glutamine_associated_state"
    else:
        label = "no_stable_component"

    verdict = {
        'gate': 'A_spatial',
        'n_patients': n_patients,
        'n_sections': len(SAMPLE_TO_PATIENT),
        'threshold': threshold,
        'n_components_pass': int(n_components_pass),
        'components_pass': df_summary[df_summary['gate_pass']]['Component'].tolist(),
        'contrast_pass': bool(contrast_pass),
        'contrast_n_D_positive': int(n_D_pos),
        'contrast_median_D': float(median_D),
        'contrast_wilcoxon_p': float(d_p),
        'gate_a_spatial_label': label,
    }

    # ── Save ──
    df_section.to_csv(OUT / 'gateA_spatial_section_rho.csv', index=False)
    df_patient.to_csv(OUT / 'gateA_spatial_patient_rho.csv', index=False)
    df_summary.to_csv(OUT / 'gateA_spatial_summary.csv', index=False)
    df_contrast.to_csv(OUT / 'gateA_spatial_contrast.csv', index=False)
    with open(OUT / 'gateA_spatial_verdict.json', 'w') as f:
        json.dump(verdict, f, indent=2)

    # ── Print ──
    print("\n" + "="*70)
    print("GATE A-SPATIAL: GREENWALD VISIUM GLUTAMINE FATE DECOMPOSITION")
    print("="*70)
    print(f"Patients: {n_patients}, Sections: {len(SAMPLE_TO_PATIENT)}")
    print(f"Direction threshold: {threshold}/{n_patients}")
    print(f"\nPer-component results:")
    print(df_summary[['Component','n_positive','n_negative','threshold',
                       'direction_pass','median_rho','wilcoxon_p','q_BH','gate_pass']].to_string(index=False))
    print(f"\nContrast (Synthesis - Catabolism):")
    print(f"  n_D_positive: {n_D_pos}/{len(D)}")
    print(f"  median_D: {median_D:.4f}")
    print(f"  Wilcoxon p: {d_p:.4f}")
    print(f"  Contrast PASS: {contrast_pass}")
    print(f"\nGate A-spatial label: {label}")
    print("="*70)

    # Single-gene summary
    print("\nSingle-gene per-patient direction counts (spatial):")
    gene_pat = df_gene.groupby(['Patient','Gene']).agg(rho=('rho','median')).reset_index()
    gene_summary = gene_pat.groupby('Gene').apply(
        lambda x: pd.Series({
            'n_positive': (x['rho'] > 0).sum(),
            'n_negative': (x['rho'] < 0).sum(),
            'median_rho': x['rho'].median(),
        })
    ).sort_values('median_rho', ascending=False)
    print(gene_summary.to_string())

if __name__ == '__main__':
    main()
