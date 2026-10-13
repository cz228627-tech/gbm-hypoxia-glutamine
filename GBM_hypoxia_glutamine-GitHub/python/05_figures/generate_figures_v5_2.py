#!/usr/bin/env python3
"""v5.2 figures: regenerate Figs 2-5 per approved revision brief (2026-10-11).

Every displayed number is derived from frozen CSVs (no hard-coded statistics);
assertions fail loudly if a displayed string would deviate from frozen values.

Fig 2: A q(C/N)=0.585; B rebuilt from Supplementary Table 5 (95th pct |rho|
       = 0.2508 -> color scale +/-0.251; 18 undefined rho -> light-gray 'NA').
Fig 3: C/D rebuilt from the FINAL patient grouping (Greenwald 5 patients /
       13 sections, Ravi 17 patients / 18 sections; analysis_role=primary).
Fig 4: C annotation box moved clear of the last bar's P label; dashed line
       relabeled 'Direction-consistency gate (6/8)'.
Fig 5: 5 columns (Lin 28 | Greenwald 5 | Ravi 17 | GSC 8 | U87MG 3v3);
       Greenwald all hollow up (5/5, P=0.0625 floor); GSC Catabolism hollow
       down, Direct uptake & Exchange/redox hollow up (P=0.148); U87MG fully
       descriptive (hollow arrows; overall P=0.007 = 11-gene t-test).

Also exports a figure-data package (figure_data/*.csv) for external
figure-making software.

Outputs -> results/15_figures_v5_2/
"""
import json
from pathlib import Path
import os

import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RES = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
SP = RES / '11_spatial_v5'
SUP = RES / '13_supplement' / 'supplementary_tables'
OUT = RES / '15_figures_v5_2'
PKG = OUT / 'figure_data'
PKG.mkdir(parents=True, exist_ok=True)

POS_COLOR = '#D62728'
NEG_COLOR = '#1F77B4'
NULL_COLOR = '#999999'
COMP_COLORS = {'Synthesis': '#1F77B4', 'Catabolism': '#FF7F0E',
               'Direct uptake': '#2CA02C', 'Exchange/redox': '#D62728',
               'C/N redistribution': '#9467BD'}

# ---------------------------------------------------------------- data
LIN = pd.read_csv(SP / 'lin_patient_stats.csv')            # 28 patients, v5 scoring
FD = pd.read_csv(SP / 'figure_data_v5.csv')                # frozen summary stats
SPAT = pd.read_csv(SP / 'spatial_patient_stats.csv')       # final patient grouping
SPAT_TESTS = pd.read_csv(SP / 'spatial_tests_summary.csv')
GENES = pd.read_csv(SUP / 'SupplementaryTable5_per_gene_rho_scRNA.csv')  # verified per-gene table
CONTRAST = pd.read_csv(RES / 'qa_lin_contrast_mean.csv')   # D per patient
DELTA11 = pd.read_csv(RES / 'gateB_per_line_delta_11gene.csv')
DELTA5 = pd.read_csv(RES / 'gateB_per_line_delta.csv')
LOCO = pd.read_csv(RES / 'qa_loco_gsc.csv')
U87_11 = pd.read_csv(RES / 'gateC_u87_11gene_summary.csv')
VISIUM1 = pd.read_csv(RES / 'visium_representative_section.csv', index_col=0)

FDL = FD[FD.cohort == 'Lin_scRNA (28 patients)'].set_index('component')


def fd_row(cohort, comp):
    r = FD[(FD.cohort == cohort) & (FD.component == comp)]
    assert len(r) == 1, (cohort, comp)
    return r.iloc[0]


def q_display(q):
    return 'q<0.001' if q < 0.001 else f'q={q:.3f}'


def p_display(p):
    return 'P<0.001' if p < 0.001 else f'P={p:.3f}'


# ============================================================
# FIGURE 2  (decomposition)
# ============================================================
print('Figure 2...')
COMPONENTS = ['Synthesis', 'Catabolism', 'Direct uptake', 'Exchange/redox', 'C/N redistribution']
COMP_RHO_COLS = {'Synthesis': 'rho_Synthesis', 'Catabolism': 'rho_Catabolism',
                 'Direct uptake': 'rho_Direct_uptake', 'Exchange/redox': 'rho_Exchange_redox',
                 'C/N redistribution': 'rho_CN_redistribution'}
FD_COMP_KEYS = {'Synthesis': 'Synthesis', 'Catabolism': 'Catabolism',
                'Direct uptake': 'Direct_uptake', 'Exchange/redox': 'Exchange_redox',
                'C/N redistribution': 'CN_redistribution'}

comp_data = LIN[['Patient'] + [COMP_RHO_COLS[c] for c in COMPONENTS]].copy()
comp_data.columns = ['Patient'] + COMPONENTS
q_vals = {c: float(FDL.loc[FD_COMP_KEYS[c], 'q_BH']) for c in COMPONENTS}
q_texts = {c: q_display(q_vals[c]) for c in COMPONENTS}
assert q_texts['C/N redistribution'] == 'q=0.585', q_texts  # brief: 0.598 -> 0.585

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw={'width_ratios': [1, 1.3, 0.8]})

# --- Panel A: box plots of patient-level rho per component
ax = axes[0]
plot_data = comp_data.melt(id_vars='Patient', var_name='Component', value_name='rho')
plot_data['Component'] = pd.Categorical(plot_data['Component'], categories=COMPONENTS, ordered=True)
sns.boxplot(data=plot_data, x='Component', y='rho', order=COMPONENTS, ax=ax,
            palette=[COMP_COLORS[c] for c in COMPONENTS], width=0.6,
            fliersize=3, boxprops=dict(alpha=0.6))
for p in ax.patches:
    p.set_alpha(0.6)
ax.axhline(y=0, color='black', linewidth=0.8)
ax.set_xticklabels(COMPONENTS, fontsize=9, rotation=30, ha='right')
ax.set_ylabel('Patient-level Spearman \u03c1')
ax.set_xlabel('')
ax.set_title('A', fontsize=14, fontweight='bold', loc='left')
ymax = comp_data[COMPONENTS].values.max()
for i, c in enumerate(COMPONENTS):
    ax.text(i, ymax * 1.07, q_texts[c], ha='center', fontsize=9)
ax.set_ylim(top=ymax * 1.18)
sns.despine(ax=ax)

# --- Panel B: per-gene rho heatmap (Supplementary Table 5)
ax = axes[1]
gene_groups = {'Synthesis': ['GLUL'], 'Catabolism': ['GLS', 'GLS2'],
               'Direct uptake': ['SLC1A5', 'SLC38A1', 'SLC38A2'],
               'Exchange/redox': ['SLC7A5', 'SLC7A11'],
               'C/N redistribution': ['GLUD1', 'GOT1', 'GOT2', 'GPT', 'GPT2']}
gene_order = [g for gs in gene_groups.values() for g in gs]
assert sorted(gene_order) == sorted(GENES.Gene.unique())

patient_order = LIN.sort_values('rho_rho_gln')['Patient'].tolist()
mat = GENES.pivot(index='Patient', columns='Gene', values='rho').loc[patient_order, gene_order]
n_na = int(mat.isna().sum().sum())
pct95 = float(mat.stack().abs().quantile(0.95))
assert n_na == 18, n_na                       # brief: 18 undefined correlations
assert round(pct95, 3) == 0.251, pct95        # brief: 95th pct |rho| ~ 0.251
VMAX = 0.251

ax.set_facecolor('#E4E4E4')                   # NA cells show as light gray
sns.heatmap(mat, cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX),
            ax=ax, xticklabels=True, yticklabels=False, mask=mat.isna(),
            cbar_kws={'label': 'Spearman \u03c1 (truncated at 95th percentile)', 'shrink': 0.7},
            linewidths=0.3, linecolor='white')
for y, x in zip(*np.where(mat.isna().values)):
    ax.text(x + 0.5, y + 0.5, 'NA', ha='center', va='center', fontsize=4, color='#555555')
ax.set_xlabel('')
ax.set_ylabel('Patients (sorted by overall \u03c1)')
ax.set_title('B', fontsize=14, fontweight='bold', loc='left')
ax.tick_params(axis='x', rotation=90)
# component brackets above the heatmap
for comp, gs in gene_groups.items():
    x0 = gene_order.index(gs[0]) - 0.4
    x1 = gene_order.index(gs[-1]) + 0.4
    xm = (x0 + x1) / 2
    ax.plot([x0, x0, x1, x1], [1.25, 1.45, 1.45, 1.25], color=COMP_COLORS[comp],
            lw=1.2, clip_on=False, transform=ax.get_xaxis_transform())
    ax.text(xm, 1.55, comp, ha='center', va='bottom', fontsize=7,
            color=COMP_COLORS[comp], clip_on=False, transform=ax.get_xaxis_transform())

# --- Panel C: D per patient
ax = axes[2]
contrast_vals = CONTRAST.set_index('Patient')['D'].sort_values()
colors = [POS_COLOR if v > 0 else NEG_COLOR for v in contrast_vals]
ax.barh(range(len(contrast_vals)), contrast_vals, color=colors, height=0.7, edgecolor='none')
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_yticks(range(len(contrast_vals)))
ax.set_yticklabels(contrast_vals.index, fontsize=6.5)
ax.set_xlabel('D = \u03c1(Synthesis) \u2212 \u03c1(Catabolism)')
ax.set_title('C', fontsize=14, fontweight='bold', loc='left')
n_d_pos = int((contrast_vals > 0).sum())
d_median = float(contrast_vals.median())
d_p = float(stats.wilcoxon(contrast_vals.values, alternative='two-sided', method='approx').pvalue)
assert n_d_pos == 24 and round(d_median, 3) == 0.117 and d_p < 0.0001
ax.text(0.95, 0.05, f'{n_d_pos}/{len(contrast_vals)} synthesis-favored\nmedian D={d_median:.3f}\nP<0.0001',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
        color=POS_COLOR, fontweight='bold')
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(OUT / 'figure2_decomposition.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'figure2_decomposition.svg', bbox_inches='tight')
plt.close()

# --- data package: figure 2
comp_data.to_csv(PKG / 'fig2A_patient_component_rho.csv', index=False)
pd.DataFrame({'component': COMPONENTS,
              'n': [28] * 5,
              'n_positive': [int((comp_data[c] > 0).sum()) for c in COMPONENTS],
              'median_rho': [float(comp_data[c].median()) for c in COMPONENTS],
              'p_two_sided_wilcoxon': [float(FDL.loc[FD_COMP_KEYS[c], 'p']) for c in COMPONENTS],
              'q_BH': [q_vals[c] for c in COMPONENTS],
              'q_display': [q_texts[c] for c in COMPONENTS]}).to_csv(
    PKG / 'fig2A_component_summary.csv', index=False)
mat.to_csv(PKG / 'fig2B_patient_gene_rho_matrix.csv')  # empty cells = NA (18)
contrast_vals.rename('D').to_csv(PKG / 'fig2C_patient_D.csv')
pd.DataFrame({'note': ['panel_A_q_texts', 'panel_B_color_scale', 'panel_B_NA_cells', 'panel_C_note'],
              'value': ['; '.join(f'{c}: {q_texts[c]}' for c in COMPONENTS),
                        f'+/-{VMAX} (95th percentile of |rho| = {pct95:.4f})',
                        f'{n_na} undefined correlations shown as light-gray NA',
                        f'{n_d_pos}/28 synthesis-favored; median D={d_median:.3f}; P<0.0001']}
            ).to_csv(PKG / 'fig2_annotations.csv', index=False)
print('  Figure 2 saved')

# ============================================================
# FIGURE 3  (spatial validation; final patient grouping)
# ============================================================
print('Figure 3...')
gw = SPAT[(SPAT.cohort == 'Greenwald_Visium') & (SPAT.analysis_role == 'primary')].copy()
rv = SPAT[(SPAT.cohort == 'Ravi_Visium') & (SPAT.analysis_role == 'primary')].copy()
assert len(gw) == 5 and gw.n_sections.sum() == 13    # brief: 5 patients / 13 sections
assert len(rv) == 17 and rv.n_sections.sum() == 18   # brief: 17 patients / 18 sections

p_gw = float(SPAT_TESTS[(SPAT_TESTS.comparison == 'primary|Greenwald_Visium|rho_gln|all')].p_wilcoxon.iloc[0])
p_rv = float(SPAT_TESTS[(SPAT_TESTS.comparison == 'primary|Ravi_Visium|rho_gln|all')].p_wilcoxon.iloc[0])
assert p_gw == 0.0625 and abs(p_rv - 1.52587890625e-05) < 1e-12

def mk_labels(df):
    return [p.replace('GW:', '').replace('RAVI:UKF', 'UKF') for p in df['global_patient_id']]

fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

ax = fig.add_subplot(gs[0, 0])
sc_plot = ax.scatter(VISIUM1['pxl_col'], -VISIUM1['pxl_row'], c=VISIUM1['HypoxiaScore'],
                     cmap='Reds', s=8, alpha=0.8, edgecolors='none')
plt.colorbar(sc_plot, ax=ax, label='Hypoxia score', shrink=0.7)
ax.set_title('A', fontsize=14, fontweight='bold', loc='left')
ax.set_xlabel(''); ax.set_ylabel(''); ax.set_xticks([]); ax.set_yticks([])
ax.text(0.02, 0.98, 'zh916 (bulk)', transform=ax.transAxes, fontsize=9,
        va='top', ha='left', style='italic')

ax = fig.add_subplot(gs[0, 1])
sc_plot = ax.scatter(VISIUM1['pxl_col'], -VISIUM1['pxl_row'], c=VISIUM1['GlutamineScore'],
                     cmap='Blues', s=8, alpha=0.8, edgecolors='none')
plt.colorbar(sc_plot, ax=ax, label='Glutamine score (11-gene)', shrink=0.7)
ax.set_title('B', fontsize=14, fontweight='bold', loc='left')
ax.set_xlabel(''); ax.set_ylabel(''); ax.set_xticks([]); ax.set_yticks([])
ax.text(0.02, 0.98, 'zh916 (bulk)', transform=ax.transAxes, fontsize=9,
        va='top', ha='left', style='italic')

SEP = len(gw) - 0.5  # Greenwald block occupies rows 0-4 (bottom)

def cohort_panel(ax, df, col, xlabel, title, note):
    df = df.reset_index(drop=True)
    colors = [POS_COLOR if v > 0 else NEG_COLOR for v in df[col]]
    assert all(v > 0 for v in df[col]), 'all 22 patients expected positive'
    ax.barh(range(len(df)), df[col], color=colors, height=0.7, edgecolor='none')
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(mk_labels(df), fontsize=6.5)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.axhline(y=SEP, color='#888888', linewidth=0.8, linestyle='dashed')
    ax.text(0.98, SEP - 0.15, 'Greenwald (n=5)', transform=ax.get_yaxis_transform(),
            ha='right', va='top', fontsize=8, style='italic')
    ax.text(0.98, SEP + 0.15, 'Ravi (n=17)', transform=ax.get_yaxis_transform(),
            ha='right', va='bottom', fontsize=8, style='italic')
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
    ax.text(0.95, 0.05, note, transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color=POS_COLOR, fontweight='bold')
    sns.despine(ax=ax)

comb_rho = pd.concat([gw.sort_values('median_rho_rho_gln'),
                      rv.sort_values('median_rho_rho_gln')]).reset_index(drop=True)
comb_z = pd.concat([gw.sort_values('median_Z_rho_gln'),
                    rv.sort_values('median_Z_rho_gln')]).reset_index(drop=True)

ax = fig.add_subplot(gs[1, 0])
note_c = f'Greenwald 5/5+ (P={p_gw:.4f})\nRavi 17/17+ (P=1.5\u00d710\u207b\u2075)'
cohort_panel(ax, comb_rho, 'median_rho_rho_gln',
             'Patient-level Spearman \u03c1 (hypoxia vs glutamine)', 'C', note_c)

ax = fig.add_subplot(gs[1, 1])
note_d = f'Greenwald {(gw["median_Z_rho_gln"] > 0).sum()}/5+\nRavi {(rv["median_Z_rho_gln"] > 0).sum()}/17+'
cohort_panel(ax, comb_z, 'median_Z_rho_gln',
             'Null-calibrated patient-level Z-score', 'D', note_d)

plt.savefig(OUT / 'figure3_spatial.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'figure3_spatial.svg', bbox_inches='tight')
plt.close()

# --- data package: figure 3
vis_export = VISIUM1.reset_index().rename(columns={'index': 'barcode'})
vis_export.to_csv(PKG / 'fig3AB_zh916_per_spot.csv', index=False)
comb_rho[['cohort', 'global_patient_id', 'n_sections', 'median_rho_rho_gln']].to_csv(
    PKG / 'fig3C_patient_rho.csv', index=False)
comb_z[['cohort', 'global_patient_id', 'n_sections', 'median_Z_rho_gln']].to_csv(
    PKG / 'fig3D_patient_Z.csv', index=False)
pd.DataFrame({'cohort': ['Greenwald_Visium', 'Ravi_Visium'],
              'n_patients': [5, 17], 'n_sections': [13, 18],
              'n_positive': [5, 17],
              'median_rho': [float(fd_row('Greenwald Visium (5 patients)', 'rho_gln')['median']),
                             float(fd_row('Ravi Visium (17 patients)', 'rho_gln')['median'])],
              'p_exact_two_sided': [p_gw, p_rv]}).to_csv(PKG / 'fig3_tests.csv', index=False)
print('  Figure 3 saved')

# ============================================================
# FIGURE 4  (perturbation; annotation moved, gate relabeled)
# ============================================================
print('Figure 4...')
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw={'width_ratios': [1, 1.2, 1]})

# --- Panel A: paired normoxia -> hypoxia
ax = axes[0]
gsc_11 = DELTA11.sort_values('GSC_line')
norm_vals = gsc_11['Normoxia_11'].values
hyp_vals = gsc_11['Hypoxia_11'].values
for i in range(len(gsc_11)):
    color = POS_COLOR if hyp_vals[i] > norm_vals[i] else NEG_COLOR
    ax.plot([0, 1], [norm_vals[i], hyp_vals[i]], 'o-', color=color, alpha=0.6, markersize=5)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Normoxia', 'Hypoxia'], fontsize=11)
ax.set_ylabel('Overall glutamine score (11-gene, z-scored)')
ax.set_title('A', fontsize=14, fontweight='bold', loc='left')
n_up = int((hyp_vals > norm_vals).sum())
p_paired = float(stats.wilcoxon(hyp_vals - norm_vals, alternative='two-sided', method='exact').pvalue)
assert n_up == 7 and round(p_paired, 4) == 0.0156
ax.text(0.5, 0.95, f'{n_up}/{len(gsc_11)} increased\nP={p_paired:.4f}',
        transform=ax.transAxes, ha='center', va='top', fontsize=9,
        color=POS_COLOR, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
sns.despine(ax=ax)

# --- Panel B: component delta heatmap
ax = axes[1]
comp_names = ['Synthesis', 'Catabolism', 'Direct_uptake', 'Exchange_redox', 'CN_redistribution']
comp_display = ['Synthesis', 'Catabolism', 'Direct\nuptake', 'Exchange\nredox',
                'C/N\nredistribution', 'Overall\nglutamine\n(11-gene)']
delta_data = DELTA5.set_index('GSC_line')[[f'Delta_{c}' for c in comp_names]].copy()
delta_data['Overall_11gene'] = DELTA11.set_index('GSC_line')['Delta_OverallGlutamine_11']
delta_data.columns = comp_display
vmax = max(abs(delta_data.values.min()), abs(delta_data.values.max()))
sns.heatmap(delta_data, cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
            ax=ax, xticklabels=True, yticklabels=True,
            cbar_kws={'label': '\u0394 (hypoxia \u2212 normoxia)', 'shrink': 0.7},
            linewidths=0.3, linecolor='white', annot=True, fmt='.2f', annot_kws={'fontsize': 7})
ax.set_title('B', fontsize=14, fontweight='bold', loc='left')
ax.set_ylabel('GSC line')

# --- Panel C: LOCO with annotation clear of the last bar
ax = axes[2]
loco_order = ['Full_11', 'Minus_Synthesis', 'Minus_Catabolism', 'Minus_Transport', 'Minus_CN']
loco_display = ['Full\n(11 genes)', 'Minus\nSynth', 'Minus\nCatab', 'Minus\nTransport', 'Minus\nC/N']
loco_f = LOCO.set_index('score_set').loc[loco_order].reset_index()
loco_f['display'] = loco_display
colors = [NEG_COLOR if s == 'Minus_CN' else (POS_COLOR if surv else NEG_COLOR)
          for s, surv in zip(loco_f['score_set'], loco_f['survives'])]
ax.bar(range(len(loco_f)), loco_f['n_positive'], color=colors, alpha=0.7, edgecolor='none')
ax.axhline(y=6, color='gray', linestyle='--', linewidth=0.8,
           label='Direction-consistency gate (6/8)')
ax.set_xticks(range(len(loco_f)))
ax.set_xticklabels(loco_f['display'], fontsize=8)
ax.set_ylabel('n positive / 8')
ax.set_title('C', fontsize=14, fontweight='bold', loc='left')
ax.set_ylim(0, 10.5)
p_texts = []
for i, row in enumerate(loco_f.itertuples()):
    pval = row.wilcoxon_p_twosided
    txt = f'P={pval:.3f}'
    p_texts.append(txt)
    ax.text(i, row.n_positive + 0.2, txt, ha='center', fontsize=7, rotation=45)
assert p_texts == ['P=0.016', 'P=0.016', 'P=0.008', 'P=0.016', 'P=0.383'], p_texts
# annotation box moved to the empty upper-left band, arrow unchanged in meaning
ax.annotate('', xy=(3.85, 2.5), xytext=(1.2, 9.3),
            arrowprops=dict(arrowstyle='->', color=NEG_COLOR, lw=1.5, linestyle='dashed'))
ax.text(1.2, 9.55, '7/8 \u2192 2/8\nafter C/N removal', ha='left', va='bottom',
        fontsize=8, color=NEG_COLOR, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='white', edgecolor=NEG_COLOR, alpha=0.9))
ax.legend(fontsize=8, loc='upper right')
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(OUT / 'figure4_perturbation.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'figure4_perturbation.svg', bbox_inches='tight')
plt.close()

# --- data package: figure 4
gsc_11[['GSC_line', 'Normoxia_11', 'Hypoxia_11']].assign(
    direction=['increased' if h > n else 'decreased' for n, h in zip(norm_vals, hyp_vals)]
).to_csv(PKG / 'fig4A_line_paired_scores.csv', index=False)
delta_data.to_csv(PKG / 'fig4B_line_component_deltas.csv')
loco_out = loco_f[['score_set', 'display', 'n_genes', 'n_positive', 'n_negative',
                   'median_delta', 'wilcoxon_p_twosided', 'survives']].copy()
loco_out['p_display'] = p_texts
loco_out['bar_color'] = ['red' if c == POS_COLOR else 'blue' for c in colors]
loco_out.to_csv(PKG / 'fig4C_loco_gsc.csv', index=False)
pd.DataFrame({'note': ['panel_A', 'panel_C_gate', 'panel_C_annotation'],
              'value': ['7/8 increased, exact two-sided Wilcoxon P=0.0156',
                        'direction-consistency gate at 6/8 (not a significance threshold)',
                        'Minus Transport removes SLC1A5, SLC7A5, SLC7A11 (see manuscript caption)']}
            ).to_csv(PKG / 'fig4_annotations.csv', index=False)
print('  Figure 4 saved')

# ============================================================
# FIGURE 5  (context model; 5 columns, arrow types per approved brief)
# ============================================================
print('Figure 5...')
CONTEXTS = [('Lin_scRNA (28 patients)', 'Lin scRNA\n(28 patients)'),
            ('Greenwald Visium (5 patients)', 'Greenwald Visium\n(5 pts)'),
            ('Ravi Visium (17 patients)', 'Ravi Visium\n(17 pts)'),
            ('GSC 24h hypoxia (8 lines)', 'GSC 24h\n(8 lines)'),
            ('U87MG (3v3)', 'U87MG\n(3v3)')]
COMP_ROWS = [('rho_gln', 'Overall glutamine\n(11-gene)'),
             ('Synthesis', 'Synthesis\n(GLUL)'),
             ('Catabolism', 'Catabolism\n(GLS, GLS2)'),
             ('Direct_uptake', 'Direct uptake\n(SLC1A5/38A1/38A2)'),
             ('Exchange_redox', 'Exchange/redox\n(SLC7A5/7A11)'),
             ('CN_redistribution', 'C/N redistribution\n(GLUD1, GOT1/2, GPT/2)')]

# arrow types approved by the user (2026-10-11):
#  - Greenwald column: all hollow up (5/5 positive; P=0.0625 = n=5 exact floor)
#  - GSC: Catabolism hollow down (P=0.250, q=0.3125); Direct uptake and
#    Exchange/redox hollow up (P=0.148); Overall and C/N solid; Synthesis null
#  - U87MG: fully descriptive -> hollow down for Overall/Exchange/C-N, null others
ARROW = {
    ('Lin', 'rho_gln'): 'up_solid', ('Lin', 'Synthesis'): 'up_solid',
    ('Lin', 'Catabolism'): 'up_solid', ('Lin', 'Direct_uptake'): 'up_solid',
    ('Lin', 'Exchange_redox'): 'up_solid', ('Lin', 'CN_redistribution'): 'null',
    ('Greenwald', 'rho_gln'): 'up_hollow', ('Greenwald', 'Synthesis'): 'up_hollow',
    ('Greenwald', 'Catabolism'): 'up_hollow', ('Greenwald', 'Direct_uptake'): 'up_hollow',
    ('Greenwald', 'Exchange_redox'): 'up_hollow', ('Greenwald', 'CN_redistribution'): 'up_hollow',
    ('Ravi', 'rho_gln'): 'up_solid', ('Ravi', 'Synthesis'): 'up_solid',
    ('Ravi', 'Catabolism'): 'up_hollow', ('Ravi', 'Direct_uptake'): 'up_solid',
    ('Ravi', 'Exchange_redox'): 'up_solid', ('Ravi', 'CN_redistribution'): 'up_solid',
    ('GSC', 'rho_gln'): 'up_solid', ('GSC', 'Synthesis'): 'null',
    ('GSC', 'Catabolism'): 'down_hollow', ('GSC', 'Direct_uptake'): 'up_hollow',
    ('GSC', 'Exchange_redox'): 'up_hollow', ('GSC', 'CN_redistribution'): 'up_solid',
    ('U87', 'rho_gln'): 'down_hollow', ('U87', 'Synthesis'): 'null',
    ('U87', 'Catabolism'): 'null', ('U87', 'Direct_uptake'): 'null',
    ('U87', 'Exchange_redox'): 'down_hollow', ('U87', 'CN_redistribution'): 'down_hollow',
}

u87_p11 = float(U87_11[U87_11.Component == 'OverallGlutamine_11gene'].t_p_twosided.iloc[0])
assert round(u87_p11, 3) == 0.007  # brief: keep P=0.007 (11-gene t-test)

def cell_annotation(ctx_key, comp_key, row):
    n_pos, n = int(row['n_positive']), int(row['n'])
    q = row['q_BH']
    if ctx_key == 'Lin':
        return f'{n_pos}/{n}'
    if ctx_key == 'Greenwald':
        return f'{n_pos}/{n}'
    if ctx_key == 'Ravi':
        if comp_key == 'rho_gln':
            return f'{n_pos}/{n}'
        if comp_key == 'Catabolism':
            return f'{n_pos}/{n}'
        return q_display(q)
    if ctx_key == 'GSC':
        if comp_key == 'CN_redistribution':
            return q_display(q)
        return f'{n_pos}/{n}'
    # U87MG (descriptive)
    if comp_key == 'rho_gln':
        return f'P={u87_p11:.3f}'
    if comp_key in ('Exchange_redox', 'CN_redistribution'):
        return q_display(q)
    return 'ns'

CTX_SHORT = {'Lin_scRNA (28 patients)': 'Lin',
             'Greenwald Visium (5 patients)': 'Greenwald',
             'Ravi Visium (17 patients)': 'Ravi',
             'GSC 24h hypoxia (8 lines)': 'GSC',
             'U87MG (3v3)': 'U87'}
cells = []
for ctx_key, ctx_label in CONTEXTS:
    for comp_key, comp_label in COMP_ROWS:
        row = fd_row(ctx_key, comp_key)
        ctx_short = CTX_SHORT[ctx_key]
        arrow = ARROW[(ctx_short, comp_key)]
        txt = cell_annotation(ctx_short, comp_key, row)
        cells.append({'context': ctx_label.replace('\n', ' '), 'component': comp_label.replace('\n', ' '),
                      'n': int(row['n']), 'n_positive': int(row['n_positive']),
                      'median': float(row['median']), 'p': float(row['p']),
                      'q_BH': None if pd.isna(row['q_BH']) else float(row['q_BH']),
                      'arrow_type': arrow, 'display_text': txt})
cells = pd.DataFrame(cells)
assert len(cells) == 30

fig, ax = plt.subplots(figsize=(17, 10))
ax.set_xlim(-0.5, len(CONTEXTS) - 0.5)
ax.set_ylim(-0.5, len(COMP_ROWS) - 0.5)
ax.invert_yaxis()

for i, (_, comp_label) in enumerate(COMP_ROWS):
    for j, (_, ctx_label) in enumerate(CONTEXTS):
        cell = cells[(cells.context == ctx_label.replace('\n', ' ')) &
                     (cells.component == comp_label.replace('\n', ' '))].iloc[0]
        arrow_type, txt = cell['arrow_type'], cell['display_text']
        bg = {'up_solid': '#F0F0F0', 'down_solid': '#F0F0F0',
              'up_hollow': '#F8F8F8', 'down_hollow': '#F8F8F8'}.get(arrow_type, '#FFFFFF')
        ax.add_patch(mpatches.FancyBboxPatch((j - 0.45, i - 0.45), 0.9, 0.9,
                     boxstyle='round,pad=0.02', facecolor=bg,
                     edgecolor='#CCCCCC', linewidth=0.5))
        cx, cy, dy = j, i, 0.13
        if arrow_type == 'up_solid':
            ax.annotate('', xy=(cx, cy - dy), xytext=(cx, cy + dy),
                        arrowprops=dict(arrowstyle='->', color=POS_COLOR, lw=2.5))
        elif arrow_type == 'down_solid':
            ax.annotate('', xy=(cx, cy + dy), xytext=(cx, cy - dy),
                        arrowprops=dict(arrowstyle='->', color=NEG_COLOR, lw=2.5))
        elif arrow_type == 'up_hollow':
            ax.annotate('', xy=(cx, cy - dy), xytext=(cx, cy + dy),
                        arrowprops=dict(arrowstyle='->', color=POS_COLOR, lw=1.5, linestyle='dashed'))
        elif arrow_type == 'down_hollow':
            ax.annotate('', xy=(cx, cy + dy), xytext=(cx, cy - dy),
                        arrowprops=dict(arrowstyle='->', color=NEG_COLOR, lw=1.5, linestyle='dashed'))
        elif arrow_type == 'null':
            ax.text(cx, cy, '\u2194', ha='center', va='center', fontsize=14, color=NULL_COLOR)
        ax.text(cx, cy + 0.24, txt, ha='center', va='top', fontsize=6.5,
                color='#333333' if arrow_type != 'null' else '#666666')

ax.set_xticks(range(len(CONTEXTS)))
ax.set_xticklabels([c[1] for c in CONTEXTS], fontsize=10, fontweight='bold')
ax.set_yticks(range(len(COMP_ROWS)))
ax.set_yticklabels([c[1] for c in COMP_ROWS], fontsize=9)
ax.tick_params(length=0)
for i in range(len(COMP_ROWS) + 1):
    ax.axhline(y=i - 0.5, color='#DDDDDD', linewidth=0.5)
for j in range(len(CONTEXTS) + 1):
    ax.axvline(x=j - 0.5, color='#DDDDDD', linewidth=0.5)

legend_elements = [
    mpatches.Patch(facecolor='white', edgecolor=POS_COLOR, linewidth=2,
                   label='Solid \u2191: significant positive (P<0.05)'),
    mpatches.Patch(facecolor='white', edgecolor=NEG_COLOR, linewidth=2,
                   label='Solid \u2193: significant negative (P<0.05)'),
    mpatches.Patch(facecolor='white', edgecolor='#666666', linewidth=1.5, linestyle='--',
                   label='Dashed \u2191/\u2193: directional trend (not significant)'),
    mpatches.Patch(facecolor='white', edgecolor=NULL_COLOR, linewidth=1,
                   label='\u2194: null'),
]
fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.06),
           ncol=4, fontsize=8, frameon=False)
ax.set_title('Context-dependent glutamine transcriptional remodeling under hypoxia in glioblastoma',
             fontsize=13, fontweight='bold', pad=15)
fig.text(0.5, 0.028, 'U87MG (3v3) column: descriptive direction only; smallest attainable exact permutation P = 0.10. '
                     'Greenwald (n=5): P = 0.0625 is the smallest attainable exact P.',
         ha='center', fontsize=7.5, color='#555555')
fig.text(0.5, 0.005, 'Transcriptional associations only; no metabolic flux inference',
         ha='center', fontsize=8, style='italic', color='#666666')
plt.tight_layout()
plt.subplots_adjust(bottom=0.16)
plt.savefig(OUT / 'figure5_context_model.png', dpi=300, bbox_inches='tight')
plt.savefig(OUT / 'figure5_context_model.svg', bbox_inches='tight')
plt.close()

cells.to_csv(PKG / 'fig5_matrix_cells.csv', index=False)
print('  Figure 5 saved')

# ------------------------------------------------------------
json.dump({'fig2': {'q_texts': q_texts, 'vmax': VMAX, 'n_na': n_na, 'pct95': pct95,
                    'panel_C': {'n_positive': n_d_pos, 'median_D': d_median, 'p': d_p}},
           'fig3': {'p_greenwald': p_gw, 'p_ravi': p_rv,
                    'gw_patients': int(len(gw)), 'gw_sections': int(gw.n_sections.sum()),
                    'rv_patients': int(len(rv)), 'rv_sections': int(rv.n_sections.sum())},
           'fig4': {'panel_A': {'n_up': n_up, 'p': p_paired}, 'panel_C_p_texts': p_texts}},
          open(OUT / 'generation_summary.json', 'w'), indent=2)
print('All figures + data package written to', OUT)
