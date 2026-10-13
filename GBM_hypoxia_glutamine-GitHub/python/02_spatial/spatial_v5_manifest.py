"""Phase B (v5): manifest-driven spatial rerun with 1000 block permutations.

Reads ONLY the frozen manifest results/10_cohort_audit/cohort_manifest_v1.csv
(md5 ffdebaf9af90cdd9182e9e10fc3d24d9). No hard-coded sample-to-patient maps.

Pipeline per section: in-tissue spots -> normalize_total(1e4) -> log1p ->
arithmetic-mean gene-set scores -> Spearman rho.
Null: KMeans spatial blocks (k = max(5, n//50), seed 42), 1000 within-block
permutations -> Z. Statistics: overall glutamine module rho, the five fate
components, and D = rho_synthesis - rho_catabolism.

Outputs -> results/11_spatial_v5/
"""
import warnings; warnings.filterwarnings('ignore')
import os
import hashlib, json
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
RESULTS_DIR = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
MANIFEST = Path(os.environ.get("GBM_COHORT_MANIFEST", RESULTS_DIR / '10_cohort_audit' / 'cohort_manifest_v1.csv')).expanduser()
MANIFEST_MD5 = 'ffdebaf9af90cdd9182e9e10fc3d24d9'
OUT = RESULTS_DIR / '11_spatial_v5'; OUT.mkdir(parents=True, exist_ok=True)
GW_DIR = Path(os.environ.get("GBM_GREENWALD_DIR", DATA_DIR / "greenwald_visium")).expanduser()
RV_DIR = Path(os.environ.get("GBM_RAVI_DIR", DATA_DIR / "ravi_visium")).expanduser()
LIN_H5AD = Path(os.environ.get("GBM_LIN_H5AD", DATA_DIR / "lin" / "gbm_integrated.h5ad")).expanduser()

N_PERM = 1000
SEED = 42

HYP = ['CA9','VEGFA','NDRG1','P4HA1','EGLN3','BNIP3','ADM','ANKRD37',
       'CDKN1A','GYS1','DDIT3','F3','JMY','LRP1','LOX','ANGPTL4']
GLN11 = ['GLS','GLS2','GLUL','GLUD1','GOT1','GOT2','GPT','GPT2','SLC1A5','SLC7A5','SLC7A11']
FATE = {
    'Synthesis':         ['GLUL'],
    'Catabolism':        ['GLS','GLS2'],
    'Direct_uptake':     ['SLC1A5','SLC38A1','SLC38A2'],
    'Exchange_redox':    ['SLC7A5','SLC7A11'],
    'CN_redistribution': ['GLUD1','GOT1','GOT2','GPT','GPT2'],
}
STATS_NAMES = ['rho_gln'] + list(FATE.keys()) + ['D']

# ---------------------------------------------------------------- manifest
m = pd.read_csv(MANIFEST, dtype=str)
md5 = hashlib.md5(MANIFEST.read_bytes()).hexdigest()
assert md5 == MANIFEST_MD5, f'manifest md5 changed: {md5}'
print(f'manifest ok, md5 {md5}, {len(m)} rows')

def cohort_rows(cohort, field):
    return m[(m.cohort == cohort) & (m[field] == 'True')]

# ---------------------------------------------------------------- scoring
def amean(X, genes, var_names):
    idx = [var_names.get_loc(g) for g in genes if g in var_names]
    if hasattr(X, 'toarray'): X = X.toarray()
    return np.asarray(X)[:, idx].mean(axis=1)

def load_section(section_id, cohort):
    if cohort.startswith('Greenwald'):
        h5 = next(GW_DIR.glob(f'*_{section_id}_filtered_feature_bc_matrix.h5'))
        posf = next(GW_DIR.glob(f'*_{section_id}_tissue_positions_list.csv.gz'))
        pos = pd.read_csv(posf, header=None, compression='gzip',
                          names=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
    else:
        h5 = RV_DIR / f'{section_id}__filtered_feature_bc_matrix.h5'
        posf = RV_DIR / f'{section_id}__tissue_positions_list.csv'
        pos = pd.read_csv(posf, header=None,
                          names=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
        if not str(pos.iloc[0, 1]).lstrip('-').isdigit():   # header row present
            pos = pos.iloc[1:].reset_index(drop=True)
    pos['in_tissue'] = pos['in_tissue'].astype(int)
    pos = pos[pos.in_tissue == 1].set_index('barcode')
    ad = sc.read_10x_h5(h5); ad.var_names_make_unique()
    ad = ad[ad.obs_names.isin(pos.index)].copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    return ad, pos

def section_stats(section_id, cohort):
    ad, pos = load_section(section_id, cohort)
    n = ad.n_obs
    vn = ad.var_names
    X = ad.X
    scores = {'rho_gln': amean(X, GLN11, vn)}
    for comp, genes in FATE.items():
        scores[comp] = amean(X, genes, vn)
    scores['D'] = None  # computed from Synthesis/Catabolism ranks
    coords = pos.loc[ad.obs_names, ['pxl_row','pxl_col']].astype(float).values

    # rank-transform once; within-block permutation of values permutes ranks identically
    r = {k: rankdata(v) for k, v in scores.items() if v is not None}
    r_h = rankdata(amean(X, HYP, vn))
    c_h = r_h - r_h.mean(); nh2 = (c_h**2).sum()

    k = min(max(5, n // 50), n)
    blocks = KMeans(n_clusters=k, random_state=SEED, n_init=10).fit_predict(coords)
    idx = [np.where(blocks == b)[0] for b in np.unique(blocks)]
    idx = [ii for ii in idx if len(ii) > 1]

    keys = [kx for kx in ['rho_gln'] + list(FATE.keys())]
    rng = np.random.RandomState(SEED)
    null = {s: np.zeros(N_PERM) for s in STATS_NAMES}
    denom = {kx: np.sqrt(((r[kx] - r[kx].mean())**2).sum()) for kx in keys}
    for i in range(N_PERM):
        rhos = {}
        for kx in keys:
            rp = r[kx].copy()
            for ii in idx:
                rp[ii] = rng.permutation(rp[ii])
            cp = rp - rp.mean()
            rhos[kx] = (c_h * cp).sum() / np.sqrt(nh2 * (cp**2).sum())
        for kx in keys:
            null[kx if kx != 'Synthesis' else 'Synthesis'][i] = rhos[kx]
        null['D'][i] = rhos['Synthesis'] - rhos['Catabolism']

    out = {'section_id': section_id, 'cohort': cohort, 'n_spots': n, 'n_blocks': k,
           'n_hyp_genes': sum(g in vn for g in HYP), 'n_gln_genes': sum(g in vn for g in GLN11)}
    obs = {}
    for kx in keys:
        rho = stats.spearmanr(r_h, r[kx])[0]
        obs[kx] = rho
        mu, sd = null[kx].mean(), null[kx].std(ddof=1)
        out[f'rho_{kx}'] = rho; out[f'Z_{kx}'] = (rho - mu) / sd
        out[f'p_{kx}'] = (np.sum(np.abs(null[kx]) >= abs(rho)) + 1) / (N_PERM + 1)
    obs['D'] = obs['Synthesis'] - obs['Catabolism']
    mu, sd = null['D'].mean(), null['D'].std(ddof=1)
    out['rho_D'] = obs['D']; out['Z_D'] = (obs['D'] - mu) / sd
    out['p_D'] = (np.sum(np.abs(null['D']) >= abs(obs['D'])) + 1) / (N_PERM + 1)
    return out

# ------------------------------------------------- rank-trick sanity check
def sanity_check():
    ad, pos = load_section('UKF243_T_ST', 'Ravi')
    vn = ad.var_names; X = ad.X
    h = amean(X, HYP, vn); q = amean(X, GLN11, vn)
    r_h, r_q = rankdata(h), rankdata(q)
    assert abs(stats.spearmanr(h, q)[0] - stats.spearmanr(r_h, r_q)[0]) < 1e-12
    print('rank-trick sanity check passed')

sanity_check()

# ---------------------------------------------------------------- run
prim_gw = cohort_rows('Greenwald_Visium', 'include_tumor_primary')
prim_rv = cohort_rows('Ravi_Visium', 'include_tumor_primary')
idhm = cohort_rows('Greenwald_Visium', 'include_exploratory_IDHm')
idhm_rv = cohort_rows('Ravi_Visium', 'include_exploratory_IDHm')

tasks = [(r.section_id, r.cohort, 'primary') for r in prim_gw.itertuples()] + \
        [(r.section_id, r.cohort, 'primary') for r in prim_rv.itertuples()] + \
        [(r.section_id, r.cohort, 'exploratory_IDHm') for r in idhm.itertuples()] + \
        [(r.section_id, r.cohort, 'exploratory_IDHm') for r in idhm_rv.itertuples()]

CACHE = OUT / 'spatial_section_stats.csv'
if CACHE.exists():
    sec_df = pd.read_csv(CACHE)
    print(f'loaded {len(sec_df)} cached section rows')
else:
    sec_rows = []
    for i, (sec, coh, role) in enumerate(tasks):
        rec = section_stats(sec, coh)
        rec['analysis_role'] = role
        man = m[m.section_id == sec].iloc[0]
        rec['global_patient_id'] = man.global_patient_id
        rec['sensitivity_exclusion_tags'] = man.sensitivity_exclusion_tags if isinstance(man.sensitivity_exclusion_tags, str) else ''
        rec['ravi_greenwald_qc_subset'] = man.ravi_greenwald_qc_subset
        sec_rows.append(rec)
        print(f'[{i+1}/{len(tasks)}] {sec} ({coh}, {role}): rho_gln={rec["rho_rho_gln"]:.3f} Z={rec["Z_rho_gln"]:.2f} D={rec["rho_D"]:.3f}')
    sec_df = pd.DataFrame(sec_rows)
    sec_df.to_csv(CACHE, index=False)

# ------------------------------------------------- patient aggregation (median)
STAT_LIST = ['rho_gln','Synthesis','Catabolism','Direct_uptake','Exchange_redox',
             'CN_redistribution','D']
agg_map = {f'median_rho_{s}': (f'rho_{s}', 'median') for s in STAT_LIST}
agg_map.update({f'median_Z_{s}': (f'Z_{s}', 'median') for s in STAT_LIST})
pat = sec_df.groupby(['analysis_role','cohort','global_patient_id']).agg(
    n_sections=('section_id','count'), **agg_map).reset_index()
pat.to_csv(OUT / 'spatial_patient_stats.csv', index=False)

# ---------------------------------------------------------------- tests
def wilcox(vals, label):
    vals = np.asarray(vals, float)
    d = {'comparison': label, 'n': len(vals), 'n_positive': int((vals > 0).sum()),
         'n_negative': int((vals < 0).sum()), 'median': float(np.median(vals))}
    nz = vals[vals != 0]
    try:
        method = 'exact' if len(nz) < 25 else 'approx'
        w = stats.wilcoxon(vals, alternative='two-sided', method=method)
        d['p_wilcoxon'] = float(w.pvalue); d['method'] = method
    except ValueError:
        d['p_wilcoxon'] = 1.0; d['method'] = 'all-zero'
    return d

test_rows = []
def patient_vals(df, col):
    return df[f'median_rho_{col}'].values

for coh in ['Greenwald_Visium', 'Ravi_Visium']:
    sub = pat[(pat.analysis_role == 'primary') & (pat.cohort == coh)]
    for stat in STAT_LIST:
        test_rows.append(wilcox(patient_vals(sub, stat), f'primary|{coh}|{stat}|all'))

# Ravi leave-one-patient-out (primary endpoint rho_gln and D)
rv = pat[(pat.analysis_role=='primary')&(pat.cohort=='Ravi_Visium')]
for stat in ['rho_gln','D']:
    for p in rv.global_patient_id:
        sub = rv[rv.global_patient_id != p]
        test_rows.append(wilcox(patient_vals(sub, stat), f'primary|Ravi_Visium|{stat}|LOPO_excl_{p}'))

# Sensitivity: Ravi restricted to Greenwald-QC subset
qc_secs = sec_df[(sec_df.analysis_role=='primary')&(sec_df.cohort=='Ravi_Visium')&
                 (sec_df.ravi_greenwald_qc_subset.astype(str).str.lower()=='true')]
qc_pats = qc_secs.groupby('global_patient_id').agg(
    n_sections=('section_id','count'),
    **{f'median_rho_{s}': (f'rho_{s}','median') for s in ['rho_gln','D']}).reset_index()
for stat in ['rho_gln','D']:
    test_rows.append(wilcox(qc_pats[f'median_rho_{stat}'].values,
                            f'sensitivity|Ravi_QCsubset|{stat}|all'))

test_df = pd.DataFrame(test_rows)
test_df.to_csv(OUT / 'spatial_tests_summary.csv', index=False)

# ---------------------------------------------------------------- Lin per-patient (deterministic)
NEO = ['AC-like','MES-like','OPC-like','NPC-like']
ad = sc.read_h5ad(LIN_H5AD)
ad.var_names = ad.var['features'].astype(str).values; ad.var_names_make_unique()
ad = ad[ad.obs['annotation'].isin(NEO)].copy()
counts = ad.obs['Patient'].value_counts()
eligible = sorted(counts[counts >= 50].index.tolist())
ad = ad[ad.obs['Patient'].isin(eligible),
        [g for g in sorted(set(HYP + GLN11 + [g for gs in FATE.values() for g in gs])) if g in ad.var_names]].copy()
vn = ad.var_names; X = ad.X
hyp = amean(X, HYP, vn)
lin_rows = []
for p_ in eligible:
    msk = (ad.obs['Patient'] == p_).values
    rec = {'Patient': p_, 'n_cells': int(msk.sum())}
    for name, gl in [('rho_gln', GLN11)] + list(FATE.items()):
        s = amean(X, gl, vn)
        rec[f'rho_{name}'], rec[f'p_{name}'] = stats.spearmanr(hyp[msk], s[msk])
    rec['rho_D'] = rec['rho_Synthesis'] - rec['rho_Catabolism']
    lin_rows.append(rec)
lin_df = pd.DataFrame(lin_rows)
lin_df.to_csv(OUT / 'lin_patient_stats.csv', index=False)
del ad

# Lin sensitivities (manifest tags)
lin_m = m[(m.cohort=='Lin_scRNA') & (m.include_tumor_primary=='True')]
tag_map = dict(zip(lin_m.source_patient_id, lin_m.sensitivity_exclusion_tags.fillna('')))
lin_tests = []
def lin_sub(drop_tag, label):
    keep = [p_ for p_ in lin_df.Patient if drop_tag not in tag_map.get(p_, '')]
    sub = lin_df[lin_df.Patient.isin(keep)]
    for stat in ['rho_gln','D']:
        lin_tests.append(wilcox(sub[f'rho_{stat}'].values, f'sensitivity|Lin|{stat}|{label}'))
lin_sub('excl_known_or_possible_recurrent', 'excl_P77_recurrent')
lin_sub('excl_equivocal_IDH', 'excl_P92_equivocal_IDH')
lin_sub('excl_BT346_hypoxic_core', 'excl_BT346_hypoxic_core')
# Lin primary (all 28)
for stat in ['rho_gln','D']:
    lin_tests.append(wilcox(lin_df[f'rho_{stat}'].values, f'primary|Lin_scRNA|{stat}|all'))
pd.DataFrame(lin_tests).to_csv(OUT / 'lin_tests_summary.csv', index=False)

# ---------------------------------------------------------------- IDHm descriptive (no P values)
idhm_pat = pat[pat.analysis_role=='exploratory_IDHm']
desc = idhm_pat[['cohort','global_patient_id','n_sections'] +
                [f'median_{c}' for c in ['rho_rho_gln','rho_D','rho_Synthesis','rho_Catabolism']]]
desc.to_csv(OUT / 'idhm_descriptive.csv', index=False)

# ---------------------------------------------------------------- figure data (per cohort x component)
def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); b = p[o] * n / np.arange(1, n + 1)
    b = np.minimum.accumulate(b[::-1])[::-1]
    r = np.empty(n); r[o] = np.clip(b, 0, 1); return r

fig_rows = []
# Lin (cell-level rho per patient; Wilcoxon on patients)
for stat in ['rho_gln'] + list(FATE.keys()):
    v = lin_df[f'rho_{stat}'].values
    t = wilcox(v, f'fig|Lin|{stat}')
    q = bh([wilcox(lin_df[f'rho_{c}'].values, '')['p_wilcoxon']
            for c in FATE.keys()])[list(FATE.keys()).index(stat)] if stat in FATE else np.nan
    fig_rows.append({'cohort': 'Lin_scRNA (28 patients)', 'component': stat, 'n': len(v),
                     'n_positive': t['n_positive'], 'median': t['median'],
                     'p': t['p_wilcoxon'], 'q_BH': q})
# Greenwald / Ravi (patient medians, from spatial_tests_summary)
for coh, label in [('Greenwald_Visium', 'Greenwald Visium (5 patients)'),
                   ('Ravi_Visium', 'Ravi Visium (17 patients)')]:
    ps = {r['comparison'].split('|')[2]: r['p_wilcoxon'] for r in test_rows
          if r['comparison'] == f'primary|{coh}|' + r['comparison'].split('|')[2] + '|all'}
    sub = pat[(pat.analysis_role == 'primary') & (pat.cohort == coh)]
    qmap = dict(zip(list(FATE.keys()), bh([ps[c] for c in FATE.keys()])))
    for stat in ['rho_gln'] + list(FATE.keys()) + ['D']:
        v = sub[f'median_rho_{stat}'].values
        fig_rows.append({'cohort': label, 'component': stat, 'n': len(v),
                         'n_positive': int((v > 0).sum()), 'median': float(np.median(v)),
                         'p': ps[stat], 'q_BH': qmap.get(stat, np.nan)})
# GSC (deltas vs normoxia; gateB_summary)
    gb = pd.read_csv(RESULTS_DIR / 'gateB_per_line_delta_11gene.csv')
    gsum = pd.read_csv(RESULTS_DIR / 'gateB_summary.csv')
name_map = {'OverallGlutamine': 'rho_gln', 'Synthesis': 'Synthesis', 'Catabolism': 'Catabolism',
            'Direct_uptake': 'Direct_uptake', 'Exchange_redox': 'Exchange_redox',
            'CN_redistribution': 'CN_redistribution'}
for _, r in gsum.iterrows():
    fig_rows.append({'cohort': 'GSC 24h hypoxia (8 lines)', 'component': name_map[r.Component],
                     'n': int(r.n_total), 'n_positive': int(r.n_positive),
                     'median': float(r.median_delta), 'p': float(r.wilcoxon_p_twosided),
                     'q_BH': float(r.q_BH_twosided) if pd.notna(r.q_BH_twosided) else np.nan})
# U87 (hypoxia vs normoxia t-test; gateC_u87_summary)
    usum = pd.read_csv(RESULTS_DIR / 'gateC_u87_summary.csv')
    usco = pd.read_csv(RESULTS_DIR / 'gateC_u87_scores.csv')
for _, r in usum.iterrows():
    comp = name_map[r.Component]
    d = usco[usco.Condition == 'hypoxia'][f'Score_{r.Component}'].values
    fig_rows.append({'cohort': 'U87MG (3v3)', 'component': comp, 'n': len(d),
                     'n_positive': int((d > 0).sum()), 'median': float(np.median(d)),
                     'p': float(r.t_p_twosided), 'q_BH': float(r.q_BH_twosided)})
pd.DataFrame(fig_rows).to_csv(OUT / 'figure_data_v5.csv', index=False)

print('\n=== HEADLINES (provisional until Phase C QC) ===')
for label, df_ in [('Greenwald n=5', pat[(pat.analysis_role=='primary')&(pat.cohort=='Greenwald_Visium')]),
                   ('Ravi n=17', rv), ('Ravi QC n=13', qc_pats)]:
    for stat in ['rho_gln','D']:
        v = df_[f'median_rho_{stat}'].values
        t = [r for r in test_rows if r['comparison'] == f'primary|{label.split()[0]}_Visium|{stat}|all'] if label != 'Ravi QC n=13' else \
            [r for r in test_rows if r['comparison'] == f'sensitivity|Ravi_QCsubset|{stat}|all']
        print(f'{label} {stat}: {int((v>0).sum())}/{len(v)} positive, median {np.median(v):.3f}, P={t[0]["p_wilcoxon"]:.4g}')
print('\nDone. Outputs in', OUT)
