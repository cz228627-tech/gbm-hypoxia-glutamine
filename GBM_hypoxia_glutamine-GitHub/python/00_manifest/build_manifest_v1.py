"""Phase A (G1-G4): build and freeze cohort_manifest_v1.csv + cross_cohort_overlap_check.csv.
Row = analysable unit (Visium section / scRNA patient / GSC line / U87 model).
Evidence strings cite source tables saved in 10_cohort_audit/source_evidence/. Nothing is imputed.
"""
import hashlib, os
import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("GBM_DATA_DIR", REPO_ROOT / "data")).expanduser()
OUT = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser() / "10_cohort_audit"
OUT.mkdir(parents=True, exist_ok=True)

COLS = ["cohort", "section_id", "sample_id", "source_patient_id", "global_patient_id", "institution",
        "source_diagnosis", "IDH_status", "IDH_evidence", "IDH_evidence_strength", "primary_recurrent",
        "age", "sex", "age_group", "evidence_tier", "disease_group",
        "include_tumor_primary", "include_perturbation_validation", "include_boundary",
        "include_exploratory_IDHm", "sensitivity_exclusion_tags", "ravi_greenwald_qc_subset",
        "cross_cohort_overlap", "flags", "notes"]
BOOL = ["include_tumor_primary", "include_perturbation_validation", "include_boundary",
        "include_exploratory_IDHm", "ravi_greenwald_qc_subset"]
rows = []


def add(**kw):
    r = {c: "" for c in COLS}
    r.update(kw)
    for b in BOOL:
        r[b] = bool(r[b]) if r[b] != "" else False
    rows.append(r)


# ============ Lin scRNA ============
lin_counts = pd.read_csv(DATA_DIR / "cohort_audit" / "lin_patient_counts.csv")
lin_s1 = pd.read_excel(DATA_DIR / "cohort_audit" / "lin_S1.xlsx", header=None)
s1 = {}
for _, r in lin_s1.iterrows():
    pid = str(r[0]).strip()
    if pid.startswith(("P", "BT", "T2", "T3", "C3")):
        s1[pid] = dict(epi=str(r[11]), sex=str(r[12]), age=str(r[13]), dx=str(r[14]), ihc=str(r[15]))
LIN_OVL = "none (Montreal/Hefei/CPTAC donors; no Freiburg, Zurich or Boston donors)"

for _, r in lin_counts.iterrows():
    pid, ds, nneo = r["Patient"], r["Dataset_id"], int(r["n_neo"])
    elig = nneo >= 50
    if ds == "PMID:32641768":
        add(cohort="Lin_scRNA", sample_id=pid, source_patient_id=pid, global_patient_id=f"MNH:{pid}",
            institution="Montreal Neurological Hospital (Couturier 2020)",
            source_diagnosis="IDH-wildtype glioblastoma (Couturier 2020 Table 1; Lin S1 'IDH- GBM')",
            IDH_status="wildtype", IDH_evidence="IDH status by DNA sequencing, certified neuropathologist (Couturier 2020 Methods/Table 1)",
            IDH_evidence_strength="confirmed", primary_recurrent="not_reported_per_case",
            age="not_reported_per_case", age_group="adult (Couturier: adult glioblastoma, cohort mean 62.3 y)",
            evidence_tier=1, disease_group="GBM", include_tumor_primary=elig,
            sensitivity_exclusion_tags=("excl_BT346_hypoxic_core" if pid == "BT346" else ""),
            cross_cohort_overlap=LIN_OVL,
            flags=("hypoxic_core_cells_likely_retained" if pid == "BT346" else ""),
            notes=f"n_neoplastic={nneo}."
                  + (" Couturier excluded necrotic-core hypoxic BT346 cells (reported 1902 tumour cells); Lin h5ad holds "
                     "5032 neoplastic BT346 cells, so these cells are likely retained." if pid == "BT346" else ""))
    elif ds == "homemade":
        if pid == "unknown":
            add(cohort="Lin_scRNA", sample_id="unknown", source_patient_id="unknown", global_patient_id="USTC:unknown",
                institution="First Affiliated Hospital of USTC", source_diagnosis="cells without patient assignment",
                IDH_status="unknown", IDH_evidence="none", IDH_evidence_strength="unreported", disease_group="GBM",
                cross_cohort_overlap=LIN_OVL, flags="excluded_unassigned", notes=f"n_neoplastic={nneo}; excluded.")
            continue
        i = s1[pid]
        idh_ihc = [x.strip() for x in i["ihc"].replace("，", ",").split(",") if "IDH" in x]
        tags, flag = [], []
        if pid == "P77":
            tags.append("excl_known_or_possible_recurrent"); flag.append("possible_recurrent")
        if pid == "P92":
            tags.append("excl_equivocal_IDH"); flag.append("IDH1_IHC_equivocal")
        if not elig:
            flag.append("excluded_lt50_neoplastic")
        add(cohort="Lin_scRNA", sample_id=pid, source_patient_id=pid, global_patient_id=f"USTC:{pid}",
            institution="First Affiliated Hospital of USTC, Hefei (Lin 2026 in-house)",
            source_diagnosis=f"{i['dx']} (Lin 2026 Supp. Table 1)", IDH_status="wildtype (source-defined)",
            IDH_evidence="IHC " + "; ".join(idh_ihc) + " (Lin S1)",
            IDH_evidence_strength=("equivocal" if pid == "P92" else "source-defined"),
            primary_recurrent=(f"possibly_recurrent (S1 epilepsy field: '{i['epi']}')" if pid == "P77"
                               else "not_reported_per_case (Lin: 'primary GBM tissue')"),
            age=i["age"], sex=i["sex"], age_group="adult", evidence_tier=1, disease_group="GBM",
            include_tumor_primary=elig, sensitivity_exclusion_tags=";".join(tags),
            cross_cohort_overlap=LIN_OVL, flags=";".join(flag), notes=f"n_neoplastic={nneo}.")
    else:
        cid = pid.split("_")[0]
        i = s1.get(cid, {})
        add(cohort="Lin_scRNA", sample_id=pid, source_patient_id=cid, global_patient_id=f"CPTAC:{cid}",
            institution="CPTAC (Wang 2021)", source_diagnosis=f"{i.get('dx','')} (Lin S1)",
            IDH_status="wildtype (source-defined)", IDH_evidence="Lin S1", IDH_evidence_strength="source-defined",
            age=i.get("age", ""), sex=i.get("sex", ""), age_group="adult", evidence_tier=1, disease_group="GBM",
            cross_cohort_overlap=LIN_OVL, flags="excluded_0_neoplastic", notes=f"snRNA; n_neoplastic={nneo}; excluded.")

# ============ Greenwald Visium (own cohort only; GSE237183) ============
GW_S1 = "Greenwald 2024 Cell Table S1 (mmc1)"
gw_pat = {"MGH258": ("F", 43, "methylated", "Massachusetts General Hospital"),
          "ZH881": ("M", 79, "partially methylated", "University Hospital Zurich"),
          "ZH916": ("M", 62, "methylated", "University Hospital Zurich"),
          "ZH1007": ("F", 50, "unmethylated", "University Hospital Zurich"),
          "ZH1019": ("F", 67, "methylated", "University Hospital Zurich")}
gw_gbm = [("mgh258", "GSM7596587", "MGH258", "MGH258", "bulk"),
          ("zh881inf", "GSM7596588", "ZH881", "ZH881_inf", "infiltrating"),
          ("zh881t1", "GSM7596589", "ZH881", "ZH881_T1", "T1 contrast-enhancing"),
          ("zh916bulk", "GSM7596590", "ZH916", "ZH916_bulk", "bulk"),
          ("zh916inf", "GSM7596591", "ZH916", "ZH916_inf", "infiltrating"),
          ("zh916t1", "GSM7596592", "ZH916", "ZH916_T1", "T1 contrast-enhancing"),
          ("zh1007inf", "GSM7596593", "ZH1007", "ZH1007_inf", "infiltrating"),
          ("zh1007nec", "GSM7596594", "ZH1007", "ZH1007_nec", "necrotic"),
          ("zh1019inf", "GSM7596595", "ZH1019", "ZH1019_inf", "infiltrating"),
          ("zh1019t1", "GSM7596596", "ZH1019", "ZH1019_T1", "T1 contrast-enhancing"),
          ("zh8811a", "GSM7596597", "ZH881", "ZH881_1A", "bulk"),
          ("zh8811b", "GSM7596598", "ZH881", "ZH881_1B", "bulk"),
          ("zh8812", "GSM7596599", "ZH881", "ZH881_2", "bulk")]
for sec, gsm, pat, s1id, region in gw_gbm:
    sex, age, mgmt, ins = gw_pat[pat]
    note = f"S1 sample {s1id}; region={region}; MGMT={mgmt}."
    flag = ""
    if pat == "ZH1019":
        note += (" S1 Patient-ID column reads 'ZH1007' for ZH1019 rows, but age (67 vs 50), MGMT (methylated vs "
                 "unmethylated) and location (left temporal vs left occipital) differ from ZH1007 -> treated as a "
                 "typo; ZH1019 is a distinct patient (consistent with GEO titles and 'four cases' in the text).")
        flag = "S1_patient_id_typo_resolved"
    add(cohort="Greenwald_Visium", section_id=sec, sample_id=gsm, source_patient_id=pat, global_patient_id=f"GW:{pat}",
        institution=ins, source_diagnosis="GBM, IDH WT, grade 4 (" + GW_S1 + ")", IDH_status="wildtype",
        IDH_evidence=GW_S1 + " 'IDH WT'; main text '13 IDH-WT GBM sections'", IDH_evidence_strength="source-defined",
        primary_recurrent="not_reported (no column in S1)", age=age, sex=sex, age_group="adult",
        evidence_tier=1, disease_group="GBM", include_tumor_primary=True,
        cross_cohort_overlap="none (Zurich/Boston donors; Greenwald's use of Ravi data is external integration, not in our files)",
        flags=flag, notes=note)

gw_idhm = [("bwh23", "GSM7596600", "BWH23", "IDH-mut oligodendroglioma", "4", "M", 32, "Brigham and Women's Hospital"),
           ("bwh24", "GSM7596601", "BWH24", "IDH-mut astrocytoma", "4", "F", 54, "Brigham and Women's Hospital"),
           ("bwh25", "GSM7596602", "BWH25", "IDH-mut astrocytoma", "4", "F", 51, "Brigham and Women's Hospital"),
           ("bwh28", "GSM7596603", "BWH28", "IDH-mut astrocytoma", "2/3", "M", 26, "Brigham and Women's Hospital"),
           ("bwh35", "GSM7596604", "BWH35", "IDH-mut oligodendroglioma", "2", "F", 36, "Brigham and Women's Hospital"),
           ("mgh259", "GSM7596605", "MGH259", "IDH-mut oligodendroglioma", "3", "F", 40, "Massachusetts General Hospital")]
for sec, gsm, pat, dx, grade, sex, age, ins in gw_idhm:
    add(cohort="Greenwald_Visium", section_id=sec, sample_id=gsm, source_patient_id=pat, global_patient_id=f"GW:{pat}",
        institution=ins, source_diagnosis=f"{dx}, grade {grade} ({GW_S1})", IDH_status="mutant",
        IDH_evidence=GW_S1, IDH_evidence_strength="source-defined", primary_recurrent="not_reported",
        age=age, sex=sex, age_group="adult", evidence_tier="", disease_group="IDH-mutant glioma",
        include_exploratory_IDHm=True, cross_cohort_overlap="none", flags="exploratory_only",
        notes="Descriptive only; no significance test. Grade per S1 (GEO grade differs for BWH23: 3).")

# ============ Ravi Visium (Zenodo 6972429) ============
RV_S1 = "Ravi 2022 Cancer Cell Table S1 (mmc1)"
ravi_s1 = {  # UKF: (age, sex, tumor, histology) from Ravi Table S1
    "275": (80, "f", "IDH-WT", "GBM"), "270": (41, "m", "IDH-MUT", "AA"), "269": (51, "m", "IDH-WT", "GBM"),
    "268": (49, "f", "IDH-MUT", "ODG"), "266": (76, "f", "IDH-WT", "GBM"), "265": (55, "m", "IDH-WT", "GBM"),
    "262": (59, "m", "IDH-WT", "GBM"), "260": (81, "m", "IDH-WT", "GBM"), "259": (75, "f", "IDH-WT", "GBM"),
    "256": (64, "f", "IDH-WT", "GBM"), "255": (81, "f", "IDH-WT", "GBM"), "251": (61, "m", "IDH-WT", "GBM"),
    "248": (44, "m", "IDH-WT", "GBM"), "243": (55, "f", "IDH-WT", "GBM"), "242": (81, "f", "IDH-WT", "GBM"),
    "241": (19, "m", "C", "Cortex"), "296": (34, "m", "IDH-WT", "GBM"), "304": (81, "m", "IDH-WT", "GBM"),
    "313": (57, "m", "IDH-WT", "GBM"), "334": (73, "f", "IDH-WT", "GBM")}
GW_RAVI_QC = {"255", "313", "243", "260", "334", "248", "259", "275", "296", "251", "266", "269", "304"}
lin_ravi = {k[1:] for k in s1 if k.startswith("T") and k[1:].isdigit()}
RV_OVL = "none (Freiburg donors; Lin reused 12 Ravi Visium sections but those are not part of our Lin scRNA set)"
ravi_dir = Path(os.environ.get("GBM_RAVI_DIR", DATA_DIR / "ravi_visium")).expanduser()
ravi_secs = sorted(f.split("__")[0] for f in os.listdir(ravi_dir) if f.endswith("filtered_feature_bc_matrix.h5"))
for sec in ravi_secs:
    num = sec.split("_")[0].replace("UKF", "")
    part = sec.split("_")[1]
    age, sex, tum, hist = ravi_s1[num]
    base = dict(cohort="Ravi_Visium", section_id=sec, sample_id=sec, source_patient_id=f"UKF{num}",
                global_patient_id=f"RAVI:UKF{num}", institution="University Clinic Freiburg", age=age, sex=sex.upper(),
                cross_cohort_overlap=RV_OVL)
    if part == "C":
        add(**base, source_diagnosis="access cortex (" + RV_S1 + ")", IDH_status="NA", IDH_evidence=RV_S1,
            IDH_evidence_strength="NA", age_group=("adult" if age >= 18 else "paediatric"),
            disease_group="non-tumor cortex", flags="not_analysed", notes="Cortex section; not analysed.")
    elif part == "IDHMutant":
        add(**base, source_diagnosis=f"IDH-mutant {'astrocytoma (AA)' if hist == 'AA' else 'oligodendroglioma (ODG)'} ({RV_S1})",
            IDH_status="mutant", IDH_evidence=RV_S1, IDH_evidence_strength="source-defined",
            primary_recurrent="not_reported", age_group="adult", evidence_tier="", disease_group="IDH-mutant glioma",
            include_exploratory_IDHm=True, flags="exploratory_only", notes="Descriptive only; grade not reported.")
    else:
        assert tum == "IDH-WT" and hist == "GBM"
        tags = [] if num in GW_RAVI_QC else ["excl_not_in_greenwald_qc_subset"]
        add(**base, source_diagnosis=f"GBM, IDH-WT ({RV_S1})", IDH_status="wildtype",
            IDH_evidence=RV_S1 + " 'IDH-WT'" + ("; also Lin 2026 S1 'IDH- GBM'" if num in lin_ravi else "")
                         + ("; also Greenwald 2024 S1 'IDH WT'" if num in GW_RAVI_QC else ""),
            IDH_evidence_strength="source-defined", primary_recurrent="not_reported (no column in S1)",
            age_group="adult", evidence_tier=1, disease_group="GBM", include_tumor_primary=True,
            sensitivity_exclusion_tags=";".join(tags), ravi_greenwald_qc_subset=(num in GW_RAVI_QC),
            flags=("" if num in lin_ravi else "not_in_Lin_S1(irrelevant_after_Ravi_S1)"),
            notes=f"region={part}.")

# ============ GSE243848 GSC ============
for g in ["G33", "G34", "G44", "G62", "G88", "G146", "G157", "G816"]:
    add(cohort="GSE243848_GSC", sample_id=g, source_patient_id=g, global_patient_id=f"GSC:{g}",
        institution="GSC panel of the Chiocca/Godlewski/Bronisz groups (same line names as Mineo 2016); GEO depositor Mossakowski Institute, Warsaw",
        source_diagnosis="patient-derived glioblastoma stem-like cells (GEO: 'tissue: glioblastoma; cell type: GSCs')",
        IDH_status="not_reported", IDH_evidence="GEO GSE243848 sample characteristics contain no IDH field",
        IDH_evidence_strength="unreported", primary_recurrent="not_reported", age="not_reported", sex="not_reported",
        age_group="not_reported", evidence_tier=2, disease_group="GBM (source-annotated)",
        include_perturbation_validation=True,
        cross_cohort_overlap="none_documented (US-derived lines; donor identities not public, so overlap cannot be formally excluded)",
        flags="IDH_not_reported_per_line", notes="Vehicle arm only, normoxia vs 1% O2 for 24 h.")

# ============ U87MG ============
add(cohort="GSE300765_U87MG", sample_id="U87MG", source_patient_id="U87MG", global_patient_id="MODEL:U87MG",
    institution="established legacy cell line", source_diagnosis="legacy glioblastoma-derived cell line",
    IDH_status="not_patient_classifiable", IDH_evidence="NA", IDH_evidence_strength="NA", evidence_tier=3,
    disease_group="model", include_boundary=True, cross_cohort_overlap="NA (model)", flags="boundary_only",
    notes="3 vs 3 replicates, 48 h 1% O2.")

m = pd.DataFrame(rows, columns=COLS)

# ---------- integrity checks ----------
tri = m[["include_tumor_primary", "include_perturbation_validation", "include_boundary"]].sum(axis=1)
assert (tri <= 1).all(), "primary logic fields not mutually exclusive"
assert not (m.include_exploratory_IDHm & (tri > 0)).any()
assert not (m.include_tumor_primary & (m.disease_group != "GBM")).any()
assert not (m.include_tumor_primary & (m.IDH_status.str.startswith("wildtype") == False)).any()
assert not (m.include_tumor_primary & (m.age_group != "adult") & (m.cohort != "Lin_scRNA")).any()
# a global_patient_id must map to exactly one cohort
assert (m.groupby("global_patient_id").cohort.nunique() == 1).all()

m.to_csv(f"{OUT}/cohort_manifest_v1.csv", index=False)

# ---------- G3 cross-cohort overlap ----------
ov = pd.DataFrame([
    ["Lin_scRNA (28)", "Ravi_Visium (17)", "Institution: Montreal/Hefei vs Freiburg. Lin reused 12 Ravi Visium sections for its own spatial work, but our Lin analysis uses only the scRNA set (Couturier BT + USTC P). No UKF ID in the scRNA obs.", "none"],
    ["Lin_scRNA (28)", "Greenwald_Visium GBM (5)", "Institution: Montreal/Hefei vs Zurich/MGH. Lin S1 does not list Greenwald samples.", "none"],
    ["Ravi_Visium (17)", "Greenwald_Visium GBM (5)", "Institution: Freiburg vs Zurich/MGH. Greenwald integrated 13 Ravi sections as an external cohort in its own analysis (Greenwald S1 'External cohort from Ravi'). GSE237183 files we use contain only ZH/MGH sections.", "none"],
    ["Couturier BT (13)", "USTC P (15)", "Different institutions (Montreal vs Hefei).", "none"],
    ["GSE243848 GSC (8)", "all tumour cohorts", "GSC lines are US-derived (Chiocca group). They do not share an institution with Freiburg, Zurich, Montreal or Hefei. MGH/BWH donors: not formally excludable because GSC donor identities are not public.", "none_documented"],
    ["within Greenwald", "zh881inf, zh881t1, zh8811a, zh8811b, zh8812", "Greenwald S1 Patient-ID = ZH881 for all five (M, 79 y, right temporal).", "same patient -> GW:ZH881"],
    ["within Greenwald", "zh1019inf/t1 vs zh1007inf/nec", "S1 Patient-ID typo ('ZH1007' on ZH1019 rows). Age, MGMT and location differ.", "distinct patients"],
    ["within Ravi", "UKF256_TC, UKF256_TI", "Same UKF number (README: number = patient).", "same patient -> RAVI:UKF256"],
], columns=["cohort_A", "cohort_B", "evidence", "conclusion"])
ov.to_csv(f"{OUT}/cross_cohort_overlap_check.csv", index=False)

md5 = hashlib.md5(open(f"{OUT}/cohort_manifest_v1.csv", "rb").read()).hexdigest()
tp = m[m.include_tumor_primary]
print(tp.groupby("cohort").agg(n_units=("global_patient_id", "size"), n_patients=("global_patient_id", "nunique")).to_string())
print("IDHm exploratory patients:", sorted(m[m.include_exploratory_IDHm].global_patient_id.unique()))
print("perturbation lines:", int(m.include_perturbation_validation.sum()), "| boundary:", int(m.include_boundary.sum()))
print("sensitivity tags:", m.sensitivity_exclusion_tags[m.sensitivity_exclusion_tags != ""].value_counts().to_dict())
print("rows:", len(m), "| md5:", md5)
