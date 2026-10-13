#!/usr/bin/env python3
"""Exact label-permutation sensitivity analysis for the U87MG (GSE300765) 3-vs-3 experiment.

With n = 3 per group there are only C(6,3) = 20 distinct group assignments, so the
smallest attainable two-sided exact permutation P value is 2/20 = 0.10. This script
enumerates all 20 splits, computes the two-sided exact permutation P for each
glutamine-fate component score, and writes u87_exact_permutation.csv.

Input : gateC_u87_scores.csv (per-sample fate scores; produced by gateC_gse300765.py)
Output: u87_exact_permutation.csv

Statistic: |mean(hypoxia) - mean(normoxia)| evaluated over every assignment of the
6 samples into two groups of 3; P = fraction of assignments with statistic >= observed.
"""
from itertools import combinations
from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = Path(os.environ.get("GBM_OUTPUT_DIR", REPO_ROOT / "results")).expanduser()
IN = Path(os.environ.get("GBM_U87_SCORES", RESULTS_DIR / "gateC_u87_scores.csv")).expanduser()
OUT = RESULTS_DIR / "11_spatial_v5" / "u87_exact_permutation.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

COMPONENTS = [
    "Score_Synthesis",
    "Score_Catabolism",
    "Score_Direct_uptake",
    "Score_Exchange_redox",
    "Score_CN_redistribution",
    "Score_OverallGlutamine",
]


def main() -> None:
    df = pd.read_csv(IN)
    hyp = df[df["Condition"] == "hypoxia"]
    nor = df[df["Condition"] == "normoxia"]
    assert len(hyp) == 3 and len(nor) == 3, "expected 3 vs 3"

    pooled = np.vstack([hyp[COMPONENTS].to_numpy(), nor[COMPONENTS].to_numpy()])
    obs_stat = np.abs(pooled[:3].mean(0) - pooled[3:].mean(0))

    rows = []
    for idx in combinations(range(6), 3):
        idx = list(idx)
        mask = np.zeros(6, dtype=bool)
        mask[idx] = True
        stat = np.abs(pooled[mask].mean(0) - pooled[~mask].mean(0))
        rows.append(stat)
    all_stats = np.array(rows)

    out = []
    for j, comp in enumerate(COMPONENTS):
        n_ge = int((all_stats[:, j] >= obs_stat[j] - 1e-12).sum())
        p_exact = n_ge / 20.0
        t, p_t = stats.ttest_ind(pooled[:3, j], pooled[3:, j], equal_var=True)
        out.append(
            {
                "component": comp.replace("Score_", ""),
                "mean_diff": pooled[:3, j].mean() - pooled[3:, j].mean(),
                "t_p": p_t,
                "exact_perm_p": p_exact,
                "n_perm_ge": n_ge,
                "n_splits": 20,
            }
        )

    res = pd.DataFrame(out)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT, index=False)
    print(res.to_string(index=False))
    print(f"\nMinimum attainable two-sided exact P with 3v3: 2/20 = 0.10")
    print(f"Written: {OUT}")


if __name__ == "__main__":
    main()
