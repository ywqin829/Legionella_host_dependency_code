"""Correlation-robust sensitivity for the LCV -> fitness-filtered-CRISPR coupling.

Motivation
----------
The headline evidence-layer coupling (LCV-associated detection vs fitness-filtered
CRISPR evidence) is reported with a gene-label permutation P and with a
multivariable logistic OR. Gene labels are not independent observations:
functionally related genes (Rab/ER-Golgi/trafficking modules) share expression
and annotation structure. A plain label permutation can therefore understate
uncertainty, and ordinary logistic standard errors can be too small.

This script adds two correlation-aware robustness analyses on top of the *frozen*
primary result. It does NOT alter candidate selection, tiers, frozen endpoints or
the published a1 coupling values (those stay untouched). It writes new outputs
only (corr_robust_logistic.tsv, corr_robust_summary.json).

Analyses
--------
(A) Cluster-robust (sandwich) inference for the *published* joint multivariable
    model (the one behind the headline "LCV OR = 9.96"): re-fits
    outcome ~ lcv + ppi + response + nstudy + chronos + common, clustered by GO
    biological-process functional module, and reports ordinary SE and
    cluster-robust SE side by side.

(B) Module-decimated matched test: from every functional module keep one random
    gene, rebuild the same confounder-matched enrichment, repeat B times. This
    removes within-module dependence by construction and reports the distribution
    of the matched fold and Fisher P.
"""
from __future__ import annotations

import collections
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
GAF = ROOT / "01_data_raw" / "annotations" / "goa_human.gaf.gz"

RNG = np.random.default_rng(42)
B_DECIM = 200


def B(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


# ---------------- Functional modules from GO biological process ----------------
def build_modules() -> dict:
    gene_terms = collections.defaultdict(set)
    with gzip.open(GAF, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[8] == "P":
                gene_terms[c[2].upper()].add(c[4])
    freq = collections.Counter()
    for v in gene_terms.values():
        for t in v:
            freq[t] += 1
    n_tot = sum(freq.values())

    def ic(t):
        return -math.log2(freq[t] / n_tot)

    # each gene -> its most specific (rarest) BP term
    return {g: max(ts, key=ic) for g, ts in gene_terms.items()}


def attach_modules(symbols: pd.Series) -> pd.Series:
    mod = symbols.astype(str).str.upper().map(module_map)
    nomod = mod.isna()
    mod.loc[nomod] = "NOMOD::" + symbols.loc[nomod].astype(str)
    return mod


module_map = build_modules()
atlas_raw = pd.read_csv(
    ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz",
    sep="\t", low_memory=False,
)


# ================= (A) cluster-robust joint model =================
# Faithful re-implementation of postfreeze_joint_evidence_model.R:
#   outcome ~ lcv + ppi + response + nstudy + chronos + common
# over the whole symbol-collapsed atlas (median per symbol; non-finite -> 0).
def build_joint_frame() -> pd.DataFrame:
    a = atlas_raw.copy()
    a["_outcome"] = a["functional_tier"].astype(str).str.startswith("F").astype(int)
    a["_lcv"] = B(a["lcv_spatial_supported"]).astype(int)
    a["_ppi"] = B(a["ppi_AB_flag"]).astype(int)
    a["_response"] = B(a["stable_cross_cohort"]).astype(int)
    a["_nstudy"] = pd.to_numeric(a["expression_n_studies"], errors="coerce")
    a["_chronos"] = pd.to_numeric(a["u937_chronos_gene_effect"], errors="coerce")
    a["_common"] = B(a["common_essential_flag"]).astype(int)
    cols = ["_outcome", "_lcv", "_ppi", "_response", "_nstudy", "_chronos", "_common"]
    g = a.groupby("human_symbol", as_index=False)[cols].median()
    g[cols] = g[cols].fillna(0.0)
    g["_module"] = attach_modules(g["human_symbol"])
    return g


def cr_joint(g: pd.DataFrame) -> list:
    terms = ["_lcv", "_ppi", "_response", "_nstudy", "_chronos", "_common"]
    labels = ["lcv", "ppi", "response", "nstudy", "chronos", "common"]
    y = g["_outcome"].astype(float).values
    X = sm.add_constant(g[terms].astype(float).values)
    groups = g["_module"].values
    model = sm.GLM(y, X, family=sm.families.Binomial())
    res_o = model.fit()
    res_c = model.fit(cov_type="cluster", cov_kwds={"groups": groups})
    names = ["const"] + labels
    rows = []
    for i, nm in enumerate(names):
        b_o, se_o = res_o.params[i], res_o.bse[i]
        b_c, se_c = res_c.params[i], res_c.bse[i]
        rows.append({
            "term": nm,
            "OR": float(np.exp(b_c)),
            "ci_lo_ordinary": float(np.exp(b_o - 1.96 * se_o)),
            "ci_hi_ordinary": float(np.exp(b_o + 1.96 * se_o)),
            "p_ordinary": float(res_o.pvalues[i]),
            "ci_lo_cluster": float(np.exp(b_c - 1.96 * se_c)),
            "ci_hi_cluster": float(np.exp(b_c + 1.96 * se_c)),
            "p_cluster": float(res_c.pvalues[i]),
        })
    return rows


# ================= (B) module-decimated matched test =================
# Universe + confounder set identical to postfreeze_a1_coupling.py.
atlas = atlas_raw.copy()
deg = (
    pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv",
                sep="\t", low_memory=False)
    [["gene", "cross_pathogen_degree", "cross_pathogen_records", "cross_pathogen_pmids"]]
    .drop_duplicates()
)
atlas["_ppi_AB"] = B(atlas["ppi_AB_flag"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fitness_ok"] = B(atlas["fitness_artifact_compatible"])
atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
     "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"]
)
atlas["_crispr_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_crispr_clean"] = atlas["_crispr_hit"] & (~atlas["_fitness_ok"])
atlas["_common_ess"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas["_nstudies"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)


def load_ui_percentile(matrix_name, ui_gsms, map_name, use_probe_map):
    mat = pd.read_csv(ROOT / "05_expression" / "matrices" / matrix_name, sep="\t", compression="gzip")
    cols = [c for c in ui_gsms if c in mat.columns]
    if not cols:
        return None
    id_col = "feature_id" if "feature_id" in mat.columns else "Geneid"
    gm = pd.read_csv(ROOT / "05_expression" / "gene_maps" / map_name, sep="\t")
    gm = gm[gm["gene_symbol"].notna()]
    if use_probe_map and "ambiguous_mapping" in gm.columns:
        gm = gm[~gm["ambiguous_mapping"].astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in gm.columns and id_col != "feature_id":
        gm = gm.rename(columns={"feature_id": id_col})
    m = mat.merge(gm, on=id_col, how="inner")
    sub = m[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    r = sub.rank(pct=True)
    return r[cols].mean(axis=1)


ui_percentiles = {
    "GSE61535": load_ui_percentile("GSE61535_log2_expression.tsv.gz",
                                   ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    "GSE135803": load_ui_percentile("GSE135803_fpkm.tsv.gz",
                                    ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    "GSE287755": load_ui_percentile("GSE287755_raw_counts.tsv.gz",
                                    ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
}
ui_stack = pd.concat(ui_percentiles.values(), axis=1, keys=list(ui_percentiles.keys()))
atlas["_ui_expr"] = atlas["human_symbol"].map(ui_stack.mean(axis=1, skipna=True))
atlas = atlas.merge(deg, left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)
atlas["_ui_expr"] = atlas["_ui_expr"].fillna(0.0)

atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
    "_ppi_AB": "max", "_lcv": "max", "_resp": "max", "_fitness_ok": "max",
    "_crispr_tested": "max", "_crispr_hit": "max", "_crispr_clean": "max",
    "_common_ess": "max", "_chronos": "median", "_nstudies": "max",
    "_ui_expr": "median", "_deg": "max", "_pmids": "max",
})
CONF = ["_deg", "_pmids", "_chronos", "_common_ess", "_ui_expr"]
atlas["_module"] = attach_modules(atlas["human_symbol"])


def universe(outcome_col):
    d = atlas.copy()
    if outcome_col.startswith("_crispr"):
        d = d[d["_crispr_tested"]]
    elif outcome_col == "_resp":
        d = d[d["_nstudies"] >= 1]
    for c in CONF:
        d = d[d[c].notna()]
    return d.reset_index(drop=True)


def matched_fold(d):
    """Reproduce the a1 confounder-matched enrichment (fold vs matched + Fisher P)."""
    d = d.reset_index(drop=True)
    X = np.column_stack([
        np.log1p(d["_deg"].clip(lower=0).replace(0, 1e-9)),
        np.log1p(d["_pmids"].clip(lower=0).replace(0, 1e-9)),
        d["_chronos"].astype(float),
        d["_common_ess"].astype(float),
        d["_ui_expr"].astype(float),
    ])
    y = d["_lcv"].astype(int).values
    if y.sum() < 5:
        return None
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(X, y)
    prop = clf.predict_proba(X)[:, 1]
    q_idx = np.where(y == 1)[0]
    ctl_idx = np.where(y == 0)[0]
    caliper = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), dtype=bool)
    positions = np.arange(len(ctl_idx))
    matched = []
    for i in q_idx[np.argsort(prop[q_idx])]:
        avail = positions[~used]
        if avail.size == 0:
            break
        avail_ctl = ctl_idx[avail]
        diff = np.abs(prop[avail_ctl] - prop[i])
        b = int(np.argmin(diff))
        if diff[b] <= caliper:
            used[avail[b]] = True
            matched.append(ctl_idx[avail[b]])
    if len(matched) < 5:
        return None
    mset = np.array(matched)
    o = d["_crispr_clean"].astype(int).values
    o_q, o_m = int(o[q_idx].sum()), int(o[mset].sum())
    rate_q, rate_m = o_q / len(q_idx), o_m / len(mset)
    _, p_fisher = stats.fisher_exact([[o_q, len(q_idx) - o_q], [o_m, len(mset) - o_m]], alternative="greater")
    return {"fold": rate_q / rate_m if rate_m > 0 else np.nan, "fisher_p": p_fisher,
            "n_query": len(q_idx), "n_matched": len(mset)}


def decimated_test(d):
    """Keep one random gene per module, recompute the matched enrichment."""
    folds, ps, nqs = [], [], []
    idx_by_mod = collections.defaultdict(list)
    for i, m in enumerate(d["_module"].values):
        idx_by_mod[m].append(i)
    mod_keys = list(idx_by_mod.keys())
    for _ in range(B_DECIM):
        keep = [int(RNG.choice(idx_by_mod[m])) for m in mod_keys]
        sub = d.iloc[sorted(keep)]
        r = matched_fold(sub)
        if r is None:
            continue
        folds.append(r["fold"]); ps.append(r["fisher_p"]); nqs.append(r["n_query"])
    return folds, ps, nqs


# ---------------- Run ----------------
g_joint = build_joint_frame()
joint_rows = cr_joint(g_joint)
n_joint = int(len(g_joint))
n_clusters_joint = int(g_joint["_module"].nunique())

d_coupling = universe("_crispr_clean")
full = matched_fold(d_coupling)
folds, ps, nqs = decimated_test(d_coupling)
n_modules = int(atlas["_module"].nunique())

df_cr = pd.DataFrame(joint_rows)
df_cr.to_csv(OUT / "corr_robust_logistic.tsv", sep="\t", index=False)

lcv_cluster = next(r for r in joint_rows if r["term"] == "lcv")
summary = {
    "A_joint_model": "outcome ~ lcv + ppi + response + nstudy + chronos + common",
    "A_n_genes": n_joint,
    "A_n_modules_in_universe": n_clusters_joint,
    "A_lcv_OR": lcv_cluster["OR"],
    "A_lcv_ci_ordinary": [lcv_cluster["ci_lo_ordinary"], lcv_cluster["ci_hi_ordinary"]],
    "A_lcv_p_ordinary": lcv_cluster["p_ordinary"],
    "A_lcv_ci_cluster": [lcv_cluster["ci_lo_cluster"], lcv_cluster["ci_hi_cluster"]],
    "A_lcv_p_cluster": lcv_cluster["p_cluster"],
    "B_n_universe": int(len(d_coupling)),
    "B_n_modules_universe": int(d_coupling["_module"].nunique()),
    "B_full_universe_matched": full,
    "B_decimation_B": B_DECIM,
    "B_decimation_n_valid": len(folds),
    "B_decimation_fold_median": float(np.median(folds)) if folds else None,
    "B_decimation_fold_p05": float(np.percentile(folds, 5)) if folds else None,
    "B_decimation_fold_p95": float(np.percentile(folds, 95)) if folds else None,
    "B_decimation_fisher_p_median": float(np.median(ps)) if ps else None,
    "B_decimation_fisher_p_p95": float(np.percentile(ps, 95)) if ps else None,
    "B_decimation_frac_fisher_p_lt_0.05": float(np.mean(np.array(ps) < 0.05)) if ps else None,
    "B_decimation_n_query_median": float(np.median(nqs)) if nqs else None,
}
(OUT / "corr_robust_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print(json.dumps(summary, indent=2, default=str))
print(df_cr.to_string(index=False))
