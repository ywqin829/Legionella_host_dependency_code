"""Species-stratified LCV -> fitness-filtered CRISPR evidence (post-freeze sensitivity).

Design is copied verbatim from the primary coupling analysis (postfreeze_a1_coupling.py):
same gene universe, same propensity model, same 1:1 nearest-neighbour matching with a
0.25*SD caliper, same 4000 label permutations, same outcome definition
(fitness-filtered CRISPR evidence = functional hit without the baseline-fitness caveat).

The ONLY thing that changes is the exposure: LCV detection is re-derived from the
species-specific LCV evidence columns instead of the pooled `lcv_spatial_supported` flag.

Purpose: test whether the LCV -> CRISPR association requires the shared human U937
cellular background in which both the human LCV proteome and the U937 CRISPR screen
were generated. This is an evidence-layer association test, not candidate validation.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)
N_PERM = 4000
N_BOOT = 2000


def B(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


# ---------------- Load (identical to A1) ----------------
atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz",
                    sep="\t", low_memory=False)
deg = (pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv",
                   sep="\t", low_memory=False)
       [["gene", "cross_pathogen_degree", "cross_pathogen_records", "cross_pathogen_pmids"]]
       .drop_duplicates())

atlas["_ppi_AB"] = B(atlas["ppi_AB_flag"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fitness_ok"] = B(atlas["fitness_artifact_compatible"])
atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
     "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"])
atlas["_crispr_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
# fitness_artifact_compatible=True means hit + baseline-essential (caveat present).
atlas["_crispr_clean"] = atlas["_crispr_hit"] & (~atlas["_fitness_ok"])
atlas["_common_ess"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas["_nstudies"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)

# species-resolved LCV exposure components
atlas["_lcv_human_wt"] = B(atlas["lcv_human_wt"])
atlas["_lcv_human_ankb"] = B(atlas["lcv_human_ankb"])
atlas["_lcv_mouse"] = B(atlas["lcv_mouse_untreated"])
atlas["_lcv_any"] = B(atlas["lcv_spatial_supported"])
atlas["_lcv_human"] = atlas["_lcv_human_wt"] | atlas["_lcv_human_ankb"]


# ---------------- Baseline UI expression percentile (identical to A1) ----------------
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
    return sub.rank(pct=True)[cols].mean(axis=1)


ui_stack = pd.concat([
    load_ui_percentile("GSE61535_log2_expression.tsv.gz",
                       ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    load_ui_percentile("GSE135803_fpkm.tsv.gz",
                       ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    load_ui_percentile("GSE287755_raw_counts.tsv.gz",
                       ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
], axis=1).mean(axis=1, skipna=True)
atlas["_ui_expr"] = atlas["human_symbol"].map(ui_stack)

atlas = atlas.merge(deg, left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)
atlas["_ui_expr"] = atlas["_ui_expr"].fillna(0.0)

# ---------------- Collapse to one row per gene (identical to A1) ----------------
source_atlas_rows = len(atlas)
atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
    "_ppi_AB": "max", "_resp": "max", "_fitness_ok": "max",
    "_crispr_tested": "max", "_crispr_hit": "max", "_crispr_clean": "max",
    "_common_ess": "max", "_chronos": "median", "_nstudies": "max",
    "_ui_expr": "median", "_deg": "max", "_pmids": "max",
    "_lcv_any": "max", "_lcv_human": "max", "_lcv_mouse": "max",
    "_lcv_human_wt": "max", "_lcv_human_ankb": "max",
})
if not atlas["human_symbol"].is_unique:
    raise RuntimeError("Species-stratified universe must contain unique human symbols")

CONF5 = ["_deg", "_pmids", "_chronos", "_common_ess", "_ui_expr"]          # A1 primary
CONF6 = CONF5 + ["_nstudies"]                                              # manuscript audited set
COVARIATE_SETS = {"a1_primary_5": CONF5, "audited_6": CONF6}


def log1p_nz(x: pd.Series) -> pd.Series:
    return np.log1p(x.replace(0, 1e-9)).astype(float)


def matched_test(conf: list[str], query_col: str, outcome_col: str, label: str) -> dict:
    """Confounder-matched enrichment test; logic copied from A1 matched_test()."""
    d = atlas.copy()
    d = d[d["_crispr_tested"]] if outcome_col.startswith("_crispr") else d
    for c in conf:
        d = d[d[c].notna()]
    d = d.reset_index(drop=True)
    d["_q"] = d[query_col].astype(bool)
    d["_o"] = d[outcome_col].astype(bool)
    nq = int(d["_q"].sum())
    if nq < 5:
        return {"label": label, "exposure": query_col, "n_query": nq,
                "n_universe": len(d), "status": "TOO_SMALL"}

    X = np.column_stack([log1p_nz(d[c]) if c in ("_deg", "_pmids") else d[c]
                         for c in conf])
    y = d["_q"].astype(int).values
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(X, y)
    prop = clf.predict_proba(X)[:, 1]

    q_idx = np.where(y == 1)[0]
    ctl_idx = np.where(y == 0)[0]
    caliper = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), dtype=bool)
    positions = np.arange(len(ctl_idx))
    matched_q, matched_c = [], []
    for i in q_idx[np.argsort(prop[q_idx])]:
        avail_pos = positions[~used]
        if avail_pos.size == 0:
            break
        avail_ctl = ctl_idx[avail_pos]
        diff = np.abs(prop[avail_ctl] - prop[i])
        best = int(np.argmin(diff))
        if diff[best] <= caliper:
            used[avail_pos[best]] = True
            matched_q.append(int(i))
            matched_c.append(int(ctl_idx[avail_pos[best]]))
    if len(matched_c) < 5:
        return {"label": label, "exposure": query_col, "n_query": nq,
                "n_universe": len(d), "n_matched": len(matched_c), "status": "MATCH_FAIL"}

    mq = np.array(matched_q)
    mset = np.array(matched_c)
    oq_all = d["_o"].values[q_idx].astype(bool)     # all queries (A1 convention)
    oq = d["_o"].values[mq].astype(bool)            # successfully matched queries
    om = d["_o"].values[mset].astype(bool)
    rate_q = float(oq_all.mean())
    rate_q_matched_queries = float(oq.mean())
    rate_m = float(om.mean())
    rate_uni = float(d["_o"].mean())
    fold = rate_q / rate_m if rate_m > 0 else np.inf

    # bootstrap CI over resampled query set and resampled matched controls
    folds = np.empty(N_BOOT)
    rds = np.empty(N_BOOT)
    for k in range(N_BOOT):
        iq = RNG.integers(0, len(oq_all), len(oq_all))
        im = RNG.integers(0, len(om), len(om))
        rq = float(oq_all[iq].mean())
        rm = float(om[im].mean())
        rds[k] = rq - rm
        folds[k] = rq / rm if rm > 0 else np.nan
    fold_ci = tuple(np.nanpercentile(folds, [2.5, 97.5]))
    rd_ci = tuple(np.percentile(rds, [2.5, 97.5]))

    table = [[int(oq_all.sum()), int(len(oq_all) - oq_all.sum())],
             [int(om.sum()), int(len(om) - om.sum())]]
    _, p_fisher = stats.fisher_exact(table, alternative="greater")

    # label-permutation p within query union matched controls (A1 convention)
    pool_o = np.concatenate([oq_all, om])
    perm_stat = np.empty(N_PERM)
    for k in range(N_PERM):
        shuf = RNG.permutation(pool_o)
        perm_stat[k] = shuf[:len(oq_all)].mean() - shuf[len(oq_all):].mean()
    obs_stat = rate_q - rate_m
    p_perm = float(((perm_stat >= obs_stat).sum() + 1) / (N_PERM + 1))

    # covariate-adjusted logistic on the full universe (identical to A1)
    Xa = np.column_stack([d["_q"].astype(float)] +
                         [log1p_nz(d[c]) if c in ("_deg", "_pmids") else d[c] for c in conf])
    try:
        model = sm.Logit(d["_o"].astype(int).values, sm.add_constant(Xa)).fit(disp=0, maxiter=200)
        or_val = float(np.exp(model.params[1]))
        ci = np.exp(model.conf_int()[1])
        p_adj = float(model.pvalues[1])
        or_ci = (float(ci[0]), float(ci[1]))
    except Exception:
        or_val, p_adj, or_ci = np.nan, np.nan, (np.nan, np.nan)

    return {
        "label": label, "exposure": query_col,
        "n_query": int(len(q_idx)), "n_matched": int(len(mset)), "n_universe": int(len(d)),
        "outcome_rate_query": rate_q, "outcome_rate_matched": rate_m,
        "outcome_rate_query_matched_only": rate_q_matched_queries,
        "outcome_rate_universe": rate_uni,
        "fold_vs_matched": fold, "fold_ci_lo": fold_ci[0], "fold_ci_hi": fold_ci[1],
        "risk_difference": rate_q - rate_m, "rd_ci_lo": rd_ci[0], "rd_ci_hi": rd_ci[1],
        "fisher_p": float(p_fisher), "perm_p": p_perm,
        "adj_logistic_OR": or_val, "adj_logistic_OR_lo": or_ci[0],
        "adj_logistic_OR_hi": or_ci[1], "adj_logistic_p": p_adj,
        "status": "OK",
    }


EXPOSURES = [
    ("_lcv_any", "LCV any species (primary exposure)"),
    ("_lcv_human", "LCV human U937 (WT or ankB)"),
    ("_lcv_mouse", "LCV mouse BMM untreated"),
    ("_lcv_human_wt", "LCV human U937 WT only"),
    ("_lcv_human_ankb", "LCV human U937 ankB only"),
]

rows = []
for set_name, conf in COVARIATE_SETS.items():
    for col, label in EXPOSURES:
        r = matched_test(conf, col, "_crispr_clean", label)
        r["covariate_set"] = set_name
        r["n_covariates"] = len(conf)
        rows.append(r)
        if r.get("status") == "OK":
            print(f"[{set_name}] {label:38s} nq={r['n_query']:5d} matched={r['n_matched']:5d} "
                  f"rate {r['outcome_rate_query']:.4f} vs {r['outcome_rate_matched']:.4f}  "
                  f"fold={r['fold_vs_matched']:.2f} [{r['fold_ci_lo']:.2f},{r['fold_ci_hi']:.2f}]  "
                  f"perm_p={r['perm_p']:.4f}  adjOR={r['adj_logistic_OR']:.2f} "
                  f"[{r['adj_logistic_OR_lo']:.2f},{r['adj_logistic_OR_hi']:.2f}] p={r['adj_logistic_p']:.4g}")
        else:
            print(f"[{set_name}] {label:38s} {r['status']}")

res = pd.DataFrame(rows)
col_order = ["covariate_set", "n_covariates", "label", "exposure", "n_query", "n_matched",
             "n_universe", "outcome_rate_query", "outcome_rate_query_matched_only",
             "outcome_rate_matched", "outcome_rate_universe",
             "fold_vs_matched", "fold_ci_lo", "fold_ci_hi", "risk_difference",
             "rd_ci_lo", "rd_ci_hi", "fisher_p", "perm_p",
             "adj_logistic_OR", "adj_logistic_OR_lo", "adj_logistic_OR_hi",
             "adj_logistic_p", "status"]
res = res[[c for c in col_order if c in res.columns]]
res.to_csv(OUT / "h1_species_stratified_lcv_crispr.tsv", sep="\t", index=False)

# ---------------- Joint model: human + mouse LCV simultaneously ----------------
joint = {}
for set_name, conf in COVARIATE_SETS.items():
    d = atlas[atlas["_crispr_tested"]].copy()
    for c in conf:
        d = d[d[c].notna()]
    X = np.column_stack([
        d["_lcv_human"].astype(float), d["_lcv_mouse"].astype(float),
    ])
    extra = np.column_stack([log1p_nz(d[c]) if c in ("_deg", "_pmids") else d[c] for c in conf])
    X = np.column_stack([X, extra])
    y = d["_crispr_clean"].astype(int).values
    try:
        m = sm.Logit(y, sm.add_constant(X)).fit(disp=0, maxiter=200)
        ci = np.exp(m.conf_int())
        # contingency between the two exposures (shared-background check)
        ct = pd.crosstab(d["_lcv_human"], d["_lcv_mouse"])
        joint[set_name] = {
            "n_universe": int(len(d)),
            "OR_lcv_human": float(np.exp(m.params[1])),
            "OR_lcv_human_lo": float(ci[1][0]), "OR_lcv_human_hi": float(ci[1][1]),
            "p_lcv_human": float(m.pvalues[1]),
            "OR_lcv_mouse": float(np.exp(m.params[2])),
            "OR_lcv_mouse_lo": float(ci[2][0]), "OR_lcv_mouse_hi": float(ci[2][1]),
            "p_lcv_mouse": float(m.pvalues[2]),
            "human_only_n": int(((d["_lcv_human"]) & (~d["_lcv_mouse"])).sum()),
            "mouse_only_n": int(((~d["_lcv_human"]) & (d["_lcv_mouse"])).sum()),
            "both_n": int((d["_lcv_human"] & d["_lcv_mouse"]).sum()),
            "neither_n": int(((~d["_lcv_human"]) & (~d["_lcv_mouse"])).sum()),
            "model": "logit(fitness-filtered CRISPR evidence ~ LCV_human + LCV_mouse + covariates)",
            "covariates": conf,
        }
    except Exception as exc:
        joint[set_name] = {"status": "MODEL_FAIL", "error": str(exc)}

jrows = []
for set_name, j in joint.items():
    if "OR_lcv_human" in j:
        jrows.append({"covariate_set": set_name, **j})
pd.DataFrame(jrows).to_csv(OUT / "h1_joint_human_mouse_lcv_model.tsv", sep="\t", index=False)

summary = {
    "analysis": "Species-stratified LCV -> fitness-filtered CRISPR evidence sensitivity",
    "design": ("Universe, propensity model, 1:1 nearest-neighbour matching (0.25*SD caliper), "
               "4000 permutations and outcome definition copied verbatim from the primary "
               "coupling analysis (postfreeze_a1_coupling.py); only the exposure changes."),
    "covariate_sets": {k: v for k, v in COVARIATE_SETS.items()},
    "n_permutations": N_PERM, "n_bootstrap": N_BOOT,
    "source_atlas_rows": int(source_atlas_rows),
    "unique_genes": int(len(atlas)),
    "joint_model": joint,
    "interpretation": ("Evidence-layer association test. It asks whether the LCV-CRISPR coupling "
                       "requires the shared human U937 cellular background, not whether any "
                       "individual candidate is validated in mouse. Species-restricted results "
                       "should not be read as candidate-level validation."),
}
(OUT / "h1_species_stratified_summary.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print("\nwrote", OUT / "h1_species_stratified_lcv_crispr.tsv")
print(json.dumps(joint, indent=2, ensure_ascii=False))
