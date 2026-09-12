from __future__ import annotations

"""Post-freeze sensitivity analyses for submission reinforcement.

This script is deliberately downstream of the frozen atlas and holdout outputs.
It does not change candidate selection, Tier assignment, holdout rules, or any
frozen primary result.  It adds:

1. per-covariate balance diagnostics and overlap-weighting sensitivity for the
   PPI(A/B) -> fitness-filtered CRISPR evidence coupling;
2. exact binomial intervals and a sparse PPI x LCV interaction audit;
3. size-matched leave-RAB1A/RAB10-out diagnostics for the 54-gene holdout panel.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.special import expit
from sklearn.linear_model import LogisticRegression


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(42)
N_PERM = 10_000


def B(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


def load_ui_percentile(matrix_name: str, ui_gsms: list[str], map_name: str, use_probe_map: bool):
    mat = pd.read_csv(ROOT / "05_expression" / "matrices" / matrix_name,
                      sep="\t", compression="gzip")
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
    merged = mat.merge(gm, on=id_col, how="inner")
    sub = merged[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    return sub.rank(pct=True)[cols].mean(axis=1)


def load_atlas() -> pd.DataFrame:
    atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz",
                        sep="\t", low_memory=False)
    deg = pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv",
                      sep="\t", low_memory=False)
    deg = deg[["gene", "cross_pathogen_degree", "cross_pathogen_pmids"]].drop_duplicates()

    atlas["_ppi"] = B(atlas["ppi_AB_flag"])
    atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
    atlas["_resp"] = B(atlas["stable_cross_cohort"])
    atlas["_fitness"] = B(atlas["fitness_artifact_compatible"])
    atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
    atlas["_clean"] = atlas["_hit"] & (~atlas["_fitness"])
    atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin([
        "F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
        "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT",
    ])
    atlas["_common"] = B(atlas["common_essential_flag"])
    atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
    atlas["_nstud"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)

    ui = pd.concat([
        load_ui_percentile("GSE61535_log2_expression.tsv.gz",
                           ["GSM1507677", "GSM1507678", "GSM1507679"],
                           "GSE61535_gene_map.tsv", True),
        load_ui_percentile("GSE135803_fpkm.tsv.gz",
                           ["GSM4030445", "GSM4030446", "GSM4030447"],
                           "GSE135803_gene_map.tsv", False),
        load_ui_percentile("GSE287755_raw_counts.tsv.gz",
                           ["GSM8750874", "GSM8750875", "GSM8750876"],
                           "GSE287755_gene_map.tsv", False),
    ], axis=1).mean(axis=1, skipna=True)
    atlas["_ui"] = atlas["human_symbol"].map(ui).fillna(0.0)
    atlas = atlas.merge(deg, left_on="human_symbol", right_on="gene", how="left")
    atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0.0)
    atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0.0)
    atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
        "_ppi": "max", "_lcv": "max", "_resp": "max", "_fitness": "max",
        "_hit": "max", "_clean": "max", "_crispr_tested": "max",
        "_common": "max", "_chronos": "median", "_nstud": "max",
        "_ui": "median", "_deg": "max", "_pmids": "max",
    })
    if not atlas["human_symbol"].is_unique:
        raise RuntimeError("Sensitivity-analysis universe must contain unique human symbols")
    return atlas


MATCH_CONF = ["_deg", "_pmids", "_chronos", "_common", "_ui", "_nstud"]
JOINT_CONF = ["_deg", "_pmids", "_chronos", "_common", "_ui"]


def log1pz(x: pd.Series) -> pd.Series:
    return np.log1p(x.replace(0, 1e-9)).astype(float)


def smd(frame: pd.DataFrame, a: np.ndarray, b: np.ndarray, col: str) -> float:
    va = frame[col].iloc[a]
    vb = frame[col].iloc[b]
    pooled = np.sqrt((va.std() ** 2 + vb.std() ** 2) / 2)
    return float((va.mean() - vb.mean()) / pooled) if pooled > 0 else 0.0


def headline_matching(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Reproduce the PPI(A/B)-to-filtered-CRISPR matching definition."""
    d = frame.loc[frame["_crispr_tested"]].copy()
    d = d.dropna(subset=MATCH_CONF + ["_clean"]).reset_index(drop=True)
    d["_q"] = d["_ppi"].astype(bool)
    d["_o"] = d["_clean"].astype(bool)
    X = np.column_stack([
        log1pz(d["_deg"]), log1pz(d["_pmids"]), d["_chronos"],
        d["_common"].astype(float), d["_ui"], d["_nstud"],
    ])
    y = d["_q"].astype(int).to_numpy()
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(X, y)
    prop = clf.predict_proba(X)[:, 1]
    q_idx = np.where(y == 1)[0]
    ctl_idx = np.where(y == 0)[0]
    caliper = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), dtype=bool)
    matched_q: list[int] = []
    matched_c: list[int] = []
    for i in q_idx[np.argsort(prop[q_idx])]:
        available = np.where(~used)[0]
        if available.size == 0:
            break
        controls = ctl_idx[available]
        diff = np.abs(prop[controls] - prop[i])
        best = int(np.argmin(diff))
        if diff[best] <= caliper:
            used[available[best]] = True
            matched_q.append(int(i))
            matched_c.append(int(controls[best]))
    mq = np.asarray(matched_q, dtype=int)
    mc = np.asarray(matched_c, dtype=int)
    balance_rows = []
    for c in MATCH_CONF:
        balance_rows.append({
            "covariate": c,
            "smd_before": smd(d, q_idx, ctl_idx, c),
            "abs_smd_before": abs(smd(d, q_idx, ctl_idx, c)),
            "smd_after": smd(d, mq, mc, c),
            "abs_smd_after": abs(smd(d, mq, mc, c)),
        })
    balance = pd.DataFrame(balance_rows)
    q_rate = float(d.loc[q_idx, "_o"].mean())
    m_rate = float(d.loc[mc, "_o"].mean())
    matched_summary = {
        "method": "propensity 1:1 nearest-neighbour, existing caliper",
        "n_universe": int(len(d)), "n_query": int(len(q_idx)),
        "n_matched": int(len(mc)), "query_rate": q_rate,
        "control_rate": m_rate, "fold_vs_matched": q_rate / m_rate if m_rate else np.inf,
        "mean_abs_smd_before": float(balance.abs_smd_before.mean()),
        "mean_abs_smd_after": float(balance.abs_smd_after.mean()),
    }

    # Overlap weighting sensitivity.  This is a marginal, alternative estimator;
    # it is reported as sensitivity and does not replace the frozen analysis.
    w = np.where(y == 1, 1.0 - prop, prop)
    treated = y == 1
    control = ~treated
    wt_q = float(np.average(d.loc[treated, "_o"], weights=w[treated]))
    wt_c = float(np.average(d.loc[control, "_o"], weights=w[control]))
    weighted_or = (wt_q / max(1e-12, 1 - wt_q)) / (wt_c / max(1e-12, 1 - wt_c))
    try:
        glm = sm.GLM(d["_o"].astype(int), sm.add_constant(d["_q"].astype(int)),
                     family=sm.families.Binomial(), freq_weights=w).fit(cov_type="HC1")
        glm_or = float(np.exp(glm.params.iloc[1]))
        glm_ci = tuple(np.exp(glm.conf_int().iloc[1]).tolist())
        glm_p = float(glm.pvalues.iloc[1])
    except Exception as exc:  # pragma: no cover - defensive for platform versions
        glm_or, glm_ci, glm_p = np.nan, (np.nan, np.nan), np.nan
        weighted_error = str(exc)
    else:
        weighted_error = ""
    ess_q = float(w[treated].sum() ** 2 / np.square(w[treated]).sum())
    ess_c = float(w[control].sum() ** 2 / np.square(w[control]).sum())
    sensitivity = pd.DataFrame([
        {**matched_summary, "method": "post-freeze propensity-match reconstruction", "estimate": matched_summary["fold_vs_matched"],
         "estimate_scale": "rate ratio", "ci_low": np.nan, "ci_high": np.nan, "p": np.nan},
        {"method": "overlap weighting", "n_universe": int(len(d)), "n_query": int(len(q_idx)),
         "n_matched": np.nan, "query_rate": wt_q, "control_rate": wt_c,
         "fold_vs_matched": wt_q / wt_c if wt_c else np.inf, "mean_abs_smd_before": np.nan,
         "mean_abs_smd_after": np.nan, "estimate": glm_or, "estimate_scale": "marginal odds ratio",
         "ci_low": glm_ci[0], "ci_high": glm_ci[1], "p": glm_p,
         "effective_n_query": ess_q, "effective_n_control": ess_c,
         "weighted_risk_difference": wt_q - wt_c, "error": weighted_error},
    ])
    return balance, sensitivity, {"matched": matched_summary, "propensity_sd": float(np.std(prop))}


def cp_interval(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if k == 0:
        lo = 0.0
    else:
        lo = float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    if k == n:
        hi = 1.0
    else:
        hi = float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def joint_audit(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = frame.loc[frame["_crispr_tested"]].dropna(subset=JOINT_CONF + ["_clean"]).copy()
    d["_o"] = d["_clean"].astype(bool)
    categories = [
        ("Universe", np.ones(len(d), dtype=bool)),
        ("PPI(A/B)", d["_ppi"].to_numpy(bool)),
        ("LCV", d["_lcv"].to_numpy(bool)),
        ("PPI(A/B)+LCV", (d["_ppi"] & d["_lcv"]).to_numpy(bool)),
    ]
    rows = []
    for label, mask in categories:
        k = int(d.loc[mask, "_o"].sum())
        n = int(mask.sum())
        lo, hi = cp_interval(k, n)
        rows.append({"group": label, "successes": k, "n": n, "rate": k / n,
                     "ci_low": lo, "ci_high": hi})
    rates = pd.DataFrame(rows)

    # Sparse interaction audit with the same confounder set as the existing A1 model.
    X = pd.DataFrame({
        "ppi": d["_ppi"].astype(float), "lcv": d["_lcv"].astype(float),
        "ppi_x_lcv": (d["_ppi"] & d["_lcv"]).astype(float),
        "log_deg": log1pz(d["_deg"]), "log_pmids": log1pz(d["_pmids"]),
        "chronos": d["_chronos"], "common_essential": d["_common"].astype(float),
        "ui_expression": d["_ui"],
    })
    X = sm.add_constant(X)
    y = d["_o"].astype(int)
    rows = []
    try:
        fit = sm.Logit(y, X).fit(disp=0, maxiter=200)
        ci = fit.conf_int()
        for term in ["ppi", "lcv", "ppi_x_lcv"]:
            rows.append({"term": term, "OR": float(np.exp(fit.params[term])),
                         "ci_low": float(np.exp(ci.loc[term, 0])),
                         "ci_high": float(np.exp(ci.loc[term, 1])),
                         "p": float(fit.pvalues[term]), "status": "OK"})
        status = "OK"
        error = ""
    except Exception as exc:
        status = "MODEL_FAIL"
        error = str(exc)
        rows.append({"term": "ppi_x_lcv", "OR": np.nan, "ci_low": np.nan,
                     "ci_high": np.nan, "p": np.nan, "status": status,
                     "error": error})
    interaction = pd.DataFrame(rows)
    return rates, interaction


def panel_leaveout() -> pd.DataFrame:
    result = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz",
                         sep="\t")
    panel = pd.read_csv(ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv",
                        sep="\t")
    d = result.loc[result["contrast_id"].eq("WT_vs_UI")].copy()
    p = np.clip(d["PValue"].to_numpy(float), 1e-300, 1 - 1e-15)
    d["z_signed"] = np.sign(d["logFC"].to_numpy(float)) * stats.norm.isf(p / 2)
    merged = d.merge(panel, on="human_symbol", how="inner")
    measured = merged.loc[merged["expected_wt_vs_uninfected"].isin(["UP", "DOWN"])].copy()
    measured["orientation"] = np.where(measured["expected_wt_vs_uninfected"].eq("UP"), 1.0, -1.0)
    measured["z_oriented"] = measured["orientation"] * measured["z_signed"]
    universe = d["z_signed"].to_numpy(float)
    scenarios = {
        "full_54_gene_panel": set(),
        "minus_RAB1A": {"RAB1A"},
        "minus_RAB10": {"RAB10"},
        "minus_RAB1A_RAB10": {"RAB1A", "RAB10"},
    }
    rows = []
    for scenario, removed in scenarios.items():
        sub = measured.loc[~measured["human_symbol"].isin(removed)]
        values = sub["z_oriented"].to_numpy(float)
        null = np.empty(N_PERM, dtype=float)
        for i in range(N_PERM):
            null[i] = float(np.mean(RNG.choice(universe, size=len(values), replace=False)))
        observed = float(np.mean(values))
        rows.append({"scenario": scenario, "removed_genes": ";".join(sorted(removed)) or "NONE",
                     "panel_n": int(len(values)), "observed_mean_oriented_z": observed,
                     "null_mean": float(np.mean(null)), "null_sd": float(np.std(null, ddof=1)),
                     "null_p95": float(np.quantile(null, 0.95)),
                     "empirical_p": float((1 + np.sum(null >= observed)) / (N_PERM + 1)),
                     "n_permutations": N_PERM,
                     "interpretation": "post-freeze descriptive influence sensitivity"})
    return pd.DataFrame(rows)


def main() -> int:
    atlas = load_atlas()
    balance, matching, matching_meta = headline_matching(atlas)
    balance.to_csv(OUT / "g_matching_balance_by_covariate.tsv", sep="\t", index=False)
    matching.to_csv(OUT / "g_matching_sensitivity.tsv", sep="\t", index=False)
    rates, interaction = joint_audit(atlas)
    rates.to_csv(OUT / "g_ppi_lcv_exact_rate_ci.tsv", sep="\t", index=False)
    interaction.to_csv(OUT / "g_ppi_lcv_interaction.tsv", sep="\t", index=False)
    leaveout = panel_leaveout()
    leaveout.to_csv(OUT / "g_panel_leaveout.tsv", sep="\t", index=False)
    summary = {
        "analysis": "postfreeze_sensitivity_reinforcement",
        "frozen_inputs_unchanged": True,
        "headline_matching": matching_meta,
        "panel_leaveout": leaveout.to_dict("records"),
        "claim_boundary": "Sensitivity analyses are descriptive and do not alter frozen candidates, tiers, holdout rules, or primary endpoint decisions.",
    }
    (OUT / "g_sensitivity_reinforcement_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
