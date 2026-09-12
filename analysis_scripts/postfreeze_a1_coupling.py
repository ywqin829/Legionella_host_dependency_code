from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)


def B(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


# ---------------- Load --------------------------------
atlas = pd.read_csv(
    ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz",
    sep="\t", low_memory=False,
)
deg = (
    pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv", sep="\t", low_memory=False)
    [["gene", "cross_pathogen_degree", "cross_pathogen_records", "cross_pathogen_pmids"]]
    .drop_duplicates()
)

# evidence tracks (gene level)
atlas["_ppi_AB"] = B(atlas["ppi_AB_flag"])
atlas["_ppi_any"] = pd.to_numeric(atlas["ppi_interaction_records"], errors="coerce").fillna(0) > 0
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fitness_ok"] = B(atlas["fitness_artifact_compatible"])
atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
     "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"]
)
atlas["_crispr_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
# NOTE: fitness_artifact_compatible=True means hit+baseline-essential (caveat present);
# Fitness-filtered CRISPR evidence = hit WITHOUT that caveat.
atlas["_crispr_clean"] = atlas["_crispr_hit"] & (~atlas["_fitness_ok"])
atlas["_common_ess"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas["_nstudies"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)

# ---------------- Baseline UI expression percentile (human) ----------
def load_ui_percentile(ds: str, matrix_name: str, ui_gsms: list[str], map_name: str, use_probe_map: bool):
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
    # rank within this dataset to a percentile (scale-free across platforms)
    r = sub.rank(pct=True)  # 0..1
    return r[cols].mean(axis=1)  # mean over UI samples

ui_percentiles = {
    "GSE61535": load_ui_percentile("GSE61535", "GSE61535_log2_expression.tsv.gz",
                                   ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    "GSE135803": load_ui_percentile("GSE135803", "GSE135803_fpkm.tsv.gz",
                                    ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    "GSE287755": load_ui_percentile("GSE287755", "GSE287755_raw_counts.tsv.gz",
                                    ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
}
ui_stack = pd.concat(ui_percentiles.values(), axis=1, keys=list(ui_percentiles.keys()))
atlas["_ui_expr"] = atlas["human_symbol"].map(ui_stack.mean(axis=1, skipna=True))
atlas["_ui_expr_measured"] = atlas["human_symbol"].isin(ui_stack.index)

# ---------------- Confounders / universe ----------------
atlas = atlas.merge(deg, left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)
atlas["_ui_expr"] = atlas["_ui_expr"].fillna(0.0)

# Expression provenance can leave more than one atlas row for a harmonized
# human symbol. All coupling analyses are gene-level, so collapse before
# defining analysis universes or fitting propensity models.
source_atlas_rows = len(atlas)
atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
    "_ppi_AB": "max", "_ppi_any": "max", "_lcv": "max", "_resp": "max",
    "_fitness_ok": "max", "_crispr_tested": "max", "_crispr_hit": "max",
    "_crispr_clean": "max", "_common_ess": "max", "_chronos": "median",
    "_nstudies": "max", "_ui_expr": "median", "_ui_expr_measured": "max",
    "_deg": "max", "_pmids": "max",
})
if not atlas["human_symbol"].is_unique:
    raise RuntimeError("Coupling universe must contain unique human symbols")

CONF = ["_deg", "_pmids", "_chronos", "_common_ess", "_ui_expr"]


def log1p_nz(x: pd.Series) -> pd.Series:
    return np.log1p(x.replace(0, 1e-9)).astype(float)


# ---------------- Matched enrichment test ----------------
def matched_test(U: pd.DataFrame, query_col: str, outcome_col: str, label: str) -> dict:
    """U: universe restricted to genes with outcome + confounders present.
    query_col/outcome_col are boolean columns. Returns confounder-matched enrichment."""
    d = U.copy()
    for c in CONF:
        d = d[d[c].notna()]
    d = d.reset_index(drop=True)
    d["_q"] = d[query_col].astype(bool)
    d["_o"] = d[outcome_col].astype(bool)
    nq = int(d["_q"].sum())
    if nq < 5:
        return {"label": label, "n_query": nq, "n_universe": len(d), "status": "TOO_SMALL"}

    X = np.column_stack([
        log1p_nz(d["_deg"]), log1p_nz(d["_pmids"]), d["_chronos"],
        d["_common_ess"].astype(float), d["_ui_expr"],
    ])
    y = d["_q"].astype(int).values
    try:
        clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
        clf.fit(X, y)
        prop = clf.predict_proba(X)[:, 1]
    except Exception as exc:  # pragma: no cover
        return {"label": label, "n_query": nq, "n_universe": len(d), "status": "MODEL_FAIL", "error": str(exc)}

    q_idx = np.where(y == 1)[0]
    ctl_idx = np.where(y == 0)[0]
    caliper = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), dtype=bool)
    positions = np.arange(len(ctl_idx))
    matched = []
    order = q_idx[np.argsort(prop[q_idx])]
    for i in order:
        avail_pos = positions[~used]
        if avail_pos.size == 0:
            break
        avail_ctl = ctl_idx[avail_pos]
        diff = np.abs(prop[avail_ctl] - prop[i])
        best = int(np.argmin(diff))
        if diff[best] <= caliper:
            used[avail_pos[best]] = True
            matched.append(ctl_idx[avail_pos[best]])
    if len(matched) < 5:
        return {"label": label, "n_query": nq, "n_universe": len(d), "n_matched": len(matched), "status": "MATCH_FAIL"}

    mset = np.array(matched)
    o_q = int(d["_o"].iloc[q_idx].sum())
    o_m = int(d["_o"].iloc[mset].sum())
    rate_q = o_q / len(q_idx)
    rate_m = o_m / len(mset)
    rate_uni = float(d["_o"].mean())
    table = [[o_q, len(q_idx) - o_q], [o_m, len(mset) - o_m]]
    odds, p_fisher = stats.fisher_exact(table, alternative="greater")

    # label-permutation p within Q ∪ matched controls
    pool_o = np.concatenate([d["_o"].iloc[q_idx].values, d["_o"].iloc[mset].values])
    n_iter = 4000
    perm_stat = np.empty(n_iter)
    for k in range(n_iter):
        shuf = RNG.permutation(pool_o)
        perm_stat[k] = shuf[:len(q_idx)].mean() - shuf[len(q_idx):].mean()
    obs_stat = rate_q - rate_m
    p_perm = float(((perm_stat >= obs_stat).sum() + 1) / (n_iter + 1))

    return {
        "label": label,
        "query": query_col, "outcome": outcome_col,
        "n_query": int(len(q_idx)), "n_matched": int(len(mset)), "n_universe": len(d),
        "outcome_rate_query": float(rate_q), "outcome_rate_matched": float(rate_m),
        "outcome_rate_universe": float(rate_uni),
        "fold_vs_matched": float(rate_q / rate_m) if rate_m > 0 else np.inf,
        "fold_vs_universe": float(rate_q / rate_uni) if rate_uni > 0 else np.inf,
        "fisher_p": float(p_fisher), "perm_p": float(p_perm),
        "status": "OK",
    }


# ---------------- Run headline couplings ----------------
def restrict(U: pd.DataFrame, outcome_col: str) -> pd.DataFrame:
    d = U.copy()
    if outcome_col.startswith("_crispr"):
        d = d[d["_crispr_tested"]]          # CRISPR-measurable universe
    elif outcome_col == "_resp":
        d = d[d["_nstudies"] >= 1]          # expression-measurable universe
    for c in CONF:
        d = d[d[c].notna()]
    return d


pairs = [
    ("_ppi_AB", "_crispr_hit", "PPI(A/B) -> CRISPR hit"),
    ("_ppi_AB", "_crispr_clean", "PPI(A/B) -> filtered CRISPR evidence"),
    ("_ppi_any", "_crispr_clean", "PPI(any) -> filtered CRISPR evidence"),
    ("_lcv", "_crispr_hit", "LCV -> CRISPR hit"),
    ("_lcv", "_crispr_clean", "LCV -> filtered CRISPR evidence"),
    ("_ppi_AB", "_lcv", "PPI(A/B) -> LCV"),
    ("_crispr_hit", "_resp", "CRISPR hit -> response"),
    ("_crispr_clean", "_resp", "filtered CRISPR evidence -> response"),
    ("_resp", "_crispr_clean", "response -> filtered CRISPR evidence"),
    ("_lcv", "_resp", "LCV -> response"),
    ("_ppi_AB", "_resp", "PPI(A/B) -> response"),
]
results = []
for q, o, label in pairs:
    U = restrict(atlas, o)
    results.append(matched_test(U, q, o, label))

# ---------------- Joint PPI+LCV increment ----------------
joint_rows = []
d = restrict(atlas, "_crispr_clean")
d = d[d["_crispr_tested"] & d["_crispr_clean"].notna()]
for c in CONF:
    d = d[d[c].notna()]
n = len(d)
ppi = d["_ppi_AB"] & d["_lcv"]           # PPI AND LCV query
j = d[ppi]
p_ppi = float(d.loc[d["_ppi_AB"], "_crispr_clean"].mean()) if d["_ppi_AB"].sum() else np.nan
p_lcv = float(d.loc[d["_lcv"], "_crispr_clean"].mean()) if d["_lcv"].sum() else np.nan
p_joint = float(d.loc[ppi, "_crispr_clean"].mean()) if ppi.sum() else np.nan
p_uni = float(d["_crispr_clean"].mean())
joint_rows.append({
    "n_universe": n,
    "n_ppi": int(d["_ppi_AB"].sum()), "n_lcv": int(d["_lcv"].sum()), "n_ppi_and_lcv": int(ppi.sum()),
    "P(fitness-filtered CRISPR evidence|PPI)": p_ppi,
    "P(fitness-filtered CRISPR evidence|LCV)": p_lcv,
    "P(fitness-filtered CRISPR evidence|PPI&LCV)": p_joint,
    "P(fitness-filtered CRISPR evidence|universe)": p_uni,
    "fitness_filtered_hits_in_ppi_and_lcv": int(d.loc[ppi, "_crispr_clean"].sum()),
})

# ---------------- Outputs ----------------
df = pd.DataFrame(results)
df.to_csv(OUT / "a1_evidence_coupling_matched.tsv", sep="\t", index=False)
pd.DataFrame(joint_rows).to_csv(OUT / "a1_joint_ppi_lcv_crispr.tsv", sep="\t", index=False)

print("=== Matched coupling results ===")
print(df.to_string(index=False))
print("\n=== Joint PPI & LCV -> filtered CRISPR evidence ===")
print(pd.DataFrame(joint_rows).to_string(index=False))

# count sanity
sanity = {
    "n_atlas": len(atlas),
    "source_atlas_rows": source_atlas_rows,
    "duplicate_atlas_rows_collapsed": source_atlas_rows - len(atlas),
    "ppi_AB": int(atlas["_ppi_AB"].sum()),
    "ppi_any": int(atlas["_ppi_any"].sum()),
    "lcv": int(atlas["_lcv"].sum()),
    "crispr_hit": int(atlas["_crispr_hit"].sum()),
    "crispr_clean": int(atlas["_crispr_clean"].sum()),
    "response": int(atlas["_resp"].sum()),
    "n_ui_expr_measured": int(atlas["_ui_expr_measured"].sum()),
}
print("\n=== Track sanity ===")
print(json.dumps(sanity, indent=2))
(OUT / "a1_track_sanity.json").write_text(json.dumps(sanity, indent=2), encoding="utf-8")
