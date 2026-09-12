from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.contingency_tables import mcnemar
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
deg = pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv", sep="\t", low_memory=False)

def B(s): return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])
atlas["_ppi"] = B(atlas["ppi_AB_flag"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fitness"] = B(atlas["fitness_artifact_compatible"])
atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_clean"] = atlas["_hit"] & (~atlas["_fitness"])
atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT", "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"])
atlas["_common"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas["_nstud"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)
atlas["_deg"] = pd.to_numeric(atlas.merge(deg[["gene", "cross_pathogen_degree", "cross_pathogen_pmids"]], left_on="human_symbol", right_on="gene", how="left")["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas.merge(deg[["gene", "cross_pathogen_pmids"]], left_on="human_symbol", right_on="gene", how="left")["cross_pathogen_pmids"], errors="coerce").fillna(0)

# baseline UI expression percentile (reuse A1 method)
def load_ui(mat, gsm, gmap, probe):
    m = pd.read_csv(ROOT / "05_expression" / "matrices" / mat, sep="\t", compression="gzip")
    idc = "feature_id" if "feature_id" in m.columns else "Geneid"
    cols = [c for c in gsm if c in m.columns]
    gm = pd.read_csv(ROOT / "05_expression" / "gene_maps" / gmap, sep="\t")
    gm = gm[gm.gene_symbol.notna()]
    if probe and "ambiguous_mapping" in gm.columns:
        gm = gm[~gm.ambiguous_mapping.astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in gm.columns and idc != "feature_id":
        gm = gm.rename(columns={"feature_id": idc})
    mg = m.merge(gm, on=idc, how="inner")
    sub = mg[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    return sub.rank(pct=True)[cols].mean(axis=1)
ui = pd.concat([
    load_ui("GSE61535_log2_expression.tsv.gz", ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    load_ui("GSE135803_fpkm.tsv.gz", ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    load_ui("GSE287755_raw_counts.tsv.gz", ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
], axis=1).mean(axis=1, skipna=True)
atlas["_ui"] = atlas.human_symbol.map(ui).fillna(0)

# Collapse harmonized aliases/provenance rows to one row per human gene before
# gene-level matching and regression.
source_atlas_rows = len(atlas)
atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
    "_ppi": "max", "_lcv": "max", "_resp": "max", "_fitness": "max",
    "_hit": "max", "_clean": "max", "_crispr_tested": "max",
    "_common": "max", "_chronos": "median", "_nstud": "max",
    "_deg": "max", "_pmids": "max", "_ui": "median",
})
if not atlas["human_symbol"].is_unique:
    raise RuntimeError("Matching audit universe must contain unique human symbols")

CONF = ["_deg", "_pmids", "_chronos", "_common", "_ui", "_nstud"]
def log1pz(x): return np.log1p(x.replace(0, 1e-9))

pairs = [
    ("_ppi", "_clean", "PPI(A/B) -> filtered CRISPR evidence"),
    ("_lcv", "_clean", "LCV -> filtered CRISPR evidence"),
    ("_ppi", "_lcv", "PPI(A/B) -> LCV"),
    ("_clean", "_resp", "filtered CRISPR evidence -> response"),
    ("_lcv", "_resp", "LCV -> response"),
    ("_ppi", "_hit", "PPI(A/B) -> CRISPR hit"),
    ("_lcv", "_hit", "LCV -> CRISPR hit"),
]

results = []
for qcol, ocol, label in pairs:
    d = atlas[atlas[ocol].notna()].copy()
    if ocol.startswith("_crispr") or ocol == "_clean" or ocol == "_hit":
        d = d[d["_crispr_tested"]]
    if ocol == "_resp":
        d = d[d["_nstud"] >= 1]
    d = d.dropna(subset=CONF + [ocol]).copy()
    d["_q"] = d[qcol].astype(bool); d["_o"] = d[ocol].astype(bool)
    nq = int(d._q.sum())
    if nq < 10: continue

    X = np.column_stack([log1pz(d["_deg"]), log1pz(d["_pmids"]), d["_chronos"], d["_common"].astype(float), d["_ui"], d["_nstud"]])
    y = d._q.astype(int).values

    # --- propensity + 1:1 nearest-neighbour matching (with caliper) ---
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(X, y)
    prop = clf.predict_proba(X)[:, 1]
    q_idx = np.where(y == 1)[0]; ctl_idx = np.where(y == 0)[0]
    cal = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), bool); positions = np.arange(len(ctl_idx)); matched_q = []; matched_c = []
    for i in q_idx[np.argsort(prop[q_idx])]:
        avail_pos = positions[~used]
        if avail_pos.size == 0: break
        ac = ctl_idx[avail_pos]; diff = np.abs(prop[ac] - prop[i]); best = int(np.argmin(diff))
        if diff[best] <= cal:
            used[avail_pos[best]] = True; matched_q.append(i); matched_c.append(ctl_idx[avail_pos[best]])
    mset = np.array(matched_c); mq = np.array(matched_q)

    # --- balance SMD before/after ---
    def smd(a, b, col):
        va = d[col].iloc[a]; vb = d[col].iloc[b]
        pooled = np.sqrt((va.std()**2 + vb.std()**2) / 2)
        return (va.mean() - vb.mean()) / pooled if pooled > 0 else 0.0
    smd_before = np.mean([abs(smd(q_idx, ctl_idx, c)) for c in CONF])
    smd_after = np.mean([abs(smd(q_idx, mset, c)) for c in CONF])

    # --- McNemar exact on matched pairs (only successfully-matched query genes) ---
    oq = d["_o"].iloc[mq].values; om = d["_o"].iloc[mset].values
    b = int((oq & ~om).sum()); c = int((~oq & om).sum())  # discordant
    mcn_p = mcnemar([[int((oq & om).sum()), b], [c, int((~oq & ~om).sum())]], exact=True).pvalue if (b + c) > 0 else np.nan

    # --- covariate-adjusted logistic regression (full universe) ---
    Xa = np.column_stack([d["_q"].astype(float)] + [log1pz(d[c]) if c in ("_deg", "_pmids") else d[c].astype(float) for c in CONF])
    try:
        model = sm.Logit(d["_o"].astype(int).values, sm.add_constant(Xa)).fit(disp=0, maxiter=100)
        or_val = float(np.exp(model.params[1])); p_adj = float(model.pvalues[1])
    except Exception:
        or_val, p_adj = np.nan, np.nan

    results.append({"coupling": label, "n_query": nq, "n_matched": len(mset),
                    "smd_before": round(smd_before, 3), "smd_after": round(smd_after, 3),
                    "mcnemar_p": round(mcn_p, 5) if mcn_p == mcn_p else "",
                    "adj_logistic_OR": round(or_val, 3) if or_val == or_val else "",
                    "adj_logistic_p": round(p_adj, 5) if p_adj == p_adj else "",
                    "discordant_pairs": f"{b}/{c}"})
    print(f"{label:30s} n_q={nq:4d} SMD {smd_before:.3f}->{smd_after:.3f}  McNemar={mcn_p:.4f}  adjOR={or_val:.2f} (p={p_adj:.4f})")

res = pd.DataFrame(results)
# BH across couplings for the two sensitivity tests
if "mcnemar_p" in res:
    mc = pd.to_numeric(res.mcnemar_p, errors="coerce")
    res["mcnemar_fdr"] = pd.Series(stats.false_discovery_control(mc, method="bh"), index=res.index) if mc.notna().any() else np.nan
if "adj_logistic_p" in res:
    lp = pd.to_numeric(res.adj_logistic_p, errors="coerce")
    res["adj_logistic_fdr"] = pd.Series(stats.false_discovery_control(lp, method="bh"), index=res.index) if lp.notna().any() else np.nan
res.to_csv(OUT / "f5_matching_audit.tsv", sep="\t", index=False)
(OUT / "f5_matching_audit_summary.json").write_text(json.dumps({
    "n_couplings": len(res),
    "source_atlas_rows": source_atlas_rows,
    "unique_genes": len(atlas),
    "duplicate_atlas_rows_collapsed": source_atlas_rows - len(atlas),
}, indent=2), encoding="utf-8")
print("\nwrote", OUT / "f5_matching_audit.tsv")
