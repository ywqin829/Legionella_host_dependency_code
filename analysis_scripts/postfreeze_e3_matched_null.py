from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)

# ---- panel ----
panel = sorted(set(pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
                   .loc[lambda d: d.contrast_id == "WT_vs_UI", "human_symbol"]))

# ---- holdout per-gene signed_z per contrast ----
h = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz", sep="\t", low_memory=False)
h["signed_z"] = np.sign(h["logFC"]) * stats.norm.ppf(np.clip(h["PValue"] / 2, 1e-300, 0.5 - 1e-300))

# ---- discovery covariates ----
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
conv = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
evz = ev.assign(az=ev.signed_z.abs())
cov = evz.groupby("human_symbol")["az"].agg(["mean", "std", "count"]).rename(columns={"mean": "resp_mean", "std": "resp_sd", "count": "n_studies"})
cov["resp_mean"] = cov["resp_mean"].fillna(0)
cov["resp_sd"] = cov["resp_sd"].fillna(cov["resp_sd"].median())

# baseline UI expression percentile (recompute: 3 human datasets, UI samples)
def load_ui(uf, mat, gsm, gmap, probe):
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
    load_ui("GSE61535", "GSE61535_log2_expression.tsv.gz", ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    load_ui("GSE135803", "GSE135803_fpkm.tsv.gz", ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    load_ui("GSE287755", "GSE287755_raw_counts.tsv.gz", ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
], axis=1).mean(axis=1, skipna=True)
cov["base_expr"] = cov.index.map(ui).fillna(0)

# ---- universe = holdout-measured genes with covariates ----
uni_genes = sorted(set(h.human_symbol) & set(cov.index) & set(cov.dropna(subset=["resp_mean", "base_expr"]).index))
cov = cov.loc[uni_genes].copy()
panel_in = [g for g in panel if g in cov.index]
print(f"universe {len(uni_genes)} genes; panel genes present {len(panel_in)}")

# ---- matching covariates (z-scored) ----
FEAT = ["resp_mean", "resp_sd", "base_expr", "n_studies"]
Z = (cov[FEAT] - cov[FEAT].mean()) / cov[FEAT].std().replace(0, 1)
panel_flag = cov.index.isin(panel_in).astype(int)

clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
clf.fit(Z.values, panel_flag)
prop = clf.predict_proba(Z.values)[:, 1]

def matched_null_sets(n_sets=5000, caliper=0.05):
    """Randomised matched null sets: for each panel gene, sample one non-panel
    gene from the propensity-score caliper pool (no replacement within a set)."""
    q = np.where(panel_flag == 1)[0]
    ctl = np.where(panel_flag == 0)[0]
    pools = {}
    for i in q:
        pool = [j for j in ctl if abs(prop[j] - prop[i]) <= caliper]
        if len(pool) < 10:  # widen caliper if too few
            pool = [j for j in ctl if abs(prop[j] - prop[i]) <= caliper * 3]
        pools[i] = np.asarray(pool)
    sets = []
    for _ in range(n_sets):
        used = set()
        chosen = []
        for i in q:
            avail = [j for j in pools[i] if j not in used]
            if not avail:
                break
            chosen.append(int(RNG.choice(avail)))
            used.add(chosen[-1])
        sets.append(chosen)
    return sets

N_SETS = 5000
null_sets = matched_null_sets(N_SETS)
# report pool sizes to confirm matching is non-degenerate
pool_sizes = []
q = np.where(panel_flag == 1)[0]
ctl = np.where(panel_flag == 0)[0]
for i in q:
    pool_sizes.append(int(sum(1 for j in ctl if abs(prop[j] - prop[i]) <= 0.05)))
print(f"matched null sets built: {len(null_sets)}; median panel-gene matched-pool size = {int(np.median(pool_sizes))} (min {min(pool_sizes)})")

multi = ["A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "E_vs_WT", "F_vs_WT",
         "FGHI_vs_WT", "I_vs_WT", "O_vs_WT", "P_vs_WT", "R_vs_WT", "S_vs_WT"]

rows = []
panel_z = h.set_index(["contrast_id", "human_symbol"])["signed_z"]
for c in multi:
    pz = panel_z.loc[(c, panel_in)].dropna()
    if len(pz) < 20:
        continue
    obs = float(pz.abs().mean())
    d = h[h.contrast_id == c].set_index("human_symbol")["signed_z"]
    null = []
    for s in null_sets:
        gz = d.loc[[cov.index[i] for i in s if cov.index[i] in d.index]].dropna()
        null.append(float(gz.abs().mean()))
    null = np.array(null)
    p_emp = float((null >= obs).mean() + 1e-6)
    rows.append({"contrast": c, "panel_mean_abs_z": round(obs, 4),
                 "null_mean": round(float(null.mean()), 4), "null_sd": round(float(null.std()), 4),
                 "null_p95": round(float(np.percentile(null, 95)), 4),
                 "empirical_p_matched_null": round(p_emp, 5)})
res = pd.DataFrame(rows)
# FDR across contrasts
res["fdr_matched_null"] = pd.Series(stats.false_discovery_control(res.empirical_p_matched_null, method="bh"), index=res.index)
res["enriched_after_fdr"] = res.fdr_matched_null < 0.05
res.to_csv(OUT / "d2_matched_null_panel_enrichment.tsv", sep="\t", index=False)

n_sig = int(res.enriched_after_fdr.sum())
summary = {
    "n_contrasts": len(res),
    "n_contrasts_still_enriched_after_matched_null_FDR": n_sig,
    "note": "Matched null sets match the 54-gene panel on discovery infection-responsiveness (mean|z|), its SD, baseline UI expression percentile and n discovery studies, via propensity nearest-neighbour.",
}
(OUT / "d2_matched_null_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(res.to_string(index=False))
print("\n" + json.dumps(summary, indent=2, ensure_ascii=False))
