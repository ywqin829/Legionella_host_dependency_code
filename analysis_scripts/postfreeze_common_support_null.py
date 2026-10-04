from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import RAB_GTPASE, ALL_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
deg = pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv", sep="\t", low_memory=False)
# Discovery/candidate null uses the truly frozen 56-gene panel.
panel = sorted(set(pd.read_csv(ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv", sep="\t")
                   .loc[:, "human_symbol"].dropna().astype(str)))

def B(s): return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])
atlas["_ppi"] = B(atlas["ppi_AB_flag"]); atlas["_apms"] = B(atlas["apms_high_confidence"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"]); atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_nstud"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)
atlas["_common"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas = atlas.merge(deg[["gene", "cross_pathogen_degree", "cross_pathogen_pmids"]], left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)

# baseline UI expression percentile (reuse A1)
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

# --- measured candidate-ascertainment features ---
atlas["_n_ev_domains"] = atlas[["_ppi", "_apms", "_hit", "_lcv", "_resp"]].sum(axis=1)
atlas["_crispr_tier"] = atlas["functional_tier"].astype(str).map(
    {"F1_GW_AND_BOTH_BATCH_CONCORDANT": 3, "F2_MULTILAYER_CONCORDANT": 2,
     "F3_SINGLE_LAYER_OR_DISCORDANT": 1, "NO_FUNCTIONAL_HIT": 0, "NO_FUNCTIONAL_RECORD": 0}).fillna(0)

# The atlas can contain more than one row per human symbol because expression
# records retain dataset-level provenance. The null samples genes, not rows, so
# collapse to one row per human symbol before estimating propensity scores.
collapse = {
    "_ppi": "max", "_apms": "max", "_hit": "max", "_lcv": "max", "_resp": "max",
    "_crispr_tier": "max", "_nstud": "max", "_ui": "median", "_deg": "max",
    "_pmids": "max", "_common": "max", "_chronos": "median",
}
d = atlas.groupby("human_symbol", as_index=False, sort=True).agg(collapse)
d["_n_ev_domains"] = d[["_ppi", "_apms", "_hit", "_lcv", "_resp"]].sum(axis=1)
d["_chronos"] = d["_chronos"].fillna(d["_chronos"].median())

FEAT = ["_n_ev_domains", "_crispr_tier", "_lcv", "_ppi", "_nstud", "_ui",
        "_deg", "_pmids", "_common", "_chronos"]
d = d.dropna(subset=FEAT).copy()
d["_in_panel"] = d.human_symbol.isin(panel).astype(int)
if not d.human_symbol.is_unique:
    raise RuntimeError("Selection-aware null universe must contain unique human symbols")
if int(d["_in_panel"].sum()) != len(panel):
    missing = sorted(set(panel) - set(d.loc[d["_in_panel"] == 1, "human_symbol"]))
    raise RuntimeError(f"Panel genes missing from null universe: {missing}")
Z = d[FEAT].astype(float)
Z = (Z - Z.mean()) / Z.std().replace(0, 1)

clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
clf.fit(Z.values, d["_in_panel"].values)
prop = clf.predict_proba(Z.values)[:, 1]
q_idx = np.where(d["_in_panel"].values == 1)[0]; ctl_idx = np.where(d["_in_panel"].values == 0)[0]
pools = {i: np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.02]) for i in q_idx}
adequate_q = [i for i,p in pools.items() if len(p) >= 20]
q_idx = np.array(adequate_q, dtype=int)
for i, p in pools.items():
    if len(p) < 20: pools[i] = np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.05])
    if len(p) < 20: pools[i] = ctl_idx

rab_set = set(RAB_GTPASE); traf_set = set(ALL_GENES)
N = 10000
counts_rab = np.empty(N); counts_traf = np.empty(N)
for k in range(N):
    used = set(); chosen = []
    for i in q_idx:
        avail = [j for j in pools[i] if j not in used]
        if not avail: break
        chosen.append(int(RNG.choice(avail))); used.add(chosen[-1])
    genes = d.iloc[chosen].human_symbol
    counts_rab[k] = len(genes[genes.isin(rab_set)])
    counts_traf[k] = len(genes[genes.isin(traf_set)])
    if (k + 1) % 2000 == 0: print(f"  perm {k+1}/{N}")

obs_rab = len(panel_set := set(panel) & rab_set)
obs_traf = len(set(panel) & traf_set)
p_rab = float(((counts_rab >= obs_rab).sum() + 1) / (N + 1))
p_traf = float(((counts_traf >= obs_traf).sum() + 1) / (N + 1))

print(f"\nobserved: RAB in panel = {obs_rab}; trafficking in panel = {obs_traf}")
print(f"selection-aware null: RAB median={np.median(counts_rab):.1f} p95={np.percentile(counts_rab,95):.0f} max={counts_rab.max():.0f}; empirical P(#RAB>={obs_rab})={p_rab:.5f}")
print(f"trafficking null: median={np.median(counts_traf):.1f} p95={np.percentile(counts_traf,95):.0f} max={counts_traf.max():.0f}; empirical P(#Traf>={obs_traf})={p_traf:.5f}")

out = {"n_pseudo_panels": N, "n_panel_original": len(panel), "n_panel_common_support": len(q_idx), "universe_n": len(d),
       "source_atlas_rows": int(len(atlas)),
       "duplicate_atlas_rows_collapsed": int(len(atlas) - len(d)),
       "obs_rab_in_panel": obs_rab, "obs_trafficking_in_panel": obs_traf,
       "null_rab_median": float(np.median(counts_rab)), "null_rab_p95": float(np.percentile(counts_rab, 95)), "null_rab_max": float(counts_rab.max()),
       "null_traf_median": float(np.median(counts_traf)), "null_traf_p95": float(np.percentile(counts_traf, 95)), "null_traf_max": float(counts_traf.max()),
       "empirical_p_rab": p_rab, "empirical_p_trafficking": p_traf,
       "empirical_p_convention": "(exceedances + 1)/(N + 1)",
       "n_enrichment_tests": 2,
       "bonferroni_adjusted_p_rab": min(1.0, 2 * p_rab),
       "bonferroni_adjusted_p_trafficking": min(1.0, 2 * p_traf),
       "matching_vars": FEAT,
       "interpretation": "Selection-aware null samples unique-gene 56-gene pseudo-panels matched on measured candidate-ascertainment features, restricted to genes with adequate common support. The RAB test uses the canonical 61-gene RAB GTPase family. Significant results indicate that the observed enrichment is not explained by the measured matching features; unmeasured ascertainment remains possible."}
(OUT / "common_support_selection_null.json").write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
np.save(OUT / "common_support_null_rab_counts.npy", counts_rab)
np.save(OUT / "common_support_null_traf_counts.npy", counts_traf)
print("\n" + json.dumps(out, indent=2, ensure_ascii=False))
