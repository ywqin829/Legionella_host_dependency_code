from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
conv = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
tier = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")

stable = set(conv.loc[conv.stable_cross_cohort.astype(str).str.upper().eq("TRUE"), "human_symbol"])
ev["fdr05"] = ev["adj.P.Val"] < 0.05
n_geo = ev[ev.fdr05].groupby("human_symbol")["dataset"].nunique()
core_ge3 = set(n_geo[n_geo >= 3].index)

rows = []
for _, r in tier.iterrows():
    g = r.human_symbol
    rows.append({
        "gene": g,
        "final_tier": r.final_tier,
        "meta_convergent_dataset_collapsed": g in stable,
        "multi_study_replicated_ge3": g in core_ge3,
        "both_meta_and_ge3": g in stable and g in core_ge3,
        "n_datasets_fdr05": int(n_geo.get(g, 0)),
        "n_human_studies": int(r.get("expression_n_human", 0)) if pd.notna(r.get("expression_n_human")) else 0,
    })
df = pd.DataFrame(rows)
df.to_csv(OUT / "d5_evidence_layers_candidates.tsv", sep="\t", index=False)

summary = {
    "meta_convergent_dataset_collapsed": len(stable),
    "multi_study_replicated_ge3": len(core_ge3),
    "both_meta_and_ge3": len(stable & core_ge3),
    "RAB1A": {"meta_convergent": "RAB1A" in stable, "multi_study_replicated": "RAB1A" in core_ge3,
              "n_datasets_fdr05": int(n_geo.get("RAB1A", 0))},
    "RAB10": {"meta_convergent": "RAB10" in stable, "multi_study_replicated": "RAB10" in core_ge3,
              "n_datasets_fdr05": int(n_geo.get("RAB10", 0))},
    "n_candidates_in_both_layers": int(((df.multi_study_replicated_ge3) & (df.meta_convergent_dataset_collapsed)).sum()),
    "n_Tier2_in_ge3_core": int(df.loc[(df.final_tier.str.startswith("TIER2")) & df.multi_study_replicated_ge3].shape[0]),
    "interpretation": f"{len(stable)} = dataset-collapsed signed meta-convergence. {len(core_ge3)} = replicated in >=3 independent datasets. {len(stable & core_ge3)} = both. RAB1A is meta-convergent but NOT >=3-study replicated; Tier 1 rests on PPI+CRISPR+LCV+holdout architecture, not multi-study expression replication.",
}
(OUT / "d5_evidence_layers_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=False))
