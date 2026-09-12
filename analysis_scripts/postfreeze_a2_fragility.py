from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)


def summarize(d: pd.DataFrame) -> dict:
    """Replicates stage5_expression_convergence.R per-gene summary."""
    direction = "UP" if (d.logFC > 0).sum() >= (d.logFC < 0).sum() else "DOWN"
    concordant = (d.logFC > 0).sum() if direction == "UP" else (d.logFC < 0).sum()
    z = d.signed_z.sum() / np.sqrt(len(d))
    return {
        "n_studies": len(d),
        "n_human": int((d.species == "HUMAN").sum()),
        "n_mouse": int((d.species == "MOUSE").sum()),
        "consensus_direction": direction,
        "direction_concordance": concordant / len(d),
        "median_logFC": float(d.logFC.median()),
        "combined_z": float(z),
        "combined_p": float(2 * stats.norm.sf(abs(z))),
        "n_study_fdr05": int((d["adj.P.Val"] < 0.05).sum()),
        "datasets": ",".join(sorted(d.dataset.unique())),
    }


def stable_cross(s: dict) -> bool:
    fdr = 2 * stats.norm.sf(abs(s["combined_z"]))
    return (s["n_studies"] >= 3 and s["n_human"] >= 2 and s["direction_concordance"] >= 0.75
            and fdr < 0.05 and s["n_study_fdr05"] >= 1)


def stable_human(s: dict) -> bool:
    fdr = 2 * stats.norm.sf(abs(s["combined_z"]))
    return (s["n_human"] >= 2 and s["direction_concordance"] >= 2 / 3 and fdr < 0.05 and s["n_study_fdr05"] >= 1)


# ---------------- Load ----------------
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
tier = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")
lcv_src = pd.read_csv(ROOT / "08_evidence_atlas" / "lcv_spatial_source_records.tsv", sep="\t")
pairs = pd.read_csv(ROOT / "03_interactome" / "legionella_human_pairs_evidence.tsv", sep="\t")

# candidate list: RAB1A + Tier2 by evidence-domain count desc, cap 16
t2 = tier[tier.final_tier.eq("TIER2_MULTIMODAL")].sort_values("domain_count", ascending=False)
cands = ["RAB1A"] + t2.human_symbol.tolist()
cands = cands[:17]

rows = []
for g in cands:
    d = ev[ev.human_symbol == g].copy()
    base = summarize(d)
    rows.append({"candidate": g, "final_tier": "Tier1" if g == "RAB1A" else "Tier2"})
    r = rows[-1]
    r["base_n_studies"] = base["n_studies"]
    r["base_n_human"] = base["n_human"]
    r["base_direction"] = base["consensus_direction"]
    r["base_combined_z"] = round(base["combined_z"], 3)
    r["base_stable_cross"] = stable_cross(base)
    r["base_datasets"] = base["datasets"]

    # leave-one-dataset-out
    flipping = []
    for ds in sorted(d.dataset.unique()):
        dd = d[d.dataset != ds]
        s = summarize(dd)
        if stable_cross(base) and not stable_cross(s):
            flipping.append(f"{ds}->UNSTABLE")
        elif not stable_cross(base) and stable_cross(s):
            flipping.append(f"{ds}->STABLE")
        else:
            flipping.append(f"{ds}:{s['consensus_direction']}/{s['n_studies']}")
    r["leave1_dataset"] = "; ".join(flipping)
    r["collapse_datasets"] = "; ".join(x.split("->")[0] for x in flipping if "->UNSTABLE" in x)

    # human-only
    hd = d[d.species == "HUMAN"]
    if len(hd) > 0:
        hs = summarize(hd)
        r["human_only_n"] = hs["n_studies"]
        r["human_only_direction"] = hs["consensus_direction"]
        r["human_only_combined_z"] = round(hs["combined_z"], 3)
        r["human_only_stable"] = stable_human(hs)
    else:
        r["human_only_n"] = 0; r["human_only_stable"] = False

    # PPI source robustness (A/B-class pmids)
    pp = pairs[pairs.gene == g]
    if len(pp):
        ab = pp[pp.evidence_class.isin(["A", "B"])] if "evidence_class" in pp else pp
        pmids = set()
        for s in ab.pmids.dropna():
            pmids.update(x.strip() for x in str(s).split(";") if x.strip())
        r["ppi_n_AB_pairs"] = int(len(ab))
        r["ppi_n_distinct_pmids"] = len(pmids)
    else:
        r["ppi_n_AB_pairs"] = 0
        r["ppi_n_distinct_pmids"] = 0

    # LCV source robustness
    ls = lcv_src[lcv_src.human_symbol == g]
    conds = set(ls.condition.astype(str))
    r["lcv_n_sources"] = int(len(conds))
    r["lcv_conditions"] = ",".join(sorted(conds))

res = pd.DataFrame(rows)
res.to_csv(OUT / "a2_evidence_fragility_matrix.tsv", sep="\t", index=False)

# aggregate summary
n_stable = int(res.base_stable_cross.sum())
n_human_stable = int(res.human_only_stable.sum())
n_collapse = int(res.collapse_datasets.astype(str).str.len().gt(0).sum())
stats_out = {
    "n_candidates": len(res),
    "n_RAB1A_stable_cross": int(res.loc[res.candidate == "RAB1A", "base_stable_cross"].iloc[0]),
    "n_RAB1A_human_stable": int(res.loc[res.candidate == "RAB1A", "human_only_stable"].iloc[0]),
    "n_Tier2_stable_cross": int(((res.candidate != "RAB1A") & res.base_stable_cross).sum()),
    "n_Tier2_human_stable": int(((res.candidate != "RAB1A") & res.human_only_stable).sum()),
    "n_Tier2_any_dataset_collapse": n_collapse,
}
(OUT / "a2_fragility_stats.json").write_text(json.dumps(stats_out, indent=2), encoding="utf-8")

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", None)
print(res.to_string(index=False))
print("\n=== summary ===")
print(json.dumps(stats_out, indent=2))
