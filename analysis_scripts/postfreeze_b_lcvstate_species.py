from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
lcv_src = pd.read_csv(ROOT / "08_evidence_atlas" / "lcv_spatial_source_records.tsv", sep="\t")
lcv_ev = pd.read_csv(ROOT / "08_evidence_atlas" / "lcv_spatial_evidence.tsv", sep="\t")
tier = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")


def B(s):
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


atlas["_crispr_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_fitness_caveat"] = B(atlas["fitness_artifact_compatible"])
atlas["_clean"] = atlas["_crispr_hit"] & (~atlas["_fitness_caveat"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])

# ---------------- B1: LCV state stratification ----------------
# per-gene LCV state from per-condition records
cond = lcv_src.groupby("human_symbol")["condition"].apply(lambda s: set(s.astype(str))).to_dict()

def lcv_state(sym):
    c = cond.get(sym, set())
    has_wt = any("WT_LCV" in x for x in c)
    has_ankb = any("ANKB" in x for x in c)
    has_mouse = any("UNTREATED_LCV_BMM" in x for x in c)
    if has_wt and has_ankb:
        base = "HUMAN_WT_AND_ANKB"
    elif has_wt:
        base = "HUMAN_WT_ONLY"
    elif has_ankb:
        base = "HUMAN_ANKB_ONLY"
    else:
        base = "NO_HUMAN_LCV"
    return base + ("_PLUS_MOUSE" if has_mouse else "")

atlas["lcv_state"] = atlas["human_symbol"].map(lcv_state)

# universe = CRISPR-tested + LCV-evaluable
U = atlas[atlas["_crispr_hit"].notna() & atlas["_lcv"].notna()].copy()
# state collapse: any human vs mouse-only vs absent
U["state_group"] = U["lcv_state"].map(lambda x: "human_LcV" if x.startswith("HUMAN") else ("mouse_only" if x.startswith("NO_HUMAN") and "PLUS_MOUSE" in x else ("no_LCV" if x.startswith("NO_HUMAN") and "PLUS_MOUSE" not in x else x)))

# clean CRISPR rate by LCV state (human vs mouse vs none)
tbl = []
for state in ["HUMAN_WT_AND_ANKB", "HUMAN_WT_ONLY", "HUMAN_ANKB_ONLY", "NO_HUMAN_LCV_PLUS_MOUSE", "NO_HUMAN_LCV"]:
    sub = U[U.lcv_state == state]
    if len(sub) < 10:
        continue
    n_clean = int(sub["_clean"].sum())
    rate = n_clean / len(sub)
    tbl.append({"lcv_state": state, "n": len(sub), "n_clean_crispr": n_clean, "clean_rate": round(rate, 4)})
B1 = pd.DataFrame(tbl)
B1.to_csv(OUT / "b1_lcv_state_x_clean_crispr.tsv", sep="\t", index=False)
print("=== B1: clean-CRISPR rate by LCV state (all CRISPR-tested genes) ===")
print(B1.to_string(index=False))

# RAB1A / RAB10 / panel trafficking genes LCV state
panel_genes = tier.human_symbol.tolist()
key = [g for g in panel_genes if g in atlas.human_symbol.values]
key_state = atlas[atlas.human_symbol.isin(key)][["human_symbol", "lcv_state", "_clean", "_lcv"]].set_index("human_symbol")
print("\n=== Key candidate LCV states ===")
for g in ["RAB1A", "RAB10", "RAB7A", "RAB8A", "RAB8B", "VAPA", "ACTR3", "SEC22B", "ARF1", "RABIF"]:
    if g in key_state.index:
        r = key_state.loc[g]
        print(f"  {g}: {r['lcv_state']}  clean_crispr={bool(r['_clean'])}")

# ---------------- B2: species-dependency score ----------------
b2 = atlas[atlas.human_symbol.isin(key)].copy()
b2 = b2[["human_symbol", "expression_n_studies", "expression_n_human", "expression_n_mouse",
         "expression_stable_cross_cohort", "lcv_state", "ppi_n_pmids", "evidence_tier"]]
# TRUE final tier from the frozen tier matrix (only RAB1A is Tier 1)
tier_map = dict(zip(tier.human_symbol, tier.final_tier))
b2["final_tier"] = b2["human_symbol"].map(tier_map).fillna(b2["evidence_tier"])
# human fraction of expression studies
b2["n_human"] = pd.to_numeric(b2["expression_n_human"], errors="coerce").fillna(0)
b2["n_mouse"] = pd.to_numeric(b2["expression_n_mouse"], errors="coerce").fillna(0)
b2["species_dependency"] = b2.apply(lambda r: "human_only" if r["n_mouse"] == 0 and r["n_human"] > 0
                                    else ("mouse_inclusive" if r["n_mouse"] > 0 and r["n_human"] > 0
                                          else ("mouse_only" if r["n_human"] == 0 and r["n_mouse"] > 0 else "none")), axis=1)
# LCV human vs mouse
b2["lcv_human"] = b2["lcv_state"].str.startswith("HUMAN")
b2["lcv_mouse"] = b2["lcv_state"].str.contains("MOUSE")
b2 = b2.sort_values(["final_tier", "human_symbol"])
b2.to_csv(OUT / "b2_species_dependency_scores.tsv", sep="\t", index=False)
print("\n=== B2: species-dependency scores (top candidates) ===")
print(b2[["human_symbol", "n_human", "n_mouse", "species_dependency", "lcv_state", "lcv_human", "lcv_mouse", "expression_stable_cross_cohort"]].head(20).to_string(index=False))

stats_out = {
    "B1_clean_crispr_by_lcv_state": B1.to_dict("records"),
    "B2_n_candidates_scored": int(len(b2)),
    "note": "lcv_state from per-condition source records; clean = CRISPR hit without fitness caveat",
}
(OUT / "b_stats.json").write_text(json.dumps(stats_out, indent=2, ensure_ascii=False), encoding="utf-8")
