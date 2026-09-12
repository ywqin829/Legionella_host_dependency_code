from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
PANEL = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
HOLDOUT = ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv"
SET_TESTS = ROOT / "09_holdout_validation" / "holdout_frozen_set_tests.tsv"
OUT = ROOT / "10_final_tiers"
CP = ROOT / "checkpoints"


def b(value) -> bool:
    return str(value).lower() == "true"


def assign(row: pd.Series) -> str:
    if row.holdout_directional_contradiction:
        return "REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION"
    if row.fitness_artifact_compatible and row.domain_count < 4:
        return "TIER3_FITNESS_CAVEAT"
    strict_tier1 = (
        row.evidence_tier == "TIER1_MULTIMODAL"
        and row.direct_evidence
        and row.functional_clean
        and row.holdout_directional_support
        and row.translatome_mechanism_support
        and not row.fitness_artifact_compatible
    )
    if strict_tier1:
        return "TIER1_EXTERNALLY_VALIDATED"
    if row.domain_count >= 3:
        return "TIER2_MULTIMODAL"
    if row.domain_count >= 2:
        return "TIER3_CONTEXT_DEPENDENT"
    return "REJECT_INSUFFICIENT_CONVERGENCE"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(ATLAS, sep="\t", low_memory=False)
    panel = pd.read_csv(PANEL, sep="\t")
    holdout = pd.read_csv(HOLDOUT, sep="\t")
    sets = pd.read_csv(SET_TESTS, sep="\t")

    cols = [
        "human_symbol", "evidence_tier", "functional_tier", "ppi_AB_flag", "apms_high_confidence",
        "fitness_artifact_compatible", "stable_cross_cohort", "any_mechanism_convergent",
        "lcv_spatial_supported", "mechanism_translation_effector_module_convergent",
        "mechanism_rna_rpf_interaction_convergent", "expression_n_human", "claim_boundary",
    ]
    x = panel[["human_symbol", "freeze_role", "expected_wt_vs_uninfected", "predefined_mechanism_route_directions"]].merge(
        atlas[[c for c in cols if c in atlas.columns]], on="human_symbol", how="left")
    for c in ["ppi_AB_flag", "apms_high_confidence", "fitness_artifact_compatible", "stable_cross_cohort", "any_mechanism_convergent", "lcv_spatial_supported", "mechanism_translation_effector_module_convergent", "mechanism_rna_rpf_interaction_convergent"]:
        if c not in x:
            x[c] = False
        x[c] = x[c].fillna(False).map(b)

    h = holdout.copy()
    h["direction_match_bool"] = h.direction_match.astype(str).str.lower().map({"true": True, "false": False})
    directional = h.loc[h.expected_direction.isin(["UP", "DOWN"])].copy()
    directional["support"] = (directional.FDR < 0.05) & directional.direction_match_bool.eq(True)
    directional["contradiction"] = (directional.FDR < 0.05) & directional.direction_match_bool.eq(False)
    dsum = directional.groupby("human_symbol").agg(
        holdout_directional_support=("support", "max"),
        holdout_directional_contradiction=("contradiction", "max"),
        holdout_directional_tests=("contrast_id", "nunique"),
    ).reset_index()
    multi = h.loc[~h.contrast_id.isin(["WT_vs_UI", "WT_vs_dotA"])].groupby("human_symbol").agg(
        holdout_multieffector_fdr05_n=("FDR", lambda v: int((v < 0.05).sum())),
        holdout_multieffector_max_abs_logfc=("logFC", lambda v: float(v.abs().max())),
    ).reset_index()
    x = x.merge(dsum, on="human_symbol", how="left").merge(multi, on="human_symbol", how="left")
    x["holdout_directional_support"] = x.holdout_directional_support.fillna(False).map(b)
    x["holdout_directional_contradiction"] = x.holdout_directional_contradiction.fillna(False).map(b)
    x["holdout_multieffector_fdr05_n"] = x.holdout_multieffector_fdr05_n.fillna(0).astype(int)
    x["holdout_multieffector_support"] = x.holdout_multieffector_fdr05_n > 0

    x["direct_evidence"] = x.ppi_AB_flag | x.apms_high_confidence
    x["functional_clean"] = x.functional_tier.isin(["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT"]) & ~x.fitness_artifact_compatible
    x["translatome_mechanism_support"] = x.mechanism_translation_effector_module_convergent | x.mechanism_rna_rpf_interaction_convergent
    domain_cols = ["direct_evidence", "functional_clean", "stable_cross_cohort", "any_mechanism_convergent", "lcv_spatial_supported", "holdout_directional_support", "holdout_multieffector_support"]
    x["domain_count"] = x[domain_cols].sum(axis=1)
    x["final_tier"] = x.apply(assign, axis=1)
    x["final_claim_boundary"] = "Computationally prioritized host dependency/response candidate; experimental causality not established"
    x = x.sort_values(["final_tier", "domain_count", "human_symbol"], ascending=[True, False, True])
    x.to_csv(OUT / "final_candidate_tier_matrix.tsv", sep="\t", index=False)

    contradictions = x.loc[x.holdout_directional_contradiction | x.fitness_artifact_compatible].copy()
    contradictions.to_csv(OUT / "candidate_contradiction_and_caveat_matrix.tsv", sep="\t", index=False)

    counts = x.final_tier.value_counts().to_dict()
    tier1 = x.loc[x.final_tier.eq("TIER1_EXTERNALLY_VALIDATED"), "human_symbol"].tolist()
    panel_passes = sets.loc[sets.passes_frozen_endpoint.astype(str).str.lower().eq("true"), ["contrast_id", "set_name"]].to_dict("records")
    stats = {
        "frozen_candidates_n": int(len(x)), "final_tier_counts": counts,
        "tier1_candidates": tier1, "directional_contradictions_n": int(x.holdout_directional_contradiction.sum()),
        "fitness_caveats_n": int(x.fitness_artifact_compatible.sum()),
        "passing_predefined_set_tests": panel_passes,
        "tier1_rule": "preliminary Tier1 + direct evidence + clean replicated CRISPR + directionally significant holdout support + translatome/mechanism support",
        "gate8_final_tiers": "PASS" if len(tier1) >= 1 and not x.final_tier.isna().any() else "FAIL",
    }
    (OUT / "stage8_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = ["# Stage 8 final tier and contradiction matrix", "",
        f"- Frozen candidates assessed: {len(x)}", f"- Final tier counts: {counts}",
        f"- Strict Tier 1 candidates: {', '.join(tier1) if tier1 else 'none'}",
        f"- Directional holdout contradictions: {stats['directional_contradictions_n']}",
        f"- Fitness caveats retained and labelled: {stats['fitness_caveats_n']}",
        f"- Gate: {stats['gate8_final_tiers']}", "",
        "The strict Tier 1 rule was not relaxed after holdout unblinding. In particular, external expression support alone could not substitute for the pre-registered translatome/lung-relevance requirement. Set-level validation did not automatically upgrade individual genes.",
        "All labels remain computational evidence states, not experimental proof of direct molecular binding, bacterial replication dependence, or therapeutic action."]
    (OUT / "STAGE8_FINAL_TIER_AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage8_final_tiers.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0 if stats["gate8_final_tiers"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
