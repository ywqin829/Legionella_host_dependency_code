from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CRISPR = ROOT / "04_genetic_dependency" / "ppi_crispr_intersection.tsv"
EXPR = ROOT / "07_convergence" / "primary_infection_convergence.tsv"
MECH = ROOT / "07_convergence" / "mechanism_convergence.tsv"
LCV = ROOT / "08_evidence_atlas" / "lcv_spatial_evidence.tsv"
TRUTH = ROOT / "00_protocol" / "KNOWN_POSITIVE_CONTROL_HOST_FACTORS.tsv"
OUT = ROOT / "08_evidence_atlas"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False}).fillna(False).astype(bool)


def tier(row: pd.Series) -> str:
    direct = bool(row.ppi_AB_flag or row.apms_high_confidence)
    functional = row.functional_tier in {"F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT"}
    functional_clean = functional and not row.fitness_artifact_compatible
    response = bool(row.stable_cross_cohort)
    mechanism = bool(row.any_mechanism_convergent)
    spatial = bool(row.lcv_spatial_supported)
    if direct and functional_clean and (response or mechanism or spatial):
        return "TIER1_MULTIMODAL"
    if direct and functional_clean:
        return "TIER2_DIRECT_FUNCTIONAL"
    if functional_clean and response and mechanism:
        return "TIER2_FUNCTIONAL_RESPONSE_MECHANISM"
    if direct and (response or mechanism or spatial):
        return "TIER2_DIRECT_ORTHOGONAL"
    if functional and row.fitness_artifact_compatible:
        return "TIER3_FITNESS_CAVEAT"
    if functional_clean or (response and mechanism) or direct or (spatial and (response or mechanism)):
        return "TIER3_SINGLE_ROUTE"
    return "BACKGROUND"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    required = [CRISPR, EXPR, MECH, LCV, TRUTH]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing integrated-atlas input(s): " + "; ".join(missing))

    base = pd.read_csv(CRISPR, sep="\t", low_memory=False)
    base["human_symbol"] = base.symbol.astype(str).str.upper()
    base_bool = ["ppi_AB_flag", "apms_high_confidence", "apms_supportive", "ppi_crispr_convergence", "fitness_artifact_compatible", "common_essential_flag", "u937_strong_dependency_flag"]
    for col in base_bool:
        base[col] = as_bool(base[col])

    expr = pd.read_csv(EXPR, sep="\t")
    expr["human_symbol"] = expr.human_symbol.astype(str).str.upper()
    expr_cols = ["human_symbol", "n_studies", "n_human", "n_mouse", "consensus_direction", "direction_concordance", "combined_fdr", "n_study_fdr05", "stable_cross_cohort", "naip_nlrc4_boundary"]
    expr = expr[expr_cols].rename(columns={c: f"expression_{c}" for c in expr_cols if c != "human_symbol"})
    expr["expression_stable_cross_cohort"] = as_bool(expr.expression_stable_cross_cohort)
    expr["expression_naip_nlrc4_boundary"] = as_bool(expr.expression_naip_nlrc4_boundary)

    mech = pd.read_csv(MECH, sep="\t")
    mech["human_symbol"] = mech.human_symbol.astype(str).str.upper()
    mech["mechanism_convergent"] = as_bool(mech.mechanism_convergent)
    route_flag = mech.pivot_table(index="human_symbol", columns="route", values="mechanism_convergent", aggfunc="max", fill_value=False)
    route_flag.columns = [f"mechanism_{x.lower()}_convergent" for x in route_flag.columns]
    route_fdr = mech.pivot_table(index="human_symbol", columns="route", values="combined_fdr", aggfunc="min")
    route_fdr.columns = [f"mechanism_{x.lower()}_fdr" for x in route_fdr.columns]
    mechanism = route_flag.join(route_fdr).reset_index()

    lcv = pd.read_csv(LCV, sep="\t")
    lcv["human_symbol"] = lcv.human_symbol.astype(str).str.upper()
    lcv_bool = ["lcv_human_wt", "lcv_human_ankb", "lcv_mouse_untreated", "lcv_spatial_supported"]
    for col in lcv_bool:
        lcv[col] = as_bool(lcv[col])

    atlas = base.merge(expr, on="human_symbol", how="outer").merge(mechanism, on="human_symbol", how="outer").merge(lcv, on="human_symbol", how="outer")
    bool_cols = [c for c in atlas.columns if c.endswith("_convergent")] + base_bool + ["expression_stable_cross_cohort", "expression_naip_nlrc4_boundary"] + lcv_bool
    for col in set(bool_cols):
        if col not in atlas:
            atlas[col] = False
        atlas[col] = atlas[col].fillna(False).astype(bool)
    atlas["stable_cross_cohort"] = atlas.expression_stable_cross_cohort
    mechanism_flags = [c for c in atlas.columns if c.startswith("mechanism_") and c.endswith("_convergent")]
    atlas["any_mechanism_convergent"] = atlas[mechanism_flags].any(axis=1)
    atlas["functional_tier"] = atlas.functional_tier.fillna("NO_FUNCTIONAL_RECORD")
    atlas["evidence_tier"] = atlas.apply(tier, axis=1)
    atlas["claim_boundary"] = np.select(
        [atlas.expression_naip_nlrc4_boundary, atlas.fitness_artifact_compatible, atlas.gw_direction.eq("KO_PROTECTS_HOST_ENRICHED_AFTER_INFECTION")],
        ["HUMAN_MOUSE_NAIP_NLRC4_BOUNDARY", "POTENTIAL_BASELINE_FITNESS_ARTIFACT", "HOST_SURVIVAL_SCREEN_ENDPOINT_NOT_REPLICATION"],
        default="STANDARD_COMPUTATIONAL_BOUNDARY")
    atlas.to_csv(OUT / "host_dependency_evidence_atlas.tsv.gz", sep="\t", index=False, compression="gzip")

    truth = pd.read_csv(TRUTH, sep="\t")
    truth["human_symbol"] = truth.gene_symbol.astype(str).str.upper()
    evidence = atlas.set_index("human_symbol")
    rows = []
    for _, item in truth.iterrows():
        gene = item.human_symbol
        if gene not in evidence.index:
            rows.append({**item.to_dict(), "recovered": False, "recovery_routes": "ABSENT_FROM_ATLAS", "evidence_tier": "ABSENT"})
            continue
        row = evidence.loc[gene]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        routes = []
        if bool(row.ppi_AB_flag): routes.append("PPI_AB")
        if bool(row.apms_high_confidence): routes.append("APMS")
        if row.functional_tier in {"F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT"}: routes.append("CRISPR_F1F2")
        if bool(row.stable_cross_cohort): routes.append("EXPRESSION")
        if bool(row.any_mechanism_convergent): routes.append("MECHANISM")
        if bool(row.lcv_human_wt): routes.append("LCV_HUMAN_WT")
        if bool(row.lcv_human_ankb): routes.append("LCV_HUMAN_ANKB")
        if bool(row.lcv_mouse_untreated): routes.append("LCV_MOUSE")
        rows.append({**item.to_dict(), "recovered": bool(routes), "recovery_routes": ";".join(routes) if routes else "NO_STRONG_ROUTE", "evidence_tier": row.evidence_tier})
    truth_out = pd.DataFrame(rows)
    truth_out.to_csv(OUT / "known_positive_control_recovery.tsv", sep="\t", index=False)

    recovery = float(truth_out.recovered.mean())
    stats = {
        "input_sha256": {p.name: sha256(p) for p in required},
        "atlas_genes": int(len(atlas)), "tier_counts": atlas.evidence_tier.value_counts().to_dict(),
        "lcv_spatial_genes": int(atlas.lcv_spatial_supported.sum()),
        "truth_set_n": int(len(truth_out)), "truth_recovered_n": int(truth_out.recovered.sum()),
        "truth_recovery_fraction": recovery, "tier_assignment_authorized": recovery >= 2 / 3,
        "holdout_accessed": False,
        "gate6_evidence_atlas": "PASS" if len(atlas) > 15000 and recovery >= 2 / 3 else "FAIL"}
    (OUT / "stage6_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = ["# Stage 6 Integrated Host Dependency Evidence Atlas", "",
        f"- Atlas genes: {stats['atlas_genes']}", f"- LCV spatial-evidence genes: {stats['lcv_spatial_genes']}",
        f"- Tier counts: {stats['tier_counts']}",
        f"- Known-positive recovery: {stats['truth_recovered_n']}/{stats['truth_set_n']} ({recovery:.1%})",
        f"- Gate: {stats['gate6_evidence_atlas']}", "",
        "Tiers are rule-based evidence states, not weighted scores. Pan-essential-compatible hits are downgraded rather than deleted. Expression-only or LCV-spatial-only genes cannot enter Tier 1. LCV detection is co-localization evidence, not proof of direct binding or functional dependency.",
        "The external holdout remains locked."]
    (OUT / "STAGE6_EVIDENCE_ATLAS.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage6c_integrated_evidence_atlas.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0 if stats["gate6_evidence_atlas"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
