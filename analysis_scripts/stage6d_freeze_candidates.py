from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
MECH = ROOT / "07_convergence" / "mechanism_convergence.tsv"
SDRF = ROOT / "02_metadata" / "E-MTAB-16655.sdrf.txt"
OUT = ROOT / "08_evidence_atlas"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(ATLAS, sep="\t", low_memory=False)
    mechanism = pd.read_csv(MECH, sep="\t")

    panel = atlas.loc[atlas.evidence_tier.str.startswith(("TIER1", "TIER2"), na=False)].copy()
    panel["freeze_role"] = panel.evidence_tier.map(lambda x: "PRIMARY" if x == "TIER1_MULTIMODAL" else "SECONDARY")
    panel["expected_wt_vs_uninfected"] = panel.expression_consensus_direction.fillna("NOT_PREDEFINED")

    convergent = mechanism.loc[mechanism.mechanism_convergent.astype(str).str.lower().eq("true")].copy()
    route_dirs = (
        convergent.groupby("human_symbol")
        .apply(lambda d: ";".join(sorted(f"{r.route}:{r.consensus_direction}" for _, r in d.iterrows())), include_groups=False)
        .rename("predefined_mechanism_route_directions")
    )
    panel = panel.merge(route_dirs, on="human_symbol", how="left")
    panel["predefined_mechanism_route_directions"] = panel.predefined_mechanism_route_directions.fillna("NONE")
    panel["holdout_interpretation_boundary"] = "Host transcriptional response; not bacterial replication, direct binding, or therapeutic causality"

    keep = [
        "human_symbol", "freeze_role", "evidence_tier", "functional_tier", "gw_direction",
        "expected_wt_vs_uninfected", "predefined_mechanism_route_directions",
        "ppi_AB_flag", "apms_high_confidence", "stable_cross_cohort",
        "any_mechanism_convergent", "lcv_spatial_supported", "fitness_artifact_compatible",
        "claim_boundary", "holdout_interpretation_boundary",
    ]
    panel = panel[keep].sort_values(["freeze_role", "evidence_tier", "human_symbol"])
    panel_path = OUT / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
    panel.to_csv(panel_path, sep="\t", index=False)

    primary = panel.loc[panel.freeze_role.eq("PRIMARY")].copy()
    primary_path = OUT / "FROZEN_PRIMARY_CANDIDATES_PRE_HOLDOUT.tsv"
    primary.to_csv(primary_path, sep="\t", index=False)

    sdrf = pd.read_csv(SDRF, sep="\t", dtype=str)
    stimulus_col = "Factor Value[stimulus]"
    sample_col = "Source Name"
    design = sdrf[[sample_col, stimulus_col]].drop_duplicates().sort_values(sample_col)
    design_path = OUT / "FROZEN_HOLDOUT_SAMPLE_DESIGN.tsv"
    design.to_csv(design_path, sep="\t", index=False)

    rules = {
        "freeze_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "holdout_accessed_for_expression_results": False,
        "primary_candidates": primary.human_symbol.tolist(),
        "candidate_panel_n": int(len(panel)),
        "primary_n": int(len(primary)),
        "secondary_n": int((panel.freeze_role == "SECONDARY").sum()),
        "fixed_contrasts": {
            "infection": "WT versus UI",
            "secretion": "WT versus dotA",
            "multi_effector": "each of A,B,C,D,E,F,FGHI,I,O,P,R,S versus WT; no pooling of distinct deletion groups",
        },
        "primary_endpoints": [
            "All four primary genes are tested without substitution.",
            "WT-versus-UI direction is compared with expected_wt_vs_uninfected; no gene is replaced after unblinding.",
            "For DOT_ICM-predefined genes, WT-versus-dotA direction is checked against the frozen mechanism direction.",
            "Across multi-effector contrasts, report gene-wise effects and 95% intervals; absence of response is not evidence of no function when redundancy is pre-registered.",
        ],
        "set_level_tests": [
            "Competitive rank enrichment of the fixed 56-gene panel against the measured gene universe for each contrast.",
            "Primary four-gene set reported separately; permutation preserves fixed set size and uses 10000 gene-label draws.",
            "Benjamini-Hochberg correction is applied across the pre-specified holdout contrasts within each endpoint family.",
        ],
        "success_rule": "At least one pre-specified primary/set-level endpoint passes FDR<0.05 with directionally interpretable effects; otherwise the atlas remains discovery-only.",
        "prohibited": [
            "No re-selection, threshold tuning, or gene substitution after opening RawCounts.csv.",
            "No interpretation of transcript changes as bacterial replication or direct effector binding.",
            "No conversion of redundancy-masked null results into negative functional evidence.",
        ],
        "inputs_sha256": {ATLAS.name: sha256(ATLAS), MECH.name: sha256(MECH), SDRF.name: sha256(SDRF)},
    }
    rules["frozen_candidate_panel_sha256"] = sha256(panel_path)
    rules["frozen_primary_panel_sha256"] = sha256(primary_path)
    rules["frozen_holdout_design_sha256"] = sha256(design_path)
    rules_path = OUT / "HOLDOUT_VALIDATION_RULES_FROZEN.json"
    rules_path.write_text(json.dumps(rules, indent=2), encoding="utf-8")

    manifest = {
        "candidate_panel": str(panel_path.relative_to(ROOT)),
        "candidate_panel_sha256": sha256(panel_path),
        "primary_panel_sha256": sha256(primary_path),
        "sample_design_sha256": sha256(design_path),
        "rules_sha256": sha256(rules_path),
        "candidate_panel_n": int(len(panel)),
        "primary_candidates": primary.human_symbol.tolist(),
        "freeze_complete": len(primary) == 4 and len(panel) >= 40,
        "holdout_expression_unlock_authorized": len(primary) == 4 and len(panel) >= 40,
    }
    manifest_path = OUT / "CANDIDATE_FREEZE_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (CP / "stage6d_candidate_freeze.done").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report = [
        "# Candidate and holdout-validation freeze", "",
        f"- Primary candidates: {', '.join(primary.human_symbol)}",
        f"- Secondary candidate panel: {manifest['candidate_panel_n'] - len(primary)} genes",
        f"- Full fixed panel: {manifest['candidate_panel_n']} genes",
        f"- Candidate-panel SHA256: `{manifest['candidate_panel_sha256']}`",
        f"- Validation-rules SHA256: `{manifest['rules_sha256']}`", "",
        "The processed holdout expression matrix had not been downloaded or inspected when this freeze was created. Sample metadata were used only to define the fixed contrasts. Candidate identities, expected directions, endpoint families, and multiplicity rules may not be changed after holdout unblinding.",
    ]
    (OUT / "STAGE6D_CANDIDATE_FREEZE.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0 if manifest["holdout_expression_unlock_authorized"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
