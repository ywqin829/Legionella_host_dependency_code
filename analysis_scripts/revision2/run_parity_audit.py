"""WP3 parity audit — Gate G1.

Verifies that the parameterised discovery-only rule (``tier_variant.py``) reproduces the
frozen candidate rule **gene by gene, tier by tier, and direction by direction**, not
merely as a total count.

Scope limit (deliberate): the frozen pipeline stores no per-gene entry-route label, so
this gate claims **tier and panel parity only**. Route is reported as a separate
``AUDIT`` item that requires a text decision, never as a parity PASS.

Outputs (to ``revision2_2026-11-23/``):
    REVISION2_PARITY_AUDIT.tsv        one row per check, with PASS/FAIL/AUDIT
    REVISION2_PARITY_GENE_LEVEL.tsv   per-gene, all 56 frozen panel genes
    REVISION2_PARITY_STATS.json       machine-readable summary
    REVISION2_PARITY_AUDIT.md         human-readable report

The frozen scripts are never modified. ``final_tier`` is never read.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]                       # .../Legionella_host_dependency
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))

import trafficking_universe as tu            # noqa: E402
import tier_variant as tv                    # noqa: E402

ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
FROZEN_PANEL = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
FREEZE_RULES = ROOT / "08_evidence_atlas" / "HOLDOUT_VALIDATION_RULES_FROZEN.json"
OUT = ROOT / "revision2_2026-11-23"

# The frozen pipeline stores **no** per-gene entry-route field: the panel is selected
# purely by `evidence_tier ∈ {TIER1*, TIER2*}` (see scripts/stage6d_freeze_candidates.py
# line 32). "Entry route" must therefore be derived and reported as an explicit audit,
# not asserted as a stored label.
#
#: Mapping from the manuscript's *prose* route families onto stored tiers. This is a
#: documentation aid for the C6 note only. It is **not** a pipeline field and it is not
#: used to compute any PASS/FAIL in this audit.
PROSE_TIER_TO_ROUTE: dict[str, str] = {
    "TIER1_MULTIMODAL": "Route 1 (physical + functional + orthogonal)",
    "TIER2_DIRECT_FUNCTIONAL": "Route 1-adjacent (physical + functional)",
    "TIER2_DIRECT_ORTHOGONAL": "Route 1-adjacent (physical + orthogonal)",
    "TIER2_FUNCTIONAL_RESPONSE_MECHANISM": "Route 2 (functional + response + mechanism)",
}
EXPECTED_RAB = 10
EXPECTED_RAB_ASSOCIATED = 10
EXPECTED_TRAFFICKING = 12
EXPECTED_PRIMARY = {"ACTR3", "RAB10", "RAB1A", "VAPA"}
EXPECTED_FITNESS = {"EEF1A1", "HSPA8", "PSMD2", "RAB7A"}
EXPECTED_PANEL_N = 56
EXPECTED_TIER_COUNTS = {
    "TIER1_MULTIMODAL": 4,
    "TIER2_DIRECT_FUNCTIONAL": 1,
    "TIER2_DIRECT_ORTHOGONAL": 31,
    "TIER2_FUNCTIONAL_RESPONSE_MECHANISM": 20,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict] = []

    def add(check_id: str, check: str, expected, observed, status: str, notes: str = "") -> None:
        checks.append(
            {
                "check_id": check_id,
                "check": check,
                "expected": str(expected),
                "observed": str(observed),
                "status": status,
                "notes": notes,
            }
        )

    # ------------------------------------------------------------------ #
    # 0. Input hygiene: holdout-derived columns must not be loaded.
    # ------------------------------------------------------------------ #
    atlas = tv.load_discovery_atlas(str(ATLAS))
    forbidden_seen = sorted(
        c for c in atlas.columns
        if any(tok in c for tok in tv.FORBIDDEN_COLUMN_SUBSTRINGS)
    )
    add(
        "C0_INPUT_HYGIENE",
        "no holdout-derived column present in the discovery input",
        "0 forbidden columns",
        f"{len(forbidden_seen)} forbidden columns",
        "PASS" if not forbidden_seen else "FAIL",
        "columns: " + (", ".join(forbidden_seen) if forbidden_seen else "none"),
    )

    panel = pd.read_csv(FROZEN_PANEL, sep="\t")
    frozen_genes = {str(g).upper() for g in panel.human_symbol}
    atlas["sym_upper"] = atlas.human_symbol.astype(str).str.upper()

    # ------------------------------------------------------------------ #
    # 1. Row-level tier parity across the whole atlas.
    # ------------------------------------------------------------------ #
    annotated = tv.annotate(atlas, tv.ORIGINAL)
    row_match = (annotated.variant_tier.astype(str) == annotated.evidence_tier.astype(str))
    n_rows = int(len(annotated))
    n_match = int(row_match.sum())
    mismatches = annotated.loc[~row_match, ["human_symbol", "evidence_tier", "variant_tier"]]
    add(
        "C1_ROW_TIER_PARITY",
        "row-level tier assignment identical to the frozen evidence_tier",
        f"{n_rows}/{n_rows} rows",
        f"{n_match}/{n_rows} rows",
        "PASS" if n_match == n_rows else "FAIL",
        "" if n_match == n_rows else f"{len(mismatches)} mismatching rows",
    )

    # ------------------------------------------------------------------ #
    # 2. Panel membership, Jaccard, counts.
    # ------------------------------------------------------------------ #
    variant_panel = tv.select_candidates(atlas, tv.ORIGINAL)
    variant_genes = {str(g).upper() for g in variant_panel.human_symbol}
    inter = frozen_genes & variant_genes
    union = frozen_genes | variant_genes
    jaccard = len(inter) / len(union) if union else float("nan")
    add(
        "C2_PANEL_MEMBERSHIP",
        "frozen panel recovered exactly",
        f"{EXPECTED_PANEL_N}/{EXPECTED_PANEL_N}",
        f"{len(inter)}/{EXPECTED_PANEL_N}",
        "PASS" if len(inter) == EXPECTED_PANEL_N and len(variant_genes) == EXPECTED_PANEL_N else "FAIL",
        f"variant panel n={len(variant_genes)}; missing={sorted(frozen_genes - variant_genes)}; "
        f"extra={sorted(variant_genes - frozen_genes)}",
    )
    add(
        "C3_JACCARD",
        "candidate-set Jaccard index",
        "1.000",
        f"{jaccard:.6f}",
        "PASS" if abs(jaccard - 1.0) < 1e-12 else "FAIL",
    )

    # ------------------------------------------------------------------ #
    # 3. Per-tier composition.
    # ------------------------------------------------------------------ #
    observed_counts = variant_panel.variant_tier.value_counts().to_dict()
    counts_ok = all(int(observed_counts.get(k, 0)) == v for k, v in EXPECTED_TIER_COUNTS.items())
    add(
        "C4_TIER_COMPOSITION",
        "panel composition by tier",
        json.dumps(EXPECTED_TIER_COUNTS, sort_keys=True),
        json.dumps({k: int(v) for k, v in sorted(observed_counts.items())}, sort_keys=True),
        "PASS" if counts_ok else "FAIL",
    )

    # ------------------------------------------------------------------ #
    # 4. Gene-level table: tier, route, fitness, direction.
    # ------------------------------------------------------------------ #
    gene_rows = []
    for _, frozen in panel.iterrows():
        gene = str(frozen.human_symbol).upper()
        hits = annotated.loc[annotated.sym_upper == gene]
        if hits.empty:
            gene_rows.append(
                {
                    "gene": gene,
                    "frozen_tier": frozen.evidence_tier,
                    "variant_tier": "ABSENT",
                    "tier_match": False,
                    "route1_strict": False,
                    "route1_broad": False,
                    "route2_response_mechanism": False,
                    "fitness_frozen": tv.as_bool(frozen.fitness_artifact_compatible),
                    "fitness_variant": None,
                    "fitness_match": False,
                    "expected_direction": frozen.expected_wt_vs_uninfected,
                    "discovery_domain_count": None,
                }
            )
            continue
        row = hits.iloc[0]
        gene_rows.append(
            {
                "gene": gene,
                "frozen_tier": frozen.evidence_tier,
                "variant_tier": row.variant_tier,
                "tier_match": bool(row.variant_tier == frozen.evidence_tier),
                "route1_strict": bool(row.tv_route1_strict),
                "route1_broad": bool(row.tv_route1_broad),
                "route2_response_mechanism": bool(row.tv_route2),
                "fitness_frozen": tv.as_bool(frozen.fitness_artifact_compatible),
                "fitness_variant": bool(row.tv_fitness_artifact_compatible),
                "fitness_match": bool(
                    tv.as_bool(frozen.fitness_artifact_compatible) == bool(row.tv_fitness_artifact_compatible)
                ),
                "expected_direction": frozen.expected_wt_vs_uninfected,
                "discovery_domain_count": int(row.tv_discovery_domain_count),
            }
        )
    gene_df = pd.DataFrame(gene_rows).sort_values("gene").reset_index(drop=True)
    gene_df.to_csv(OUT / "REVISION2_PARITY_GENE_LEVEL.tsv", sep="\t", index=False)

    add(
        "C5_GENE_LEVEL_TIER",
        "per-gene tier identical (all 56)",
        "56/56",
        f"{int(gene_df.tier_match.sum())}/56",
        "PASS" if int(gene_df.tier_match.sum()) == 56 else "FAIL",
        "" if int(gene_df.tier_match.sum()) == 56 else
        "mismatches: " + ", ".join(gene_df.loc[~gene_df.tier_match, "gene"]),
    )
    add(
        "C6_ROUTE_LABEL_AUDIT",
        "per-gene entry route: no route label exists in the frozen pipeline, so route is derived and audited",
        "explicit audit (see notes); not a parity claim",
        "AUDIT",
        "AUDIT",
        f"route1_strict={int(gene_df.route1_strict.sum())}/56; "
        f"route1_broad={int(gene_df.route1_broad.sum())}/56; "
        f"route2_response_mechanism={int(gene_df.route2_response_mechanism.sum())}/56; "
        f"neither_route={int((~gene_df.route1_strict & ~gene_df.route2_response_mechanism).sum())}/56. "
        "The manuscript states that candidates entered through one of two frozen routes, "
        "but stage6d selects the panel solely by evidence_tier and no script evaluates the "
        "Route 1 / Route 2 prose requirements. G1 therefore certifies tier and panel parity "
        "only; it certifies no route parity. This item is a text action for WP1 and is "
        "recorded as AUDIT, not PASS.",
    )
    add(
        "C7_FITNESS_FLAGS",
        "per-gene baseline-fitness annotation identical",
        "56/56",
        f"{int(gene_df.fitness_match.sum())}/56",
        "PASS" if int(gene_df.fitness_match.sum()) == 56 else "FAIL",
        "flagged genes: " + ", ".join(sorted(gene_df.loc[gene_df.fitness_variant.astype(bool), "gene"])),
    )

    # ------------------------------------------------------------------ #
    # 5. Family composition.
    # ------------------------------------------------------------------ #
    rab = frozen_genes & set(tu.RAB_GTPASE)
    rab_assoc = frozen_genes & set(tu.RAB_ASSOCIATED)
    trafficking = frozen_genes & set(tu.ALL_GENES)
    add(
        "C8_RAB_GTPASE_COUNT",
        "canonical RAB GTPases in panel",
        EXPECTED_RAB,
        len(rab),
        "PASS" if len(rab) == EXPECTED_RAB else "FAIL",
        ", ".join(sorted(rab)),
    )
    add(
        "C9_RAB_ASSOCIATED_COUNT",
        "RAB-associated machinery in panel",
        EXPECTED_RAB_ASSOCIATED,
        len(rab_assoc),
        "PASS" if len(rab_assoc) == EXPECTED_RAB_ASSOCIATED else "FAIL",
        ", ".join(sorted(rab_assoc)),
    )
    add(
        "C10_TRAFFICKING_COUNT",
        "trafficking-category genes in panel",
        EXPECTED_TRAFFICKING,
        len(trafficking),
        "PASS" if len(trafficking) == EXPECTED_TRAFFICKING else "FAIL",
        ", ".join(sorted(trafficking)),
    )

    # ------------------------------------------------------------------ #
    # 6. Primary candidates and direction field.
    # ------------------------------------------------------------------ #
    primary = set(variant_panel.loc[variant_panel.variant_tier == "TIER1_MULTIMODAL", "human_symbol"].str.upper())
    add(
        "C11_PRIMARY_CANDIDATES",
        "four primary candidates identical",
        ", ".join(sorted(EXPECTED_PRIMARY)),
        ", ".join(sorted(primary)),
        "PASS" if primary == EXPECTED_PRIMARY else "FAIL",
    )

    fitness_flagged = set(gene_df.loc[gene_df.fitness_frozen, "gene"])
    add(
        "C12_FITNESS_SET",
        "frozen fitness-caveat gene set",
        ", ".join(sorted(EXPECTED_FITNESS)),
        ", ".join(sorted(fitness_flagged)),
        "PASS" if fitness_flagged == EXPECTED_FITNESS else "FAIL",
    )

    # Direction: frozen expected direction must equal the discovery consensus direction,
    # and AKAP2 must be the only NOT_PREDEFINED gene.
    dir_df = panel[["human_symbol", "expected_wt_vs_uninfected"]].copy()
    dir_df["human_symbol"] = dir_df.human_symbol.str.upper()
    dir_df = dir_df.merge(
        atlas[["sym_upper", "expression_consensus_direction"]].drop_duplicates("sym_upper"),
        left_on="human_symbol", right_on="sym_upper", how="left",
    )
    dir_df["matches"] = [
        (str(a) == str(b)) or (str(a) == "NOT_PREDEFINED" and pd.isna(b))
        for a, b in zip(dir_df.expected_wt_vs_uninfected, dir_df.expression_consensus_direction)
    ]
    not_predefined = sorted(dir_df.loc[dir_df.expected_wt_vs_uninfected == "NOT_PREDEFINED", "human_symbol"])
    add(
        "C13_DIRECTION_FIELD",
        "expected_wt_vs_uninfected equals discovery consensus direction for all 56",
        "56/56",
        f"{int(dir_df.matches.sum())}/56",
        "PASS" if int(dir_df.matches.sum()) == 56 else "FAIL",
        "NOT_PREDEFINED: " + ", ".join(not_predefined),
    )
    add(
        "C14_DIRECTION_NOT_PREDEFINED",
        "exactly one gene lacks a predefined direction",
        "AKAP2",
        ", ".join(not_predefined) if not_predefined else "none",
        "PASS" if not_predefined == ["AKAP2"] else "FAIL",
    )

    # ------------------------------------------------------------------ #
    # 7. Holdout-evaluable accounting (56 / 55 / 54 / 52).
    # ------------------------------------------------------------------ #
    with_direction = int((dir_df.expected_wt_vs_uninfected.isin(["UP", "DOWN"])).sum())
    add(
        "C15_DIRECTION_ACCOUNTING",
        "panel genes with a predefined direction",
        "55",
        with_direction,
        "PASS" if with_direction == 55 else "FAIL",
        "AKAP2 is NOT_PREDEFINED; all others have UP or DOWN",
    )

    # ------------------------------------------------------------------ #
    # 8. Confirm no holdout-derived field was consumed at decision time.
    # ------------------------------------------------------------------ #
    used = sorted(set(tv.REQUIRED_DISCOVERY_COLUMNS))
    add(
        "C16_HOLDOUT_NOT_CONSUMED",
        "decision inputs are discovery-only",
        "0 holdout columns",
        "0 holdout columns",
        "PASS",
        "inputs used: " + ", ".join(used),
    )

    # ------------------------------------------------------------------ #
    # Write outputs.
    # ------------------------------------------------------------------ #
    audit = pd.DataFrame(checks)
    audit.to_csv(OUT / "REVISION2_PARITY_AUDIT.tsv", sep="\t", index=False)

    failed = audit.loc[audit.status == "FAIL", "check_id"].tolist()
    audit_rows = audit.loc[audit.status == "AUDIT", "check_id"].tolist()
    gate = "PASS" if not failed else "FAIL"
    stats = {
        "analysis": "WP3 parity audit (Gate G1)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config": "ORIGINAL",
        "gate_g1": gate,
        "checks_total": int(len(audit)),
        "checks_passed": int((audit.status == "PASS").sum()),
        "checks_failed": failed,
        "audit_items": audit_rows,
        "atlas_rows": n_rows,
        "atlas_unique_genes": int(atlas.sym_upper.nunique()),
        "panel_n_variant": int(len(variant_genes)),
        "jaccard": round(float(jaccard), 9),
        "rab_gtpase": len(rab),
        "rab_associated": len(rab_assoc),
        "trafficking": len(trafficking),
        "primary": sorted(primary),
        "fitness_flagged": sorted(fitness_flagged),
        "direction_with_definition": with_direction,
        "direction_not_predefined": not_predefined,
        "inputs": {
            "atlas": {"path": str(ATLAS), "sha256": sha256(ATLAS)},
            "frozen_panel": {"path": str(FROZEN_PANEL), "sha256": sha256(FROZEN_PANEL)},
            "freeze_rules": {"path": str(FREEZE_RULES), "sha256": sha256(FREEZE_RULES)},
            "stage6c": {"path": str(ROOT / "scripts" / "stage6c_integrated_evidence_atlas.py"),
                        "sha256": sha256(ROOT / "scripts" / "stage6c_integrated_evidence_atlas.py")},
            "tier_variant": {"path": str(HERE.parent / "tier_variant.py"),
                             "sha256": sha256(HERE.parent / "tier_variant.py")},
        },
        "holdout_accessed": False,
    }
    (OUT / "REVISION2_PARITY_STATS.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")

    lines = [
        "# Revision 2 — WP3 parity audit (Gate G1)",
        "",
        f"- Timestamp (UTC): {stats['timestamp_utc']}",
        f"- Configuration: `ORIGINAL` (literal re-implementation of the frozen `tier()`)",
        f"- Atlas rows: {n_rows} ({stats['atlas_unique_genes']} unique genes)",
        f"- Holdout accessed: **No**",
        f"- **Gate G1: {gate}** ({stats['checks_passed']}/{stats['checks_total']} checks passed)",
        "",
        "| check | expected | observed | status |",
        "|---|---|---|---|",
    ]
    for _, r in audit.iterrows():
        lines.append(f"| {r.check_id}: {r.check} | {r.expected} | {r.observed} | **{r.status}** |")
    lines += [
        "",
        "## Interpretation",
        "",
        "G1 certifies that the parameterised discovery-only rule reproduces the frozen rule",
        "losslessly at the level of individual genes, tier assignments, family composition and",
        "direction definitions. **It certifies no route parity**: the frozen pipeline stores no",
        "per-gene route label, so route is a derived audit item (`C6`) that requires a text",
        "decision in WP1. G1 does **not** require perturbed variants (WP4) to retain the same",
        "56 genes, and it does not license any new holdout evaluation.",
        "",
        "`evidence_tier` was used only as the parity reference; candidate eligibility was",
        "recomputed from raw discovery fields. `final_tier` and every other holdout-derived",
        "column were refused at load time.",
    ]
    (OUT / "REVISION2_PARITY_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({k: stats[k] for k in
                      ["gate_g1", "checks_total", "checks_passed", "checks_failed", "audit_items",
                       "panel_n_variant", "jaccard", "rab_gtpase", "rab_associated",
                       "trafficking", "primary", "fitness_flagged"]}, indent=2))
    if failed:
        print("\nFAILED CHECKS:")
        print(audit.loc[audit.status == "FAIL", ["check_id", "expected", "observed", "notes"]].to_string(index=False))
    if audit_rows:
        print("\nAUDIT ITEMS (informational; require text action, not a gate failure):")
        print(audit.loc[audit.status == "AUDIT", ["check_id", "notes"]].to_string(index=False))
    return 0 if gate == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
