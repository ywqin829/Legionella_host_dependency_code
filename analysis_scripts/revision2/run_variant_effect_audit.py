"""WP4 variant-effect audit — Gate G2.

Measures what each declared rule variant actually does, and fails if the measured effect
class differs from the class declared *before* measurement in ``tier_variant.py``.

This gate exists to prevent two specific failure modes:

1. **Post-hoc re-description.** A variant whose panel is unchanged may not be sold as
   panel-robustness evidence. ``VARIANT_CLASS`` fixes the class up front, and
   ``DECLARED_EFFECT`` fixes the expected direction of panel movement up front.
2. **Undeclared no-ops.** A variant that changes nothing carries no information. Any
   variant with zero relabelled rows *and* an unchanged panel must either be declared
   ``noop`` or fail the gate.

Gate G2 checks
--------------
G2a  every variant's measured class equals ``VARIANT_CLASS[name]``
G2b  every variant's panel direction equals ``DECLARED_EFFECT[name]['panel']``
G2c  every variant's row-change presence equals ``DECLARED_EFFECT[name]['rows']``
G2d  no undeclared no-op exists
G2e  the removed ``A8_T1_require_three_domains`` is confirmed a mathematical no-op
G2f  the domain-count series is monotone non-increasing in panel size
G2g  ``ORIGINAL`` recovers the frozen 56-gene panel exactly

Outputs (to ``revision2_2026-11-23/``):
    REVISION2_VARIANT_EFFECTS.tsv          one row per variant
    REVISION2_VARIANT_GENE_RETENTION.tsv   one row per original panel gene
    REVISION2_VARIANT_EFFECTS.json         machine-readable summary
    REVISION2_VARIANT_EFFECTS.md           human-readable report

Discovery-only throughout: ``final_tier`` and every holdout column are refused at load.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))

import trafficking_universe as tu            # noqa: E402
import tier_variant as tv                    # noqa: E402

ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
FROZEN_PANEL = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
FREEZE_RULES = ROOT / "08_evidence_atlas" / "HOLDOUT_VALIDATION_RULES_FROZEN.json"
OUT = ROOT / "revision2_2026-11-23"

EXPECTED_PANEL_N = 56
EXPECTED_RAB = 10
EXPECTED_RAB_ASSOCIATED = 10
EXPECTED_TRAFFICKING = 12

RAB_GTPASE = {g.upper() for g in tu.RAB_GTPASE}
RAB_ASSOCIATED = {g.upper() for g in tu.RAB_ASSOCIATED}
TRAFFICKING = {g.upper() for g in tu.ALL_GENES}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def panel_direction(original: set[str], candidate: set[str]) -> str:
    added = candidate - original
    removed = original - candidate
    if added and not removed:
        return "expand"
    if removed and not added:
        return "contract"
    if not added and not removed:
        return "unchanged"
    return "mixed"


def measured_class(panel_changed: bool, rows_changed: int) -> str:
    if not panel_changed and rows_changed == 0:
        return "noop"
    if panel_changed:
        return "panel_effect"
    return "tier_only"


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else float("nan")


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
    # Load and reproduce the frozen panel with ORIGINAL.
    # ------------------------------------------------------------------ #
    atlas = tv.load_discovery_atlas(str(ATLAS))
    atlas["sym_upper"] = atlas.human_symbol.astype(str).str.upper()

    frozen = pd.read_csv(FROZEN_PANEL, sep="\t")
    frozen_genes = {str(g).upper() for g in frozen.human_symbol}

    original_annotated = tv.annotate(atlas, tv.ORIGINAL)
    original_tiers = original_annotated.variant_tier.astype(str).to_numpy()
    original_genes = tv.candidate_genes(atlas, tv.ORIGINAL)
    original_panel = tv.select_candidates(atlas, tv.ORIGINAL)

    add(
        "G2g_ORIGINAL_RECOVERS_PANEL",
        "ORIGINAL reproduces the frozen pre-holdout panel exactly",
        f"{EXPECTED_PANEL_N} genes, Jaccard 1.000000",
        f"{len(original_genes)} genes, Jaccard {jaccard(frozen_genes, original_genes):.6f}",
        "PASS" if original_genes == frozen_genes else "FAIL",
        f"missing={sorted(frozen_genes - original_genes)}; extra={sorted(original_genes - frozen_genes)}",
    )

    # ------------------------------------------------------------------ #
    # Measure every declared variant.
    # ------------------------------------------------------------------ #
    rows: list[dict] = []
    gene_panels: dict[str, set[str]] = {}

    for name, cfg in tv.ALL_CONFIGS.items():
        if name == "ORIGINAL":
            continue
        annotated = tv.annotate(atlas, cfg)
        tiers = annotated.variant_tier.astype(str).to_numpy()
        rows_changed = int((tiers != original_tiers).sum())

        panel = tv.select_candidates(atlas, cfg)
        genes = {str(g).upper() for g in panel.human_symbol}
        gene_panels[name] = genes

        added = sorted(genes - original_genes)
        removed = sorted(original_genes - genes)
        direction = panel_direction(original_genes, genes)
        klass_measured = measured_class(genes != original_genes, rows_changed)

        tier_counts = panel.variant_tier.value_counts().to_dict()
        n_tier1 = int(tier_counts.get(tv.TIER1_MULTIMODAL, 0))
        n_tier2 = int(sum(v for k, v in tier_counts.items() if str(k).startswith("TIER2")))

        declared = tv.DECLARED_EFFECT[name]
        declared_class = tv.VARIANT_CLASS[name]

        rows.append(
            {
                "variant": name,
                "changed_requirement": cfg.changed_requirement,
                "declared_class": declared_class,
                "measured_class": klass_measured,
                "class_match": klass_measured == declared_class,
                "declared_panel_direction": declared["panel"],
                "measured_panel_direction": direction,
                "direction_match": direction == declared["panel"],
                "declared_rows": declared["mode"],
                "measured_rows_changed": rows_changed,
                "rows_match": (
                    (rows_changed > 0) if declared["mode"] == "tier_relabel" else (rows_changed == 0)
                ),
                "panel_n": len(genes),
                "jaccard_vs_original": round(jaccard(original_genes, genes), 6),
                "added_n": len(added),
                "removed_n": len(removed),
                "added_genes": ",".join(added),
                "removed_genes": ",".join(removed),
                "tier1_n": n_tier1,
                "tier2_n": n_tier2,
                "primary_genes": ",".join(
                    sorted(panel.loc[panel.variant_tier == tv.TIER1_MULTIMODAL, "human_symbol"].str.upper())
                ),
                "rab_gtpase_retained": len(genes & RAB_GTPASE),
                "rab_associated_retained": len(genes & RAB_ASSOCIATED),
                "trafficking_retained": len(genes & TRAFFICKING),
                "rab_gtpase_genes_retained": ",".join(sorted(genes & RAB_GTPASE)),
                "notes": cfg.notes,
            }
        )

    effects = pd.DataFrame(rows)
    # Order: panel-effect, tier-only, then the removed no-op.
    class_order = {"panel_effect": 0, "tier_only": 1, "noop": 2}
    effects["_sort"] = effects.declared_class.map(class_order)
    effects = effects.sort_values(["_sort", "variant"]).drop(columns="_sort").reset_index(drop=True)
    effects.to_csv(OUT / "REVISION2_VARIANT_EFFECTS.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ #
    # Measure the removed no-op explicitly, from its retained config.
    # ------------------------------------------------------------------ #
    noop_rows: list[dict] = []
    for name, cfg in tv.REMOVED_NOOP.items():
        tiers = tv.annotate(atlas, cfg).variant_tier.astype(str).to_numpy()
        changed = int((tiers != original_tiers).sum())
        genes = tv.candidate_genes(atlas, cfg)
        noop_rows.append(
            {
                "variant": name,
                "rows_changed": changed,
                "panel_n": len(genes),
                "panel_identical": genes == original_genes,
            }
        )
    noop_df = pd.DataFrame(noop_rows)

    # ------------------------------------------------------------------ #
    # Per-gene retention across panel-effect variants.
    #
    # Two blocks, never pooled: the domain-count series is a deliberately extreme
    # structural series (>=5 domains admits 1 gene by design) and pooling it with local
    # single-requirement perturbations would let one extreme variant drive the summary.
    #
    # Identical panels are merged before counting. A7 and B1 produce the same 36-gene
    # set, so counting both rows would give that one observation two votes and inflate
    # every "retained in >= k" tally. Each *distinct* gene set contributes one vote.
    # ------------------------------------------------------------------ #
    LOCAL_PANEL_EFFECT = [
        "A1_admit_F3",
        "A2_fitness_does_not_block",
        "A4_T2_orthogonal_no_spatial",
        "A6_T2_orthogonal_require_response_and_mechanism",
        "A7_drop_T2_functional_response_branch",
        "B1_direct_only",
    ]
    DOMAIN_PANEL_EFFECT = [
        "D3_panel_min_discovery_domains",
        "D4_panel_min_discovery_domains",
        "D5_panel_min_discovery_domains",
    ]
    assert set(LOCAL_PANEL_EFFECT) | set(DOMAIN_PANEL_EFFECT) == set(tv.PANEL_EFFECT_VARIANTS)

    def merge_block(names: list[str]) -> list[tuple[str, set[str]]]:
        """Collapse identical gene sets to one observation each."""
        seen: dict[frozenset[str], list[str]] = {}
        for name in names:
            seen.setdefault(frozenset(gene_panels[name]), []).append(name)
        merged = [(" = ".join(sorted(v)), set(k)) for k, v in seen.items()]
        merged.sort(key=lambda item: item[0])
        return merged

    local_obs = merge_block(LOCAL_PANEL_EFFECT)
    domain_obs = merge_block(DOMAIN_PANEL_EFFECT)

    retention = []
    for gene in sorted(original_genes):
        kept_local = [label for label, panel in local_obs if gene in panel]
        kept_domain = [label for label, panel in domain_obs if gene in panel]
        retention.append(
            {
                "gene": gene,
                "n_local_panel_effect": len(local_obs),
                "n_retained_local": len(kept_local),
                "retained_in_all_local": len(kept_local) == len(local_obs),
                "variants_retaining_local": ";".join(kept_local),
                "variants_dropping_local": ";".join(
                    label for label, panel in local_obs if gene not in panel
                ),
                "n_domain_panel_effect": len(domain_obs),
                "n_retained_domain": len(kept_domain),
                "retained_in_all_domain": len(kept_domain) == len(domain_obs),
                "rab_gtpase": gene in RAB_GTPASE,
                "rab_associated": gene in RAB_ASSOCIATED,
                "trafficking": gene in TRAFFICKING,
            }
        )
    ret_df = pd.DataFrame(retention).sort_values(
        ["n_retained_local", "n_retained_domain", "gene"], ascending=[True, True, True]
    ).reset_index(drop=True)
    ret_df.to_csv(OUT / "REVISION2_VARIANT_GENE_RETENTION.tsv", sep="\t", index=False)

    # Single source of truth for retention, reused by the report text and by stats.
    retention_stats = {
        "n_local_variant_rows": len(LOCAL_PANEL_EFFECT),
        "n_local_observations": len(local_obs),
        "n_domain_observations": len(domain_obs),
        "retained_in_all_local": int(ret_df.retained_in_all_local.sum()),
        "retained_in_all_domain": int(ret_df.retained_in_all_domain.sum()),
        "retained_in_no_local_variant": int((ret_df.n_retained_local == 0).sum()),
        "retained_in_at_least_local_minus_one": int(
            (ret_df.n_retained_local >= len(local_obs) - 1).sum()
        ),
        "local_retention_distribution": {
            str(k): int((ret_df.n_retained_local == k).sum())
            for k in range(len(local_obs) + 1)
        },
        "rab_retention_local_min": int(ret_df.loc[ret_df.rab_gtpase, "n_retained_local"].min()),
        "rab_retention_local_max": int(ret_df.loc[ret_df.rab_gtpase, "n_retained_local"].max()),
        "trafficking_retention_local_min": int(ret_df.loc[ret_df.trafficking, "n_retained_local"].min()),
        "trafficking_retention_local_max": int(ret_df.loc[ret_df.trafficking, "n_retained_local"].max()),
    }

    # ------------------------------------------------------------------ #
    # Gate checks.
    # ------------------------------------------------------------------ #
    bad_class = effects.loc[~effects.class_match, "variant"].tolist()
    add(
        "G2a_DECLARED_CLASS",
        "measured effect class equals declared effect class",
        "all match",
        f"{int(effects.class_match.sum())}/{len(effects)} match",
        "PASS" if not bad_class else "FAIL",
        "" if not bad_class else "mismatched: " + ", ".join(bad_class),
    )

    bad_dir = effects.loc[~effects.direction_match, "variant"].tolist()
    add(
        "G2b_PANEL_DIRECTION",
        "measured panel direction equals declared panel direction",
        "all match",
        f"{int(effects.direction_match.sum())}/{len(effects)} match",
        "PASS" if not bad_dir else "FAIL",
        "" if not bad_dir else "mismatched: " + ", ".join(bad_dir),
    )

    bad_rows = effects.loc[~effects.rows_match, "variant"].tolist()
    add(
        "G2c_ROW_CHANGES",
        "relabelling behaviour equals declaration (tier relabel vs admission filter)",
        "all match",
        f"{int(effects.rows_match.sum())}/{len(effects)} match",
        "PASS" if not bad_rows else "FAIL",
        "" if not bad_rows else "mismatched: " + ", ".join(bad_rows),
    )

    undeclared_noop = effects.loc[
        (effects.measured_class == "noop") & (effects.declared_class != "noop"), "variant"
    ].tolist()
    add(
        "G2d_NO_UNDECLARED_NOOP",
        "no measured no-op survived into the formal variant set",
        "0 undeclared no-ops",
        f"{len(undeclared_noop)} undeclared no-ops",
        "PASS" if not undeclared_noop else "FAIL",
        "" if not undeclared_noop else "undeclared: " + ", ".join(undeclared_noop),
    )

    noop_ok = bool(
        len(noop_df) == 1
        and int(noop_df.iloc[0].rows_changed) == 0
        and bool(noop_df.iloc[0].panel_identical)
    )
    add(
        "G2e_REMOVED_A8_3_IS_NOOP",
        "removed A8>=3 variant is a mathematical no-op",
        "0 rows changed, panel identical",
        f"{int(noop_df.iloc[0].rows_changed)} rows changed, panel identical={bool(noop_df.iloc[0].panel_identical)}"
        if len(noop_df) == 1 else "missing",
        "PASS" if noop_ok else "FAIL",
        "recorded as removed from formal WP4 after detection",
    )

    dom = effects.loc[effects.variant.str.startswith("D"), ["variant", "panel_n"]]
    dom = dom.assign(cut=dom.variant.str.extract(r"D(\d)")[0].astype(int)).sort_values("cut")
    monotone = dom.panel_n.is_monotonic_decreasing
    add(
        "G2f_DOMAIN_SERIES_MONOTONE",
        "domain-count series is monotone non-increasing in panel size",
        "panel_n(>=3) >= panel_n(>=4) >= panel_n(>=5)",
        " -> ".join(f">={int(c)}:{int(n)}" for c, n in zip(dom.cut, dom.panel_n)),
        "PASS" if monotone else "FAIL",
        "domain-count sensitivity analysis, not a symmetric threshold perturbation",
    )

    # ------------------------------------------------------------------ #
    # A7 vs B1: same number is not the same set.
    # ------------------------------------------------------------------ #
    a7 = gene_panels["A7_drop_T2_functional_response_branch"]
    b1 = gene_panels["B1_direct_only"]
    a7_b1_identical = a7 == b1
    a7_b1_symmetric_difference = sorted(a7 ^ b1)
    add(
        "G2h_A7_B1_SET_EQUALITY",
        "A7 and B1 panel gene sets compared as sets, not counts",
        "explicit comparison (see notes)",
        "identical" if a7_b1_identical else f"differ in {len(a7_b1_symmetric_difference)} genes",
        "PASS",
        (
            "A7 and B1 produce the same gene set; both rows are retained for the rule record "
            "but they count as ONE robustness observation and must be merged in figures."
            if a7_b1_identical
            else "A7 and B1 differ; symmetric difference: " + ", ".join(a7_b1_symmetric_difference)
        ),
    )

    # ------------------------------------------------------------------ #
    # Pairwise set identity among panel-effect variants.
    #
    # Two different rule changes can land on the same gene set. If they do, they are one
    # robustness observation, not two, and figures must merge them.
    # ------------------------------------------------------------------ #
    groups: dict[frozenset[str], list[str]] = {}
    for name in tv.PANEL_EFFECT_VARIANTS:
        groups.setdefault(frozenset(gene_panels[name]), []).append(name)
    identical_groups = [sorted(v) for k, v in groups.items() if len(v) > 1]
    identical_groups.sort()
    add(
        "G2i_PAIRWISE_SET_IDENTITY",
        "panel-effect variants grouped by exact gene-set identity",
        "explicit grouping (see notes)",
        f"{len(identical_groups)} non-singleton group(s)",
        "PASS",
        (
            "no two panel-effect variants share a gene set"
            if not identical_groups
            else "identical sets: " + " | ".join(" = ".join(g) for g in identical_groups)
            + ". Count each group once as robustness evidence; keep every rule row in the audit table."
        ),
    )

    # ------------------------------------------------------------------ #
    # Retention must be counted per distinct panel, not per rule row.
    # ------------------------------------------------------------------ #
    n_local_rows = len(LOCAL_PANEL_EFFECT)
    n_local_obs = len(local_obs)
    add(
        "G2j_RETENTION_DEDUPLICATED",
        "retention votes counted per distinct gene set, not per rule row",
        f"{n_local_rows} rows -> {n_local_obs} distinct observations "
        f"(1 merged group: A7 = B1)",
        f"{n_local_obs} distinct local observations",
        "PASS" if n_local_obs == n_local_rows - len(identical_groups)
        and n_local_obs > 0 else "FAIL",
        "merged groups: " + " | ".join(label for label, _ in local_obs if "=" in label),
    )

    # ------------------------------------------------------------------ #
    # Retention report text, built once and checked before it is written.
    #
    # Every "... retention: a-b of N" denominator must be the number of *distinct local
    # panels*, never the number of trafficking (or RAB) genes inside a single panel.
    # An earlier draft printed "of 12" because it reused the panel composition constant
    # EXPECTED_TRAFFICKING as if it were an observation count. The emitted text is
    # therefore parsed back and asserted rather than trusted.
    # ------------------------------------------------------------------ #
    retention_section = [
        "",
        "## Gene retention across panel-effect variants",
        "",
        f"Local single-requirement perturbations ({len(local_obs)} distinct observations "
        f"from {len(LOCAL_PANEL_EFFECT)} rule rows; A7 and B1 merged):",
        f"- Panel genes retained by all {len(local_obs)}: "
        f"{retention_stats['retained_in_all_local']}/{len(original_genes)}",
        f"- Retained by at least {len(local_obs) - 1} of {len(local_obs)}: "
        f"{retention_stats['retained_in_at_least_local_minus_one']}/{len(original_genes)}",
        f"- Retained by none: {retention_stats['retained_in_no_local_variant']}",
        f"- Retention distribution (0..{len(local_obs)}): "
        + ", ".join(
            f"{k}:{v}"
            for k, v in retention_stats["local_retention_distribution"].items()
        ),
        f"- RAB GTPase retention: {retention_stats['rab_retention_local_min']}–"
        f"{retention_stats['rab_retention_local_max']} of "
        f"{len(local_obs)} distinct local panels",
        f"- Trafficking retention: {retention_stats['trafficking_retention_local_min']}–"
        f"{retention_stats['trafficking_retention_local_max']} of "
        f"{len(local_obs)} distinct local panels",
        "",
        f"Structural domain-count series ({len(domain_obs)}):",
        f"- Panel genes retained at every cut-off: "
        f"{retention_stats['retained_in_all_domain']}/{len(original_genes)}",
        "",
        "A7 vs B1 gene sets: "
        + ("**identical** — reported as one robustness observation, both rows kept in the audit table."
           if a7_b1_identical else "different — reported separately."),
        "",
        "Exact gene-set identity among all panel-effect variants:",
        (
            "- none — every panel-effect variant yields a distinct gene set"
            if not identical_groups
            else "\n".join("- " + " = ".join(g) for g in identical_groups)
        ),
        "",
    ]

    retention_denominators = [
        int(re.search(r"\d+", line.rsplit(" of ", 1)[-1]).group())
        for line in retention_section
        if " retention: " in line
    ]
    add(
        "G2k_RETENTION_DENOMINATOR",
        "retention report denominators equal the distinct-panel count",
        f"exactly 2 denominators, both == {len(local_obs)}",
        f"{len(retention_denominators)} denominator(s): {retention_denominators}",
        "PASS" if len(retention_denominators) == 2
        and all(d == len(local_obs) for d in retention_denominators) else "FAIL",
        "guards the earlier 'of 12' error, which used the panel trafficking count as a "
        "variant count",
    )

    # ------------------------------------------------------------------ #
    # Write summary.
    # ------------------------------------------------------------------ #
    audit = pd.DataFrame(checks)
    audit.to_csv(OUT / "REVISION2_VARIANT_EFFECT_AUDIT.tsv", sep="\t", index=False)
    failed = audit.loc[audit.status == "FAIL", "check_id"].tolist()
    gate = "PASS" if not failed else "FAIL"

    stats = {
        "analysis": "WP4 variant-effect audit (Gate G2)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "gate_g2": gate,
        "checks_total": int(len(audit)),
        "checks_passed": int((audit.status == "PASS").sum()),
        "checks_failed": failed,
        "original_panel_n": len(original_genes),
        "panel_effect_variants": list(tv.PANEL_EFFECT_VARIANTS),
        "tier_only_variants": list(tv.TIER_ONLY_VARIANTS),
        "removed_noop": noop_df.to_dict(orient="records"),
        "domain_series": {
            str(int(c)): {"panel_n": int(n), "rab_gtpase": int(r), "trafficking": int(t)}
            for c, n, r, t in zip(
                dom.cut,
                dom.panel_n,
                effects.set_index("variant").loc[dom.variant, "rab_gtpase_retained"],
                effects.set_index("variant").loc[dom.variant, "trafficking_retained"],
            )
        },
        "a7_b1_identical": a7_b1_identical,
        "identical_panel_effect_groups": identical_groups,
        "gene_retention": retention_stats,
        "input_hashes": {
            "atlas": {"path": str(ATLAS), "sha256": sha256(ATLAS)},
            "frozen_panel": {"path": str(FROZEN_PANEL), "sha256": sha256(FROZEN_PANEL)},
            "freeze_rules": {"path": str(FREEZE_RULES), "sha256": sha256(FREEZE_RULES)},
            "tier_variant": {"path": str(HERE.parent / "tier_variant.py"),
                             "sha256": sha256(HERE.parent / "tier_variant.py")},
        },
        "holdout_accessed": False,
    }
    (OUT / "REVISION2_VARIANT_EFFECTS.json").write_text(
        json.dumps(stats, indent=2, default=str), encoding="utf-8"
    )

    # ------------------------------------------------------------------ #
    # Markdown report.
    # ------------------------------------------------------------------ #
    lines = [
        "# Revision 2 — WP4 variant-effect audit (Gate G2)",
        "",
        f"- Timestamp (UTC): {stats['timestamp_utc']}",
        f"- Holdout accessed: **No**",
        f"- **Gate G2: {gate}** ({stats['checks_passed']}/{stats['checks_total']} checks passed)",
        f"- Original panel: {len(original_genes)} genes",
        "",
        "## Gate checks",
        "",
        "| check | expected | observed | status |",
        "|---|---|---|---|",
    ]
    for _, r in audit.iterrows():
        lines.append(f"| {r.check_id}: {r.check} | {r.expected} | {r.observed} | **{r.status}** |")

    def table(subset: pd.DataFrame, title: str) -> None:
        lines.extend(
            [
                "",
                f"## {title}",
                "",
                "| variant | panel n | Jaccard | added | removed | Tier1 | Tier2 | RAB | RAB-assoc | trafficking | rows relabelled |",
                "|---|---|---|---|---|---|---|---|---|---|---|",
            ]
        )
        for _, r in subset.iterrows():
            lines.append(
                f"| {r.variant} | {r.panel_n} | {r.jaccard_vs_original:.3f} | {r.added_n} | "
                f"{r.removed_n} | {r.tier1_n} | {r.tier2_n} | {r.rab_gtpase_retained} | "
                f"{r.rab_associated_retained} | {r.trafficking_retained} | {r.measured_rows_changed} |"
            )

    table(effects.loc[effects.declared_class == "panel_effect"], "Panel-effect variants")
    table(effects.loc[effects.declared_class == "tier_only"], "Tier-only variants")

    lines.extend(
        [
            "",
            "## Domain-count sensitivity series",
            "",
            "Domain-count sensitivity analysis, **not** a symmetric threshold perturbation.",
            "",
            "| min discovery domains | panel n | RAB GTPases | trafficking |",
            "|---|---|---|---|",
        ]
    )
    for cut, payload in stats["domain_series"].items():
        lines.append(f"| >= {cut} | {payload['panel_n']} | {payload['rab_gtpase']} | {payload['trafficking']} |")

    lines.extend(
        [
            "",
            "## Removed no-op",
            "",
            "| variant | rows relabelled | panel n | panel identical |",
            "|---|---|---|---|",
        ]
    )
    for r in noop_df.to_dict(orient="records"):
        lines.append(
            f"| {r['variant']} | {r['rows_changed']} | {r['panel_n']} | {r['panel_identical']} |"
        )

    lines.extend(retention_section)
    lines.extend(
        [
            "## Interpretation boundary",
            "",
            "Tier-only variants are **not** panel-robustness evidence. They show that tier and",
            "primary assignment can move while membership is fixed. No count of variants that",
            "\"kept 56 genes\" is reported as a robustness score.",
        ]
    )
    (OUT / "REVISION2_VARIANT_EFFECTS.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps(
        {
            "gate_g2": gate,
            "checks_passed": f"{stats['checks_passed']}/{stats['checks_total']}",
            "failed": failed,
            "domain_series": stats["domain_series"],
            "a7_b1_identical": a7_b1_identical,
            "retained_all_local": stats["gene_retention"]["retained_in_all_local"],
            "n_local_observations": stats["gene_retention"]["n_local_observations"],
            "retained_at_least_local_minus_one":
                stats["gene_retention"]["retained_in_at_least_local_minus_one"],
            "retained_in_no_local_variant": stats["gene_retention"]["retained_in_no_local_variant"],
        },
        indent=2,
    ))
    if failed:
        print("\nFAILED CHECKS:")
        print(audit.loc[audit.status == "FAIL", ["check_id", "expected", "observed", "notes"]].to_string(index=False))
    return 0 if gate == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
