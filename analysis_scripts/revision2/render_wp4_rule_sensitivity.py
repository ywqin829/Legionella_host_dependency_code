"""Render the discovery-only WP4 candidate-rule sensitivity package.

This script consumes only the frozen discovery panel, discovery atlas, trafficking
gene-set definitions, and the already-audited Gate G2 outputs. It does not read or
evaluate the transcriptional holdout.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revision2_2026-11-23"
VARIANT_TSV = OUT / "REVISION2_VARIANT_EFFECTS.tsv"
VARIANT_JSON = OUT / "REVISION2_VARIANT_EFFECTS.json"
RETENTION_TSV = OUT / "REVISION2_VARIANT_GENE_RETENTION.tsv"
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
FROZEN = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
TRAFFIC_CODE = ROOT / "scripts" / "trafficking_universe.py"

TABLE_OUT = OUT / "REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv"
REPORT_OUT = OUT / "REVISION2_WP4_RULE_SENSITIVITY_REPORT.md"
SVG_OUT = OUT / "REVISION2_WP4_RULE_SENSITIVITY_FIGURE.svg"
PNG_OUT = OUT / "REVISION2_WP4_RULE_SENSITIVITY_FIGURE.png"
STATS_OUT = OUT / "REVISION2_WP4_RULE_SENSITIVITY_STATS.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fisher_row(query_n: int, member_n: int, universe_n: int, set_n: int) -> tuple[float, float, float]:
    """Return fold, one-sided Fisher P and two-family Bonferroni P."""
    if not (0 <= member_n <= query_n <= universe_n and member_n <= set_n <= universe_n):
        raise ValueError((query_n, member_n, universe_n, set_n))
    table = [
        [member_n, query_n - member_n],
        [set_n - member_n, universe_n - query_n - set_n + member_n],
    ]
    _, p_value = stats.fisher_exact(table, alternative="greater")
    fold = (member_n / query_n) / (set_n / universe_n)
    return float(fold), float(p_value), float(min(1.0, 2.0 * p_value))


def shortened(variant: str) -> str:
    return {
        "ORIGINAL": "Original",
        "A1_admit_F3": "A1",
        "A2_fitness_does_not_block": "A2",
        "A4_T2_orthogonal_no_spatial": "A4",
        "A6_T2_orthogonal_require_response_and_mechanism": "A6",
        "A7_drop_T2_functional_response_branch": "A7 / B1",
        "B1_direct_only": "B1 (same panel as A7)",
        "D3_panel_min_discovery_domains": "D3 (≥3 domains)",
        "D4_panel_min_discovery_domains": "D4 (≥4 domains)",
        "D5_panel_min_discovery_domains": "D5 (≥5 domains)",
        "A3_T1_orthogonal_response_only": "A3",
        "A5_T1_require_response_and_mechanism": "A5",
        "A8_T1_require_four_domains": "A8 (Tier 1 ≥4 domains)",
        "B2_drop_T3_single_route": "B2",
    }[variant]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gate = json.loads(VARIANT_JSON.read_text(encoding="utf-8"))
    if gate["gate_g2"] != "PASS" or gate["checks_passed"] != gate["checks_total"]:
        raise RuntimeError("WP4 rendering requires a fully passing Gate G2")
    if gate.get("holdout_accessed") is not False:
        raise RuntimeError("Gate G2 does not certify discovery-only execution")

    variants = pd.read_csv(VARIANT_TSV, sep="\t", keep_default_na=False)
    retention = pd.read_csv(RETENTION_TSV, sep="\t")
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["human_symbol"], low_memory=False)
    frozen = pd.read_csv(FROZEN, sep="\t")

    sys.path.insert(0, str(ROOT / "scripts"))
    from trafficking_universe import ALL_GENES, GENE_TO_FAMILY, RAB_ASSOCIATED, RAB_GTPASE  # noqa: PLC0415

    universe = set(atlas["human_symbol"].dropna().astype(str))
    rab = universe & set(RAB_GTPASE)
    rab_assoc = universe & set(RAB_ASSOCIATED)
    trafficking = universe & set(GENE_TO_FAMILY)
    if (len(universe), len(rab), len(rab_assoc), len(ALL_GENES), len(trafficking)) != (26761, 61, 77, 276, 275):
        raise AssertionError("Background universe or gene-set sizes changed")

    frozen_genes = set(frozen["human_symbol"].dropna().astype(str))
    original = {
        "variant": "ORIGINAL",
        "changed_requirement": "Frozen candidate rule",
        "declared_class": "reference",
        "measured_class": "reference",
        "panel_n": len(frozen_genes),
        "jaccard_vs_original": 1.0,
        "added_n": 0,
        "removed_n": 0,
        "added_genes": "",
        "removed_genes": "",
        "tier1_n": int((frozen["evidence_tier"] == "TIER1_MULTIMODAL").sum()),
        "tier2_n": int(frozen["evidence_tier"].str.startswith("TIER2_").sum()),
        "primary_genes": ",".join(sorted(frozen.loc[frozen["freeze_role"] == "PRIMARY", "human_symbol"])),
        "rab_gtpase_retained": len(frozen_genes & rab),
        "rab_associated_retained": len(frozen_genes & rab_assoc),
        "trafficking_retained": len(frozen_genes & trafficking),
        "formal_reporting": "reference",
        "duplicate_panel_group": "",
    }
    if (original["panel_n"], original["tier1_n"], original["tier2_n"]) != (56, 4, 52):
        raise AssertionError("Frozen panel composition changed")

    keep_columns = list(original)
    rows = [original]
    for _, source in variants.iterrows():
        row = {column: source.get(column, "") for column in keep_columns}
        row["formal_reporting"] = "included"
        row["duplicate_panel_group"] = ""
        if source["variant"] == "B1_direct_only":
            row["formal_reporting"] = "audit row; merged with A7 in summaries/figure"
            row["duplicate_panel_group"] = "A7=B1"
        elif source["variant"] == "A7_drop_T2_functional_response_branch":
            row["duplicate_panel_group"] = "A7=B1"
        elif source["measured_class"] == "tier_only":
            row["formal_reporting"] = "tier-only; panel membership fixed by construction"
        rows.append(row)

    table = pd.DataFrame(rows)
    numeric = [
        "panel_n", "jaccard_vs_original", "added_n", "removed_n", "tier1_n", "tier2_n",
        "rab_gtpase_retained", "rab_associated_retained", "trafficking_retained",
    ]
    for column in numeric:
        table[column] = pd.to_numeric(table[column])

    for prefix, count_column, set_size in [
        ("rab", "rab_gtpase_retained", len(rab)),
        ("trafficking", "trafficking_retained", len(trafficking)),
    ]:
        computed = table.apply(
            lambda row: fisher_row(
                int(row["panel_n"]), int(row[count_column]), len(universe), set_size
            ),
            axis=1,
        )
        table[f"{prefix}_fold"] = [value[0] for value in computed]
        table[f"{prefix}_fisher_p"] = [value[1] for value in computed]
        table[f"{prefix}_bonferroni_p_two_families"] = [value[2] for value in computed]

    table.insert(1, "display_label", table["variant"].map(shortened))
    table.to_csv(TABLE_OUT, sep="\t", index=False, float_format="%.10g")

    # The formal figure counts A7 and B1 once because their selected gene sets are identical.
    figure_ids = [
        "ORIGINAL", "A1_admit_F3", "A2_fitness_does_not_block",
        "A4_T2_orthogonal_no_spatial", "A6_T2_orthogonal_require_response_and_mechanism",
        "A7_drop_T2_functional_response_branch", "D3_panel_min_discovery_domains",
        "D4_panel_min_discovery_domains", "D5_panel_min_discovery_domains",
    ]
    figure = table.set_index("variant").loc[figure_ids].reset_index()
    labels = figure["display_label"].tolist()
    colors = ["#606060"] + ["#377eb8"] * 5 + ["#e68613"] * 3
    y = np.arange(len(figure))

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8})
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.6), gridspec_kw={"width_ratios": [1.2, 1.45, 0.9]})

    ax = axes[0]
    ax.barh(y, figure["panel_n"], color=colors, edgecolor="white")
    for pos, row in figure.iterrows():
        ax.text(row["panel_n"] + 1.2, pos, f"n={int(row['panel_n'])}; J={row['jaccard_vs_original']:.2f}", va="center", fontsize=7)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, max(100, int(figure["panel_n"].max()) + 14))
    ax.set_xlabel("Selected genes")
    ax.set_title("A  Panel membership")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    rab_pct = 100 * figure["rab_gtpase_retained"] / figure["panel_n"]
    tr_pct = 100 * figure["trafficking_retained"] / figure["panel_n"]
    ax.scatter(rab_pct, y - 0.12, marker="o", s=30, color="#4daf4a", label="RAB GTPase")
    ax.scatter(tr_pct, y + 0.12, marker="s", s=30, color="#984ea3", label="Trafficking")
    for pos, row in figure.iterrows():
        ax.text(rab_pct.iloc[pos] + 1.1, pos - 0.12, str(int(row["rab_gtpase_retained"])), va="center", fontsize=6.5)
        ax.text(tr_pct.iloc[pos] + 1.1, pos + 0.12, str(int(row["trafficking_retained"])), va="center", fontsize=6.5)
    ax.set_yticks(y, [])
    ax.tick_params(axis="y", left=False)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of selected panel (%)\n(labels show gene counts)")
    ax.set_title("B  Family signal")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)

    ax = axes[2]
    distribution = gate["gene_retention"]["local_retention_distribution"]
    retained = np.array([3, 4, 5])
    counts = np.array([distribution[str(value)] for value in retained])
    ax.bar(retained, counts, color="#377eb8", width=0.7)
    for x_value, count in zip(retained, counts):
        ax.text(x_value, count + 0.7, str(count), ha="center", fontsize=8)
    ax.set_xticks(retained)
    ax.set_ylim(0, max(counts) + 5)
    ax.set_xlabel("Local panels retaining each\noriginal gene (of 5)")
    ax.set_ylabel("Original genes")
    ax.set_title("C  Per-gene retention")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Discovery-only candidate-rule sensitivity", fontsize=11, y=0.995)
    fig.text(0.01, 0.01, "Blue: local one-requirement variants; orange: structural domain-count series. A7 and B1 are one identical panel.", fontsize=7)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(SVG_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=600, bbox_inches="tight")
    plt.close(fig)

    local_ids = [
        "A1_admit_F3", "A2_fitness_does_not_block", "A4_T2_orthogonal_no_spatial",
        "A6_T2_orthogonal_require_response_and_mechanism",
        "A7_drop_T2_functional_response_branch",
    ]
    local = table.set_index("variant").loc[local_ids]
    domain = table.set_index("variant").loc[
        ["D3_panel_min_discovery_domains", "D4_panel_min_discovery_domains", "D5_panel_min_discovery_domains"]
    ]
    all_family_significant = bool(
        (table.loc[table["measured_class"].isin(["panel_effect", "reference"]),
                   ["rab_bonferroni_p_two_families", "trafficking_bonferroni_p_two_families"]] < 0.05).all().all()
    )

    report = f"""# WP4 discovery-only candidate-rule sensitivity report

## Scope and gate

- Gate G2: **PASS ({gate['checks_passed']}/{gate['checks_total']})**.
- Transcriptional holdout accessed: **no**.
- Frozen reference panel: **56 genes** (Tier 1 = 4; Tier 2 = 52).
- A7 and B1 selected the same 36 genes and are counted once in summaries and the figure.
- A8 at a three-domain Tier 1 threshold was a mathematical no-op and is excluded from formal results.

## Membership sensitivity

Across five distinct local one-requirement panels, panel size ranged from **{int(local.panel_n.min())} to {int(local.panel_n.max())}** and Jaccard similarity to the frozen panel ranged from **{local.jaccard_vs_original.min():.3f} to {local.jaccard_vs_original.max():.3f}**. Of the 56 frozen genes, **13** were retained in all five panels, **40** in at least four, **16** in exactly three, and **none** in fewer than three.

The structural discovery-domain series was more restrictive: the ≥3, ≥4 and ≥5 thresholds retained **{int(domain.iloc[0].panel_n)}, {int(domain.iloc[1].panel_n)} and {int(domain.iloc[2].panel_n)}** genes. The ≥5 result is a declared boundary case with one admitted gene, not a local robustness perturbation.

## Family signal

RAB GTPase counts across the local panels were **{int(local.rab_gtpase_retained.min())}-{int(local.rab_gtpase_retained.max())}** and trafficking counts were **{int(local.trafficking_retained.min())}-{int(local.trafficking_retained.max())}**. The structural domain series retained RAB/trafficking counts of **5/7, 2/3 and 1/1**. The curated trafficking list contains **{len(ALL_GENES)} genes**, of which **{len(trafficking)}** are represented in the 26,761-gene atlas and therefore enter the enrichment denominator. Using this atlas background and the same one-sided Fisher test as the frozen analysis, both prespecified family signals remained significant after correction for the two families within each formally reported membership-changing panel, including the declared one-gene boundary row: **{'yes' if all_family_significant else 'no'}**. These P values are descriptive sensitivity checks on reused gene sets, not independent discovery tests, and no across-variant multiplicity claim is made. Fold enrichment is reported in the table but is not treated as a stability criterion because it is mechanically affected by panel size.

## Tier-only sensitivity

A3, A5 and A8 (Tier 1 ≥4 domains) kept the 56-gene panel fixed but reduced Tier 1 membership from 4 to 1. B2 also kept panel membership fixed and retained four Tier 1 genes, while changing 1,707 row-level tier assignments. These are tier-assignment results, not evidence of panel-membership robustness.

## Interpretation fixed before manuscript wording

The panel was stable to several local rule modifications but sensitive to structural assumptions concerning cross-domain evidence requirements. Tier assignments were additionally sensitive to some rule changes that did not alter panel membership.

## Files and reproducibility

- Reader-facing statistics: `{TABLE_OUT}`
- Editable figure: `{SVG_OUT}`
- 600-dpi preview: `{PNG_OUT}`
- Machine-readable audit: `{STATS_OUT}`
- Per-gene source: `{RETENTION_TSV}`
"""
    REPORT_OUT.write_text(report, encoding="utf-8")

    inputs = [VARIANT_TSV, VARIANT_JSON, RETENTION_TSV, ATLAS, FROZEN, TRAFFIC_CODE]
    outputs = [TABLE_OUT, REPORT_OUT, SVG_OUT, PNG_OUT]
    audit = {
        "analysis": "WP4 discovery-only candidate-rule sensitivity formal report",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "gate_g2": gate["gate_g2"],
        "holdout_accessed": False,
        "universe_n": len(universe),
        "gene_set_sizes_in_atlas": {
            "rab_gtpase": len(rab),
            "rab_associated": len(rab_assoc),
            "trafficking": len(trafficking),
        },
        "curated_trafficking_gene_list_n": len(ALL_GENES),
        "a7_b1_counted_once": True,
        "formal_membership_panels_n": int((table["measured_class"].isin(["panel_effect", "reference"])).sum() - 1),
        "all_formal_membership_panels_rab_and_trafficking_bonferroni_p_lt_0_05": all_family_significant,
        "local_retention": gate["gene_retention"],
        "input_hashes": {str(path.resolve()): sha256(path) for path in inputs},
        "output_hashes": {str(path.resolve()): sha256(path) for path in outputs},
    }
    STATS_OUT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "outputs": [str(path) for path in outputs + [STATS_OUT]]}, indent=2))


if __name__ == "__main__":
    main()
