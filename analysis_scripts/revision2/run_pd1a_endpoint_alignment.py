"""Run the frozen PD1A endpoint–modality alignment diagnostic.

PD1A is a post hoc explanatory analysis. It constructs four new panels from
discovery-only fields, evaluates them on the already-used WP6 holdout endpoint,
and refuses interpretation if any archived WP5/WP6 object has changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))

from build_wp5_comparator_panels import (  # noqa: E402
    ATLAS_COLUMNS,
    collapse_atlas,
    load_discovery_atlas,
    randomized_rank,
)

OUT = ROOT / "revision2_2026-11-23"
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
HOLDOUT = ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz"
SPEC = OUT / "REVISION2_PD1A_EXECUTION_SPEC.md"
WP5_MANIFEST = OUT / "REVISION2_WP5_FREEZE_MANIFEST.json"
WP6_MANIFEST = OUT / "REVISION2_WP6_MANIFEST.json"
WP5_PANELS = OUT / "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz"
WP6_RESULTS = OUT / "REVISION2_WP6_PANEL_RESULTS.tsv.gz"
WP6_SUMMARY = OUT / "REVISION2_WP6_STRATEGY_SUMMARY.tsv"

PANELS_OUT = OUT / "REVISION2_PD1A_PANELS.tsv.gz"
DRAW_OUT = OUT / "REVISION2_PD1A_DRAW_RESULTS.tsv.gz"
SUMMARY_OUT = OUT / "REVISION2_PD1A_STRATEGY_SUMMARY.tsv"
CONTRAST_OUT = OUT / "REVISION2_PD1A_CONTRASTS.tsv"
AUDIT_OUT = OUT / "REVISION2_PD1A_AUDIT.tsv"
STATS_OUT = OUT / "REVISION2_PD1A_STATS.json"
REPORT_OUT = OUT / "REVISION2_PD1A_REPORT.md"
MANIFEST_OUT = OUT / "REVISION2_PD1A_MANIFEST.json"

SPEC_HASH = "46492e13bd4a0dac8d553f899d65c67b0ef33cbfd70a26f772d4d7d91cd1d9f4"
PANEL_SIZE = 56
N_DRAWS = 1000
TOL = 1e-12

STRATEGIES = {
    "D1_COUNT4_NO_RESPONSE": {
        "label": "Four-domain count excluding response",
        "code": 81,
        "universe": "all",
        "key": "count4_no_response",
    },
    "D2_RESPONSE_ONLY": {
        "label": "Response domain only",
        "code": 82,
        "universe": "all",
        "key": "response_only",
    },
    "D3_COUNT4_NO_RESPONSE_EXPR": {
        "label": "Four-domain count excluding response | expression eligible",
        "code": 83,
        "universe": "expr",
        "key": "count4_no_response",
    },
    "D4_RESPONSE_ONLY_EXPR": {
        "label": "Response domain only | expression eligible",
        "code": 84,
        "universe": "expr",
        "key": "response_only",
    },
}

REFERENCE_MAP = {
    "D1_COUNT4_NO_RESPONSE": "S4_COUNT",
    "D2_RESPONSE_ONLY": "S4_COUNT",
    "D3_COUNT4_NO_RESPONSE_EXPR": "S7_COUNT_EXPR",
    "D4_RESPONSE_ONLY_EXPR": "S7_COUNT_EXPR",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def validate_archives() -> dict[str, str]:
    """Verify archived scientific inputs/outputs named by WP5/WP6 manifests.

    The project lock is an append-only governance document and was deliberately
    updated after WP5. Its historical WP5 hash is retained in that manifest,
    but its current hash is provenance rather than an equality requirement.
    """
    if sha256(SPEC) != SPEC_HASH:
        raise RuntimeError("PD1A execution specification hash mismatch")
    checked: dict[str, str] = {rel(SPEC): SPEC_HASH}
    for manifest_path in (WP5_MANIFEST, WP6_MANIFEST):
        manifest = read_json(manifest_path)
        checked[rel(manifest_path)] = sha256(manifest_path)
        for section in ("inputs", "outputs"):
            for name, expected in dict(manifest.get(section, {})).items():
                path = ROOT / str(name)
                if not path.is_file():
                    raise FileNotFoundError(path)
                actual = sha256(path)
                if path == ROOT / "13_postfreeze_interpretation" / "REVISION2_POSTREVIEW_ANALYSIS_LOCK.md":
                    checked[rel(path)] = actual
                    continue
                if actual.lower() != str(expected).lower():
                    raise RuntimeError(f"archived {section[:-1]} hash mismatch: {path}")
                checked[rel(path)] = actual
    return checked


def construct_discovery_panels() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    raw = load_discovery_atlas(str(ATLAS), ATLAS_COLUMNS)
    universe, collapse_stats = collapse_atlas(raw)
    universe["count4_no_response"] = (
        universe[["direct", "functional_clean", "mechanism", "spatial"]]
        .astype(int)
        .sum(axis=1)
    )
    universe["response_only"] = universe["response"].astype(int)
    expected_count5 = universe["count4_no_response"] + universe["response_only"]
    if not expected_count5.eq(universe["discovery_domain_count"]).all():
        raise RuntimeError("five-domain decomposition failed")

    expr = universe.loc[universe.expression_n_studies.fillna(0).ge(1)].copy()
    rows: list[dict[str, object]] = []
    for strategy_id, spec in STRATEGIES.items():
        pool = universe if spec["universe"] == "all" else expr
        for draw in range(1, N_DRAWS + 1):
            selected = randomized_rank(pool, [str(spec["key"])], int(spec["code"]), draw)
            for rank, item in enumerate(selected.itertuples(index=False), start=1):
                rows.append(
                    {
                        "strategy_id": strategy_id,
                        "strategy_label": spec["label"],
                        "draw": draw,
                        "rank": rank,
                        "gene": item.gene,
                        "expected_direction": item.expected_direction,
                        "direction_status": (
                            "EVALUABLE" if item.expected_direction in {"UP", "DOWN"} else "NOT_EVALUABLE"
                        ),
                        "expression_n_studies": item.expression_n_studies,
                        "direct": int(item.direct),
                        "functional_clean": int(item.functional_clean),
                        "response": int(item.response),
                        "mechanism": int(item.mechanism),
                        "spatial": int(item.spatial),
                        "count4_no_response": int(item.count4_no_response),
                        "response_only": int(item.response_only),
                        "discovery_domain_count": int(item.discovery_domain_count),
                    }
                )
    return pd.DataFrame(rows), universe, collapse_stats


def signed_z(frame: pd.DataFrame) -> np.ndarray:
    p = np.clip(frame.PValue.to_numpy(float), 1e-300, 1 - 1e-15)
    return np.sign(frame.logFC.to_numpy(float)) * norm.isf(p / 2)


def load_holdout_z() -> tuple[dict[str, float], int]:
    h = pd.read_csv(HOLDOUT, sep="\t", low_memory=False)
    h = h.loc[h.contrast_id.eq("WT_vs_UI")].copy()
    h["gene"] = h.human_symbol.fillna("").astype(str).str.strip().str.upper()
    h = h.loc[h.gene.ne("")].copy()
    if h.gene.duplicated().any():
        raise RuntimeError("duplicate WT_vs_UI gene symbols in holdout")
    h["signed_z"] = signed_z(h)
    h = h.loc[np.isfinite(h.signed_z)].copy()
    return h.set_index("gene").signed_z.to_dict(), int(len(h))


def evaluate(
    panels: pd.DataFrame, zmap: dict[str, float], archived_panels: pd.DataFrame
) -> pd.DataFrame:
    frozen = set(archived_panels.loc[archived_panels.strategy_id.eq("S1_FROZEN"), "gene"])
    archived_sets = {
        (sid, int(draw)): set(d.gene.astype(str))
        for (sid, draw), d in archived_panels.loc[
            archived_panels.strategy_id.isin({"S4_COUNT", "S7_COUNT_EXPR"})
        ].groupby(["strategy_id", "draw"], sort=False)
    }
    rows: list[dict[str, object]] = []
    for (sid, draw), d in panels.groupby(["strategy_id", "draw"], sort=False):
        direction = d.expected_direction.map({"UP": 1.0, "DOWN": -1.0})
        present = d.gene.map(zmap).notna()
        valid = direction.notna() & present
        raw = d.loc[valid, "gene"].map(zmap).to_numpy(float)
        oriented = direction.loc[valid].to_numpy(float) * raw
        if len(oriented) == 0 or not np.isfinite(oriented).all():
            raise RuntimeError(f"non-finite or empty result for {sid}/{draw}")
        genes = set(d.gene.astype(str))
        reference_id = REFERENCE_MAP[str(sid)]
        reference = archived_sets[(reference_id, int(draw))]
        rows.append(
            {
                "strategy_id": sid,
                "strategy_label": d.strategy_label.iloc[0],
                "draw": int(draw),
                "selected_n": int(len(d)),
                "evaluable_n": int(len(oriented)),
                "not_evaluable_direction_n": int(direction.isna().sum()),
                "holdout_missing_n": int((direction.notna() & ~present).sum()),
                "statistic": float(np.mean(oriented)),
                "direction_concordance": float(np.mean(oriented > 0)),
                "overlap_frozen_n": int(len(genes & frozen)),
                "reference_strategy_id": reference_id,
                "overlap_reference_same_draw_n": int(len(genes & reference)),
            }
        )
    return pd.DataFrame(rows)


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    fields = [
        "evaluable_n",
        "statistic",
        "direction_concordance",
        "overlap_frozen_n",
        "overlap_reference_same_draw_n",
    ]
    for sid, d in results.groupby("strategy_id", sort=False):
        row: dict[str, object] = {
            "strategy_id": sid,
            "strategy_label": d.strategy_label.iloc[0],
            "draws": int(len(d)),
            "selected_n": int(d.selected_n.iloc[0]),
        }
        for field in fields:
            values = d[field].to_numpy(float)
            row[f"{field}_median"] = float(np.median(values))
            row[f"{field}_p025"] = float(np.quantile(values, 0.025))
            row[f"{field}_p975"] = float(np.quantile(values, 0.975))
            row[f"{field}_min"] = float(np.min(values))
            row[f"{field}_max"] = float(np.max(values))
        rows.append(row)
    return pd.DataFrame(rows)


def build_contrasts(summary: pd.DataFrame, archived_summary: pd.DataFrame) -> pd.DataFrame:
    new = summary.set_index("strategy_id").statistic_median.to_dict()
    old = archived_summary.set_index("strategy_id").statistic_median.to_dict()
    definitions = [
        ("C1_COUNT5_MINUS_COUNT4", "S4_COUNT", "D1_COUNT4_NO_RESPONSE"),
        ("C2_COUNT5_MINUS_RESPONSE", "S4_COUNT", "D2_RESPONSE_ONLY"),
        ("C3_COUNT4_MINUS_FROZEN", "D1_COUNT4_NO_RESPONSE", "S1_FROZEN"),
        ("C4_RESPONSE_MINUS_FROZEN", "D2_RESPONSE_ONLY", "S1_FROZEN"),
        ("C5_COUNT5_EXPR_MINUS_COUNT4_EXPR", "S7_COUNT_EXPR", "D3_COUNT4_NO_RESPONSE_EXPR"),
        ("C6_COUNT5_EXPR_MINUS_RESPONSE_EXPR", "S7_COUNT_EXPR", "D4_RESPONSE_ONLY_EXPR"),
        ("C7_COUNT4_EXPR_MINUS_FROZEN", "D3_COUNT4_NO_RESPONSE_EXPR", "S1_FROZEN"),
        ("C8_RESPONSE_EXPR_MINUS_FROZEN", "D4_RESPONSE_ONLY_EXPR", "S1_FROZEN"),
    ]

    def value(strategy: str) -> float:
        return float(new[strategy] if strategy in new else old[strategy])

    return pd.DataFrame(
        [
            {
                "contrast_id": cid,
                "left_strategy": left,
                "right_strategy": right,
                "left_statistic_median": value(left),
                "right_statistic_median": value(right),
                "median_difference": value(left) - value(right),
                "inference_class": "post_hoc_descriptive",
            }
            for cid, left, right in definitions
        ]
    )


def audit(
    panels: pd.DataFrame,
    universe: pd.DataFrame,
    results: pd.DataFrame,
    archived_summary: pd.DataFrame,
    collapse_stats: dict[str, int],
    wp5_manifest: dict[str, object],
    archive_hashes_before: dict[str, str],
) -> pd.DataFrame:
    checks: list[tuple[str, bool, str]] = []

    def add(check_id: str, passed: bool, detail: str) -> None:
        checks.append((check_id, bool(passed), detail))

    add("P01_ARCHIVE_HASHES", True, f"{len(archive_hashes_before)} archived/spec hashes verified")
    add(
        "P02_COLLAPSE_REPRODUCED",
        collapse_stats == dict(wp5_manifest["collapse_stats"]),
        json.dumps(collapse_stats, sort_keys=True),
    )
    draw_sizes = panels.groupby(["strategy_id", "draw"]).size()
    draw_unique = panels.groupby(["strategy_id", "draw"]).gene.nunique()
    add("P03_PANEL_SIZE", bool(draw_sizes.eq(PANEL_SIZE).all()), f"draws={len(draw_sizes)}")
    add("P04_PANEL_UNIQUE", bool(draw_unique.eq(PANEL_SIZE).all()), "56 unique genes per draw")
    draw_counts = panels.groupby("strategy_id").draw.nunique().to_dict()
    add("P05_DRAW_COUNTS", all(v == N_DRAWS for v in draw_counts.values()) and len(draw_counts) == 4, str(draw_counts))
    count4 = universe[["direct", "functional_clean", "mechanism", "spatial"]].astype(int).sum(axis=1)
    add("P06_COUNT4_DEFINITION", bool(count4.eq(universe.count4_no_response).all()), "response excluded")
    add("P07_RESPONSE_DEFINITION", bool(universe.response.astype(int).eq(universe.response_only).all()), "response indicator only")
    forbidden = {"logFC", "PValue", "FDR", "signed_z", "oriented_z", "contrast_id"}
    add("P08_NO_HOLDOUT_IN_PANELS", not bool(forbidden & set(panels.columns)), str(sorted(forbidden & set(panels.columns))))
    allowed_directions = {"UP", "DOWN", "NOT_EVALUABLE"}
    add("P09_DIRECTION_LABELS", set(panels.expected_direction.astype(str)).issubset(allowed_directions), str(sorted(set(panels.expected_direction))))
    old = archived_summary.set_index("strategy_id").statistic_median
    expected = {"S1_FROZEN": 3.818344814213004, "S4_COUNT": 5.423932581137072, "S7_COUNT_EXPR": 5.4371622307997285}
    old_match = all(abs(float(old[k]) - v) <= TOL for k, v in expected.items())
    add("P10_ARCHIVED_G3_RESULTS", old_match, str({k: float(old[k]) for k in expected}))
    finite_cols = ["evaluable_n", "statistic", "direction_concordance", "overlap_frozen_n", "overlap_reference_same_draw_n"]
    finite = np.isfinite(results[finite_cols].to_numpy(float)).all()
    add("P11_FINITE_RESULTS", bool(finite), f"rows={len(results)}")
    add("P12_RESULT_STRUCTURE", len(results) == 4 * N_DRAWS and results.groupby("strategy_id").draw.nunique().eq(N_DRAWS).all(), "4,000 draw results")
    archive_hashes_after = validate_archives()
    add("P13_ARCHIVES_UNCHANGED", archive_hashes_before == archive_hashes_after, "before/after hashes identical")
    return pd.DataFrame(checks, columns=["check_id", "pass", "detail"]).assign(
        status=lambda x: np.where(x["pass"], "PASS", "FAIL")
    )[["check_id", "status", "detail"]]


def preflight() -> None:
    hashes = validate_archives()
    raw = load_discovery_atlas(str(ATLAS), ATLAS_COLUMNS)
    universe, collapse_stats = collapse_atlas(raw)
    wp5 = read_json(WP5_MANIFEST)
    if collapse_stats != dict(wp5["collapse_stats"]):
        raise RuntimeError("collapsed atlas no longer reproduces WP5")
    print(
        json.dumps(
            {
                "status": "PASS",
                "formal_holdout_statistics_calculated": False,
                "verified_hashes": len(hashes),
                "discovery_genes": len(universe),
                "collapse_stats": collapse_stats,
            },
            indent=2,
        )
    )


def execute() -> None:
    outputs = [PANELS_OUT, DRAW_OUT, SUMMARY_OUT, CONTRAST_OUT, AUDIT_OUT, STATS_OUT, REPORT_OUT, MANIFEST_OUT]
    existing = [str(p) for p in outputs if p.exists()]
    if existing:
        raise FileExistsError(f"PD1A refuses to overwrite existing outputs: {existing}")
    started = datetime.now(timezone.utc).isoformat()
    archive_hashes = validate_archives()
    panels, universe, collapse_stats = construct_discovery_panels()
    archived_panels = pd.read_csv(WP5_PANELS, sep="\t", low_memory=False)
    archived_summary = pd.read_csv(WP6_SUMMARY, sep="\t")
    zmap, holdout_n = load_holdout_z()
    results = evaluate(panels, zmap, archived_panels)
    summary = summarize(results)
    contrasts = build_contrasts(summary, archived_summary)
    wp5 = read_json(WP5_MANIFEST)
    audit_table = audit(panels, universe, results, archived_summary, collapse_stats, wp5, archive_hashes)
    if audit_table.status.ne("PASS").any():
        raise RuntimeError("PD1A audit failed: " + ", ".join(audit_table.loc[audit_table.status.ne("PASS"), "check_id"]))

    panels.to_csv(PANELS_OUT, sep="\t", index=False, compression="gzip")
    results.to_csv(DRAW_OUT, sep="\t", index=False, compression="gzip")
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    contrasts.to_csv(CONTRAST_OUT, sep="\t", index=False)
    audit_table.to_csv(AUDIT_OUT, sep="\t", index=False)

    stat_map = summary.set_index("strategy_id").statistic_median.to_dict()
    stats = {
        "stage": "PD1A",
        "classification": "post_hoc_explanatory",
        "gate": "PASS",
        "checks_passed": int(audit_table.status.eq("PASS").sum()),
        "checks_total": int(len(audit_table)),
        "holdout_contrast": "WT_vs_UI",
        "holdout_finite_genes": holdout_n,
        "panel_size": PANEL_SIZE,
        "draws_per_strategy": N_DRAWS,
        "strategy_statistic_medians": {k: float(v) for k, v in stat_map.items()},
        "contrasts": {
            r.contrast_id: float(r.median_difference) for r in contrasts.itertuples(index=False)
        },
        "interpretation_boundary": (
            "Descriptive diagnostic on the already-used holdout; not an independent validation "
            "or a direct superiority test."
        ),
    }
    write_json(STATS_OUT, stats)

    report_lines = [
        "# PD1A endpoint–modality alignment diagnostic",
        "",
        "**Status:** PASS  ",
        "**Classification:** post hoc explanatory diagnostic on the already-used PMA–THP-1 holdout.",
        "",
        "## Strategy results",
        "",
        "| Strategy | Median evaluable n | Median statistic | 2.5th–97.5th percentile | Median concordance |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in summary.itertuples(index=False):
        report_lines.append(
            f"| {r.strategy_id} | {r.evaluable_n_median:.1f} | {r.statistic_median:.4f} | "
            f"{r.statistic_p025:.4f}–{r.statistic_p975:.4f} | {r.direction_concordance_median:.3f} |"
        )
    report_lines.extend(
        [
            "",
            "## Prespecified descriptive contrasts",
            "",
            "| Contrast | Median difference |",
            "|---|---:|",
        ]
    )
    for r in contrasts.itertuples(index=False):
        report_lines.append(f"| {r.contrast_id} | {r.median_difference:.4f} |")
    report_lines.extend(
        [
            "",
            "## Inference boundary",
            "",
            "These results diagnose whether the response domain contributes to the five-domain count result. "
            "They do not constitute independent validation, a direct test of strategy superiority, or evidence "
            "of functional dependency. WP5 and WP6 archives were hash-verified before and after execution.",
            "",
        ]
    )
    REPORT_OUT.write_text("\n".join(report_lines), encoding="utf-8")

    output_hashes = {
        rel(p): sha256(p)
        for p in (PANELS_OUT, DRAW_OUT, SUMMARY_OUT, CONTRAST_OUT, AUDIT_OUT, STATS_OUT, REPORT_OUT)
    }
    manifest = {
        "stage": "PD1A",
        "status": "COMPLETE",
        "classification": "post_hoc_explanatory",
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": archive_hashes,
        "script": {rel(Path(__file__).resolve()): sha256(Path(__file__).resolve())},
        "outputs": output_hashes,
        "gate": "PASS",
        "checks_passed": int(len(audit_table)),
        "checks_total": int(len(audit_table)),
    }
    write_json(MANIFEST_OUT, manifest)
    print(json.dumps(stats, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
    else:
        execute()


if __name__ == "__main__":
    main()
