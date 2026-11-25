"""Build and freeze reviewer-requested discovery-only comparator panels (WP5).

This script never reads holdout data. It constructs four main strategies and
three expression-eligibility sensitivities, freezes stochastic tie handling,
and emits a machine-auditable G2b record.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))

from tier_variant import ORIGINAL, annotate, as_bool, load_discovery_atlas  # noqa: E402

ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
FROZEN = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
RAW_CRISPR = ROOT / "01_data_raw" / "crispr" / "Table_S1.xlsx"
LOCK = ROOT / "13_postfreeze_interpretation" / "REVISION2_POSTREVIEW_ANALYSIS_LOCK.md"
OUT = ROOT / "revision2_2026-11-23"

PANEL_OUT = OUT / "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz"
SUMMARY_OUT = OUT / "REVISION2_WP5_COMPARATOR_SUMMARY.tsv"
SCALE_OUT = OUT / "REVISION2_WP5_CRISPR_SCALE_AUDIT.tsv"
REPORT_OUT = OUT / "REVISION2_WP5_COMPARATOR_REPORT.md"
AUDIT_OUT = OUT / "REVISION2_WP5_G2B_AUDIT.tsv"
MANIFEST_OUT = OUT / "REVISION2_WP5_FREEZE_MANIFEST.json"

BASE_SEED = 42
N_DRAWS = 1000
PANEL_SIZE = 56
TOL = 1e-12

ATLAS_COLUMNS = (
    "human_symbol",
    "gw_signed_confidence",
    "gw_pvalue",
    "gw_effect",
    "lcv_human_wt",
    "lcv_human_ankb",
    "lcv_mouse_untreated",
    "lcv_spatial_supported",
    "expression_n_studies",
    "expression_consensus_direction",
    "functional_tier",
    "ppi_AB_flag",
    "apms_high_confidence",
    "fitness_artifact_compatible",
    "stable_cross_cohort",
    "any_mechanism_convergent",
)

STRATEGIES = {
    "S1_FROZEN": {"code": 1, "label": "Frozen two-route panel", "draws": 1, "universe": "all"},
    "S2_CRISPR": {"code": 2, "label": "CRISPR-only", "draws": 1, "universe": "all"},
    "S3_LCV": {"code": 3, "label": "LCV-only", "draws": N_DRAWS, "universe": "all"},
    "S4_COUNT": {"code": 4, "label": "Unweighted discovery-domain count", "draws": N_DRAWS, "universe": "all"},
    "S5_CRISPR_EXPR": {"code": 5, "label": "CRISPR-only | expression eligible", "draws": 1, "universe": "expr"},
    "S6_LCV_EXPR": {"code": 6, "label": "LCV-only | expression eligible", "draws": N_DRAWS, "universe": "expr"},
    "S7_COUNT_EXPR": {"code": 7, "label": "Unweighted count | expression eligible", "draws": N_DRAWS, "universe": "expr"},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite_unique(series: pd.Series, name: str, gene: str) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna().to_numpy(float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return math.nan
    if float(np.max(vals) - np.min(vals)) > TOL:
        raise ValueError(f"conflicting {name} values for {gene}: {sorted(set(vals))}")
    return float(vals[0])


def collapse_atlas(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = annotate(raw, ORIGINAL).copy()
    frame["gene"] = frame["human_symbol"].fillna("").astype(str).str.strip().str.upper()
    frame = frame.loc[frame.gene.ne("")].copy()
    sizes = frame.groupby("gene", sort=True).size()
    rows: list[dict[str, object]] = []
    direction_conflicts = 0
    for gene, group in frame.groupby("gene", sort=True):
        directions = sorted(
            {
                str(x).strip().upper()
                for x in group["expression_consensus_direction"].dropna()
                if str(x).strip().upper() in {"UP", "DOWN"}
            }
        )
        if len(directions) == 1:
            direction = directions[0]
        elif len(directions) == 0:
            direction = "NOT_EVALUABLE"
        else:
            direction = "NOT_EVALUABLE"
            direction_conflicts += 1

        human = any(as_bool(x) for x in group["lcv_human_wt"]) or any(
            as_bool(x) for x in group["lcv_human_ankb"]
        )
        mouse = any(as_bool(x) for x in group["lcv_mouse_untreated"])
        lcv_wt = any(as_bool(x) for x in group["lcv_human_wt"])
        lcv_ankb = any(as_bool(x) for x in group["lcv_human_ankb"])
        direct = any(as_bool(x) for x in group["tv_direct_evidence"])
        functional_clean = any(as_bool(x) for x in group["tv_functional_clean"])
        response = any(as_bool(x) for x in group["stable_cross_cohort"])
        mechanism = any(as_bool(x) for x in group["any_mechanism_convergent"])
        spatial = any(as_bool(x) for x in group["lcv_spatial_supported"])
        rows.append(
            {
                "gene": gene,
                "gw_signed_confidence": finite_unique(group["gw_signed_confidence"], "gw_signed_confidence", gene),
                "gw_pvalue": finite_unique(group["gw_pvalue"], "gw_pvalue", gene),
                "gw_effect": finite_unique(group["gw_effect"], "gw_effect", gene),
                "lcv_human_wt": lcv_wt,
                "lcv_human_ankb": lcv_ankb,
                "lcv_mouse_untreated": mouse,
                "lcv_study_count": int(human) + int(mouse),
                "lcv_condition_count": int(lcv_wt) + int(lcv_ankb) + int(mouse),
                "direct": direct,
                "functional_clean": functional_clean,
                "response": response,
                "mechanism": mechanism,
                "spatial": spatial,
                "discovery_domain_count": sum((direct, functional_clean, response, mechanism, spatial)),
                "expression_n_studies": float(pd.to_numeric(group["expression_n_studies"], errors="coerce").max()),
                "expected_direction": direction,
            }
        )
    out = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    stats = {
        "atlas_rows_nonempty_symbol": int(len(frame)),
        "unique_genes": int(len(out)),
        "duplicated_symbols": int((sizes > 1).sum()),
        "duplicated_rows": int(sizes[sizes > 1].sum()),
        "direction_conflicts_set_not_evaluable": int(direction_conflicts),
    }
    return out, stats


def crispr_scale_audit() -> tuple[pd.DataFrame, str]:
    screen = pd.read_excel(RAW_CRISPR, sheet_name="Genome-wide screen")
    signed = pd.to_numeric(screen["Signed Confidence Score"], errors="coerce")
    unsigned = pd.to_numeric(screen["Confidence Score"], errors="coerce")
    valid = signed.notna() & unsigned.notna() & np.isfinite(signed) & np.isfinite(unsigned)
    diff = (signed[valid].abs() - unsigned[valid]).abs()
    identity_pass = bool(len(diff) > 0 and float(diff.max()) <= TOL)
    key = "abs(gw_signed_confidence)" if identity_pass else "-log10(gw_pvalue)"
    pos = signed[signed > 0].abs()
    neg = signed[signed < 0].abs()
    metrics = [
        ("source_rows", len(screen), "count"),
        ("finite_identity_pairs", int(valid.sum()), "count"),
        ("max_abs_identity_difference", float(diff.max()) if len(diff) else math.nan, "absolute difference"),
        ("identity_tolerance", TOL, "absolute tolerance"),
        ("identity_pass", identity_pass, "boolean"),
        ("positive_score_n", int((signed > 0).sum()), "count"),
        ("negative_score_n", int((signed < 0).sum()), "count"),
        ("zero_score_n", int((signed == 0).sum()), "count"),
        ("positive_abs_median", float(pos.median()), "score magnitude"),
        ("negative_abs_median", float(neg.median()), "score magnitude"),
        ("positive_abs_p90", float(pos.quantile(0.9)), "score magnitude"),
        ("negative_abs_p90", float(neg.quantile(0.9)), "score magnitude"),
        ("selected_key", key, "locked decision"),
    ]
    return pd.DataFrame(metrics, columns=["metric", "value", "unit_or_role"]), key


def deterministic_rank(frame: pd.DataFrame, key: pd.Series) -> pd.DataFrame:
    ranked = frame.assign(_primary=pd.to_numeric(key, errors="coerce").fillna(-np.inf))
    return ranked.sort_values(["_primary", "gene"], ascending=[False, True], kind="mergesort").head(PANEL_SIZE)


def randomized_rank(frame: pd.DataFrame, keys: list[str], strategy_code: int, draw: int) -> pd.DataFrame:
    ranked = frame.copy()
    for key in keys:
        ranked[key] = pd.to_numeric(ranked[key], errors="coerce").fillna(-np.inf)
    ranked = ranked.sort_values(keys, ascending=[False] * len(keys), kind="mergesort")
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, strategy_code, draw]))
    blocks = []
    grouper: str | list[str] = keys[0] if len(keys) == 1 else keys
    for _, block in ranked.groupby(grouper, sort=False, dropna=False):
        blocks.append(block.iloc[rng.permutation(len(block))])
    return pd.concat(blocks, ignore_index=True).head(PANEL_SIZE)


def panel_rows(strategy_id: str, label: str, draw: int, selected: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for rank, item in enumerate(selected.itertuples(index=False), start=1):
        rows.append(
            {
                "strategy_id": strategy_id,
                "strategy_label": label,
                "draw": draw,
                "rank": rank,
                "gene": item.gene,
                "expected_direction": item.expected_direction,
                "direction_status": "EVALUABLE" if item.expected_direction in {"UP", "DOWN"} else "NOT_EVALUABLE",
                "expression_n_studies": item.expression_n_studies,
                "gw_abs_signed_confidence": abs(item.gw_signed_confidence) if np.isfinite(item.gw_signed_confidence) else math.nan,
                "lcv_study_count": item.lcv_study_count,
                "lcv_condition_count": item.lcv_condition_count,
                "discovery_domain_count": item.discovery_domain_count,
            }
        )
    return rows


def construct_panels(universe: pd.DataFrame, selected_key: str, frozen: pd.DataFrame) -> pd.DataFrame:
    all_rows: list[dict[str, object]] = []
    expr_universe = universe.loc[universe.expression_n_studies.fillna(0).ge(1)].copy()
    universe_map = universe.set_index("gene")

    frozen_genes = frozen["human_symbol"].astype(str).str.strip().str.upper().tolist()
    missing = sorted(set(frozen_genes) - set(universe_map.index))
    if missing:
        raise ValueError(f"frozen genes missing from discovery universe: {missing}")
    frozen_selected = universe_map.loc[frozen_genes].reset_index()
    frozen_direction = (
        frozen.assign(gene=frozen["human_symbol"].astype(str).str.strip().str.upper())
        .set_index("gene")["expected_wt_vs_uninfected"]
        .astype(str)
        .str.upper()
        .to_dict()
    )
    frozen_selected["expected_direction"] = frozen_selected["gene"].map(frozen_direction)
    all_rows.extend(panel_rows("S1_FROZEN", STRATEGIES["S1_FROZEN"]["label"], 0, frozen_selected))

    for strategy_id in ("S2_CRISPR", "S5_CRISPR_EXPR"):
        pool = universe if STRATEGIES[strategy_id]["universe"] == "all" else expr_universe
        if selected_key == "abs(gw_signed_confidence)":
            key = pool["gw_signed_confidence"].abs()
        else:
            p = pd.to_numeric(pool["gw_pvalue"], errors="coerce")
            key = pd.Series(np.where((p > 0) & np.isfinite(p), -np.log10(p), np.nan), index=pool.index)
        selected = deterministic_rank(pool, key)
        all_rows.extend(panel_rows(strategy_id, STRATEGIES[strategy_id]["label"], 0, selected))

    for strategy_id in ("S3_LCV", "S4_COUNT", "S6_LCV_EXPR", "S7_COUNT_EXPR"):
        spec = STRATEGIES[strategy_id]
        pool = universe if spec["universe"] == "all" else expr_universe
        keys = ["lcv_study_count", "lcv_condition_count"] if "LCV" in strategy_id else ["discovery_domain_count"]
        for draw in range(1, N_DRAWS + 1):
            selected = randomized_rank(pool, keys, int(spec["code"]), draw)
            all_rows.extend(panel_rows(strategy_id, str(spec["label"]), draw, selected))
    return pd.DataFrame(all_rows)


def summarize_panels(panels: pd.DataFrame, frozen_genes: set[str], universe_sizes: dict[str, int]) -> pd.DataFrame:
    rows = []
    for strategy_id, strategy in panels.groupby("strategy_id", sort=False):
        per_draw = []
        hashes = set()
        for draw, panel in strategy.groupby("draw", sort=True):
            genes = set(panel.gene)
            hashes.add(hashlib.sha256("\n".join(sorted(genes)).encode()).hexdigest())
            per_draw.append(
                {
                    "draw": draw,
                    "size": len(panel),
                    "evaluable_n": int(panel.direction_status.eq("EVALUABLE").sum()),
                    "overlap_frozen": len(genes & frozen_genes),
                }
            )
        d = pd.DataFrame(per_draw)
        universe_label = str(STRATEGIES[strategy_id]["universe"])
        rows.append(
            {
                "strategy_id": strategy_id,
                "strategy_label": STRATEGIES[strategy_id]["label"],
                "universe": universe_label,
                "universe_n": universe_sizes[universe_label],
                "draws": len(d),
                "unique_panels": len(hashes),
                "panel_size_min": int(d["size"].min()),
                "panel_size_max": int(d["size"].max()),
                "evaluable_n_median": float(d["evaluable_n"].median()),
                "evaluable_n_p05": float(d["evaluable_n"].quantile(0.05)),
                "evaluable_n_p95": float(d["evaluable_n"].quantile(0.95)),
                "overlap_frozen_median": float(d["overlap_frozen"].median()),
                "overlap_frozen_p05": float(d["overlap_frozen"].quantile(0.05)),
                "overlap_frozen_p95": float(d["overlap_frozen"].quantile(0.95)),
            }
        )
    return pd.DataFrame(rows)


def add_check(checks: list[dict[str, object]], check: str, passed: bool, detail: str) -> None:
    checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = load_discovery_atlas(str(ATLAS), columns=ATLAS_COLUMNS)
    universe, collapse_stats = collapse_atlas(raw)
    scale, selected_key = crispr_scale_audit()
    frozen = pd.read_csv(FROZEN, sep="\t", low_memory=False)
    panels = construct_panels(universe, selected_key, frozen)
    frozen_genes = set(frozen.human_symbol.astype(str).str.strip().str.upper())
    universe_sizes = {
        "all": int(len(universe)),
        "expr": int(universe.expression_n_studies.fillna(0).ge(1).sum()),
    }
    summary = summarize_panels(panels, frozen_genes, universe_sizes)

    scale.to_csv(SCALE_OUT, sep="\t", index=False)
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    with gzip.GzipFile(filename=str(PANEL_OUT), mode="wb", mtime=0) as handle:
        handle.write(panels.to_csv(sep="\t", index=False).encode("utf-8"))

    checks: list[dict[str, object]] = []
    group_sizes = panels.groupby(["strategy_id", "draw"]).size()
    unique_sizes = panels.groupby(["strategy_id", "draw"]).gene.nunique()
    actual_draws = panels.groupby("strategy_id").draw.nunique().to_dict()
    expected_draws = {k: int(v["draws"]) for k, v in STRATEGIES.items()}
    add_check(checks, "G2b01_DISCOVERY_ONLY_COLUMNS", not any("holdout" in c.lower() or "final_tier" in c.lower() for c in raw.columns), f"loaded {len(raw.columns)} declared discovery columns")
    add_check(checks, "G2b02_GENE_UNIVERSE_UNIQUE", universe.gene.is_unique, f"{len(universe)} unique upper-case symbols")
    add_check(checks, "G2b03_DUPLICATE_POLICY", collapse_stats["duplicated_symbols"] == 171 and collapse_stats["direction_conflicts_set_not_evaluable"] == 65, json.dumps(collapse_stats, sort_keys=True))
    add_check(checks, "G2b04_CRISPR_IDENTITY", str(scale.loc[scale.metric.eq("identity_pass"), "value"].iloc[0]).lower() == "true", f"selected_key={selected_key}")
    add_check(checks, "G2b05_ALL_PANELS_SIZE_56", bool(group_sizes.eq(PANEL_SIZE).all()), f"min={group_sizes.min()}, max={group_sizes.max()}")
    add_check(checks, "G2b06_WITHIN_PANEL_UNIQUE", bool(unique_sizes.eq(PANEL_SIZE).all()), f"min={unique_sizes.min()}, max={unique_sizes.max()}")
    add_check(checks, "G2b07_DRAW_COUNTS", actual_draws == expected_draws, f"actual={actual_draws}")
    allowed_directions = set(panels.expected_direction.unique()) <= {"UP", "DOWN", "NOT_EVALUABLE", "NOT_PREDEFINED"}
    add_check(checks, "G2b08_DIRECTION_VOCABULARY", allowed_directions, f"values={sorted(panels.expected_direction.unique())}")
    frozen_panel = panels.loc[panels.strategy_id.eq("S1_FROZEN")]
    add_check(checks, "G2b09_FROZEN_PANEL_EXACT", set(frozen_panel.gene) == frozen_genes, f"observed={len(set(frozen_panel.gene))}, expected={len(frozen_genes)}")
    frozen_expected = frozen.set_index(frozen.human_symbol.astype(str).str.upper()).expected_wt_vs_uninfected.astype(str).str.upper().to_dict()
    direction_match = all(row.expected_direction == frozen_expected[row.gene] for row in frozen_panel.itertuples())
    add_check(checks, "G2b10_FROZEN_DIRECTIONS_EXACT", direction_match, f"matched={sum(row.expected_direction == frozen_expected[row.gene] for row in frozen_panel.itertuples())}/56")
    expr_genes = set(universe.loc[universe.expression_n_studies.fillna(0).ge(1), "gene"])
    expr_rows = panels.loc[panels.strategy_id.isin({"S5_CRISPR_EXPR", "S6_LCV_EXPR", "S7_COUNT_EXPR"})]
    add_check(checks, "G2b11_EXPR_UNIVERSE_ENFORCED", set(expr_rows.gene) <= expr_genes, f"eligible_universe={len(expr_genes)}")
    add_check(checks, "G2b12_BASE_SEED", BASE_SEED == 42, f"base_seed={BASE_SEED}")
    add_check(checks, "G2b13_NO_MOUSE_LFQ_FIELD", not any("lfq" in c.lower() or "abundance" in c.lower() for c in raw.columns), "mouse LFQ/abundance not loaded")
    reproducibility_cases = [
        ("S3_LCV", 1),
        ("S4_COUNT", N_DRAWS),
        ("S6_LCV_EXPR", 1),
        ("S7_COUNT_EXPR", N_DRAWS),
    ]
    reproducibility_matches = []
    for strategy_id, draw in reproducibility_cases:
        spec = STRATEGIES[strategy_id]
        pool = universe if spec["universe"] == "all" else universe.loc[universe.expression_n_studies.fillna(0).ge(1)].copy()
        keys = ["lcv_study_count", "lcv_condition_count"] if "LCV" in strategy_id else ["discovery_domain_count"]
        regenerated = randomized_rank(pool, keys, int(spec["code"]), draw).gene.tolist()
        archived = panels.loc[
            panels.strategy_id.eq(strategy_id) & panels.draw.eq(draw)
        ].sort_values("rank").gene.tolist()
        reproducibility_matches.append(regenerated == archived)
    add_check(
        checks,
        "G2b14_STOCHASTIC_REPRODUCIBILITY",
        all(reproducibility_matches),
        f"representative_draws={len(reproducibility_cases)}, matches={sum(reproducibility_matches)}",
    )
    all_pass = all(c["status"] == "PASS" for c in checks)
    audit = pd.DataFrame(checks)
    audit.to_csv(AUDIT_OUT, sep="\t", index=False)

    report = [
        "# Revision 2 WP5 comparator-panel freeze",
        "",
        f"- Gate G2b: **{'PASS' if all_pass else 'FAIL'}** ({sum(c['status'] == 'PASS' for c in checks)}/{len(checks)} checks)",
        f"- Holdout accessed: **false**",
        f"- Common discovery universe: {universe_sizes['all']:,} unique genes",
        f"- Expression-eligible universe: {universe_sizes['expr']:,} genes",
        f"- CRISPR ranking key: `{selected_key}`",
        f"- CRISPR source identity: |Signed Confidence Score| = Confidence Score for {int(scale.loc[scale.metric.eq('finite_identity_pairs'), 'value'].iloc[0]):,} finite rows",
        f"- Atlas case-collapse audit: {collapse_stats['duplicated_symbols']} duplicated symbols; {collapse_stats['direction_conflicts_set_not_evaluable']} opposing directions set to NOT_EVALUABLE",
        f"- Random tie draws: {N_DRAWS:,}; base seed: {BASE_SEED}",
        "",
        "## Frozen strategies",
        "",
        summary.to_markdown(index=False),
        "",
        "## Interpretation boundary",
        "",
        "These are discovery-only comparator panels. No holdout statistic has been calculated. Stochastic panels are frozen as complete gene lists for later one-time WP6 evaluation; they may not be regenerated or filtered after holdout access.",
    ]
    REPORT_OUT.write_text("\n".join(report) + "\n", encoding="utf-8")

    outputs = [PANEL_OUT, SUMMARY_OUT, SCALE_OUT, REPORT_OUT, AUDIT_OUT]
    manifest = {
        "stage": "WP5_G2b",
        "gate": "PASS" if all_pass else "FAIL",
        "holdout_accessed": False,
        "base_seed": BASE_SEED,
        "n_stochastic_draws": N_DRAWS,
        "panel_size": PANEL_SIZE,
        "selected_crispr_key": selected_key,
        "universe_sizes": universe_sizes,
        "collapse_stats": collapse_stats,
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (ATLAS, FROZEN, RAW_CRISPR, LOCK)},
        "script": {str(Path(__file__).relative_to(ROOT)): sha256(Path(__file__))},
        "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
        "strategy_draws": actual_draws,
        "checks_passed": sum(c["status"] == "PASS" for c in checks),
        "checks_total": len(checks),
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"gate": manifest["gate"], "checks": f"{manifest['checks_passed']}/{manifest['checks_total']}", "panels_rows": len(panels), "selected_key": selected_key, "outputs": [str(p) for p in outputs + [MANIFEST_OUT]]}, indent=2))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
