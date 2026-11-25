"""Execute the frozen WP6/G3 holdout comparator evaluation exactly once.

``--preflight`` validates immutable inputs without calculating holdout panel
statistics. ``--execute`` performs the complete formal evaluation and refuses to
overwrite a completed run.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import math
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revision2_2026-11-23"
PANELS = OUT / "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz"
HOLDOUT = ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz"
PRIOR_SET_TESTS = ROOT / "09_holdout_validation" / "holdout_frozen_set_tests.tsv"
STAGE7B = ROOT / "scripts" / "stage7b_holdout_validation.py"
RANDOM_NULL = ROOT / "13_postfreeze_interpretation" / "f6b_directional_null_100k.npy"
GAF = ROOT / "01_data_raw" / "annotations" / "goa_human.gaf.gz"
TRAFFICKING_SCRIPT = ROOT / "scripts" / "trafficking_universe.py"
EXEC_SPEC = OUT / "REVISION2_WP6_EXECUTION_SPEC.md"

PANEL_RESULTS = OUT / "REVISION2_WP6_PANEL_RESULTS.tsv.gz"
STRATEGY_SUMMARY = OUT / "REVISION2_WP6_STRATEGY_SUMMARY.tsv"
MAIN_CONTRASTS = OUT / "REVISION2_WP6_MAIN_CONTRASTS.tsv"
EQUAL_N = OUT / "REVISION2_WP6_EQUAL_N_SENSITIVITY.tsv"
MODULE_BOOT = OUT / "REVISION2_WP6_MODULE_BOOTSTRAP.tsv.gz"
G3_AUDIT = OUT / "REVISION2_WP6_G3_AUDIT.tsv"
STATS = OUT / "REVISION2_WP6_STATS.json"
REPORT = OUT / "REVISION2_WP6_REPORT.md"
MANIFEST = OUT / "REVISION2_WP6_MANIFEST.json"
RUN_MARKER = OUT / "REVISION2_WP6_RUN_MARKER.json"
FAILED_ATTEMPT = OUT / "REVISION2_WP6_ATTEMPT1_FAILED.json"

EXPECTED_HASHES = {
    PANELS: "9c070d7e8b3824bf8cfd5dbd3eb9f641d697ecb9c79106126b57fc70cc9d4d75",
    HOLDOUT: "1817f35940d9b5c700cfcadb90ea3a6833217b10482e2c3370c31ecdd826fbff",
    STAGE7B: "730f37cb6e6a117e65187f4987a94ddc2978107714921ec57b97e11905d27d4f",
    RANDOM_NULL: "1ceb832a742357a7731014d3622419116ffaa9120adbf86e993a4a027af34387",
    GAF: "db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30",
    EXEC_SPEC: "484d95d4a5922f6e6ee37d7a646b44f53cc495475366ada982a7ef681416553c",
}

N_PERM = 10_000
N_EQUAL = 2_000
N_BOOT = 2_000
EXPECTED_PANEL_ROWS = 224_168
EXPECTED_DRAWS = {
    "S1_FROZEN": 1, "S2_CRISPR": 1, "S3_LCV": 1000,
    "S4_COUNT": 1000, "S5_CRISPR_EXPR": 1,
    "S6_LCV_EXPR": 1000, "S7_COUNT_EXPR": 1000,
}
STRATEGY_CODES = {s: i + 1 for i, s in enumerate(EXPECTED_DRAWS)}
MAIN_COMPARATORS = ["S2_CRISPR", "S3_LCV", "S4_COUNT"]
TOL = 1e-12


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def signed_z(frame: pd.DataFrame) -> np.ndarray:
    """Verbatim numerical definition from the frozen stage-7B script."""
    p = np.clip(frame.PValue.to_numpy(float), 1e-300, 1 - 1e-15)
    return np.sign(frame.logFC.to_numpy(float)) * norm.isf(p / 2)


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def input_hashes() -> dict[str, str]:
    return {str(p.relative_to(ROOT)): sha256(p) for p in EXPECTED_HASHES}


def validate_hashes() -> None:
    for path, expected in EXPECTED_HASHES.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual.lower() != expected:
            raise RuntimeError(f"immutable input hash mismatch: {path}: {actual} != {expected}")


def load_panels() -> pd.DataFrame:
    p = pd.read_csv(PANELS, sep="\t")
    required = {"strategy_id", "strategy_label", "draw", "rank", "gene", "expected_direction"}
    if not required.issubset(p.columns):
        raise RuntimeError(f"missing panel columns: {sorted(required - set(p.columns))}")
    if len(p) != EXPECTED_PANEL_ROWS:
        raise RuntimeError(f"panel rows {len(p)} != {EXPECTED_PANEL_ROWS}")
    if set(p.strategy_id) != set(EXPECTED_DRAWS):
        raise RuntimeError("strategy identifiers differ from frozen specification")
    for sid, n_draws in EXPECTED_DRAWS.items():
        d = p.loc[p.strategy_id.eq(sid)]
        if d.draw.nunique() != n_draws or not d.groupby("draw").size().eq(56).all():
            raise RuntimeError(f"invalid draw structure for {sid}")
        if d.groupby("draw").gene.nunique().ne(56).any():
            raise RuntimeError(f"duplicate genes within a {sid} draw")
    if not set(p.expected_direction.astype(str)).issubset({"UP", "DOWN", "NOT_EVALUABLE", "NOT_PREDEFINED"}):
        raise RuntimeError("unexpected direction label")
    return p


def preflight() -> dict[str, object]:
    validate_hashes()
    p = load_panels()
    # Header/schema checks only: no formal panel statistic is calculated here.
    h = pd.read_csv(HOLDOUT, sep="\t", nrows=5)
    expected_h = {"contrast_id", "human_symbol", "logFC", "PValue"}
    if not expected_h.issubset(h.columns):
        raise RuntimeError(f"holdout schema is missing {sorted(expected_h - set(h.columns))}")
    null = np.load(RANDOM_NULL, mmap_mode="r")
    if null.shape != (100_000,) or not np.isfinite(null).all():
        raise RuntimeError("archived matched-random null is not the expected finite 100,000-vector")
    out = {
        "status": "PASS", "formal_holdout_statistics_calculated": False,
        "panel_rows": int(len(p)), "strategies": int(p.strategy_id.nunique()),
        "draw_panels": int(p.groupby(["strategy_id", "draw"]).ngroups),
        "input_hashes": input_hashes(),
    }
    print(json.dumps(out, indent=2))
    return out


def load_trafficking_sets() -> tuple[set[str], set[str], set[str]]:
    import importlib.util
    spec = importlib.util.spec_from_file_location("wp6_trafficking", TRAFFICKING_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not import trafficking universe")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return set(mod.RAB_GTPASE), set(mod.RAB_ASSOCIATED), set(mod.ALL_GENES)


def build_modules() -> dict[str, str]:
    """Reproduce the ST15c rarest-GO-BP assignment convention."""
    gene_terms: dict[str, set[str]] = collections.defaultdict(set)
    with gzip.open(GAF, "rt", errors="ignore") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[8] == "P":
                gene_terms[c[2].upper()].add(c[4])
    freq: collections.Counter[str] = collections.Counter()
    for terms in gene_terms.values():
        freq.update(terms)
    n_tot = sum(freq.values())

    def ic(term: str) -> float:
        return -math.log2(freq[term] / n_tot)

    return {gene: max(terms, key=ic) for gene, terms in gene_terms.items()}


def panel_permutation_null(universe: np.ndarray, n: int) -> np.ndarray:
    """Frozen empirical_greater loop with an order-independent fixed generator."""
    rng = np.random.default_rng(np.random.SeedSequence([42, 600, n]))
    null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        sample = rng.choice(universe, size=n, replace=False)
        null[i] = np.mean(sample)
    return null


def module_descriptors(modules: list[str]) -> dict[str, float | int]:
    sizes = pd.Series(modules).value_counts().to_numpy(int)
    return {
        "go_module_n": int(len(sizes)),
        "go_module_size_median": float(np.median(sizes)),
        "go_module_size_max": int(np.max(sizes)),
        "go_singleton_module_fraction": float(np.mean(sizes == 1)),
    }


def evaluate_panels(
    panels: pd.DataFrame, zmap: dict[str, float], universe: np.ndarray,
    module_map: dict[str, str], rab: set[str], rab_assoc: set[str], traffic: set[str],
) -> tuple[pd.DataFrame, dict[tuple[str, int], tuple[np.ndarray, np.ndarray, list[str]]]]:
    frozen_genes = set(panels.loc[panels.strategy_id.eq("S1_FROZEN"), "gene"])
    null_cache: dict[int, np.ndarray] = {}
    value_cache: dict[tuple[str, int], tuple[np.ndarray, np.ndarray, list[str]]] = {}
    rows: list[dict[str, object]] = []
    for (sid, draw), d in panels.groupby(["strategy_id", "draw"], sort=False):
        directions = d.expected_direction.map({"UP": 1.0, "DOWN": -1.0})
        present = d.gene.map(zmap).notna()
        valid = directions.notna() & present
        genes = d.loc[valid, "gene"].astype(str).tolist()
        signs = directions.loc[valid].to_numpy(float)
        raw_z = d.loc[valid, "gene"].map(zmap).to_numpy(float)
        oriented = signs * raw_z
        n = len(oriented)
        if n == 0:
            raise RuntimeError(f"no evaluable genes: {sid}/{draw}")
        if n not in null_cache:
            null_cache[n] = panel_permutation_null(universe, n)
        statistic = float(np.mean(oriented))
        perm_p = float((1 + np.sum(null_cache[n] >= statistic)) / (N_PERM + 1))
        modules = [module_map.get(g, f"NOMOD::{g}") for g in genes]
        desc = module_descriptors(modules)
        all_genes = set(d.gene.astype(str))
        value_cache[(sid, int(draw))] = (oriented, np.array(genes, dtype=object), modules)
        rows.append({
            "strategy_id": sid, "strategy_label": str(d.strategy_label.iloc[0]),
            "draw": int(draw), "selected_n": int(len(d)), "evaluable_n": n,
            "not_evaluable_direction_n": int(directions.isna().sum()),
            "holdout_missing_n": int((directions.notna() & ~present).sum()),
            "statistic": statistic, "permutation_p": perm_p,
            "direction_concordance": float(np.mean(oriented > 0)),
            "overlap_frozen_n": int(len(all_genes & frozen_genes)),
            "rab_gtpase_n": int(len(all_genes & rab)),
            "rab_associated_n": int(len(all_genes & rab_assoc)),
            "trafficking_n": int(len(all_genes & traffic)), **desc,
        })
    return pd.DataFrame(rows), value_cache


def q025(x: pd.Series | np.ndarray) -> float:
    return float(np.quantile(np.asarray(x, float), 0.025))


def q975(x: pd.Series | np.ndarray) -> float:
    return float(np.quantile(np.asarray(x, float), 0.975))


def summarize_strategies(panel_results: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "evaluable_n", "statistic", "permutation_p", "direction_concordance",
        "overlap_frozen_n", "rab_gtpase_n", "rab_associated_n", "trafficking_n",
        "go_module_n", "go_module_size_median", "go_module_size_max",
        "go_singleton_module_fraction",
    ]
    rows = []
    for sid, d in panel_results.groupby("strategy_id", sort=False):
        row: dict[str, object] = {
            "strategy_id": sid, "strategy_label": d.strategy_label.iloc[0],
            "draws": len(d), "selected_n": 56,
            "strategy_estimand": "single frozen panel" if len(d) == 1 else "median of 1,000 frozen tie-resolution draws",
        }
        for field in fields:
            row[f"{field}_median"] = float(d[field].median())
            row[f"{field}_p025"] = q025(d[field])
            row[f"{field}_p975"] = q975(d[field])
            row[f"{field}_min"] = float(d[field].min())
            row[f"{field}_max"] = float(d[field].max())
        rows.append(row)
    return pd.DataFrame(rows)


def equal_n_sensitivity(
    panel_results: pd.DataFrame,
    cache: dict[tuple[str, int], tuple[np.ndarray, np.ndarray, list[str]]],
) -> tuple[pd.DataFrame, int]:
    target_n = int(panel_results.evaluable_n.min())
    rows = []
    for sid in EXPECTED_DRAWS:
        keys = [k for k in cache if k[0] == sid]
        rng = np.random.default_rng(np.random.SeedSequence([42, 610, STRATEGY_CODES[sid]]))
        vals = np.empty(N_EQUAL, dtype=float)
        for i in range(N_EQUAL):
            key = keys[int(rng.integers(len(keys)))] if len(keys) > 1 else keys[0]
            oriented = cache[key][0]
            idx = rng.choice(len(oriented), size=target_n, replace=False)
            vals[i] = float(np.mean(oriented[idx]))
        rows.append({
            "strategy_id": sid, "target_n": target_n, "repetitions": N_EQUAL,
            "statistic_median": float(np.median(vals)),
            "statistic_p025": q025(vals), "statistic_p975": q975(vals),
            "sampling_rule": "draw-uniform then gene subsample" if len(keys) > 1 else "gene subsample",
            "backfill": False,
        })
    return pd.DataFrame(rows), target_n


def weighted_module_mean(values: np.ndarray, modules: list[str], multiplicity: collections.Counter[str]) -> float:
    weights = np.array([multiplicity[m] for m in modules], dtype=float)
    if weights.sum() == 0:
        return math.nan
    return float(np.average(values, weights=weights))


def module_bootstrap(
    cache: dict[tuple[str, int], tuple[np.ndarray, np.ndarray, list[str]]],
) -> tuple[pd.DataFrame, dict[str, dict[str, float | int]]]:
    frozen = cache[("S1_FROZEN", 0)]
    all_rows = []
    summaries: dict[str, dict[str, float | int]] = {}
    for comparator in MAIN_COMPARATORS:
        keys = [k for k in cache if k[0] == comparator]
        rng = np.random.default_rng(np.random.SeedSequence([42, 620, STRATEGY_CODES[comparator]]))
        deltas = []
        skipped = 0
        union_counts = []
        for rep in range(N_BOOT):
            key = keys[int(rng.integers(len(keys)))] if len(keys) > 1 else keys[0]
            comp = cache[key]
            union = sorted(set(frozen[2]) | set(comp[2]))
            union_counts.append(len(union))
            sampled = rng.choice(np.array(union, dtype=object), size=len(union), replace=True)
            mult = collections.Counter(sampled.tolist())
            tf = weighted_module_mean(frozen[0], frozen[2], mult)
            tc = weighted_module_mean(comp[0], comp[2], mult)
            delta = tf - tc
            if not np.isfinite(delta):
                skipped += 1
                continue
            deltas.append(delta)
            all_rows.append({"comparator_id": comparator, "replicate": rep + 1, "delta": delta, "union_module_n": len(union), "comparator_draw": int(key[1])})
        if len(deltas) != N_BOOT:
            raise RuntimeError(f"module bootstrap produced {skipped} non-finite replicates for {comparator}")
        summaries[comparator] = {
            "module_bootstrap_repetitions": len(deltas),
            "module_delta_median": float(np.median(deltas)),
            "module_delta_p025": q025(deltas), "module_delta_p975": q975(deltas),
            "union_module_n_min": int(min(union_counts)), "union_module_n_max": int(max(union_counts)),
        }
    return pd.DataFrame(all_rows), summaries


def main_contrasts(
    strategy_summary: pd.DataFrame, boot_summary: dict[str, dict[str, float | int]],
) -> pd.DataFrame:
    stats = strategy_summary.set_index("strategy_id").statistic_median
    frozen = float(stats["S1_FROZEN"])
    random_stats = np.load(RANDOM_NULL).astype(float)
    reference_delta = frozen - random_stats
    rows = []
    for comparator in MAIN_COMPARATORS:
        comp = float(stats[comparator])
        delta = frozen - comp
        p = float((1 + np.sum(reference_delta <= delta)) / (len(reference_delta) + 1))
        rows.append({
            "frozen_id": "S1_FROZEN", "comparator_id": comparator,
            "frozen_statistic": frozen, "comparator_statistic": comp,
            "observed_delta": delta, "random_calibration_p": p,
            **boot_summary[comparator],
            "inference_boundary": "matched-random calibration; not a direct superiority test",
        })
    out = pd.DataFrame(rows)
    out["random_calibration_p_holm"] = multipletests(out.random_calibration_p, method="holm")[1]
    return out


def audit_rows(
    panels: pd.DataFrame, panel_results: pd.DataFrame, summary: pd.DataFrame,
    equal_n: pd.DataFrame, boot: pd.DataFrame, contrasts: pd.DataFrame,
    zmap: dict[str, float], target_n: int,
) -> pd.DataFrame:
    checks: list[tuple[str, bool, str]] = []
    checks.append(("G3a_INPUT_HASHES", all(sha256(p).lower() == h for p, h in EXPECTED_HASHES.items()), "all immutable input hashes match"))
    checks.append(("G3b_PANEL_ARCHIVE", len(panels) == EXPECTED_PANEL_ROWS, f"rows={len(panels)}"))
    checks.append(("G3c_ALL_PANELS", len(panel_results) == 4003, f"evaluated={len(panel_results)}"))
    checks.append(("G3d_PANEL_SIZE", panel_results.selected_n.eq(56).all(), "all selected_n=56"))
    checks.append(("G3e_STRATEGIES", len(summary) == 7, f"strategies={len(summary)}"))
    frozen = panel_results.loc[panel_results.strategy_id.eq("S1_FROZEN")].iloc[0]
    checks.append(("G3f_FROZEN_EVALUABLE", int(frozen.evaluable_n) == 54, f"n={int(frozen.evaluable_n)}"))
    checks.append(("G3g_FROZEN_STATISTIC", abs(float(frozen.statistic) - 3.818344814213004) < TOL, f"T={float(frozen.statistic):.15g}"))
    # Explicit replay of the once-oriented frozen statistic.
    fp = panels.loc[panels.strategy_id.eq("S1_FROZEN")].copy()
    fp["d"] = fp.expected_direction.map({"UP": 1.0, "DOWN": -1.0})
    fp["z"] = fp.gene.map(zmap)
    replay = float((fp.d * fp.z).dropna().mean())
    checks.append(("G3h_ORIENTATION_ONCE", abs(replay - float(frozen.statistic)) < TOL, f"replay={replay:.15g}"))
    checks.append(("G3i_PERMUTATIONS", panel_results.permutation_p.between(1 / 10001, 1).all(), f"B={N_PERM} per panel"))
    checks.append(("G3j_EQUAL_N", len(equal_n) == 7 and equal_n["repetitions"].eq(N_EQUAL).all() and equal_n["backfill"].eq(False).all(), f"target_n={target_n}; B={N_EQUAL}; no backfill"))
    checks.append(("G3k_BOOTSTRAP", len(boot) == 3 * N_BOOT and boot.groupby("comparator_id").size().eq(N_BOOT).all(), f"rows={len(boot)}"))
    checks.append(("G3l_MAIN_CONTRASTS", list(contrasts.comparator_id) == MAIN_COMPARATORS, ",".join(contrasts.comparator_id)))
    checks.append(("G3m_HOLM", contrasts.random_calibration_p_holm.between(0, 1).all(), "Holm applied to exactly 3 calibration P values"))
    checks.append(("G3n_MODULES", panel_results.go_module_n.gt(0).all() and panel_results.go_singleton_module_fraction.between(0, 1).all(), "all panel module descriptors valid"))
    checks.append(("G3o_STOCHASTIC_ESTIMAND", all(int(summary.set_index("strategy_id").loc[s, "draws"]) == 1000 for s in ["S3_LCV", "S4_COUNT", "S6_LCV_EXPR", "S7_COUNT_EXPR"]), "four stochastic strategies use 1,000 draws"))
    return pd.DataFrame([{"check_id": cid, "status": "PASS" if ok else "FAIL", "detail": detail} for cid, ok, detail in checks])


def execute(resume_technical_failure: bool = False) -> None:
    if RUN_MARKER.exists():
        marker = json.loads(RUN_MARKER.read_text(encoding="utf-8"))
        formal_outputs = [PANEL_RESULTS, STRATEGY_SUMMARY, MAIN_CONTRASTS, EQUAL_N, MODULE_BOOT, G3_AUDIT, STATS, REPORT, MANIFEST]
        recoverable = (
            resume_technical_failure
            and marker.get("status") == "FAILED"
            and marker.get("error") == "'function' object has no attribute 'eq'"
            and not any(p.exists() for p in formal_outputs)
            and not FAILED_ATTEMPT.exists()
        )
        if not recoverable:
            raise RuntimeError(f"formal WP6 run marker already exists with status={marker.get('status')}; refusing overwrite")
        # Preserve the unmodified failure record before allowing a result-unseen,
        # output-free recovery from the audit-column name collision.
        write_json(FAILED_ATTEMPT, marker)
    validate_hashes()
    panels = load_panels()
    start = datetime.now(timezone.utc).isoformat()
    write_json(RUN_MARKER, {"status": "RUNNING", "started_utc": start, "input_hashes": input_hashes()})
    try:
        holdout = pd.read_csv(HOLDOUT, sep="\t")
        wt = holdout.loc[holdout.contrast_id.eq("WT_vs_UI")].copy()
        if len(wt) != 10846 or wt.human_symbol.nunique() != 10846:
            raise RuntimeError(f"unexpected WT_vs_UI shape: {wt.shape}, unique={wt.human_symbol.nunique()}")
        wt["z"] = signed_z(wt)
        if not np.isfinite(wt.z).all():
            raise RuntimeError("non-finite signed z in WT_vs_UI universe")
        zmap = dict(zip(wt.human_symbol.astype(str), wt.z.astype(float), strict=True))
        universe = wt.z.to_numpy(float)
        module_map = build_modules()
        rab, rab_assoc, traffic = load_trafficking_sets()
        panel_results, cache = evaluate_panels(panels, zmap, universe, module_map, rab, rab_assoc, traffic)
        summary = summarize_strategies(panel_results)
        equal_n, target_n = equal_n_sensitivity(panel_results, cache)
        boot, boot_summary = module_bootstrap(cache)
        contrasts = main_contrasts(summary, boot_summary)
        audit = audit_rows(panels, panel_results, summary, equal_n, boot, contrasts, zmap, target_n)
        if audit.status.ne("PASS").any():
            failed = audit.loc[audit.status.ne("PASS")].to_dict("records")
            raise RuntimeError(f"G3 audit failed before output release: {failed}")

        panel_results.to_csv(PANEL_RESULTS, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 9})
        summary.to_csv(STRATEGY_SUMMARY, sep="\t", index=False)
        contrasts.to_csv(MAIN_CONTRASTS, sep="\t", index=False)
        equal_n.to_csv(EQUAL_N, sep="\t", index=False)
        boot.to_csv(MODULE_BOOT, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 9})
        audit.to_csv(G3_AUDIT, sep="\t", index=False)

        sidx = summary.set_index("strategy_id")
        result_stats = {
            "gate": "G3", "status": "PASS", "checks_pass": int(audit.status.eq("PASS").sum()),
            "checks_fail": 0, "panel_draws_evaluated": int(len(panel_results)),
            "strategies": int(len(summary)), "selected_n_each": 56,
            "equal_n_target": target_n, "permutations_per_panel": N_PERM,
            "equal_n_repetitions": N_EQUAL, "module_bootstrap_repetitions_per_contrast": N_BOOT,
            "frozen_evaluable_n": int(sidx.loc["S1_FROZEN", "evaluable_n_median"]),
            "frozen_statistic": float(sidx.loc["S1_FROZEN", "statistic_median"]),
            "main_contrasts": contrasts.to_dict("records"),
            "claim_boundary": "cross-context transcriptional corroboration only; no functional-dependency validation or direct superiority test",
        }
        write_json(STATS, result_stats)
        report_lines = [
            "# WP6/G3 one-time holdout comparison", "",
            f"- G3: **PASS** ({len(audit)}/{len(audit)} checks).",
            f"- Evaluated all {len(panel_results):,} frozen panels across seven strategies; every panel selected 56 genes.",
            f"- Equal-n sensitivity target: {target_n}; no back-filling.",
            f"- Frozen panel: {int(sidx.loc['S1_FROZEN', 'evaluable_n_median'])} evaluable genes; statistic {float(sidx.loc['S1_FROZEN', 'statistic_median']):.6f}.",
            "", "## Strategy statistics", "",
            "Strategy-level values for stochastic comparators are medians across 1,000 frozen tie-resolution draws.", "",
            "| Strategy | Evaluable n (median) | Statistic (median) | 2.5%-97.5% across draws | Direction concordance |",
            "|---|---:|---:|---:|---:|",
        ]
        for _, r in summary.iterrows():
            report_lines.append(f"| {r.strategy_id} | {r.evaluable_n_median:.0f} | {r.statistic_median:.4f} | {r.statistic_p025:.4f} to {r.statistic_p975:.4f} | {r.direction_concordance_median:.3f} |")
        report_lines += ["", "## Frozen-minus-comparator contrasts", "",
            "Bootstrap intervals quantify sensitivity to GO-module composition. Calibration P values compare each comparator with the archived matched-random reference and are not direct tests of frozen-panel superiority.", "",
            "| Comparator | Observed delta | Module bootstrap 2.5%-97.5% | Random calibration P | Holm P |",
            "|---|---:|---:|---:|---:|",
        ]
        for _, r in contrasts.iterrows():
            report_lines.append(f"| {r.comparator_id} | {r.observed_delta:.4f} | {r.module_delta_p025:.4f} to {r.module_delta_p975:.4f} | {r.random_calibration_p:.6g} | {r.random_calibration_p_holm:.6g} |")
        report_lines += ["", "## Interpretation boundary", "",
            "This is the single post-review comparative use of the frozen WT-vs-uninfected holdout. It assesses directional transcriptional corroboration, not Legionella-specific functional dependency, direct binding, or therapeutic causality. No panel, rank, direction, endpoint, seed, or inferential rule was changed after result access.", ""]
        REPORT.write_text("\n".join(report_lines), encoding="utf-8")

        output_paths = [PANEL_RESULTS, STRATEGY_SUMMARY, MAIN_CONTRASTS, EQUAL_N, MODULE_BOOT, G3_AUDIT, STATS, REPORT]
        manifest = {
            "status": "COMPLETE", "started_utc": start,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "input_hashes": input_hashes(),
            "script": {str(Path(__file__).resolve().relative_to(ROOT)): sha256(Path(__file__).resolve())},
            "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in output_paths},
        }
        write_json(MANIFEST, manifest)
        manifest["outputs"][str(MANIFEST.relative_to(ROOT))] = sha256(MANIFEST)
        # Manifest does not self-hash internally; marker records its final external hash.
        marker = {
            "status": "COMPLETE", "started_utc": start,
            "completed_utc": manifest["completed_utc"],
            "g3": "PASS", "manifest_sha256": sha256(MANIFEST),
            "completed_formal_output_runs": 1,
            "technical_attempts": 2 if FAILED_ATTEMPT.exists() else 1,
            "attempt1_failed_before_output_release": FAILED_ATTEMPT.exists(),
        }
        write_json(RUN_MARKER, marker)
        print(json.dumps(result_stats, indent=2))
    except Exception as exc:
        write_json(RUN_MARKER, {
            "status": "FAILED", "started_utc": start,
            "failed_utc": datetime.now(timezone.utc).isoformat(),
            "error": str(exc), "traceback": traceback.format_exc(),
        })
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--execute", action="store_true")
    mode.add_argument("--resume-technical-failure", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
    else:
        execute(resume_technical_failure=args.resume_technical_failure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
