"""Construct outcome-blind entropy-balanced opportunity references for PD1C.

This Stage-1 executable intentionally has no holdout or WP6 endpoint path.  It
uses only frozen comparator panels and prespecified opportunity covariates.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))
from build_wp5_comparator_panels import ATLAS_COLUMNS, collapse_atlas, load_discovery_atlas  # noqa: E402

OUT = ROOT / "revision2_2026-11-23"
SPEC = OUT / "REVISION2_PD1C_STAGE1_SPEC.md"
SPEC_AMENDMENT = OUT / "REVISION2_PD1C_STAGE1_SPEC_AMENDMENT_01.md"
PANELS = OUT / "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz"
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
GAF = ROOT / "01_data_raw" / "annotations" / "goa_human.gaf.gz"
OBO = ROOT / "01_data_raw" / "annotations" / "go-basic.obo"
DEGREE = ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv"
MATRIX_DIR = ROOT / "05_expression" / "matrices"
MAP_DIR = ROOT / "05_expression" / "gene_maps"

WEIGHTS_OUT = OUT / "REVISION2_PD1C_STAGE1_WEIGHTS.tsv.gz"
BALANCE_OUT = OUT / "REVISION2_PD1C_STAGE1_BALANCE.tsv"
DIAGNOSTICS_OUT = OUT / "REVISION2_PD1C_STAGE1_DIAGNOSTICS.tsv"
AUDIT_OUT = OUT / "REVISION2_PD1C_STAGE1_AUDIT.tsv"
STATS_OUT = OUT / "REVISION2_PD1C_STAGE1_STATS.json"
REPORT_OUT = OUT / "REVISION2_PD1C_STAGE1_REPORT.md"
MANIFEST_OUT = OUT / "REVISION2_PD1C_STAGE1_MANIFEST.json"
RUN_MARKER = OUT / "REVISION2_PD1C_STAGE1_RUN_MARKER.json"

EXPECTED_HASHES = {
    SPEC: "632aef1404d375d8c5814066c5edaf093d6805fda946db56950bd2da7a3a62d0",
    SPEC_AMENDMENT: "a32ac478f1fc7043a49aee6e80a904d41e156cc42e56e612b435dd2210624a35",
    PANELS: "9c070d7e8b3824bf8cfd5dbd3eb9f641d697ecb9c79106126b57fc70cc9d4d75",
    ATLAS: "9f1ccdc73894da2e1e4ee251db296193df63cd9493b3e16383913e5a3dcb3da5",
    GAF: "db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30",
    OBO: "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc",
    DEGREE: "3ec0a0296adc970753ca8a6fdd6195fc7eb48ff74b0abca2d5857ebd9c337676",
    MATRIX_DIR / "GSE61535_log2_expression.tsv.gz": "dc4892602fe379fc5539b6ef54ae59b6e152296187143aa369ca8601b56a8bed",
    MATRIX_DIR / "GSE135803_fpkm.tsv.gz": "269ba7b5978335a98efae5853a94ebab5423f92b6f12b7eba719e81f3c787ebc",
    MATRIX_DIR / "GSE287755_raw_counts.tsv.gz": "f177240d44285daee83aada53cb8c23dea2366ea2c9f881a58adca069919f8cf",
    MAP_DIR / "GSE61535_gene_map.tsv": "26abf7d237eed33380488e12c11a1b980c3a437160ba9d6e839439ee0c4addeb",
    MAP_DIR / "GSE135803_gene_map.tsv": "9d07b1f12dddac0bd6efbb84b60d66c4350d904afe2465259c33a66616bd55fc",
    MAP_DIR / "GSE287755_gene_map.tsv": "aefe372ca6a66027f81f016e0f8efb3bdf7c922f9bbc34260271fd4961d2e276",
}

PANEL_SIZE = 56
MEMBRANE_ROOT = "GO:0016020"
IMBO_ROOT = "GO:0043231"
STRATEGIES = ["S1_FROZEN", "S2_CRISPR", "S3_LCV", "S4_COUNT",
              "S5_CRISPR_EXPR", "S6_LCV_EXPR", "S7_COUNT_EXPR"]
STRATEGY_UNIVERSE = {s: ("expr" if s.endswith("_EXPR") else "all") for s in STRATEGIES}
LINEAR = [
    "log1p_go_terms_n", "membrane_flag", "imbo_flag", "baseline_ui_percentile",
    "expression_n_studies", "log1p_cross_pathogen_degree", "log1p_cross_pathogen_pmids",
]
CONTINUOUS = [
    "log1p_go_terms_n", "baseline_ui_percentile", "expression_n_studies",
    "log1p_cross_pathogen_degree", "log1p_cross_pathogen_pmids",
]
PROHIBITED_PATH_TOKENS = ("holdout", "wp6")
PROHIBITED_OUTCOME_COLUMNS = {"raw_signed_z", "oriented_z", "pvalue", "logfc"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def validate_inputs() -> dict[str, str]:
    checked: dict[str, str] = {}
    for path, expected in EXPECTED_HASHES.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual.lower() != expected:
            raise RuntimeError(f"immutable input hash mismatch: {rel(path)}: {actual} != {expected}")
        checked[rel(path)] = actual
    return checked


def audit_no_outcome_input(universe: pd.DataFrame) -> tuple[bool, str]:
    path_hits = sorted(
        rel(path) for path in EXPECTED_HASHES
        if any(token in rel(path).lower() for token in PROHIBITED_PATH_TOKENS)
    )
    column_hits = sorted(PROHIBITED_OUTCOME_COLUMNS & {str(c).lower() for c in universe.columns})
    ok = not path_hits and not column_hits
    detail = f"prohibited_paths={path_hits}; prohibited_columns={column_hits}"
    return ok, detail


def load_obo(path: Path) -> dict[str, list[str]]:
    parents: dict[str, list[str]] = {}
    obsolete: set[str] = set()
    cur: str | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]" or line.startswith("["):
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:].strip(); parents.setdefault(cur, [])
            elif cur and line.startswith("is_a: GO:"):
                parents[cur].append(line[6:].split("!")[0].strip())
            elif cur and line.startswith("relationship: part_of GO:"):
                parents[cur].append(line[len("relationship: part_of "):].split("!")[0].strip())
            elif cur and line.startswith("is_obsolete: true"):
                obsolete.add(cur)
    for term in obsolete:
        parents.pop(term, None)
    return parents


def descendants(parents: dict[str, list[str]], root: str) -> set[str]:
    children: dict[str, list[str]] = {}
    for child, pars in parents.items():
        for parent in pars:
            children.setdefault(parent, []).append(child)
    out, stack = {root}, [root]
    while stack:
        for child in children.get(stack.pop(), []):
            if child not in out:
                out.add(child); stack.append(child)
    return out


def load_gaf(path: Path) -> dict[str, set[str]]:
    terms: dict[str, set[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or (fields[3] and "NOT" in fields[3].split("|")):
                continue
            gene, go_id = fields[2].strip().upper(), fields[4].strip()
            if gene and go_id.startswith("GO:"):
                terms.setdefault(gene, set()).add(go_id)
    return terms


def load_ui(mat: str, gsm: list[str], gmap: str, probe: bool) -> pd.Series:
    matrix = pd.read_csv(MATRIX_DIR / mat, sep="\t", compression="gzip")
    id_col = "feature_id" if "feature_id" in matrix.columns else "Geneid"
    cols = [c for c in gsm if c in matrix.columns]
    mapping = pd.read_csv(MAP_DIR / gmap, sep="\t")
    mapping = mapping.loc[mapping.gene_symbol.notna()].copy()
    if probe and "ambiguous_mapping" in mapping.columns:
        mapping = mapping.loc[~mapping.ambiguous_mapping.astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in mapping.columns and id_col != "feature_id":
        mapping = mapping.rename(columns={"feature_id": id_col})
    merged = matrix.merge(mapping, on=id_col, how="inner")
    sub = merged[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    out = sub.rank(pct=True)[cols].mean(axis=1)
    out.index = out.index.astype(str).str.upper()
    return out


def build_covariate_universe() -> tuple[pd.DataFrame, dict[str, int]]:
    raw = load_discovery_atlas(str(ATLAS), ATLAS_COLUMNS)
    universe, collapse_stats = collapse_atlas(raw)
    universe["gene"] = universe.gene.astype(str).str.upper()
    parents = load_obo(OBO)
    membrane, imbo = descendants(parents, MEMBRANE_ROOT), descendants(parents, IMBO_ROOT)
    terms = load_gaf(GAF)
    universe["go_terms_n"] = universe.gene.map(lambda g: len(terms.get(g, set()))).astype(float)
    universe["membrane_flag"] = universe.gene.map(lambda g: bool(terms.get(g, set()) & membrane)).astype(float)
    universe["imbo_flag"] = universe.gene.map(lambda g: bool(terms.get(g, set()) & imbo)).astype(float)
    ui = pd.concat([
        load_ui("GSE61535_log2_expression.tsv.gz", ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
        load_ui("GSE135803_fpkm.tsv.gz", ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
        load_ui("GSE287755_raw_counts.tsv.gz", ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
    ], axis=1).mean(axis=1, skipna=True)
    universe["baseline_ui_percentile"] = universe.gene.map(ui).fillna(0.0)
    degree = pd.read_csv(DEGREE, sep="\t", low_memory=False)
    degree["gene"] = degree.gene.astype(str).str.upper()
    degree = degree.groupby("gene", as_index=False)[["cross_pathogen_degree", "cross_pathogen_pmids"]].max()
    universe = universe.merge(degree, on="gene", how="left")
    for col in ["cross_pathogen_degree", "cross_pathogen_pmids", "expression_n_studies"]:
        universe[col] = pd.to_numeric(universe[col], errors="coerce").fillna(0)
    universe["log1p_go_terms_n"] = np.log1p(universe.go_terms_n)
    universe["log1p_cross_pathogen_degree"] = np.log1p(universe.cross_pathogen_degree)
    universe["log1p_cross_pathogen_pmids"] = np.log1p(universe.cross_pathogen_pmids)
    if universe.gene.duplicated().any() or not np.isfinite(universe[LINEAR].to_numpy(float)).all():
        raise RuntimeError("invalid unique/covariate universe")
    stats = dict(collapse_stats)
    stats.update({"membrane_descendants_including_root": len(membrane),
                  "imbo_descendants_including_root": len(imbo), "universe_n": len(universe)})
    return universe, stats


def validate_panels(panels: pd.DataFrame, universe: pd.DataFrame) -> None:
    required = {"strategy_id", "draw", "gene", "expected_direction"}
    if not required <= set(panels.columns):
        raise RuntimeError(f"missing panel columns: {sorted(required - set(panels.columns))}")
    known = set(universe.gene)
    for sid, group in panels.groupby("strategy_id"):
        if sid not in STRATEGIES or group.groupby("draw").size().ne(PANEL_SIZE).any():
            raise RuntimeError(f"invalid frozen panels for {sid}")
        if group.groupby("draw").gene.nunique().ne(PANEL_SIZE).any():
            raise RuntimeError(f"duplicate frozen panel genes for {sid}")
        missing = set(group.gene.astype(str).str.upper()) - known
        if missing:
            raise RuntimeError(f"unmapped panel genes for {sid}: {sorted(missing)[:10]}")


def weighted_mean_var(x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.sum(x * w[:, None], axis=0)
    var = np.sum(((x - mean) ** 2) * w[:, None], axis=0)
    return mean, var


def balance_strategy(sid: str, universe: pd.DataFrame, panels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    base = universe.loc[universe.expression_n_studies.ge(1)].copy() if STRATEGY_UNIVERSE[sid] == "expr" else universe.copy()
    base = base.loc[base.expected_direction.isin(["UP", "DOWN"])].copy().reset_index(drop=True)
    psub = panels.loc[panels.strategy_id.eq(sid)].copy()
    psub["gene"] = psub.gene.astype(str).str.upper()
    n_draws = int(psub.draw.nunique())
    counts = psub.groupby("gene").size()
    base["selection_slots"] = base.gene.map(counts).fillna(0).astype(float)
    query = base.loc[base.selection_slots.gt(0)].copy()
    controls = base.loc[base.selection_slots.eq(0)].copy()
    eligible_genes = set(base.gene)
    eligible_panel_slots = int(psub.gene.isin(eligible_genes).sum())
    if not math.isclose(float(query.selection_slots.sum()), eligible_panel_slots):
        raise RuntimeError(f"query-slot mismatch for {sid}")
    q_weight = query.selection_slots.to_numpy(float)
    q_weight /= q_weight.sum()

    control_raw = controls[LINEAR].to_numpy(float)
    query_raw = query[LINEAR].to_numpy(float)
    centre = control_raw.mean(axis=0)
    scale = control_raw.std(axis=0, ddof=0)
    scale[scale == 0] = 1.0
    control_z = (control_raw - centre) / scale
    query_z = (query_raw - centre) / scale
    continuous_idx = [LINEAR.index(c) for c in CONTINUOUS]
    control_sq0 = control_z[:, continuous_idx] ** 2
    query_sq0 = query_z[:, continuous_idx] ** 2
    sq_centre = control_sq0.mean(axis=0)
    sq_scale = control_sq0.std(axis=0, ddof=0)
    sq_scale[sq_scale == 0] = 1.0
    control_sq = (control_sq0 - sq_centre) / sq_scale
    query_sq = (query_sq0 - sq_centre) / sq_scale
    features = np.column_stack([control_z, control_sq])
    q_features = np.column_stack([query_z, query_sq])
    target = np.sum(q_features * q_weight[:, None], axis=0)

    def objective(lam: np.ndarray) -> tuple[float, np.ndarray]:
        scores = features @ lam
        value = float(logsumexp(scores) - math.log(len(scores)) - np.dot(target, lam))
        w = softmax(scores)
        grad = np.sum(features * w[:, None], axis=0) - target
        return value, grad

    result = minimize(lambda x: objective(x), np.zeros(features.shape[1]), jac=True,
                      method="L-BFGS-B", options={"maxiter": 3000, "ftol": 1e-15, "gtol": 1e-10, "maxls": 100})
    weights = softmax(features @ result.x)
    moment_error = np.sum(features * weights[:, None], axis=0) - target
    max_moment_error = float(np.max(np.abs(moment_error)))
    ess = float(1.0 / np.sum(weights ** 2))

    q_mean, q_var = weighted_mean_var(query_raw, q_weight)
    c_mean, c_var = weighted_mean_var(control_raw, weights)
    denom = np.sqrt((q_var + c_var) / 2)
    denom[denom == 0] = 1.0
    linear_smd = (q_mean - c_mean) / denom

    feature_names = LINEAR + [f"square_{c}" for c in CONTINUOUS]
    balance_rows = []
    weighted_feature_mean = np.sum(features * weights[:, None], axis=0)
    for i, name in enumerate(feature_names):
        kind = "linear" if i < len(LINEAR) else "squared_moment"
        difference = float(target[i] - weighted_feature_mean[i])
        balance_rows.append({
            "strategy_id": sid, "feature": name, "feature_class": kind,
            "target_standardized_mean": float(target[i]),
            "weighted_control_standardized_mean": float(weighted_feature_mean[i]),
            "standardized_mean_difference": float(linear_smd[i]) if kind == "linear" else difference,
            "absolute_difference": abs(float(linear_smd[i]) if kind == "linear" else difference),
        })
    weights_df = pd.DataFrame({
        "strategy_id": sid, "gene": controls.gene.to_numpy(),
        "expected_direction": controls.expected_direction.to_numpy(), "weight": weights,
    })
    diagnostic = {
        "strategy_id": sid, "universe_type": STRATEGY_UNIVERSE[sid], "panel_draws": n_draws,
        "archived_panel_slots": int(n_draws * PANEL_SIZE),
        "direction_eligible_panel_slots": eligible_panel_slots,
        "query_unique_genes": int(len(query)), "query_slots": int(query.selection_slots.sum()),
        "control_pool_n": int(len(controls)), "optimizer_success": bool(result.success),
        "optimizer_status": int(result.status), "optimizer_message": str(result.message),
        "optimizer_iterations": int(result.nit), "max_abs_moment_error": max_moment_error,
        "max_abs_linear_smd": float(np.max(np.abs(linear_smd))),
        "max_abs_squared_difference": float(np.max(np.abs(moment_error[len(LINEAR):]))),
        "ess": ess, "ess_fraction": ess / len(controls), "max_weight": float(np.max(weights)),
        "weight_sum": float(np.sum(weights)), "min_weight": float(np.min(weights)),
    }
    return weights_df, pd.DataFrame(balance_rows), diagnostic


def add_audit(rows: list[dict[str, str]], check: str, ok: bool, detail: str) -> None:
    rows.append({"check": check, "status": "PASS" if ok else "FAIL", "detail": detail})


def execute() -> None:
    outputs = [WEIGHTS_OUT, BALANCE_OUT, DIAGNOSTICS_OUT, AUDIT_OUT, STATS_OUT,
               REPORT_OUT, MANIFEST_OUT, RUN_MARKER]
    existing = [str(p) for p in outputs if p.exists()]
    if existing:
        raise RuntimeError(f"refusing to overwrite PD1C Stage-1 outputs: {existing}")
    checked = validate_inputs()
    universe, universe_stats = build_covariate_universe()
    panels = pd.read_csv(PANELS, sep="\t")
    panels["gene"] = panels.gene.astype(str).str.upper()
    validate_panels(panels, universe)
    all_weights, all_balance, diagnostics = [], [], []
    for sid in STRATEGIES:
        print(f"PD1C Stage 1 balancing {sid} ...", flush=True)
        w, b, d = balance_strategy(sid, universe, panels)
        all_weights.append(w); all_balance.append(b); diagnostics.append(d)
    weights = pd.concat(all_weights, ignore_index=True)
    balance = pd.concat(all_balance, ignore_index=True)
    diag = pd.DataFrame(diagnostics)
    weights.to_csv(WEIGHTS_OUT, sep="\t", index=False, compression="gzip")
    balance.to_csv(BALANCE_OUT, sep="\t", index=False)
    diag.to_csv(DIAGNOSTICS_OUT, sep="\t", index=False)

    audits: list[dict[str, str]] = []
    add_audit(audits, "C01_INPUT_HASHES", len(checked) == len(EXPECTED_HASHES), f"{len(checked)} files verified")
    add_audit(audits, "C02_PANEL_ROWS", len(panels) == 224168, f"rows={len(panels)}")
    no_outcome_ok, no_outcome_detail = audit_no_outcome_input(universe)
    add_audit(audits, "C03_NO_OUTCOME_INPUT", no_outcome_ok, no_outcome_detail)
    for row in diag.itertuples(index=False):
        sid = row.strategy_id
        add_audit(audits, f"C04_OPTIMIZER_{sid}", bool(row.optimizer_success) and row.max_abs_moment_error <= 1e-7,
                  f"success={row.optimizer_success}; max_error={row.max_abs_moment_error:.3g}")
        add_audit(audits, f"C05_LINEAR_BALANCE_{sid}", row.max_abs_linear_smd <= .10,
                  f"max_abs_smd={row.max_abs_linear_smd:.6g}")
        add_audit(audits, f"C06_SQUARED_BALANCE_{sid}", row.max_abs_squared_difference <= .10,
                  f"max_abs_difference={row.max_abs_squared_difference:.3g}")
        add_audit(audits, f"C07_ESS_{sid}", row.ess >= 560 and row.ess_fraction >= .02,
                  f"ess={row.ess:.3f}; fraction={row.ess_fraction:.6f}")
        add_audit(audits, f"C08_MAX_WEIGHT_{sid}", row.max_weight <= .01,
                  f"max_weight={row.max_weight:.8f}")
        sid_w = weights.loc[weights.strategy_id.eq(sid), "weight"].to_numpy(float)
        valid = np.isfinite(sid_w).all() and (sid_w >= 0).all() and abs(sid_w.sum() - 1) <= 1e-10
        add_audit(audits, f"C09_WEIGHT_VALIDITY_{sid}", bool(valid),
                  f"n={len(sid_w)}; sum={sid_w.sum():.12g}; min={sid_w.min():.3g}")
    audit = pd.DataFrame(audits)
    audit.to_csv(AUDIT_OUT, sep="\t", index=False)
    status = "PASS" if audit.status.eq("PASS").all() else "FAIL"
    stats = {
        "analysis": "PD1C Stage 1 outcome-blind entropy balancing", "status": status,
        "stage2_permitted": status == "PASS", "outcome_loaded": False,
        "covariates": LINEAR, "squared_terms": CONTINUOUS,
        "thresholds": {"max_abs_linear_smd": .10, "max_abs_squared_difference": .10,
                       "max_abs_moment_error": 1e-7, "min_ess": 560,
                       "min_ess_fraction": .02, "max_weight": .01},
        "universe_stats": universe_stats,
        "strategy_diagnostics": diagnostics,
        "interpretive_boundary": "Stage 1 contains no holdout outcome and supports no manuscript endpoint conclusion.",
    }
    write_json(STATS_OUT, stats)
    lines = [
        "# PD1C Stage 1 outcome-blind opportunity balancing", "",
        f"- Gate: **{status}**", "- Holdout outcome loaded: **No**",
        f"- Conditional Stage 2 permitted: **{'Yes' if status == 'PASS' else 'No'}**", "",
        "## Strategy diagnostics", "", diag.to_markdown(index=False), "",
        "## Executable audit", "", audit.to_markdown(index=False), "",
        "## Interpretation boundary", "",
        "This stage only constructs and audits covariate-balanced control weights. It contains no infection-transcription endpoint and provides no biological or comparative-performance result.", "",
    ]
    REPORT_OUT.write_text("\n".join(lines), encoding="utf-8")
    write_json(RUN_MARKER, {"status": status, "stage2_permitted": status == "PASS",
                            "completed_utc": datetime.now(timezone.utc).isoformat(),
                            "spec_sha256": sha256(SPEC), "outcome_loaded": False})
    output_hashes = {rel(p): sha256(p) for p in [WEIGHTS_OUT, BALANCE_OUT, DIAGNOSTICS_OUT,
                                                  AUDIT_OUT, STATS_OUT, REPORT_OUT, RUN_MARKER]}
    write_json(MANIFEST_OUT, {
        "analysis": "PD1C_STAGE1", "status": status,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": {rel(Path(__file__)): sha256(Path(__file__))},
        "inputs": checked, "outputs": output_hashes,
    })
    print(json.dumps({"status": status, "stage2_permitted": status == "PASS",
                      "outputs": output_hashes}, indent=2), flush=True)


def preflight() -> None:
    checked = validate_inputs()
    universe, universe_stats = build_covariate_universe()
    panels = pd.read_csv(PANELS, sep="\t")
    panels["gene"] = panels.gene.astype(str).str.upper()
    validate_panels(panels, universe)
    print(json.dumps({"status": "PASS", "formal_weights_calculated": False,
                      "outcome_loaded": False, "panel_rows": len(panels),
                      "panel_draws": panels.groupby(["strategy_id", "draw"]).ngroups,
                      "universe_stats": universe_stats, "verified_inputs": checked}, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    preflight() if args.preflight else execute()


if __name__ == "__main__":
    main()
