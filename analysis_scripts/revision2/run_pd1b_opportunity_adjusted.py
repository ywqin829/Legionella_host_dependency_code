"""Run frozen PD1B annotation/study-opportunity matched comparator audit.

This is a post hoc diagnostic on the already-used WT-versus-uninfected
transcriptional holdout.  It estimates how each frozen comparator strategy
performs relative to panels matched on prespecified annotation, study-intensity
and baseline-expression covariates.  It cannot restore a general superiority
claim and does not modify G3/WP6.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "revision2"))
from build_wp5_comparator_panels import ATLAS_COLUMNS, collapse_atlas, load_discovery_atlas  # noqa: E402

OUT = ROOT / "revision2_2026-11-23"
SPEC = OUT / "REVISION2_PD1B_EXECUTION_SPEC.md"
PANELS = OUT / "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz"
WP6_RESULTS = OUT / "REVISION2_WP6_PANEL_RESULTS.tsv.gz"
WP5_MANIFEST = OUT / "REVISION2_WP5_FREEZE_MANIFEST.json"
WP6_MANIFEST = OUT / "REVISION2_WP6_MANIFEST.json"
ATLAS = ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz"
HOLDOUT = ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz"
GAF = ROOT / "01_data_raw" / "annotations" / "goa_human.gaf.gz"
OBO = ROOT / "01_data_raw" / "annotations" / "go-basic.obo"
DEGREE = ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv"
MATRIX_DIR = ROOT / "05_expression" / "matrices"
MAP_DIR = ROOT / "05_expression" / "gene_maps"

MATCHED_OUT = OUT / "REVISION2_PD1B_PANEL_MATCHED_RESULTS.tsv.gz"
SUMMARY_OUT = OUT / "REVISION2_PD1B_STRATEGY_SUMMARY.tsv"
BALANCE_OUT = OUT / "REVISION2_PD1B_BALANCE.tsv"
STRATA_OUT = OUT / "REVISION2_PD1B_DETECTABILITY_STRATA.tsv"
PERCENTILE_OUT = OUT / "REVISION2_PD1B_FROZEN_PERCENTILES.tsv"
AUDIT_OUT = OUT / "REVISION2_PD1B_AUDIT.tsv"
STATS_OUT = OUT / "REVISION2_PD1B_STATS.json"
REPORT_OUT = OUT / "REVISION2_PD1B_REPORT.md"
MANIFEST_OUT = OUT / "REVISION2_PD1B_MANIFEST.json"
RUN_MARKER = OUT / "REVISION2_PD1B_RUN_MARKER.json"

EXPECTED_HASHES = {
    SPEC: "d786579be757a3ffd3e3bc5f02172aabe60e50103bd0dfaaeee6946851afaf11",
    PANELS: "9c070d7e8b3824bf8cfd5dbd3eb9f641d697ecb9c79106126b57fc70cc9d4d75",
    WP6_RESULTS: "2a603bd481a86a34cbc0232e1bf18bc7408809e8d750b8090522bc371acfd058",
    ATLAS: "9f1ccdc73894da2e1e4ee251db296193df63cd9493b3e16383913e5a3dcb3da5",
    HOLDOUT: "1817f35940d9b5c700cfcadb90ea3a6833217b10482e2c3370c31ecdd826fbff",
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

N_REFERENCE = 100
PANEL_SIZE = 56
MEMBRANE_ROOT = "GO:0016020"
IMBO_ROOT = "GO:0043231"
STRATEGY_CODE = {
    "S1_FROZEN": 1, "S2_CRISPR": 2, "S3_LCV": 3, "S4_COUNT": 4,
    "S5_CRISPR_EXPR": 5, "S6_LCV_EXPR": 6, "S7_COUNT_EXPR": 7,
}
STRATEGY_UNIVERSE = {s: ("expr" if s.endswith("_EXPR") else "all") for s in STRATEGY_CODE}
COVARIATES = [
    "log1p_go_terms_n", "membrane_flag", "imbo_flag", "baseline_ui_percentile",
    "expression_n_studies", "log1p_cross_pathogen_degree", "log1p_cross_pathogen_pmids",
]


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
    for manifest_path in (WP5_MANIFEST, WP6_MANIFEST):
        if not manifest_path.is_file():
            raise FileNotFoundError(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        checked[rel(manifest_path)] = sha256(manifest_path)
        for section in ("inputs", "outputs"):
            for name, expected in dict(manifest.get(section, {})).items():
                path = ROOT / str(name)
                if not path.is_file():
                    raise FileNotFoundError(path)
                actual = sha256(path)
                # The append-only lock is expected to change after each governed phase.
                if path.name != "REVISION2_POSTREVIEW_ANALYSIS_LOCK.md" and actual.lower() != str(expected).lower():
                    raise RuntimeError(f"archived {section[:-1]} hash mismatch: {rel(path)}")
                checked[rel(path)] = actual
    return checked


def load_obo(path: Path) -> dict[str, list[str]]:
    parents: dict[str, list[str]] = {}
    obsolete: set[str] = set()
    cur: str | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("["):
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:].strip()
                parents.setdefault(cur, [])
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
    out = {root}
    stack = [root]
    while stack:
        term = stack.pop()
        for child in children.get(term, []):
            if child not in out:
                out.add(child)
                stack.append(child)
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


def build_universe() -> tuple[pd.DataFrame, dict[str, int]]:
    raw = load_discovery_atlas(str(ATLAS), ATLAS_COLUMNS)
    universe, collapse_stats = collapse_atlas(raw)
    universe["gene"] = universe.gene.astype(str).str.upper()

    parents = load_obo(OBO)
    membrane = descendants(parents, MEMBRANE_ROOT)
    imbo = descendants(parents, IMBO_ROOT)
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
    universe["cross_pathogen_degree"] = pd.to_numeric(universe.cross_pathogen_degree, errors="coerce").fillna(0)
    universe["cross_pathogen_pmids"] = pd.to_numeric(universe.cross_pathogen_pmids, errors="coerce").fillna(0)
    universe["expression_n_studies"] = pd.to_numeric(universe.expression_n_studies, errors="coerce").fillna(0)
    universe["log1p_go_terms_n"] = np.log1p(universe.go_terms_n)
    universe["log1p_cross_pathogen_degree"] = np.log1p(universe.cross_pathogen_degree)
    universe["log1p_cross_pathogen_pmids"] = np.log1p(universe.cross_pathogen_pmids)

    holdout = pd.read_csv(HOLDOUT, sep="\t", low_memory=False)
    wt = holdout.loc[holdout.contrast_id.eq("WT_vs_UI")].copy()
    wt["gene"] = wt.human_symbol.astype(str).str.upper()
    if wt.gene.duplicated().any():
        raise RuntimeError("duplicate WT_vs_UI genes")
    p = np.clip(pd.to_numeric(wt.PValue, errors="coerce").to_numpy(float), 1e-300, 1 - 1e-15)
    wt["raw_signed_z"] = np.sign(pd.to_numeric(wt.logFC, errors="coerce").to_numpy(float)) * norm.isf(p / 2)
    universe = universe.merge(wt[["gene", "raw_signed_z"]], on="gene", how="left")
    universe["direction_sign"] = universe.expected_direction.map({"UP": 1.0, "DOWN": -1.0})
    universe["oriented_z"] = universe.raw_signed_z * universe.direction_sign
    if universe.gene.duplicated().any() or not np.isfinite(universe[COVARIATES].to_numpy(float)).all():
        raise RuntimeError("invalid unique/covariate universe")
    ontology = {"membrane_descendants_including_root": len(membrane), "imbo_descendants_including_root": len(imbo)}
    ontology.update(collapse_stats)
    return universe, ontology


def validate_panels(panels: pd.DataFrame) -> None:
    required = {"strategy_id", "draw", "gene", "expected_direction"}
    if not required <= set(panels.columns):
        raise RuntimeError(f"missing panel columns: {sorted(required - set(panels.columns))}")
    for sid, group in panels.groupby("strategy_id"):
        if sid not in STRATEGY_CODE or group.groupby("draw").size().ne(PANEL_SIZE).any():
            raise RuntimeError(f"invalid frozen panels for {sid}")
        if group.groupby("draw").gene.nunique().ne(PANEL_SIZE).any():
            raise RuntimeError(f"duplicate frozen panel genes for {sid}")


def fit_propensity(frame: pd.DataFrame, panels: pd.DataFrame, sid: str) -> tuple[np.ndarray, pd.Series, int]:
    sub = panels.loc[panels.strategy_id.eq(sid)].copy()
    n_draws = int(sub.draw.nunique())
    counts = sub.groupby(sub.gene.astype(str).str.upper()).size()
    freq = frame.gene.map(counts).fillna(0.0) / n_draws
    x = frame[COVARIATES].astype(float)
    x = (x - x.mean()) / x.std(ddof=0).replace(0, 1)
    model = sm.GLM(freq.to_numpy(float), sm.add_constant(x, has_constant="add"),
                   family=sm.families.Binomial(), freq_weights=np.repeat(n_draws, len(frame)))
    result = model.fit(maxiter=200, disp=0)
    prop = np.clip(result.predict(sm.add_constant(x, has_constant="add")), 1e-9, 1 - 1e-9)
    if not np.isfinite(prop).all():
        raise RuntimeError(f"nonfinite propensity for {sid}")
    return np.asarray(prop, float), freq, n_draws


def smd(a_sum: np.ndarray, a_sumsq: np.ndarray, a_n: int,
        b_sum: np.ndarray, b_sumsq: np.ndarray, b_n: int) -> np.ndarray:
    am, bm = a_sum / a_n, b_sum / b_n
    av = np.maximum((a_sumsq - a_sum * a_sum / a_n) / max(a_n - 1, 1), 0)
    bv = np.maximum((b_sumsq - b_sum * b_sum / b_n) / max(b_n - 1, 1), 0)
    den = np.sqrt((av + bv) / 2)
    den[den == 0] = 1.0
    return (am - bm) / den


def evaluate_strategy(sid: str, universe: pd.DataFrame, panels: pd.DataFrame,
                      writer: csv.DictWriter) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    base = universe.loc[universe.expression_n_studies.ge(1)].copy() if STRATEGY_UNIVERSE[sid] == "expr" else universe.copy()
    frame = base.loc[base.expected_direction.isin(["UP", "DOWN"]) & np.isfinite(base.oriented_z)].copy().reset_index(drop=True)
    prop, selection_freq, n_draws = fit_propensity(frame, panels, sid)
    frame["propensity"] = prop
    frame["selection_frequency"] = selection_freq.to_numpy(float)
    controls = np.where(frame.selection_frequency.eq(0).to_numpy())[0]
    if len(controls) < PANEL_SIZE:
        raise RuntimeError(f"insufficient zero-frequency controls for {sid}")
    gene_to_idx = pd.Series(frame.index.to_numpy(), index=frame.gene).to_dict()
    panel_sub = panels.loc[panels.strategy_id.eq(sid)].copy()
    panel_sub["gene"] = panel_sub.gene.astype(str).str.upper()

    selected_genes = sorted(set(panel_sub.gene) & set(gene_to_idx))
    pools: dict[str, tuple[np.ndarray, str]] = {}
    for gene in selected_genes:
        idx = int(gene_to_idx[gene])
        delta = np.abs(prop[controls] - prop[idx])
        pool = controls[delta <= 0.02]
        label = "0.02"
        if len(pool) < 20:
            pool = controls[delta <= 0.05]
            label = "0.05"
        if len(pool) < 20:
            order = np.argsort(delta, kind="mergesort")[:min(200, len(controls))]
            pool = controls[order]
            label = "nearest200"
        pools[gene] = (np.asarray(pool, int), label)

    values = frame[COVARIATES].to_numpy(float)
    oriented = frame.oriented_z.to_numpy(float)
    q_sum = np.zeros(len(COVARIATES)); q_sumsq = np.zeros(len(COVARIATES)); q_n = 0
    m_sum = np.zeros(len(COVARIATES)); m_sumsq = np.zeros(len(COVARIATES)); m_n = 0
    control_vals = values[controls]
    c_sum = control_vals.sum(axis=0); c_sumsq = (control_vals ** 2).sum(axis=0)
    summaries: list[dict[str, object]] = []
    relax_counts = {"0.02": 0, "0.05": 0, "nearest200": 0}

    for panel_draw, panel in panel_sub.groupby("draw", sort=True):
        query_genes = [g for g in panel.gene if g in gene_to_idx]
        query_idx = np.array([gene_to_idx[g] for g in query_genes], dtype=int)
        observed = float(oriented[query_idx].mean())
        qv = values[query_idx]
        q_sum += qv.sum(axis=0); q_sumsq += (qv ** 2).sum(axis=0); q_n += len(query_idx)
        ordering = sorted(query_genes, key=lambda g: (len(pools[g][0]), g))
        refs = np.empty(N_REFERENCE, dtype=float)
        for reference_draw in range(N_REFERENCE):
            rng = np.random.default_rng(np.random.SeedSequence([42, 611, STRATEGY_CODE[sid], int(panel_draw), reference_draw]))
            used: set[int] = set()
            chosen: list[int] = []
            for gene in ordering:
                pool, label = pools[gene]
                available = pool[~np.isin(pool, np.fromiter(used, dtype=int) if used else np.array([], dtype=int))]
                if len(available) == 0:
                    raise RuntimeError(f"unmatched slot {sid}/{panel_draw}/{reference_draw}/{gene}")
                pick = int(rng.choice(available))
                used.add(pick); chosen.append(pick); relax_counts[label] += 1
            if len(chosen) != len(query_idx) or len(set(chosen)) != len(chosen):
                raise RuntimeError("matched reference size/uniqueness failure")
            chosen_idx = np.asarray(chosen, int)
            statistic = float(oriented[chosen_idx].mean())
            refs[reference_draw] = statistic
            mv = values[chosen_idx]
            m_sum += mv.sum(axis=0); m_sumsq += (mv ** 2).sum(axis=0); m_n += len(chosen_idx)
            writer.writerow({
                "strategy_id": sid, "panel_draw": int(panel_draw), "reference_draw": reference_draw,
                "selected_n": len(query_idx), "observed_statistic": f"{observed:.12g}",
                "matched_statistic": f"{statistic:.12g}", "opportunity_excess": f"{observed-statistic:.12g}",
            })
        summaries.append({
            "strategy_id": sid, "panel_draw": int(panel_draw), "selected_n": len(query_idx),
            "observed_statistic": observed, "matched_median": float(np.median(refs)),
            "matched_p025": float(np.quantile(refs, .025)), "matched_p975": float(np.quantile(refs, .975)),
            "opportunity_excess": observed - float(np.median(refs)),
        })

    before = smd(q_sum, q_sumsq, q_n, c_sum, c_sumsq, len(controls))
    after = smd(q_sum, q_sumsq, q_n, m_sum, m_sumsq, m_n)
    balance = [{"strategy_id": sid, "covariate": cov, "smd_before": float(b), "smd_after": float(a),
                "abs_smd_after": float(abs(a))} for cov, b, a in zip(COVARIATES, before, after)]
    detail = {
        "strategy_id": sid, "universe_n": int(len(base)), "evaluable_universe_n": int(len(frame)),
        "zero_frequency_controls_n": int(len(controls)), "panel_draws": n_draws,
        "query_slots": q_n, "matched_slots": m_n, "max_abs_smd_after": float(np.max(np.abs(after))),
        "pool_slot_counts": relax_counts,
    }
    return summaries, balance, detail


def detectability_strata(universe: pd.DataFrame, panels: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sid, psub in panels.groupby("strategy_id", sort=True):
        base = universe.loc[universe.expression_n_studies.ge(1)].copy() if STRATEGY_UNIVERSE[sid] == "expr" else universe.copy()
        base["detectability_quartile"] = pd.qcut(base.baseline_ui_percentile.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
        lookup = base.set_index("gene")
        psub = psub.copy(); psub["gene"] = psub.gene.astype(str).str.upper()
        draw_rows: list[dict[str, object]] = []
        for draw, panel in psub.groupby("draw"):
            mapped = lookup.reindex(panel.gene)
            for quartile in ["Q1", "Q2", "Q3", "Q4"]:
                q = mapped.loc[mapped.detectability_quartile.astype(str).eq(quartile)]
                finite = q.loc[np.isfinite(q.oriented_z)]
                draw_rows.append({"draw": draw, "detectability_quartile": quartile, "selected_n": len(q),
                                  "evaluable_n": len(finite), "mean_oriented_z": finite.oriented_z.mean()})
        d = pd.DataFrame(draw_rows)
        for quartile, group in d.groupby("detectability_quartile", sort=True):
            rows.append({
                "strategy_id": sid, "detectability_quartile": quartile,
                "selected_n_median": float(group.selected_n.median()),
                "evaluable_n_median": float(group.evaluable_n.median()),
                "evaluable_rate_median": float((group.evaluable_n / group.selected_n.replace(0, np.nan)).median()),
                "mean_oriented_z_median": float(group.mean_oriented_z.median()),
                "mean_oriented_z_p025": float(group.mean_oriented_z.quantile(.025)),
                "mean_oriented_z_p975": float(group.mean_oriented_z.quantile(.975)),
            })
    return pd.DataFrame(rows)


def frozen_percentiles(wp6: pd.DataFrame) -> pd.DataFrame:
    frozen = float(wp6.loc[wp6.strategy_id.eq("S1_FROZEN"), "statistic"].iloc[0])
    rows = []
    for sid in ("S3_LCV", "S4_COUNT", "S6_LCV_EXPR", "S7_COUNT_EXPR"):
        vals = pd.to_numeric(wp6.loc[wp6.strategy_id.eq(sid), "statistic"], errors="coerce").dropna().to_numpy()
        rows.append({"comparator_strategy": sid, "frozen_statistic": frozen, "comparator_draws": len(vals),
                     "frozen_percentile": float(((vals <= frozen).sum() + 1) / (len(vals) + 1)),
                     "comparator_median": float(np.median(vals))})
    return pd.DataFrame(rows)


def add_audit(rows: list[dict[str, object]], check: str, status: str, detail: str) -> None:
    rows.append({"check": check, "status": status, "detail": detail})


def execute() -> None:
    outputs = [MATCHED_OUT, SUMMARY_OUT, BALANCE_OUT, STRATA_OUT, PERCENTILE_OUT,
               AUDIT_OUT, STATS_OUT, REPORT_OUT, MANIFEST_OUT, RUN_MARKER]
    existing = [str(p) for p in outputs if p.exists()]
    if existing:
        raise RuntimeError(f"refusing to overwrite PD1B outputs: {existing}")
    input_hash = validate_inputs()
    panels = pd.read_csv(PANELS, sep="\t")
    validate_panels(panels)
    universe, ontology_stats = build_universe()
    wp6 = pd.read_csv(WP6_RESULTS, sep="\t")

    all_summaries: list[dict[str, object]] = []
    all_balance: list[dict[str, object]] = []
    details: dict[str, object] = {}
    with gzip.open(MATCHED_OUT, "wt", encoding="utf-8", newline="") as handle:
        fields = ["strategy_id", "panel_draw", "reference_draw", "selected_n", "observed_statistic", "matched_statistic", "opportunity_excess"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sid in STRATEGY_CODE:
            print(f"PD1B matching {sid} ...", flush=True)
            summaries, balance, detail = evaluate_strategy(sid, universe, panels, writer)
            all_summaries.extend(summaries); all_balance.extend(balance); details[sid] = detail

    panel_summary = pd.DataFrame(all_summaries)
    strategy_summary = panel_summary.groupby("strategy_id", sort=True).agg(
        panel_draws=("panel_draw", "nunique"), selected_n_median=("selected_n", "median"),
        observed_median=("observed_statistic", "median"), observed_p025=("observed_statistic", lambda x: x.quantile(.025)),
        observed_p975=("observed_statistic", lambda x: x.quantile(.975)), matched_median=("matched_median", "median"),
        matched_p025=("matched_median", lambda x: x.quantile(.025)), matched_p975=("matched_median", lambda x: x.quantile(.975)),
        opportunity_excess_median=("opportunity_excess", "median"),
        opportunity_excess_p025=("opportunity_excess", lambda x: x.quantile(.025)),
        opportunity_excess_p975=("opportunity_excess", lambda x: x.quantile(.975)),
    ).reset_index()
    strategy_summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    balance = pd.DataFrame(all_balance)
    balance.to_csv(BALANCE_OUT, sep="\t", index=False)
    strata = detectability_strata(universe, panels)
    strata.to_csv(STRATA_OUT, sep="\t", index=False)
    percentiles = frozen_percentiles(wp6)
    percentiles.to_csv(PERCENTILE_OUT, sep="\t", index=False)

    audits: list[dict[str, object]] = []
    add_audit(audits, "P01_INPUT_HASHES", "PASS", f"{len(input_hash)} files verified")
    add_audit(audits, "P02_PANEL_ROWS", "PASS" if len(panels) == 224168 else "FAIL", str(len(panels)))
    expected_refs = int(panels.groupby(["strategy_id", "draw"]).ngroups * N_REFERENCE)
    actual_refs = sum(int(d["panel_draws"]) * N_REFERENCE for d in details.values())
    add_audit(audits, "P03_REFERENCE_ROWS", "PASS" if actual_refs == expected_refs else "FAIL", f"{actual_refs}/{expected_refs}")
    add_audit(audits, "P04_FINITE_RESULTS", "PASS" if np.isfinite(panel_summary.select_dtypes("number")).all().all() else "FAIL", f"panel summaries={len(panel_summary)}")
    for sid, detail in details.items():
        max_smd = float(detail["max_abs_smd_after"])
        status = "PASS" if max_smd <= .10 else ("AUDIT" if max_smd <= .20 else "FAIL")
        add_audit(audits, f"P05_BALANCE_{sid}", status, f"max_abs_smd_after={max_smd:.6f}")
    # Archived WP6 statistics must be exactly reproduced by the selected evaluable genes.
    mismatches = []
    for row in panel_summary.itertuples(index=False):
        archived = wp6.loc[wp6.strategy_id.eq(row.strategy_id) & wp6.draw.eq(row.panel_draw), "statistic"]
        if len(archived) != 1 or not math.isclose(float(archived.iloc[0]), float(row.observed_statistic), abs_tol=1e-12, rel_tol=0):
            mismatches.append(f"{row.strategy_id}:{row.panel_draw}")
    add_audit(audits, "P06_WP6_STATISTIC_REPRODUCTION", "PASS" if not mismatches else "FAIL", f"mismatches={len(mismatches)}")
    forbidden = {"discovery_domain_count", "direct", "functional_clean", "response", "mechanism", "spatial",
                 "gw_signed_confidence", "lcv_study_count", "ppi_AB_flag", "raw_signed_z", "oriented_z"}
    add_audit(audits, "P07_NO_CIRCULAR_COVARIATES", "PASS" if not forbidden.intersection(COVARIATES) else "FAIL", str(COVARIATES))
    audit = pd.DataFrame(audits)
    audit.to_csv(AUDIT_OUT, sep="\t", index=False)
    overall = "FAIL" if audit.status.eq("FAIL").any() else ("PASS_WITH_AUDIT" if audit.status.eq("AUDIT").any() else "PASS")

    count_excess = float(strategy_summary.loc[strategy_summary.strategy_id.eq("S4_COUNT"), "opportunity_excess_median"].iloc[0])
    frozen_excess = float(strategy_summary.loc[strategy_summary.strategy_id.eq("S1_FROZEN"), "opportunity_excess_median"].iloc[0])
    stats = {
        "analysis": "PD1B annotation/study-opportunity matched comparator audit",
        "status": overall, "post_hoc_same_holdout": True, "n_reference_per_panel": N_REFERENCE,
        "covariates": COVARIATES, "strategy_details": details, "ontology_stats": ontology_stats,
        "count_minus_frozen_opportunity_excess": count_excess - frozen_excess,
        "interpretive_boundary": "Post hoc explanatory diagnostic on the already-used holdout; it does not provide independent validation or restore a general superiority claim.",
    }
    write_json(STATS_OUT, stats)

    lines = [
        "# PD1B annotation/study-opportunity matched comparator audit", "",
        f"- Gate: **{overall}**", f"- Reference panels per archived panel: {N_REFERENCE}",
        "- This was a post hoc diagnostic on the already-used WT-versus-uninfected transcriptional holdout.",
        "- Matching used only prespecified GO opportunity, broad membrane/organelle annotation, baseline expression, discovery expression-study count, and cross-pathogen degree/publication count.",
        "- Discovery-domain count, its component domains, CRISPR/LCV/PPI metrics and the holdout outcome were prohibited matching covariates.", "",
        "## Strategy-level results", "", strategy_summary.to_markdown(index=False), "",
        "## Balance", "", balance.to_markdown(index=False), "",
        "## Frozen-panel percentile in stochastic comparator distributions", "", percentiles.to_markdown(index=False), "",
        "## Interpretation boundary", "",
        "These results estimate attenuation after adjustment for measured annotation and study opportunity. Residual differences may reflect unmeasured ascertainment, selection-rule structure, or biology. They are not an independent benchmark and cannot establish general prioritisation superiority.", "",
    ]
    REPORT_OUT.write_text("\n".join(lines), encoding="utf-8")
    write_json(RUN_MARKER, {"status": overall, "completed_utc": datetime.now(timezone.utc).isoformat(), "spec_sha256": sha256(SPEC)})
    output_hashes = {rel(p): sha256(p) for p in [MATCHED_OUT, SUMMARY_OUT, BALANCE_OUT, STRATA_OUT, PERCENTILE_OUT, AUDIT_OUT, STATS_OUT, REPORT_OUT, RUN_MARKER]}
    write_json(MANIFEST_OUT, {
        "analysis": "PD1B", "status": overall, "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": {rel(Path(__file__)): sha256(Path(__file__))}, "inputs": input_hash, "outputs": output_hashes,
    })
    print(json.dumps({"status": overall, "outputs": output_hashes}, indent=2), flush=True)


def preflight() -> None:
    checked = validate_inputs()
    panels = pd.read_csv(PANELS, sep="\t")
    validate_panels(panels)
    headers = pd.read_csv(HOLDOUT, sep="\t", nrows=3)
    if not {"contrast_id", "human_symbol", "logFC", "PValue"} <= set(headers.columns):
        raise RuntimeError("holdout schema mismatch")
    print(json.dumps({"status": "PASS", "formal_statistics_calculated": False,
                      "panel_rows": len(panels), "panel_draws": panels.groupby(["strategy_id", "draw"]).ngroups,
                      "verified_inputs": checked}, indent=2), flush=True)


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
