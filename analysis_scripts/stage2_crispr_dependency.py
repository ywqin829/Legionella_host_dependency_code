from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_data_raw" / "crispr"
DEP_RAW = ROOT / "01_data_raw" / "depmap"
INTERACTOME = ROOT / "03_interactome" / "legionella_host_target_summary.tsv"
OUT = ROOT / "04_genetic_dependency"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def clean_symbol(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.upper().replace({"": pd.NA, "NAN": pd.NA})


def direction(effect: float | int | None) -> str:
    if pd.isna(effect):
        return "UNRESOLVED"
    if float(effect) > 0:
        return "KO_PROTECTS_HOST_ENRICHED_AFTER_INFECTION"
    if float(effect) < 0:
        return "KO_SENSITIZES_HOST_DEPLETED_AFTER_INFECTION"
    return "NO_DIRECTION"


def sign_value(value: float | int | None) -> int:
    if pd.isna(value) or float(value) == 0:
        return 0
    return 1 if float(value) > 0 else -1


def load_optional_common_essentials() -> tuple[set[str], dict]:
    """Load a user-supplied official DepMap common-essential list if present.

    Automated download may be blocked by human verification. We never infer or
    fabricate the list from portal HTML.
    """
    candidates = [
        DEP_RAW / "CRISPRInferredCommonEssentials.csv",
        DEP_RAW / "AchillesCommonEssentialControls.csv",
    ]
    for path in candidates:
        if not path.exists() or path.stat().st_size < 100:
            continue
        try:
            df = pd.read_csv(path, header=None)
            symbols = set(
                df.iloc[:, 0]
                .astype(str)
                .str.extract(r"^\s*([^\s,(]+)", expand=False)
                .str.upper()
                .dropna()
            )
            if len(symbols) >= 100:
                return symbols, {
                    "status": "AVAILABLE",
                    "source_file": str(path),
                    "sha256": sha256(path),
                    "n_symbols": len(symbols),
                }
        except Exception as exc:  # fail-soft external annotation
            return set(), {"status": "UNRESOLVED", "reason": f"{type(exc).__name__}: {exc}"}
    return set(), {
        "status": "PENDING_EXTERNAL_FILE",
        "reason": "Official DepMap download endpoint required human verification; no list was inferred from HTML.",
    }


def load_optional_u937_effect() -> tuple[dict[str, float], dict]:
    path = DEP_RAW / "CRISPRGeneEffect.csv"
    if not path.exists() or path.stat().st_size < 1000:
        return {}, {"status": "PENDING_EXTERNAL_FILE", "reason": "CRISPRGeneEffect.csv not locally available"}
    try:
        df = pd.read_csv(path, index_col=0, low_memory=False)
        u937_ids = [x for x in df.index.astype(str) if x in {"ACH-000406", "U937", "U-937"}]
        if not u937_ids:
            return {}, {"status": "UNRESOLVED", "reason": "U937 model row not found"}
        row = df.loc[u937_ids[0]]
        values: dict[str, float] = {}
        for col, val in row.items():
            symbol = str(col).split(" (")[0].strip().upper()
            v = pd.to_numeric(pd.Series([val]), errors="coerce").iloc[0]
            if pd.notna(v):
                values[symbol] = float(v)
        return values, {
            "status": "AVAILABLE",
            "source_file": str(path),
            "sha256": sha256(path),
            "model_id": u937_ids[0],
            "sign_convention": "negative raw Chronos score indicates stronger dependency",
        }
    except Exception as exc:
        return {}, {"status": "UNRESOLVED", "reason": f"{type(exc).__name__}: {exc}"}


def prep_batch(path: Path, sheet: str, prefix: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet)
    df = df.rename(columns={"#Gene ID": "gene_id", "Symbol": "symbol"})
    df["symbol"] = clean_symbol(df["symbol"])
    df[f"{prefix}_effect"] = numeric(df["Combo Effect"])
    df[f"{prefix}_signed_confidence"] = numeric(df["Signed Combo Confidence Score"])
    df[f"{prefix}_hit95"] = numeric(df["Hits (95% CI)"]).fillna(0).eq(1)
    df[f"{prefix}_direction"] = df[f"{prefix}_effect"].map(direction)
    keep = ["gene_id", "symbol", f"{prefix}_effect", f"{prefix}_signed_confidence", f"{prefix}_hit95", f"{prefix}_direction"]
    return df[keep].dropna(subset=["symbol"]).drop_duplicates("symbol")


def functional_tier(row: pd.Series) -> str:
    gw = bool(row.get("gw_hit10", False))
    u = bool(row.get("undiff_hit95", False))
    d = bool(row.get("diff_hit95", False))
    signs = [x for x in [row.get("gw_sign", 0), row.get("undiff_sign", 0), row.get("diff_sign", 0)] if x != 0]
    concordant = bool(signs) and len(set(signs)) == 1
    if gw and u and d and concordant:
        return "F1_GW_AND_BOTH_BATCH_CONCORDANT"
    if concordant and ((gw and (u or d)) or (u and d)):
        return "F2_MULTILAYER_CONCORDANT"
    if gw or u or d:
        return "F3_SINGLE_LAYER_OR_DISCORDANT"
    return "NO_FUNCTIONAL_HIT"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    required = [RAW / f"Table_S{i}.xlsx" for i in range(1, 6)]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing CRISPR supplements: {missing}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        gw = pd.read_excel(RAW / "Table_S1.xlsx", sheet_name="Genome-wide screen")
        undiff = prep_batch(RAW / "Table_S3.xlsx", "Undifferentiated U937 batch", "undiff")
        diff = prep_batch(RAW / "Table_S3.xlsx", "Differentiated U937 batch", "diff")
        apms = pd.read_excel(RAW / "Table_S4.xlsx", sheet_name="RABIF AP-MS SAINT analysis")

    gw = gw.rename(columns={"Gene ID": "gene_id", "Symbol": "symbol", "Gene Info": "gene_info"})
    gw["symbol"] = clean_symbol(gw["symbol"])
    gw["gw_effect"] = numeric(gw["Effect"])
    gw["gw_signed_confidence"] = numeric(gw["Signed Confidence Score"])
    gw["gw_pvalue"] = numeric(gw["P value"])
    for level in (5, 10, 20):
        gw[f"gw_hit{level}"] = numeric(gw[f"{level}% FDR"]).fillna(0).eq(1)
    gw["gw_direction"] = gw["gw_effect"].map(direction)
    gw["gw_sign"] = gw["gw_effect"].map(sign_value)
    gw_keep = [
        "gene_id", "symbol", "gene_info", "gw_effect", "gw_signed_confidence", "gw_pvalue",
        "gw_hit5", "gw_hit10", "gw_hit20", "gw_direction", "gw_sign",
    ]
    gw = gw[gw_keep].dropna(subset=["symbol"]).drop_duplicates("symbol")

    merged = gw.merge(undiff, on="symbol", how="outer", suffixes=("", "_undiff_id"))
    merged = merged.merge(diff, on="symbol", how="outer", suffixes=("", "_diff_id"))
    for col in ["gw_hit5", "gw_hit10", "gw_hit20", "undiff_hit95", "diff_hit95"]:
        merged[col] = merged[col].astype("boolean").fillna(False).astype(bool)
    merged["undiff_sign"] = merged["undiff_effect"].map(sign_value)
    merged["diff_sign"] = merged["diff_effect"].map(sign_value)
    merged["all_nonzero_directions_concordant"] = merged.apply(
        lambda r: len({x for x in [r.gw_sign, r.undiff_sign, r.diff_sign] if x != 0}) <= 1,
        axis=1,
    )
    merged["functional_tier"] = merged.apply(functional_tier, axis=1)

    common_essential, common_meta = load_optional_common_essentials()
    u937_effect, u937_meta = load_optional_u937_effect()
    if common_meta["status"] == "AVAILABLE":
        merged["common_essential_flag"] = merged.symbol.isin(common_essential)
    else:
        merged["common_essential_flag"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")
    if u937_meta["status"] == "AVAILABLE":
        merged["u937_chronos_gene_effect"] = merged.symbol.map(u937_effect)
        merged["u937_strong_dependency_flag"] = merged["u937_chronos_gene_effect"].le(-0.5)
    else:
        merged["u937_chronos_gene_effect"] = np.nan
        merged["u937_strong_dependency_flag"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")

    merged["fitness_artifact_compatible"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")
    available_fitness = merged["common_essential_flag"].notna() | merged["u937_strong_dependency_flag"].notna()
    depleted = merged["gw_sign"].lt(0)
    essential = merged["common_essential_flag"].fillna(False) | merged["u937_strong_dependency_flag"].fillna(False)
    merged.loc[available_fitness, "fitness_artifact_compatible"] = (depleted & essential)[available_fitness]
    caveat_true = merged["fitness_artifact_compatible"].eq(True).fillna(False).to_numpy(dtype=bool)
    caveat_false = merged["fitness_artifact_compatible"].eq(False).fillna(False).to_numpy(dtype=bool)
    merged["baseline_fitness_caveat"] = np.select(
        [
            caveat_true,
            caveat_false,
        ],
        ["POTENTIAL_BASELINE_FITNESS_ARTIFACT", "NO_FLAG_FROM_AVAILABLE_DEPMAP_ANNOTATION"],
        default="UNRESOLVED_EXTERNAL_DEPMAP_DATA",
    )

    apms = apms.rename(
        columns={
            "preyProtein": "prey_uniprot",
            "preyGene": "symbol",
            "Fold change vs. control (FCB)": "fold_change_vs_control",
            "SAINT score": "saint_score",
        }
    )
    apms["symbol"] = clean_symbol(apms["symbol"])
    apms["saint_score"] = numeric(apms["saint_score"])
    apms["fold_change_vs_control"] = numeric(apms["fold_change_vs_control"])
    apms["apms_high_confidence"] = apms.saint_score.ge(0.90) & apms.fold_change_vs_control.ge(3.0)
    apms["apms_supportive"] = apms.saint_score.ge(0.80) & apms.fold_change_vs_control.ge(3.0)

    ppi = pd.read_csv(INTERACTOME, sep="\t")
    ppi["symbol"] = clean_symbol(ppi["gene"])
    ppi = ppi.rename(
        columns={
            "human_accession": "ppi_human_accession",
            "best_class": "ppi_best_class",
            "bacterial_partners": "ppi_bacterial_partners",
            "interaction_records": "ppi_interaction_records",
            "n_pmids": "ppi_n_pmids",
        }
    )
    ppi_keep = [
        "symbol", "ppi_human_accession", "ppi_best_class", "ppi_bacterial_partners",
        "ppi_interaction_records", "ppi_n_pmids",
    ]
    integrated = merged.merge(ppi[ppi_keep].drop_duplicates("symbol"), on="symbol", how="left")
    apms_gene = apms.groupby("symbol", as_index=False).agg(
        apms_max_saint=("saint_score", "max"),
        apms_max_fc=("fold_change_vs_control", "max"),
        apms_high_confidence=("apms_high_confidence", "max"),
        apms_supportive=("apms_supportive", "max"),
    )
    integrated = integrated.merge(apms_gene, on="symbol", how="left")
    integrated["ppi_AB_flag"] = integrated.ppi_best_class.isin(["A", "B"])
    integrated["ppi_crispr_convergence"] = integrated.ppi_AB_flag & integrated.functional_tier.ne("NO_FUNCTIONAL_HIT")

    gw.sort_values(["gw_hit10", "gw_signed_confidence"], ascending=[False, False]).to_csv(
        OUT / "genomewide_crispr.tsv", sep="\t", index=False
    )
    undiff.to_csv(OUT / "batch_undifferentiated.tsv", sep="\t", index=False)
    diff.to_csv(OUT / "batch_differentiated.tsv", sep="\t", index=False)
    merged.to_csv(OUT / "crispr_reproducibility.tsv", sep="\t", index=False)
    apms.to_csv(OUT / "rabif_apms.tsv", sep="\t", index=False)
    integrated.to_csv(OUT / "ppi_crispr_intersection.tsv", sep="\t", index=False)

    stats = {
        "input_sha256": {p.name: sha256(p) for p in required},
        "genomewide_genes": int(len(gw)),
        "gw_hits_5pct_fdr": int(gw.gw_hit5.sum()),
        "gw_hits_10pct_fdr": int(gw.gw_hit10.sum()),
        "gw_hits_20pct_fdr": int(gw.gw_hit20.sum()),
        "gw_10pct_direction_counts": gw.loc[gw.gw_hit10, "gw_direction"].value_counts().to_dict(),
        "undifferentiated_batch_hits": int(undiff.undiff_hit95.sum()),
        "differentiated_batch_hits": int(diff.diff_hit95.sum()),
        "functional_tier_counts": merged.functional_tier.value_counts().to_dict(),
        "rabif_apms_high_confidence": int(apms.apms_high_confidence.sum()),
        "rabif_apms_supportive": int(apms.apms_supportive.sum()),
        "ppi_AB_crispr_convergence": int(integrated.ppi_crispr_convergence.sum()),
        "common_essential_annotation": common_meta,
        "u937_gene_effect_annotation": u937_meta,
        "gate2_main_screen": "PASS",
        "gate2b_fitness_caveat": "PASS" if common_meta["status"] == "AVAILABLE" else "PENDING_EXTERNAL_FILE",
        "screen_endpoint_warning": "Enrichment/depletion after infection is not equivalent to bacterial entry, replication, or host-death mechanism.",
    }
    (OUT / "stage2_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")

    report = [
        "# Stage 2 Functional Dependency Audit",
        "",
        "## Main screen",
        "",
        f"- Genome-wide genes: {stats['genomewide_genes']}",
        f"- Hits at 5/10/20% FDR: {stats['gw_hits_5pct_fdr']} / {stats['gw_hits_10pct_fdr']} / {stats['gw_hits_20pct_fdr']}",
        f"- 10% FDR directions: {stats['gw_10pct_direction_counts']}",
        f"- Undifferentiated/differentiated batch hits: {stats['undifferentiated_batch_hits']} / {stats['differentiated_batch_hits']}",
        f"- Functional tier counts: {stats['functional_tier_counts']}",
        "",
        "## Orthogonal RABIF AP-MS",
        "",
        f"- High-confidence (SAINT >= 0.90 and FC >= 3): {stats['rabif_apms_high_confidence']}",
        f"- Supportive (SAINT >= 0.80 and FC >= 3): {stats['rabif_apms_supportive']}",
        f"- Class A/B PPI plus any CRISPR hit convergence: {stats['ppi_AB_crispr_convergence']}",
        "",
        "## Fitness caveat",
        "",
        f"- Common-essential annotation: {common_meta['status']}",
        f"- U937 Chronos annotation: {u937_meta['status']}",
        "- Missing DepMap files do not invalidate the Legionella screen; they prevent infection-specificity claims until annotated.",
        "- Negative raw Chronos score means stronger dependency when the optional file becomes available.",
        "",
        "## Claim boundary",
        "",
        "A positive infection-screen effect means the knockout was enriched after infection; a negative effect means depletion. Neither alone identifies whether the mechanism is entry, intracellular replication, general fitness, or host-cell death.",
    ]
    (OUT / "STAGE2_FUNCTIONAL_DEPENDENCY.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage2.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
