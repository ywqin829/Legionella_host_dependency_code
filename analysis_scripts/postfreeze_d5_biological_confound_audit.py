from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)
N_MATCHED_NULL = 5000
MATCH_POOL = 50


def signed_z(frame: pd.DataFrame) -> np.ndarray:
    p = np.clip(frame["PValue"].to_numpy(float), 1e-300, 1 - 1e-15)
    return np.sign(frame["logFC"].to_numpy(float)) * stats.norm.isf(p / 2)


def parse_condition(stimulus: str) -> str:
    if stimulus.startswith("Uninfected"):
        return "UI"
    if stimulus.startswith("WT "):
        return "WT"
    match = re.match(r"Delta-([A-Za-z]+)", stimulus)
    if not match:
        raise ValueError(f"Unrecognized stimulus: {stimulus}")
    label = match.group(1)
    return "dotA" if label.lower() == "dota" else label


def audit_sdrf() -> pd.DataFrame:
    path = ROOT / "02_metadata" / "E-MTAB-16655.sdrf.txt"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    source_i = header.index("Source Name")
    stimulus_i = header.index("Characteristics[stimulus]")
    # Each biological sample occurs three times in SDRF (R1/I1/I2 files).
    samples = {}
    for row in rows[1:]:
        samples.setdefault(row[source_i], row[stimulus_i])
    records = []
    for source, stimulus in sorted(samples.items()):
        records.append(
            {
                "source_name": source,
                "official_condition": parse_condition(stimulus),
                "official_stimulus": stimulus,
            }
        )
    sample_df = pd.DataFrame(records)
    counts = (
        sample_df.groupby("official_condition", as_index=False)
        .agg(n_biological_samples=("source_name", "nunique"), official_stimuli=("official_stimulus", lambda x: " | ".join(sorted(x))))
    )
    order = ["UI", "WT", "dotA", "A", "B", "C", "D", "E", "F", "FGHI", "I", "O", "P", "R", "S"]
    counts["order"] = counts["official_condition"].map({x: i for i, x in enumerate(order)})
    counts = counts.sort_values("order").drop(columns="order")
    counts["label_source"] = "E-MTAB-16655 SDRF official sample-level stimulus"
    counts["gene_content_source"] = np.where(
        counts.official_condition.isin(["A", "B", "C", "D", "E", "F", "I"]),
        "external strain-region source plus 130b annotation; separate from SDRF label evidence",
        np.where(
            counts.official_condition.eq("FGHI"),
            "composite label is official; F/H/I content mapped externally; G content pending",
            np.where(
                counts.official_condition.isin(["O", "P", "R", "S"]),
                "official label present; deletion content pending primary strain table",
                "not applicable",
            ),
        ),
    )
    counts["mapping_status"] = np.where(
        counts.official_condition.isin(["O", "P", "R", "S"]),
        "LABEL_CONFIRMED_CONTENT_PENDING",
        np.where(counts.official_condition.eq("FGHI"), "LABEL_CONFIRMED_CONTENT_PARTIAL", "LABEL_CONFIRMED"),
    )
    counts.to_csv(OUT / "d5_holdout_official_condition_audit.tsv", sep="\t", index=False)
    sample_df.to_csv(OUT / "d5_holdout_official_sample_manifest.tsv", sep="\t", index=False)
    return counts


def build_matching_features() -> tuple[pd.DataFrame, list[str]]:
    norm = pd.read_csv(ROOT / "01_data_raw" / "holdout" / "E-MTAB-16655" / "FilteredNormLogCPM.csv")
    raw = pd.read_csv(ROOT / "01_data_raw" / "holdout" / "E-MTAB-16655" / "RawCounts.csv")
    sample_cols = [c for c in norm.columns if re.match(r"^(UI|WT|dotA|A|B|C|D|E|F|I|FGHI|O|P|R|S)\d+$", c)]
    raw_sample_cols = [c for c in raw.columns if c != "Geneid"]
    norm["gene_length"] = (norm["end_position"] - norm["start_position"] + 1).clip(lower=1)
    norm["baseline_ui_mean_logcpm"] = norm[["UI2", "UI3", "UI4"]].mean(axis=1)
    norm["across_sample_sd_logcpm"] = norm[sample_cols].std(axis=1)
    norm["infection_abs_logfc"] = (norm[["WT2", "WT3", "WT4"]].mean(axis=1) - norm[["UI2", "UI3", "UI4"]].mean(axis=1)).abs()
    raw["mean_raw_count"] = raw[raw_sample_cols].mean(axis=1)
    features = norm.merge(raw[["Geneid", "mean_raw_count"]], on="Geneid", how="left")
    features = features.rename(columns={"external_gene_name": "human_symbol"})
    features = features.dropna(subset=["human_symbol"]).copy()
    features["human_symbol"] = features.human_symbol.astype(str).str.upper()
    features["log_gene_length"] = np.log1p(features.gene_length)
    features["log_mean_raw_count"] = np.log1p(features.mean_raw_count.fillna(0))
    features = features.sort_values("mean_raw_count", ascending=False).drop_duplicates("human_symbol")
    match_cols = [
        "baseline_ui_mean_logcpm",
        "across_sample_sd_logcpm",
        "infection_abs_logfc",
        "log_gene_length",
        "log_mean_raw_count",
    ]
    return features, match_cols


def robust_scale(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    scaled = frame[columns].astype(float).copy()
    for column in columns:
        median = scaled[column].median()
        mad = np.median(np.abs(scaled[column] - median))
        scale = 1.4826 * mad
        if not np.isfinite(scale) or scale <= 0:
            scale = scaled[column].std()
        scaled[column] = (scaled[column] - median) / (scale if scale > 0 else 1.0)
    return scaled


def generate_matched_sets(features: pd.DataFrame, panel: list[str], match_cols: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    features = features.set_index("human_symbol", drop=False)
    panel = [g for g in panel if g in features.index]
    universe = features.loc[~features.index.isin(panel)].copy()
    all_frame = pd.concat([features.loc[panel], universe], axis=0)
    scaled = robust_scale(all_frame, match_cols)
    pools: dict[str, list[str]] = {}
    pool_rows = []
    for gene in panel:
        query = features.loc[gene]
        same_type = universe[universe.gene_biotype.eq(query.gene_biotype)]
        candidates = same_type if len(same_type) >= MATCH_POOL else universe
        distance = ((scaled.loc[candidates.index, match_cols] - scaled.loc[gene, match_cols]) ** 2).sum(axis=1) ** 0.5
        nearest = distance.nsmallest(max(MATCH_POOL * 4, MATCH_POOL))
        pools[gene] = nearest.index.tolist()
        for rank, (candidate, dist) in enumerate(nearest.iloc[:10].items(), start=1):
            pool_rows.append({"panel_gene": gene, "matched_gene": candidate, "rank": rank, "distance": float(dist)})
    matched = np.empty((N_MATCHED_NULL, len(panel)), dtype=object)
    for iteration in range(N_MATCHED_NULL):
        used: set[str] = set()
        for j in RNG.permutation(len(panel)):
            gene = panel[j]
            choices = [x for x in pools[gene][:MATCH_POOL] if x not in used]
            if not choices:
                choices = [x for x in pools[gene] if x not in used]
            if not choices:
                raise RuntimeError(f"No unused match available for {gene}")
            chosen = str(RNG.choice(choices))
            matched[iteration, j] = chosen
            used.add(chosen)
    pool_df = pd.DataFrame(pool_rows)
    pool_df.to_csv(OUT / "d5_matched_null_nearest_genes.tsv", sep="\t", index=False)
    return matched, features.loc[panel].reset_index(drop=True)


def matched_null_audit() -> tuple[pd.DataFrame, dict]:
    holdout = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz", sep="\t", low_memory=False)
    frozen = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
    panel = sorted(frozen.loc[frozen.contrast_id.eq("WT_vs_UI"), "human_symbol"].astype(str).str.upper().unique())
    features, match_cols = build_matching_features()
    available = set(holdout.human_symbol.astype(str).str.upper())
    panel = [g for g in panel if g in available and g in set(features.human_symbol)]
    matched_sets, panel_features = generate_matched_sets(features, panel, match_cols)
    panel_features.to_csv(OUT / "d5_frozen_panel_matching_features.tsv", sep="\t", index=False)

    multi = ["A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "E_vs_WT", "F_vs_WT", "FGHI_vs_WT", "I_vs_WT", "O_vs_WT", "P_vs_WT", "R_vs_WT", "S_vs_WT"]
    passing = {"A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "FGHI_vs_WT", "S_vs_WT"}
    rows = []
    global_values = []
    panel_values = []
    for contrast in multi:
        data = holdout.loc[holdout.contrast_id.eq(contrast)].copy()
        data["human_symbol"] = data.human_symbol.astype(str).str.upper()
        data["signed_z"] = signed_z(data)
        zmap = data.drop_duplicates("human_symbol").set_index("human_symbol").signed_z
        observed = float(zmap.reindex(panel).abs().mean())
        null = np.array([zmap.reindex(row).abs().mean() for row in matched_sets], dtype=float)
        null = null[np.isfinite(null)]
        empirical_p = float((1 + np.sum(null >= observed)) / (len(null) + 1))
        global_mean = float(data.signed_z.abs().mean())
        global_values.append(global_mean)
        panel_values.append(observed)
        rows.append(
            {
                "contrast": contrast,
                "passes_frozen": contrast in passing,
                "panel_n": len(panel),
                "panel_mean_abs_z": observed,
                "global_mean_abs_z": global_mean,
                "matched_null_mean": float(np.mean(null)),
                "matched_null_sd": float(np.std(null, ddof=1)),
                "matched_null_p95": float(np.quantile(null, 0.95)),
                "matched_empirical_p": empirical_p,
            }
        )
    result = pd.DataFrame(rows)
    result["matched_empirical_fdr"] = multipletests(result.matched_empirical_p, method="fdr_bh")[1]
    result.to_csv(OUT / "d5_global_response_amplitude_matched_null.tsv", sep="\t", index=False)
    rho, rho_p = stats.spearmanr(panel_values, global_values)
    summary = {
        "analysis_name": "global response-amplitude confounding sensitivity analysis",
        "explicit_nonclaim": "No bacterial burden was measured in E-MTAB-16655; global host-transcriptional disturbance is not a bacterial-load measurement.",
        "n_official_conditions": 15,
        "n_multi_effector_conditions": 12,
        "n_panel_genes_matched": len(panel),
        "n_matched_null_sets": N_MATCHED_NULL,
        "matching_variables": match_cols,
        "panel_vs_global_spearman": {"rho": float(rho), "p": float(rho_p)},
        "contrasts_panel_above_matched_null_p05": int((result.matched_empirical_p < 0.05).sum()),
        "contrasts_panel_above_matched_null_fdr05": int((result.matched_empirical_fdr < 0.05).sum()),
        "replication_evidence_status": "Originating investigators report no gross intracellular-replication defect in this mutant library, but quantitative strain-level burden at the 18 h RNA-seq time point is not present in the deposited dataset and remains a residual uncertainty.",
    }
    (OUT / "d5_biological_confound_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return result, summary


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    conditions = audit_sdrf()
    matched, summary = matched_null_audit()
    print(conditions[["official_condition", "n_biological_samples", "mapping_status"]].to_string(index=False))
    print(matched[["contrast", "panel_mean_abs_z", "matched_null_mean", "matched_empirical_p", "matched_empirical_fdr"]].to_string(index=False))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
