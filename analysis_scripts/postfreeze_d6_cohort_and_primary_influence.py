from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
KEY = ["RAB1A", "RAB10", "ACTR3", "VAPA"]


def stouffer(values: pd.Series) -> float:
    x = values.dropna().to_numpy(float)
    return float(x.sum() / np.sqrt(len(x))) if len(x) else np.nan


def direction_summary(frame: pd.DataFrame) -> tuple[str, float]:
    signs = np.sign(frame.logFC.dropna().to_numpy(float))
    if not len(signs):
        return "NA", np.nan
    up = float(np.mean(signs > 0))
    down = float(np.mean(signs < 0))
    return ("UP", up) if up >= down else ("DOWN", down)


def cohort_layers() -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence = pd.read_csv(
        ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
        sep="\t",
        compression=None,
        low_memory=False,
    )
    convergence = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
    evidence["study_fdr05"] = evidence["adj.P.Val"] < 0.05
    counts = (
        evidence.groupby("human_symbol", as_index=False)
        .agg(
            n_studies_observed=("dataset", "nunique"),
            n_human_observed=("species", lambda x: int((x == "HUMAN").sum())),
            n_study_fdr05_observed=("study_fdr05", "sum"),
        )
    )
    hit_counts = evidence.loc[evidence.study_fdr05].groupby("human_symbol")["dataset"].nunique()
    human_counts = evidence.loc[evidence.species.eq("HUMAN")].groupby("human_symbol")["dataset"].nunique()
    counts["n_study_fdr05_observed"] = counts.human_symbol.map(hit_counts).fillna(0).astype(int)
    counts["n_human_observed"] = counts.human_symbol.map(human_counts).fillna(0).astype(int)
    layers = convergence.merge(counts, on="human_symbol", how="left")
    layers["meta_convergent_response"] = layers.stable_cross_cohort.astype(bool)
    layers["multistudy_replicated_core"] = layers.n_study_fdr05_observed.ge(3)
    layers["recommended_label"] = np.select(
        [layers.multistudy_replicated_core, layers.meta_convergent_response],
        ["MULTISTUDY_REPLICATED_CORE_GE3", "META_CONVERGENT_RESPONSE"],
        default="NOT_IN_EITHER_LAYER",
    )
    layers.to_csv(OUT / "d6_expression_evidence_layers.tsv", sep="\t", index=False)

    key_rows = []
    for gene in KEY:
        gene_ev = evidence.loc[evidence.human_symbol.eq(gene)].copy()
        base = layers.loc[layers.human_symbol.eq(gene)]
        if base.empty:
            continue
        base = base.iloc[0]
        human = gene_ev.loc[gene_ev.species.eq("HUMAN")]
        human_direction, human_concordance = direction_summary(human)
        leave_one = []
        for dataset in sorted(gene_ev.dataset.unique()):
            remaining = gene_ev.loc[~gene_ev.dataset.eq(dataset)]
            direction, concordance = direction_summary(remaining)
            leave_one.append(
                f"drop {dataset}: z={stouffer(remaining.signed_z):.3f},dir={direction},conc={concordance:.2f},n={remaining.dataset.nunique()}"
            )
        key_rows.append(
            {
                "candidate": gene,
                "meta_convergent_response": bool(base.meta_convergent_response),
                "multistudy_replicated_core_ge3": bool(base.multistudy_replicated_core),
                "n_studies": int(base.n_studies),
                "n_human": int(base.n_human),
                "n_study_fdr05": int(base.n_study_fdr05),
                "combined_z_all": float(base.combined_z),
                "direction_concordance_all": float(base.direction_concordance),
                "human_only_n": int(human.dataset.nunique()),
                "human_only_combined_z": stouffer(human.signed_z),
                "human_only_direction": human_direction,
                "human_only_direction_concordance": human_concordance,
                "leave_one_dataset_out": " | ".join(leave_one),
            }
        )
    key = pd.DataFrame(key_rows)
    key.to_csv(OUT / "d6_key_candidate_cohort_fragility.tsv", sep="\t", index=False)
    return layers, key


def primary_influence() -> pd.DataFrame:
    individual = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
    data = individual.loc[(individual.contrast_id.eq("WT_vs_UI")) & individual.human_symbol.isin(KEY)].copy()
    p = np.clip(data.PValue.to_numpy(float), 1e-300, 1 - 1e-15)
    data["direction_oriented_z"] = np.where(
        data.expected_direction.eq("UP"),
        1.0,
        -1.0,
    ) * np.sign(data.logFC.to_numpy(float)) * stats.norm.isf(p / 2)
    all_mean = float(data.direction_oriented_z.mean())
    rows = []
    for _, row in data.sort_values("human_symbol").iterrows():
        remaining = data.loc[~data.human_symbol.eq(row.human_symbol), "direction_oriented_z"]
        leave_mean = float(remaining.mean())
        rows.append(
            {
                "candidate": row.human_symbol,
                "logFC": float(row.logFC),
                "FDR": float(row.FDR),
                "direction_oriented_z": float(row.direction_oriented_z),
                "all_four_mean_z": all_mean,
                "leave_candidate_out_mean_z": leave_mean,
                "candidate_above_group_mean": bool(row.direction_oriented_z > all_mean),
                "interpretation": "post-freeze descriptive influence only; does not alter the prespecified four-gene test",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "d6_primary4_influence_decomposition.tsv", sep="\t", index=False)
    return result


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    layers, key = cohort_layers()
    influence = primary_influence()
    summary = {
        "meta_convergent_response_n": int(layers.meta_convergent_response.sum()),
        "multistudy_replicated_core_ge3_n_all_genes": int(layers.multistudy_replicated_core.sum()),
        "multistudy_replicated_core_ge3_n_within_meta_convergent": int(
            (layers.meta_convergent_response & layers.multistudy_replicated_core).sum()
        ),
        "key_candidate_summary": key.to_dict("records"),
        "primary4_influence_statement": "The prespecified four-gene set remains a valid frozen test, but its directional signal is heterogeneous and is driven predominantly by RAB1A and RAB10.",
    }
    (OUT / "d6_cohort_and_primary_influence_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "key_candidate_summary"}, indent=2))
    print(key.to_string(index=False))
    print(influence.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
