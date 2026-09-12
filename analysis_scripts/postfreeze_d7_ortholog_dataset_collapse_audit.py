from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"


def collapse_gene_dataset(group: pd.DataFrame) -> pd.Series:
    # A dataset is one independent study unit. Multiple mouse genes mapped to the
    # same human symbol are collapsed, never counted as independent studies.
    signed_z = float(group.signed_z.mean())
    logfc = float(group.logFC.median())
    return pd.Series(
        {
            "species": group.species.iloc[0],
            "signed_z": signed_z,
            "logFC": logfc,
            "adj.P.Val_min_descriptive": float(group["adj.P.Val"].min()),
            "study_fdr_hit_any": bool((group["adj.P.Val"] < 0.05).any()),
            "n_original_rows_collapsed": int(len(group)),
            "original_symbols": ";".join(sorted(set(group.original_symbol.astype(str)))),
        }
    )


def summarize_gene(group: pd.DataFrame) -> pd.Series:
    direction = "UP" if (group.logFC > 0).sum() >= (group.logFC < 0).sum() else "DOWN"
    concordant = int((group.logFC > 0).sum()) if direction == "UP" else int((group.logFC < 0).sum())
    z = float(group.signed_z.sum() / np.sqrt(len(group)))
    return pd.Series(
        {
            "n_studies": int(group.dataset.nunique()),
            "n_human": int(group.loc[group.species.eq("HUMAN"), "dataset"].nunique()),
            "n_mouse": int(group.loc[group.species.eq("MOUSE"), "dataset"].nunique()),
            "consensus_direction": direction,
            "direction_concordance": float(concordant / len(group)),
            "median_logFC": float(group.logFC.median()),
            "combined_z": z,
            "combined_p": float(2 * (0.5 * np.math.erfc(abs(z) / np.sqrt(2)))) if hasattr(np, "math") else np.nan,
            "n_study_fdr05": int(group.study_fdr_hit_any.sum()),
            "datasets": ";".join(sorted(group.dataset.unique())),
            "n_rows_collapsed_total": int(group.n_original_rows_collapsed.sum()),
            "n_dataset_units_with_multiple_mappings": int((group.n_original_rows_collapsed > 1).sum()),
        }
    )


def main() -> int:
    evidence = pd.read_csv(
        ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
        sep="\t",
        compression=None,
        low_memory=False,
    )
    original = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
    collapsed = (
        evidence.groupby(["human_symbol", "dataset"], sort=False, as_index=False)
        .apply(collapse_gene_dataset, include_groups=False)
        .reset_index(drop=True)
    )
    collapsed.to_csv(OUT / "d7_gene_dataset_collapsed_evidence.tsv", sep="\t", index=False)

    summary = collapsed.groupby("human_symbol", sort=False).apply(summarize_gene, include_groups=False).reset_index()
    # scipy-free two-sided normal tail for compatibility.
    from scipy.stats import norm

    summary["combined_p"] = 2 * norm.sf(summary.combined_z.abs())
    summary["combined_fdr"] = multipletests(summary.combined_p, method="fdr_bh")[1]
    summary["stable_cross_cohort_corrected"] = (
        summary.n_studies.ge(3)
        & summary.n_human.ge(2)
        & summary.direction_concordance.ge(0.75)
        & summary.combined_fdr.lt(0.05)
        & summary.n_study_fdr05.ge(1)
    )
    summary.to_csv(OUT / "d7_primary_infection_convergence_dataset_collapsed.tsv", sep="\t", index=False)

    compare = original[["human_symbol", "stable_cross_cohort", "combined_z", "combined_fdr", "n_studies", "n_study_fdr05"]].merge(
        summary[
            [
                "human_symbol",
                "stable_cross_cohort_corrected",
                "combined_z",
                "combined_fdr",
                "n_studies",
                "n_study_fdr05",
                "n_dataset_units_with_multiple_mappings",
            ]
        ],
        on="human_symbol",
        how="outer",
        suffixes=("_original", "_corrected"),
    )
    compare["stable_status_changed"] = compare.stable_cross_cohort.fillna(False).astype(bool) != compare.stable_cross_cohort_corrected.fillna(False).astype(bool)
    compare.to_csv(OUT / "d7_original_vs_dataset_collapsed_convergence.tsv", sep="\t", index=False)

    tiers = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")
    candidate_compare = tiers[["human_symbol", "final_tier", "stable_cross_cohort"]].merge(
        summary[["human_symbol", "stable_cross_cohort_corrected", "combined_z", "combined_fdr", "n_studies", "n_study_fdr05"]],
        on="human_symbol",
        how="left",
    )
    candidate_compare["expression_stability_changed"] = candidate_compare.stable_cross_cohort.fillna(False).astype(bool) != candidate_compare.stable_cross_cohort_corrected.fillna(False).astype(bool)
    candidate_compare.to_csv(OUT / "d7_frozen_candidate_dataset_collapse_sensitivity.tsv", sep="\t", index=False)

    report = {
        "duplicate_human_gene_dataset_pairs": int((evidence.groupby(["human_symbol", "dataset"]).size() > 1).sum()),
        "extra_rows_that_were_not_independent_studies": int((evidence.groupby(["human_symbol", "dataset"]).size() - 1).clip(lower=0).sum()),
        "original_stable_n": int(original.stable_cross_cohort.astype(bool).sum()),
        "dataset_collapsed_stable_n": int(summary.stable_cross_cohort_corrected.sum()),
        "stable_status_changed_n": int(compare.stable_status_changed.sum()),
        "frozen_candidate_expression_flag_changed_n": int(candidate_compare.expression_stability_changed.sum()),
        "frozen_candidate_changes": candidate_compare.loc[candidate_compare.expression_stability_changed, ["human_symbol", "final_tier", "stable_cross_cohort", "stable_cross_cohort_corrected"]].to_dict("records"),
        "aggregation_rule": "One human_symbol x dataset unit. Multiple original mouse symbols mapped to one human symbol are represented by mean signed Z and median logFC; the dataset contributes once.",
        "claim_boundary": "This is a post-freeze bug/sensitivity audit. It does not optimize or replace the frozen holdout; any changed expression flags must be transparently reported.",
    }
    (OUT / "d7_dataset_collapse_audit_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
