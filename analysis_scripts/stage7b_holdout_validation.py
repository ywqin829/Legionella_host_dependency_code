from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz"
PANEL = ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
RULES = ROOT / "08_evidence_atlas" / "HOLDOUT_VALIDATION_RULES_FROZEN.json"
OUT = ROOT / "09_holdout_validation"
CP = ROOT / "checkpoints"
RNG = np.random.default_rng(42)
N_PERM = 10_000


def empirical_greater(values: np.ndarray, size: int, observed: float, absolute: bool = False) -> float:
    null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        sample = RNG.choice(values, size=size, replace=False)
        null[i] = np.mean(np.abs(sample)) if absolute else np.mean(sample)
    return float((1 + np.sum(null >= observed)) / (N_PERM + 1))


def signed_z(frame: pd.DataFrame) -> np.ndarray:
    p = np.clip(frame.PValue.to_numpy(float), 1e-300, 1 - 1e-15)
    return np.sign(frame.logFC.to_numpy(float)) * norm.isf(p / 2)


def route_direction(text: str, route: str) -> str | None:
    for item in str(text).split(";"):
        if item.startswith(route + ":"):
            return item.split(":", 1)[1]
    return None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    result = pd.read_csv(RESULTS, sep="\t")
    panel = pd.read_csv(PANEL, sep="\t")
    rules = json.loads(RULES.read_text(encoding="utf-8"))
    primary = panel.loc[panel.freeze_role.eq("PRIMARY")].copy()
    all_genes = set(result.human_symbol)
    panel = panel.loc[panel.human_symbol.isin(all_genes)].copy()
    primary = primary.loc[primary.human_symbol.isin(all_genes)].copy()

    individual_rows = []
    for contrast, d in result.groupby("contrast_id", sort=False):
        hit = d.merge(panel, on="human_symbol", how="inner")
        for _, r in hit.iterrows():
            expected = None
            if contrast == "WT_vs_UI" and r.expected_wt_vs_uninfected in {"UP", "DOWN"}:
                expected = r.expected_wt_vs_uninfected
            elif contrast == "WT_vs_dotA":
                expected = route_direction(r.predefined_mechanism_route_directions, "DOT_ICM")
            observed = "UP" if r.logFC > 0 else "DOWN" if r.logFC < 0 else "ZERO"
            individual_rows.append({
                "contrast_id": contrast, "human_symbol": r.human_symbol, "freeze_role": r.freeze_role,
                "evidence_tier": r.evidence_tier, "logFC": r.logFC, "PValue": r.PValue, "FDR": r.FDR,
                "expected_direction": expected or "NOT_PREDEFINED", "observed_direction": observed,
                "direction_match": (observed == expected) if expected else np.nan,
                "fitness_artifact_compatible": r.fitness_artifact_compatible,
            })
    individual = pd.DataFrame(individual_rows)
    individual.to_csv(OUT / "holdout_frozen_candidate_results.tsv", sep="\t", index=False)

    set_rows = []
    for contrast, d in result.groupby("contrast_id", sort=False):
        d = d.copy(); d["z"] = signed_z(d)
        zmap = d.set_index("human_symbol").z
        universe = d.z.to_numpy(float)
        available_panel = panel.loc[panel.human_symbol.isin(zmap.index)]

        if contrast == "WT_vs_UI":
            for set_name, genes in [("PRIMARY4_DIRECTIONAL", primary), ("PANEL_DIRECTIONAL", available_panel)]:
                genes = genes.loc[genes.expected_wt_vs_uninfected.isin(["UP", "DOWN"])]
                oriented = np.array([(1 if r.expected_wt_vs_uninfected == "UP" else -1) * zmap[r.human_symbol] for _, r in genes.iterrows()])
                observed = float(np.mean(oriented))
                p = empirical_greater(universe, len(oriented), observed, absolute=False)
                set_rows.append({"endpoint_family": "INFECTION_DIRECTIONAL", "contrast_id": contrast, "set_name": set_name, "set_n": len(oriented), "score": observed, "empirical_p": p, "score_definition": "mean direction-oriented signed Z"})
        elif contrast == "WT_vs_dotA":
            dot = available_panel.copy()
            dot["expected"] = dot.predefined_mechanism_route_directions.map(lambda x: route_direction(x, "DOT_ICM"))
            dot = dot.loc[dot.expected.isin(["UP", "DOWN"])]
            oriented = np.array([(1 if r.expected == "UP" else -1) * zmap[r.human_symbol] for _, r in dot.iterrows()])
            observed = float(np.mean(oriented))
            p = empirical_greater(universe, len(oriented), observed, absolute=False)
            set_rows.append({"endpoint_family": "DOT_ICM_DIRECTIONAL", "contrast_id": contrast, "set_name": "DOT_ICM_PREDEFINED", "set_n": len(oriented), "score": observed, "empirical_p": p, "score_definition": "mean direction-oriented signed Z"})
        else:
            vals = zmap.loc[available_panel.human_symbol].to_numpy(float)
            observed = float(np.mean(np.abs(vals)))
            p = empirical_greater(universe, len(vals), observed, absolute=True)
            set_rows.append({"endpoint_family": "MULTI_EFFECTOR_TWO_SIDED", "contrast_id": contrast, "set_name": "FROZEN_PANEL56", "set_n": len(vals), "score": observed, "empirical_p": p, "score_definition": "mean absolute signed Z"})

    sets = pd.DataFrame(set_rows)
    sets["empirical_fdr_within_family"] = np.nan
    for family, idx in sets.groupby("endpoint_family").groups.items():
        sets.loc[idx, "empirical_fdr_within_family"] = multipletests(sets.loc[idx, "empirical_p"], method="fdr_bh")[1]
    sets["passes_frozen_endpoint"] = sets.empirical_fdr_within_family < 0.05
    sets.to_csv(OUT / "holdout_frozen_set_tests.tsv", sep="\t", index=False)

    ui_primary = individual.loc[(individual.contrast_id == "WT_vs_UI") & (individual.freeze_role == "PRIMARY")]
    primary_direction_matches = int(ui_primary.direction_match.fillna(False).sum())
    primary_genomewide_hits = int((ui_primary.FDR < 0.05).sum())
    passed = sets.loc[sets.passes_frozen_endpoint]
    status = "PASS" if len(passed) else "FAIL_DISCOVERY_ONLY"
    stats = {
        "frozen_panel_n": int(rules["candidate_panel_n"]), "panel_measured_n": int(len(panel)),
        "primary_frozen_n": int(rules["primary_n"]), "primary_measured_n": int(len(primary)),
        "permutations_per_set_test": N_PERM, "set_tests_n": int(len(sets)),
        "set_tests_pass_fdr05_n": int(len(passed)), "passing_set_tests": passed[["contrast_id", "set_name"]].to_dict("records"),
        "primary_wt_ui_direction_matches": primary_direction_matches,
        "primary_wt_ui_genomewide_fdr05": primary_genomewide_hits,
        "holdout_validation_status": status,
        "claim_boundary": "Validation concerns host transcriptional response only; it does not establish bacterial replication dependency or direct effector binding",
    }
    (OUT / "stage7b_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = ["# Stage 7B frozen external-holdout validation", "",
        f"- Fixed candidates measured: {stats['panel_measured_n']}/{stats['frozen_panel_n']}",
        f"- Primary candidates measured: {stats['primary_measured_n']}/{stats['primary_frozen_n']}",
        f"- Primary WT-vs-uninfected direction matches: {primary_direction_matches}/{len(ui_primary)}",
        f"- Primary genes with genome-wide FDR < 0.05 in WT-vs-uninfected: {primary_genomewide_hits}/{len(ui_primary)}",
        f"- Fixed set tests passing within-family FDR < 0.05: {len(passed)}/{len(sets)}",
        f"- Holdout status: **{status}**", "",
        "Candidate identities, directions, contrasts, and multiplicity rules were frozen before the expression matrices were downloaded. Multi-effector deletion groups were tested separately. Null results for single or grouped perturbations are not interpreted as proof of no function where effector redundancy is a registered boundary.",
        "This validation concerns host transcriptional responses and cannot establish bacterial replication dependence, direct molecular binding, or therapeutic causality."]
    (OUT / "STAGE7B_FROZEN_HOLDOUT_VALIDATION.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage7b_holdout_validation.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
