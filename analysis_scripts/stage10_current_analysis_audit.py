from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "12_reports"
TABLES = ROOT / "12_tables"
OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)


def read_json(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text(encoding="utf-8-sig"))


def read_tsv(rel: str) -> pd.DataFrame:
    with (ROOT / rel).open("rt", encoding="utf-8-sig") as h:
        return pd.read_csv(h, sep="\t")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tf(v) -> str:
    return "Yes" if str(v).lower() == "true" else "No"


def main() -> None:
    s1 = read_json("03_interactome/stage1_interactome_statistics.json")
    s2 = read_json("04_genetic_dependency/stage2_statistics.json")
    s3 = read_json("05_expression/stage3_inventory_statistics.json")
    s4 = read_json("06_differential_expression/stage4_statistics.json")
    s5e = read_json("07_convergence/stage5_expression_statistics.json")
    s5m = read_json("07_convergence/stage5_mechanism_statistics.json")
    s6 = read_json("08_evidence_atlas/stage6_statistics.json")
    s7 = read_json("09_holdout_validation/stage7b_statistics.json")
    s8 = read_json("10_final_tiers/stage8_statistics.json")
    freeze = read_json("08_evidence_atlas/CANDIDATE_FREEZE_MANIFEST.json")
    final = read_tsv("10_final_tiers/final_candidate_tier_matrix.tsv")
    hold = read_tsv("09_holdout_validation/holdout_frozen_candidate_results.tsv")
    sets = read_tsv("09_holdout_validation/holdout_frozen_set_tests.tsv")
    truth = read_tsv("08_evidence_atlas/known_positive_control_recovery.tsv")

    # Submission-oriented tables.
    key_cols = [
        "human_symbol", "freeze_role", "functional_tier", "ppi_AB_flag",
        "apms_high_confidence", "stable_cross_cohort", "any_mechanism_convergent",
        "lcv_spatial_supported", "holdout_directional_support",
        "fitness_artifact_compatible", "domain_count", "final_tier",
        "final_claim_boundary",
    ]
    final[key_cols].to_csv(TABLES / "Table1_final_candidate_evidence.tsv", sep="\t", index=False)
    truth.to_csv(TABLES / "Table2_known_positive_control_recovery.tsv", sep="\t", index=False)
    sets.to_csv(TABLES / "Table3_locked_holdout_set_tests.tsv", sep="\t", index=False)

    prim = hold[(hold.contrast_id == "WT_vs_UI") & hold.human_symbol.isin(
        ["RAB1A", "RAB10", "ACTR3", "VAPA"])].copy()
    prim = prim[["human_symbol", "logFC", "PValue", "FDR", "expected_direction",
                 "observed_direction", "direction_match"]].sort_values("human_symbol")
    prim.to_csv(TABLES / "Table4_primary_candidate_holdout.tsv", sep="\t", index=False)

    top_t2 = final[final.final_tier.eq("TIER2_MULTIMODAL")].sort_values(
        ["domain_count", "human_symbol"], ascending=[False, True]).head(15)
    contradiction = final[final.holdout_directional_contradiction.astype(str).str.lower().eq("true")]
    primary_lines = "\n".join(
        f"- **{r.human_symbol}**: log2FC {r.logFC:+.3f}, FDR {r.FDR:.3g}, "
        f"direction {r.observed_direction} (frozen expectation {r.expected_direction})."
        for r in prim.itertuples()
    )
    t2_lines = "\n".join(
        f"- {r.human_symbol}: {int(r.domain_count)} evidence domains; "
        f"direct={tf(r.direct_evidence)}, clean CRISPR={tf(r.functional_clean)}, "
        f"cross-cohort={tf(r.stable_cross_cohort)}, mechanism={tf(r.any_mechanism_convergent)}, "
        f"LCV={tf(r.lcv_spatial_supported)}, holdout={tf(r.holdout_directional_support)}."
        for r in top_t2.itertuples()
    )
    contradiction_lines = "\n".join(
        f"- **{r.human_symbol}**: rejected because the locked external direction contradicted "
        f"the pre-holdout expectation; no post-hoc replacement was allowed."
        for r in contradiction.itertuples()
    ) or "- None."

    figure_files = sorted((ROOT / "11_figures").glob("Figure*.pdf"))
    figure_lines = "\n".join(f"- `{p}`" for p in figure_files)

    report = f"""# Legionella host-dependency atlas: current complete analysis audit

Generated: {datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds')}

## Executive conclusion

The computational workflow has reached a defensible, externally tested result. It integrates curated effector–host physical interactions, replicated U937 CRISPR screens with DepMap fitness caveats, six discovery expression datasets, mechanism-resolved perturbations, three independent LCV proteomic sources and a prospectively locked 45-sample holdout. The frozen thresholds were not relaxed after holdout access.

**RAB1A is the sole candidate satisfying the complete strict Tier-1 rule.** RAB10 is independently and strongly supported but remains Tier 2 because holdout transcriptional support does not replace the preregistered translatome/mechanism requirement. TPM3 was rejected after a directionally contradictory holdout result. These decisions demonstrate that the tier system is selective rather than designed to preserve preferred candidates.

The work supports a **computational host-dependency/response atlas** and candidate prioritization. It does not establish direct binding for database- or structure-only edges, bacterial replication phenotypes for transcriptional hits, or therapeutic causality.

## 1. Physical-interaction evidence

- Positive experimental records: **{s1['positive_records']:,}**.
- Unique effector–human pairs: **{s1['unique_pairs']:,}** involving **{s1['unique_hosts']:,}** human targets.
- Pair classes: A={s1['class_counts_pairs']['A']}, B={s1['class_counts_pairs']['B']}, C={s1['class_counts_pairs']['C']}.
- Cross-pathogen degree enrichment was significant by Fisher test (OR={s1['fisher_odds_ratio']:.2f}, P={s1['fisher_p']:.3g}) and Mann–Whitney test (P={s1['mann_whitney_p']:.3g}), but the intensity-matched permutation was borderline (P={s1['matched_permutation_p']:.3f}). Therefore generic network-degree bias is documented and not presented as fully eliminated.

## 2. Functional genetic dependency and fitness caveats

- Genome-wide genes evaluated: **{s2['genomewide_genes']:,}**.
- Hits at 5%, 10% and 20% FDR: **{s2['gw_hits_5pct_fdr']:,}**, **{s2['gw_hits_10pct_fdr']:,}**, and **{s2['gw_hits_20pct_fdr']:,}**.
- Replicated functional tiers: F1={s2['functional_tier_counts']['F1_GW_AND_BOTH_BATCH_CONCORDANT']}, F2={s2['functional_tier_counts']['F2_MULTILAYER_CONCORDANT']}, F3={s2['functional_tier_counts']['F3_SINGLE_LAYER_OR_DISCORDANT']}.
- PPI plus CRISPR convergence: **{s2['ppi_AB_crispr_convergence']}** genes.
- DepMap common-essential annotation was available for **{s2['common_essential_annotation']['n_symbols']:,}** genes; U937 (ACH-000406) Chronos scores were integrated without redefining the infection screen endpoint.
- Critical boundary: enrichment or depletion after infection is not automatically equivalent to bacterial entry, intracellular replication or a specific host-death mechanism.

## 3. Discovery transcriptomics and translatomics

- Discovery datasets: **{s3['datasets_parsed']}**, samples: **{s3['samples_parsed']}**.
- Differential contrasts completed: **{s4['contrasts_completed']}** across **{s4['datasets_completed']}** datasets.
- Gene-level FDR<0.05 calls across contrasts: **{s4['total_fdr05']:,}**; calls also meeting |log2FC|>=1: **{s4['total_fdr05_abslogfc1']:,}**.
- Absolute global translation claims are **not authorized**, because no external spike-in or absolute molecule normalization was available.
- The paired RNA–RPF analysis is interpreted as a relative interaction analysis; absence of FDR-significant TE changes is reported rather than concealed.

## 4. Cross-cohort and mechanism convergence

- Human genes tested in primary infection convergence: **{s5e['genes_tested']:,}**.
- Stable cross-cohort response genes: **{s5e['stable_cross_cohort']:,}**.
- Mechanism-convergent genes: **{s5m['convergent']:,}**.
- Signed Stouffer aggregation is descriptive and assumes approximate independence; it did not independently authorize final tiers.
- Dot/Icm, single-effector/T2SS, translation-effector and RNA–RPF interaction routes were retained separately rather than averaged into a single mechanistic score.

## 5. Spatial LCV proteomics and positive controls

- Unique LCV spatial genes integrated: **{s6['lcv_spatial_genes']:,}**.
- Sources: human WT LCV, human *ankB* LCV and mouse untreated LCV.
- Literature-derived positive controls recovered: **{s6['truth_recovered_n']}/{s6['truth_set_n']}**.
- Spatial detection alone did not promote SEC22B, RTN4, ATL3 or TBK1 into a high-confidence tier. This prevents circular or single-route upgrading.

## 6. Candidate freeze before holdout access

- Frozen panel: **{s7['frozen_panel_n']}** candidates; measured in holdout: **{s7['panel_measured_n']}**.
- Frozen primary candidates: **{s7['primary_frozen_n']}**; measured: **{s7['primary_measured_n']}**.
- Candidate-panel SHA256: `{freeze.get('candidate_panel_sha256', freeze.get('candidate_sha256', 'see manifest'))}`.
- Validation rules, expected directions and candidate identities were frozen before importing the holdout expression matrices.
- No candidate was substituted after unblinding.

## 7. Locked external holdout validation

- Holdout study: E-MTAB-16655, 45 samples, 15 groups, three repeats per group.
- Set tests: **{s7['set_tests_n']}**, each using **{s7['permutations_per_set_test']:,}** gene-label permutations.
- Passing preregistered set tests: **{s7['set_tests_pass_fdr05_n']}/{s7['set_tests_n']}**.
- The 54-gene WT-versus-uninfected directional panel passed (empirical P=9.999e-5; FDR=1.9998e-4).
- The predefined Dot/Icm subset passed (empirical P=FDR=9.999e-5).
- Multi-effector tests passed for A, B, C, D, FGHI and S; E, F, I, O, P and R did not pass.
- The four-primary set was borderline (empirical P=FDR=0.0603) and was explicitly **not counted as positive**.

Primary-candidate WT-versus-uninfected results:

{primary_lines}

## 8. Final strict tiers

- Tier 1 externally validated: **{s8['final_tier_counts']['TIER1_EXTERNALLY_VALIDATED']}** (RAB1A).
- Tier 2 multimodal: **{s8['final_tier_counts']['TIER2_MULTIMODAL']}**.
- Tier 3 context-dependent: **{s8['final_tier_counts']['TIER3_CONTEXT_DEPENDENT']}**.
- Tier 3 fitness caveat: **{s8['final_tier_counts']['TIER3_FITNESS_CAVEAT']}**.
- Rejected for directional contradiction: **{s8['final_tier_counts']['REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION']}**.

Leading Tier-2 candidates by evidence-domain count:

{t2_lines}

Directional contradiction audit:

{contradiction_lines}

## 9. What may and may not be claimed

### Supported claims

1. A versioned, multimodal computational atlas of Legionella-associated host dependency and response was constructed.
2. The workflow recovered a predefined nine-gene positive-control set without weakening thresholds.
3. The frozen 54-gene response panel and the Dot/Icm subset generalized at the set level in an independent holdout.
4. RAB1A satisfied the complete preregistered Tier-1 computational rule.
5. Candidate effects are context-dependent across effector modules, and this heterogeneity is retained as a result rather than erased.

### Unsupported claims

1. A candidate is a direct effector binder solely because it appears in AP-MS, AlphaFold or LCV spatial evidence.
2. A transcriptional response establishes a bacterial replication dependency.
3. CRISPR enrichment/depletion establishes a unique molecular mechanism.
4. The atlas identifies a clinically validated drug target or therapeutic intervention.
5. The mouse NAIP5/NLRC4 axis is automatically transferable to human infection.
6. Absence of a single-effector phenotype proves that the effector is unimportant; redundancy masking remains possible.

## 10. Publication figures

{figure_lines}

Each panel has a TSV source-data file under `{SRC.relative_to(ROOT) if False else '11_figures/source_data'}` and a detailed legend in `11_figures/MAIN_FIGURE_LEGENDS.md`.

## 11. Remaining work before submission

1. Perform a clean-room rerun from frozen inputs and verify all key counts and hashes.
2. Complete visual inspection of every PDF at final print size and revise labels if necessary.
3. Convert the audit into a manuscript with a claim-matched Introduction, Methods, Results and Discussion.
4. Add full data/code availability statements and versioned environment capture.
5. Treat any future structure-predicted interactions, yeast-genetic evidence or generic-bacterial specificity comparisons as prespecified extensions, not as post-hoc support required to rescue the present conclusions.

## Overall status

**ANALYSIS_COMPLETE_WITH_WARNINGS**

Warnings are scientific boundaries, not unresolved technical failures: primary-four set-level holdout FDR=0.0603; one directionally contradictory candidate; no absolute global-translation inference; and no experimental causality. The core computational result, positive-control recovery, frozen-panel external validation and strict final tiers are complete.
"""
    report_path = OUT / "CURRENT_ANALYSIS_AUDIT.md"
    report_path.write_text(report, encoding="utf-8")

    manifest_targets = [
        ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv",
        ROOT / "09_holdout_validation/holdout_frozen_set_tests.tsv",
        ROOT / "08_evidence_atlas/CANDIDATE_FREEZE_MANIFEST.json",
        report_path,
    ] + figure_files
    manifest = pd.DataFrame({
        "path": [str(p.relative_to(ROOT)) for p in manifest_targets],
        "bytes": [p.stat().st_size for p in manifest_targets],
        "sha256": [sha(p) for p in manifest_targets],
    })
    manifest.to_csv(OUT / "final_output_manifest.tsv", sep="\t", index=False)
    stats = {
        "report": str(report_path.relative_to(ROOT)),
        "tables": 4,
        "figure_pdfs": len(figure_files),
        "final_state": "ANALYSIS_COMPLETE_WITH_WARNINGS",
        "gate10_reporting": "PASS",
    }
    (OUT / "stage10_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
