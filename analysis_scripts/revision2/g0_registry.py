"""Declarative claim registry for Gate G0 (WP1 claim ledger).

Every quantitative statement in the revised manuscript is expressed here as one or more
*machine-checkable* assertions that point at an archived output file, an exact column or a
named derived calculation. Nothing in this file is a judgement: a claim's status is
computed by ``run_claim_ledger.py`` from assertion outcomes, never written by hand.

Assertion kinds
---------------
``json``        value at a dotted path in a JSON file
``tsv_cell``    one column of the unique row selected by ``where``
``tsv_row``     several columns of the unique row selected by ``where``
``tsv_count``   number of rows matching ``where``
``tsv_nunique`` number of distinct values of ``column``
``text``        literal token(s) that must occur in the pinned manuscript
``path``        file must exist
``sha``         file SHA-256 must equal the pinned value
``calc``        named derived quantity from ``CALCULATIONS``

Status vocabulary (computed, never assigned by hand)
----------------------------------------------------
AUTO_PASS      every assertion passed
FAIL           at least one assertion did not match
STALE          a pinned file hash no longer matches (the analysis moved under the claim)
MANUAL_PASS    no machine assertion is possible; a recorded manual check exists
NOT_TESTABLE   no machine assertion and no manual check

A claim may carry ``text_action`` when the number is correct but the sentence needs work;
that is logged separately and never silently converted into a pass.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

MANUSCRIPT = "revision2_2026-11-23/REVISION2_MANUSCRIPT_WORKING_DRAFT.md"

#: The manuscript snapshot this ledger audits. G0 is invalid for any other snapshot.
MANUSCRIPT_SHA256 = "e6c3bf070604723921d4aebd69d9cb6a3d700f5316cd6f83674b4073d046002a"

#: Frozen inputs. A mismatch against any of these makes the claim STALE.
PINNED_INPUTS: dict[str, str] = {
    "08_evidence_atlas/host_dependency_evidence_atlas.tsv.gz":
        "9f1ccdc73894da2e1e4ee251db296193df63cd9493b3e16383913e5a3dcb3da5",
    "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv":
        "16cdb3bef2601c3c9f8d278034bf85da1b04f4e1c7a5789d915c7a6d6b0085d0",
    "08_evidence_atlas/HOLDOUT_VALIDATION_RULES_FROZEN.json":
        "4c3839ddd09a56b30ce78aeff5f622ae6f1f8acb68a7c1836b49960542eaa79f",
    "08_evidence_atlas/CANDIDATE_FREEZE_MANIFEST.json":
        "2bc6e63742f083c89d2ad1365c881f7265acbb5dc35322d4a00f9c4729d340db",
    "scripts/stage6c_integrated_evidence_atlas.py":
        "ffc4b852c5e441f6d98c55cbe299a91784eadfcd421d0ebb13357e0458ae58e5",
    "scripts/stage6d_freeze_candidates.py":
        "d2d6cd8a0e04cad8ae3c34bfffbe74d7081d1f310e78dcfffef47f628f5aca06",
    "scripts/stage7b_holdout_validation.py":
        "730f37cb6e6a117e65187f4987a94ddc2978107714921ec57b97e11905d27d4f",
    "scripts/stage8_final_tier_matrix.py":
        "f2f2ea96db6e3a15124a13b9dc194c1bff181db38ddbdd7f3f067970d72e6b2f",
    "scripts/trafficking_universe.py":
        "54b7e259c60665c1991fbe51302ee99c2fe88dd82ca8f4ebe62087619918d4ef",
    "revision2_2026-11-23/REVISION2_WP4_RULE_SENSITIVITY_STATS.json":
        "8303398945c330fa64fe4926a5052f0cc31f08f05586736978af7acb7dd33e97",
    "revision2_2026-11-23/REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv":
        "495644656beaaf95541bfbf7ef51421beb9e49a7dff5f2983e94a97d1ec99704",
    MANUSCRIPT:
        MANUSCRIPT_SHA256,
}


# --------------------------------------------------------------------------- #
# Assertion constructors
# --------------------------------------------------------------------------- #

def jx(file: str, path: str, expect, rel: float = 1e-3) -> dict:
    return {"kind": "json", "file": file, "path": path, "expect": expect, "rel": rel}


def tj(file: str, where: dict, column: str, expect, rel: float = 1e-3) -> dict:
    return {"kind": "tsv_cell", "file": file, "where": where, "column": column,
            "expect": expect, "rel": rel}


def tr(file: str, where: dict, expect: dict, rel: float = 1e-3) -> dict:
    return {"kind": "tsv_row", "file": file, "where": where, "expect": expect, "rel": rel}


def tc(file: str, expect: int, where: dict | None = None) -> dict:
    return {"kind": "tsv_count", "file": file, "where": where, "expect": expect}


def tn(file: str, column: str, expect: int, where: dict | None = None) -> dict:
    return {"kind": "tsv_nunique", "file": file, "column": column, "expect": expect,
            "where": where}


def tx(*tokens: str) -> dict:
    return {"kind": "text", "tokens": list(tokens)}


def px(file: str) -> dict:
    return {"kind": "path", "file": file}


def sh(file: str, expect: str) -> dict:
    return {"kind": "sha", "file": file, "expect": expect}


def cx(name: str, expect, rel: float = 1e-3) -> dict:
    return {"kind": "calc", "name": name, "expect": expect, "rel": rel}


# --------------------------------------------------------------------------- #
# Derived calculations
# --------------------------------------------------------------------------- #

def _atlas() -> pd.DataFrame:
    return pd.read_csv(ROOT / "08_evidence_atlas/host_dependency_evidence_atlas.tsv.gz",
                       sep="\t", low_memory=False)


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "1.0", "yes"})


def _atlas_unique_symbols() -> int:
    return int(_atlas().human_symbol.astype(str).str.upper().nunique())


def _atlas_rows() -> int:
    return int(len(_atlas()))


def _panel_genes() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv",
                        sep="\t")
    return int(frame.human_symbol.astype(str).str.upper().nunique())


def _panel_predefined_direction() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv",
                        sep="\t")
    return int(frame.expected_wt_vs_uninfected.isin(["UP", "DOWN"]).sum())


def _panel_not_predefined() -> str:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv",
                        sep="\t")
    return ",".join(sorted(frame.loc[
        frame.expected_wt_vs_uninfected == "NOT_PREDEFINED", "human_symbol"
    ].astype(str)))


def _ppi_lcv_intersection() -> int:
    """Atlas-level PPI(A/B) ∧ LCV intersection = 19.

    **Not used by any claim.** The manuscript's denominator of 18 is the *coupling
    universe* (16,667 genes with both layers evaluable), which the archived
    ``g_ppi_lcv_exact_rate_ci.tsv`` row reports directly. This function is retained only
    to document the two-population distinction; asserting 18 from it would be wrong.
    """
    atlas = _atlas()
    agg = atlas.assign(
        sym=atlas.human_symbol.astype(str).str.upper(),
        ppi=_truthy(atlas.ppi_AB_flag),
        lcv=_truthy(atlas.lcv_spatial_supported),
    ).groupby("sym")[["ppi", "lcv"]].max()
    return int((agg.ppi & agg.lcv).sum())


def _holdout_genes() -> int:
    frame = pd.read_csv(ROOT / "09_holdout_validation/holdout_all_contrasts.tsv.gz",
                        sep="\t", usecols=["human_symbol"])
    return int(frame.human_symbol.astype(str).nunique())


def _rab1a_bh_rank() -> int:
    frame = pd.read_csv(ROOT / "09_holdout_validation/holdout_all_contrasts.tsv.gz",
                        sep="\t")
    sub = frame.loc[frame.contrast_id == "WT_vs_UI", ["human_symbol", "FDR"]]
    sub = sub.dropna().sort_values("FDR", kind="mergesort").reset_index(drop=True)
    hits = sub.index[sub.human_symbol.astype(str).str.upper() == "RAB1A"]
    return int(hits[0]) + 1 if len(hits) else -1


def _lcv_condition_sources() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/lcv_spatial_evidence.tsv", sep="\t")
    tokens: set[str] = set()
    for value in frame.lcv_sources.dropna().astype(str):
        tokens.update(t.strip() for t in value.split(";") if t.strip())
    return len(tokens)


def _lcv_studies() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/lcv_spatial_evidence.tsv", sep="\t")
    tokens: set[str] = set()
    for value in frame.lcv_sources.dropna().astype(str):
        for t in value.split(";"):
            if t.strip():
                tokens.add(t.split(":")[0].strip())
    return len(tokens)


def _known_positive_recovered() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/known_positive_lcv_recovery.tsv", sep="\t")
    return int(_truthy(frame.lcv_spatial_supported).sum())


def _known_positive_total() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/known_positive_lcv_recovery.tsv", sep="\t")
    return int(frame.human_symbol.astype(str).nunique())


def _known_positive_not_recovered() -> str:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/known_positive_lcv_recovery.tsv", sep="\t")
    missing = frame.loc[~_truthy(frame.lcv_spatial_supported), "human_symbol"].astype(str)
    return ",".join(sorted(missing))


def _ppi_class_count() -> int:
    stats = json.loads((ROOT / "03_interactome/stage1_interactome_statistics.json")
                       .read_text(encoding="utf-8"))
    return len(stats["class_counts_pairs"])


def _holdout_design_rows() -> int:
    frame = pd.read_csv(ROOT / "08_evidence_atlas/FROZEN_HOLDOUT_SAMPLE_DESIGN.tsv", sep="\t")
    return int(len(frame))


def _global_bh_passes() -> int:
    frame = pd.read_csv(
        ROOT / "13_postfreeze_interpretation/SupplementaryTable_holdout_BH_comparison.tsv",
        sep="\t")
    return int(_truthy(frame.pass_global_bh).sum())


def _holdout_set_test_passes() -> int:
    frame = pd.read_csv(ROOT / "09_holdout_validation/holdout_frozen_set_tests.tsv", sep="\t")
    return int(_truthy(frame.passes_frozen_endpoint).sum())


def _holdout_set_tests() -> int:
    frame = pd.read_csv(ROOT / "09_holdout_validation/holdout_frozen_set_tests.tsv", sep="\t")
    return int(len(frame))


def _a4_member_count(label_suffix: str) -> int:
    frame = pd.read_csv(ROOT / "13_postfreeze_interpretation/a4_panel_enrichment.tsv",
                        sep="\t")
    row = frame.loc[frame.label == label_suffix]
    if len(row) != 1:
        raise ValueError(f"a4 label not unique: {label_suffix}")
    return int(row.iloc[0].members_in_query)


def _a4_background_members(label_suffix: str) -> int:
    frame = pd.read_csv(ROOT / "13_postfreeze_interpretation/a4_panel_enrichment.tsv",
                        sep="\t")
    row = frame.loc[frame.label == label_suffix]
    return int(row.iloc[0].members_in_background)


def _trafficking_module():
    import importlib
    sys_path = str(ROOT / "scripts")
    import sys as _sys
    if sys_path not in _sys.path:
        _sys.path.insert(0, sys_path)
    return importlib.import_module("trafficking_universe")


def _rab_gtpase_size() -> int:
    return len(set(_trafficking_module().RAB_GTPASE))


def _rab_associated_size() -> int:
    return len(set(_trafficking_module().RAB_ASSOCIATED))


def _trafficking_size() -> int:
    return len(set(_trafficking_module().ALL_GENES))


def _trafficking_families() -> int:
    return len(_trafficking_module().TRAFFICKING)


def _coupling_permutations() -> int:
    """Read ``n_iter`` from ``scripts/postfreeze_a1_coupling.py`` by AST.

    The coupling P values depend on this count, but the frozen script is never imported
    or executed here: importing it would run the whole analysis. Parsing the assignment
    keeps the check machine-verifiable while leaving the frozen file untouched.
    """
    import ast

    source = (ROOT / "scripts/postfreeze_a1_coupling.py").read_text(encoding="utf-8")
    tree = ast.parse(source, filename="postfreeze_a1_coupling.py")
    values: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "n_iter":
                    value = node.value
                    if isinstance(value, ast.Constant) and isinstance(value.value, int):
                        values.append(int(value.value))
                    else:
                        raise ValueError("n_iter is assigned a non-constant expression")
    if not values:
        raise ValueError("no integer n_iter assignment found in postfreeze_a1_coupling.py")
    if len(set(values)) != 1:
        raise ValueError(f"inconsistent n_iter assignments: {values}")
    return values[0]


CALCULATIONS = {
    "atlas_unique_symbols": _atlas_unique_symbols,
    "atlas_rows": _atlas_rows,
    "panel_genes": _panel_genes,
    "panel_predefined_direction": _panel_predefined_direction,
    "panel_not_predefined": _panel_not_predefined,
    # `ppi_lcv_intersection` is deliberately NOT registered: the atlas-level value (19)
    # is not the denominator the manuscript uses (18, coupling universe).
    "holdout_genes": _holdout_genes,
    "rab1a_bh_rank": _rab1a_bh_rank,
    "lcv_condition_sources": _lcv_condition_sources,
    "lcv_studies": _lcv_studies,
    "known_positive_recovered": _known_positive_recovered,
    "known_positive_total": _known_positive_total,
    "known_positive_not_recovered": _known_positive_not_recovered,
    "ppi_class_count": _ppi_class_count,
    "holdout_design_rows": _holdout_design_rows,
    "global_bh_passes": _global_bh_passes,
    "holdout_set_test_passes": _holdout_set_test_passes,
    "holdout_set_tests": _holdout_set_tests,
    "a4_rab_members": lambda: _a4_member_count("RAB_GTPase_family | panel | atlas"),
    "a4_trafficking_members": lambda: _a4_member_count("trafficking_universe | panel | atlas"),
    "a4_trafficking_background": lambda: _a4_background_members("trafficking_universe | panel | atlas"),
    "rab_gtpase_size": _rab_gtpase_size,
    "rab_associated_size": _rab_associated_size,
    "trafficking_universe_size": _trafficking_size,
    "trafficking_families": _trafficking_families,
    "coupling_permutations": _coupling_permutations,
}


# --------------------------------------------------------------------------- #
# Claim specifications
# --------------------------------------------------------------------------- #

S1 = "03_interactome/stage1_interactome_statistics.json"
S2 = "04_genetic_dependency/stage2_statistics.json"
S3 = "02_metadata/stage3_design_statistics.json"
S5 = "07_convergence/stage5_expression_statistics.json"
S6 = "08_evidence_atlas/stage6_statistics.json"
S6B = "08_evidence_atlas/stage6b_statistics.json"
S7B = "09_holdout_validation/stage7b_statistics.json"
S8 = "10_final_tiers/stage8_statistics.json"
A1 = "13_postfreeze_interpretation/a1_evidence_coupling_matched.tsv"
F5 = "13_postfreeze_interpretation/f5_matching_audit.tsv"
A4 = "13_postfreeze_interpretation/a4_panel_enrichment.tsv"
GS = "13_postfreeze_interpretation/g_sensitivity_reinforcement_summary.json"
GM = "13_postfreeze_interpretation/g_matching_sensitivity.tsv"
GI = "13_postfreeze_interpretation/g_ppi_lcv_interaction.tsv"
GR = "13_postfreeze_interpretation/g_ppi_lcv_exact_rate_ci.tsv"
JM = "13_postfreeze_interpretation/joint_multivariable_evidence_model.tsv"
CR = "13_postfreeze_interpretation/corr_robust_logistic.tsv"
CRS = "13_postfreeze_interpretation/corr_robust_summary.json"
H1 = "13_postfreeze_interpretation/h1_species_stratified_summary.json"
P0B = "13_postfreeze_interpretation/p0b_detectability_adjusted_summary.json"
F7 = "13_postfreeze_interpretation/f7_selection_aware_enrichment.json"
F6B = "13_postfreeze_interpretation/f6b_directional_100k.json"
CAM = "13_postfreeze_interpretation/g_camera_gene_set_sensitivity.json"
LO = "13_postfreeze_interpretation/g_panel_leaveout.tsv"
D5 = "13_postfreeze_interpretation/d5_evidence_layers_summary.json"
DM = "13_postfreeze_interpretation/depmap_coessentiality_stats.json"
D5B = "13_postfreeze_interpretation/d5_biological_confound_summary.json"
HC = "09_holdout_validation/holdout_frozen_candidate_results.tsv"
HS = "09_holdout_validation/holdout_frozen_set_tests.tsv"
BH = "13_postfreeze_interpretation/SupplementaryTable_holdout_BH_comparison.tsv"
ATLAS = "08_evidence_atlas/host_dependency_evidence_atlas.tsv.gz"
PANEL = "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv"
TRAFFIC = "scripts/trafficking_universe.py"
PARITY = "revision2_2026-11-23/REVISION2_PARITY_GENE_LEVEL.tsv"
WP4_STATS = "revision2_2026-11-23/REVISION2_WP4_RULE_SENSITIVITY_STATS.json"
WP4_TABLE = "revision2_2026-11-23/REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv"

ORIG = "original pre-holdout freeze"
POST = "original post-freeze"


def spec(claim_id, claim_text, estimate, source_file, generating_script, analysis_status,
         location, assertions, manual_note="", text_action="") -> dict:
    return {
        "claim_id": claim_id,
        "claim_text": claim_text,
        "estimate_ci_p": estimate,
        "source_file": source_file,
        "generating_script": generating_script,
        "analysis_status": analysis_status,
        "manuscript_location": location,
        "assertions": assertions,
        "manual_note": manual_note,
        "text_action": text_action,
    }


SPECS: list[dict] = [
    # ---------------- data assembly ----------------
    spec("WP1-001", "Discovery expression cohorts", "6 cohorts, 104 samples",
         f"{S3}; {S5}", "stage3/stage5", ORIG, "Results L41",
         [jx(S3, "n_datasets", 6), jx(S3, "n_samples", 104),
          # stage5 counts *datasets* entering the convergence step (6). stage3's
          # `primary_contrasts` is 11 and is a contrast-level quantity, not a dataset count,
          # so it must not be used to support the "six expression datasets" sentence.
          jx(S5, "primary_contrasts", 6),
          tx("six expression datasets representing 104 samples")]),

    spec("WP1-002", "Curated physical-interaction records",
         "152 records -> 79 pairs -> 54 host proteins",
         S1, "stage1_interactome.py", ORIG, "Results L41, L43",
         [jx(S1, "positive_records", 152), jx(S1, "unique_pairs", 79),
          jx(S1, "unique_hosts", 54), px(S1),
          tx("152 curated physical-interaction records")]),

    spec("WP1-003", "LCV-associated proteomic condition-level sources",
         "3 condition-level sources from 2 studies",
         "08_evidence_atlas/lcv_spatial_evidence.tsv", "stage6c_integrated_evidence_atlas.py",
         ORIG, "Results L41; Methods L113",
         [cx("lcv_condition_sources", 3), cx("lcv_studies", 2),
          jx(S6, "lcv_spatial_genes", 3122),
          tx("three LCV-associated proteomic condition-level sources from two studies")]),

    spec("WP1-004", "Holdout sample size", "45 samples",
         "08_evidence_atlas/FROZEN_HOLDOUT_SAMPLE_DESIGN.tsv",
         "stage6d_freeze_candidates.py", ORIG, "Results L41; Methods L113",
         [cx("holdout_design_rows", 45), px("08_evidence_atlas/FROZEN_HOLDOUT_SAMPLE_DESIGN.tsv"),
          tx("locked 45-sample PMA-THP-1 dataset")]),

    spec("WP1-005", "Enrichment universe", "26,761 unique genes", ATLAS,
         "stage6c_integrated_evidence_atlas.py", ORIG, "Results L41; Methods L113, L151",
         [cx("atlas_unique_symbols", 26761), cx("atlas_rows", 26932),
          jx("13_postfreeze_interpretation/a1_track_sanity.json",
             "duplicate_atlas_rows_collapsed", 171),
          tx("26,761 unique human genes")]),

    spec("WP1-006", "PPI class composition", "79 pairs; A=27, B=47, C=5; 54 hosts",
         "03_interactome/legionella_human_pairs_evidence.tsv", "stage1_interactome.py",
         ORIG, "Results L43",
         [jx(S1, "class_counts_pairs", {"A": 27, "B": 47, "C": 5}),
          jx(S1, "class_counts_hosts", {"A": 20, "B": 30, "C": 4}),
          cx("ppi_class_count", 3),
          tx("79 host-pathogen pairs involving 54 host proteins")]),

    # ---------------- cross-modal coupling ----------------
    spec("WP1-007", "Cross-pathogen degree enrichment before matching",
         "OR=3.42, P=0.00175", S1, "stage1_interactome.py", ORIG, "Results L43",
         [jx(S1, "fisher_odds_ratio", 3.42, rel=2e-3), jx(S1, "fisher_p", 0.00175, rel=5e-3),
          tx("OR=3.42, P=0.00175")]),

    spec("WP1-008", "Cross-pathogen degree enrichment after matching", "P=0.065", S1,
         "stage1_interactome.py", ORIG, "Results L43",
         [jx(S1, "matched_permutation_p", 0.065, rel=5e-3), tx("P=0.065")]),

    spec("WP1-009", "Fitness-filtered CRISPR evidence genes", "595 genes",
         "13_postfreeze_interpretation/a1_track_sanity.json", "postfreeze_a1_coupling.py",
         POST, "Results L45",
         [jx("13_postfreeze_interpretation/a1_track_sanity.json", "crispr_clean", 595),
          px(S2),
          tx("595 genes were termed fitness-filtered CRISPR evidence genes")]),

    spec("WP1-010", "Meta-convergent expression genes", "2,601 genes", S5,
         "stage5 (expression convergence)", ORIG, "Results L45",
         [jx(S5, "stable_cross_cohort", 2601), tx("2,601 meta-convergent genes")]),

    spec("WP1-011", ">=3-dataset replication", "490 genes; 432 overlap meta-convergent",
         D5, "postfreeze_d5_*", POST, "Results L45",
         [jx(D5, "multi_study_replicated_ge3", 490), jx(D5, "both_meta_and_ge3", 432),
          tx("490 genes met the independent >=3-dataset replication criterion"),
          tx("432 genes overlapped the meta-convergent set")]),

    spec("WP1-012", "LCV positive controls recovered", "8 of 9 by LCV; RABIF by AP-MS",
         "08_evidence_atlas/known_positive_lcv_recovery.tsv",
         "stage6c_integrated_evidence_atlas.py", ORIG, "Results L45",
         [cx("known_positive_recovered", 8), cx("known_positive_total", 9),
          cx("known_positive_not_recovered", "RABIF"),
          tx("recovered eight of the nine prespecified positive-control genes")]),

    spec("WP1-013", "PPI(A/B) -> CRISPR matched enrichment", "2.5-fold, permutation P=0.070",
         A1, "postfreeze_a1_coupling.py", POST, "Results L51; Abstract",
         [tr(A1, {"label": "PPI(A/B) -> filtered CRISPR evidence"},
             {"fold_vs_matched": 2.5, "perm_p": 0.0705}),
          tx("2.5-fold enrichment (permutation P=0.070")]),

    spec("WP1-014", "PPI(A/B) -> CRISPR covariate-adjusted logistic", "OR=3.1, P=0.004",
         F5, "postfreeze_f5_matching_audit.py", POST, "Results L51",
         [tr(F5, {"coupling": "PPI(A/B) -> filtered CRISPR evidence"},
             {"adj_logistic_OR": 3.063, "adj_logistic_p": 0.00442, "n_query": 43}),
          tx("OR=3.1 (P=0.004)")]),

    spec("WP1-015", "Propensity-match reconstruction rate ratio and balance",
         "RR=3.33; mean |SMD|=0.114", GS, "postfreeze_g_sensitivity(*)", POST, "Results L51",
         [jx(GS, "headline_matching.matched.fold_vs_matched", 3.3333),
          jx(GS, "headline_matching.matched.mean_abs_smd_after", 0.11395),
          tx("rate ratio of 3.33"), tx("mean absolute SMD=0.114")]),

    spec("WP1-016", "Overlap-weighted PPI sensitivity",
         "OR=2.53 (95% CI 0.56-11.51), P=0.230", GM, "postfreeze_g_sensitivity(*)", POST,
         "Results L51; Fig. 5D",
         [tj(GM, {"method": "overlap weighting"}, "estimate", 2.53, rel=5e-3),
          tj(GM, {"method": "overlap weighting"}, "ci_low", 0.555971, rel=5e-3),
          tj(GM, {"method": "overlap weighting"}, "ci_high", 11.513077, rel=5e-3),
          tx("OR=2.53 (95% CI 0.56-11.51, P=0.230")]),

    spec("WP1-017", "LCV -> CRISPR matched enrichment", "5.3-fold, P=8.9e-52", A1,
         "postfreeze_a1_coupling.py", POST, "Results L53; Abstract",
         [tr(A1, {"label": "LCV -> filtered CRISPR evidence"},
             {"fold_vs_matched": 5.275, "fisher_p": 8.85e-52}),
          tx("fold enrichment 5.3, P=8.9x10-52")]),

    spec("WP1-018", "LCV -> CRISPR adjusted OR", "OR=5.5", F5,
         "postfreeze_g_ppi_lcv_interaction(*)", POST, "Results L53",
         [tj(F5, {"coupling": "LCV -> filtered CRISPR evidence"}, "adj_logistic_OR", 5.453),
          tx("adjusted OR=5.5")]),

    spec("WP1-019", "PPI/LVC conjunction rate", "7/18 = 38.9% (95% CI 17.3-64.3%)", GR,
         "postfreeze_g_ppi_lcv_interaction(*)", POST, "Results L53; Fig. 5B",
         [tr(GR, {"group": "PPI(A/B)+LCV"},
             {"successes": 7, "n": 18, "rate": 0.388889,
              "ci_low": 0.172986, "ci_high": 0.642549}),
          # Denominator note: this 18 is the *coupling universe* intersection
          # (16,667 genes with both layers evaluated), not the atlas-level intersection
          # of ppi_AB_flag and lcv_spatial_supported, which is 19. The archived
          # exact-rate row above is the authoritative source for 7/18, so no atlas-level
          # `ppi_lcv_intersection` calculation is asserted here.
          jx("13_postfreeze_interpretation/h1_species_stratified_summary.json",
             "joint_model.audited_6.n_universe", 16667),
          tx("7 of 18 genes (38.9%; exact 95% CI 17.3-64.3%")]),

    spec("WP1-020", "PPI x LCV interaction", "OR=0.46 (95% CI 0.09-2.26), P=0.336", GI,
         "postfreeze_g_ppi_lcv_interaction(*)", POST, "Results L53; Fig. 5 legend",
         [tr(GI, {"term": "ppi_x_lcv"},
             {"OR": 0.4568, "ci_low": 0.0925, "ci_high": 2.2556, "p": 0.3362}),
          tx("OR=0.46, 95% CI 0.09-2.26, P=0.336")]),

    spec("WP1-021", "PPI(A/B) -> LCV matched enrichment", "fold=0.78, P=0.90", A1,
         "postfreeze_a1_coupling.py", POST, "Results L55",
         [tr(A1, {"label": "PPI(A/B) -> LCV"}, {"fold_vs_matched": 0.7826, "perm_p": 0.9025}),
          tx("fold enrichment 0.78, P=0.90")]),

    spec("WP1-022", "CRISPR -> response and LCV -> response adjusted and matched estimates",
         "adjusted OR=0.743/0.731; adjusted P=0.00739/0.00000 archived; McNemar P=0.38278/0.1012",
         F5, "postfreeze_f5_matching_audit.py", POST, "Results L55",
         [tr(F5, {"coupling": "filtered CRISPR evidence -> response"},
             {"adj_logistic_OR": 0.743, "adj_logistic_p": 0.00739,
              "mcnemar_p": 0.38278, "n_query": 538}),
          tr(F5, {"coupling": "LCV -> response"},
             {"adj_logistic_OR": 0.731, "adj_logistic_p": 0.0,
              "mcnemar_p": 0.1012, "n_query": 2874}),
          tx("CRISPR→response OR=0.743, P=0.00739"),
          tx("LCV→response OR=0.731, P rounded to 0.00000 at the archived five-decimal precision"),
          tx("matched-pair McNemar tests were null (P=0.38278 and P=0.1012"),
          tx("estimator discordance also precludes interpreting them as robust depletion")]),

    spec("WP1-023", "Joint multivariable model, LCV term",
         "OR=9.96 (95% CI 8.34-11.9), P=1.1e-140", JM,
         "postfreeze_a1_coupling.py (joint model)", POST, "Results L57; Fig. 5E",
         [tr(JM, {"term": "lcv"},
             {"estimate": 9.9553, "conf.low": 8.3444, "conf.high": 11.9232,
              "p.value": 1.092e-140}),
          tx("OR=9.96, 95% CI 8.34-11.9, P=1.1x10-140")]),

    spec("WP1-024", "Joint multivariable model, PPI term", "OR=6.01, P=4.6e-6", JM,
         "postfreeze_a1_coupling.py (joint model)", POST, "Results L57",
         [tr(JM, {"term": "ppi"}, {"estimate": 6.0143, "p.value": 4.552e-6}),
          tx("OR=6.01, P=4.6x10-6")]),

    spec("WP1-025", "Joint multivariable model, response term", "OR=0.93, P=0.49", JM,
         "postfreeze_a1_coupling.py (joint model)", POST, "Results L57; Discussion L89",
         [tr(JM, {"term": "response"}, {"estimate": 0.9334, "p.value": 0.4866}),
          tx("OR=0.93, P=0.49")]),

    spec("WP1-026", "Module-cluster-robust joint model",
         "OR=9.96 (95% CI 8.12-12.21), P=6.4e-108", CR,
         "postfreeze_corr_robust(*)", POST, "Results L57; Fig. 5E",
         [tr(CR, {"term": "lcv"},
             {"ci_lo_cluster": 8.1173, "ci_hi_cluster": 12.2095, "p_cluster": 6.4338e-108}),
          tx("module-cluster-robust OR=9.96, 95% CI 8.12-12.21, P=6.4x10-108")]),

    spec("WP1-027", "Species-stratified joint model, mouse term",
         "OR=6.06 (95% CI 4.89-7.51), P=4.7e-61", H1, "postfreeze_h1_species(*)", POST,
         "Results L61; Fig. 5G",
         [jx(H1, "joint_model.audited_6.OR_lcv_mouse", 6.0594),
          jx(H1, "joint_model.audited_6.OR_lcv_mouse_lo", 4.8911),
          jx(H1, "joint_model.audited_6.OR_lcv_mouse_hi", 7.5067),
          jx(H1, "joint_model.audited_6.p_lcv_mouse", 4.6602e-61),
          tx("OR=6.06, 95% CI 4.89-7.51, P=4.7x10-61")]),

    spec("WP1-028", "Species-stratified joint model, human term",
         "OR=1.33 (95% CI 1.07-1.66), P=0.0097", H1, "postfreeze_h1_species(*)", POST,
         "Results L61; Abstract",
         [jx(H1, "joint_model.audited_6.OR_lcv_human", 1.33267),
          jx(H1, "joint_model.audited_6.OR_lcv_human_lo", 1.07197),
          jx(H1, "joint_model.audited_6.OR_lcv_human_hi", 1.65676),
          jx(H1, "joint_model.audited_6.p_lcv_human", 0.0097167),
          tx("OR=1.33, 95% CI 1.07-1.66, P=0.0097")]),

    spec("WP1-029", "LCV species overlap", "1,282 mouse-only; 806 human-only; 790 both", H1,
         "postfreeze_h1_species(*)", POST, "Results L61",
         [jx(H1, "joint_model.audited_6.mouse_only_n", 1282),
          jx(H1, "joint_model.audited_6.human_only_n", 806),
          jx(H1, "joint_model.audited_6.both_n", 790),
          tx("1,282 genes were detected only in the mouse layer and 806 only in the human layer, "
             "with 790 in both")]),

    spec("WP1-030", "High-confidence mouse LCV exposure",
         "fold 4.9 -> 6.8; P=1.8e-46; coverage median 38.6%", P0B,
         "postfreeze_p0b_abundance(*)", POST, "Results L61; Fig. 5G",
         [jx(P0B, "full_mouse_exposure.fold", 4.8806),
          jx(P0B, "highconf_mouse_exposure.fold", 6.8108),
          jx(P0B, "highconf_mouse_exposure.fisher_p", 1.7599e-46),
          jx(P0B, "mouse_seq_cov_median_among_detected", 38.6),
          tx("matched enrichment fold 4.9 to 6.8, P=1.8x10-46")]),

    spec("WP1-031", "Abundance dose-response", "1.2% -> 29.5%; trend OR=2.26",
         P0B, "postfreeze_p0b_abundance(*)", POST, "Results L61; Discussion L99",
         [jx(P0B, "dose_response_levels.0.crispr_clean_rate", 0.012127),
          jx(P0B, "dose_response_levels.3.crispr_clean_rate", 0.295224),
          jx(P0B, "dose_response_trend_OR_per_level", 2.2555),
          tx("rate rises monotonically from 1.2% to 29.5%")]),

    # ---------------- panel enrichment ----------------
    spec("WP1-032", "RAB GTPase panel enrichment", "10/61; 78.3-fold, P=5.7e-17", A4,
         "postfreeze_a4_panel(*)", POST, "Results L65; Abstract",
         [cx("a4_rab_members", 10),
          cx("rab_gtpase_size", 61),
          tr(A4, {"label": "RAB_GTPase_family | panel | atlas"},
             {"fold": 78.3402, "fisher_p": 5.7205e-17}),
          tx("78.3-fold enrichment, P=5.7x10-17")]),

    spec("WP1-033", "RAB-associated panel enrichment", "10/77; 62.1-fold, P=6.8e-16", A4,
         "postfreeze_a4_panel(*)", POST, "Results L65",
         [cx("rab_associated_size", 77),
          tr(A4, {"label": "RAB_associated_machinery | panel | atlas"},
             {"members_in_query": 10, "fold": 62.0617, "fisher_p": 6.7867e-16}),
          tx("62.1-fold enrichment (P=6.8x10-16)")]),

    spec("WP1-034", "Trafficking panel enrichment", "12/275; 20.9-fold, P=4.1e-13", A4,
         "postfreeze_a4_panel(*)", POST, "Results L65; Abstract",
         [cx("a4_trafficking_members", 12), cx("a4_trafficking_background", 275),
          tr(A4, {"label": "trafficking_universe | panel | atlas"},
             {"fold": 20.8527, "fisher_p": 4.0748e-13}),
          tx("20.9-fold enrichment, P=4.1x10-13")]),

    spec("WP1-035", "Selection-aware null", "0 of 10,000 pseudo-panels; Bonferroni P~2e-4",
         F7, "postfreeze_f7_selection_null(*)", POST, "Results L65; Supp Fig. 6B",
         [jx(F7, "n_pseudo_panels", 10000), jx(F7, "n_enrichment_tests", 2),
          jx(F7, "null_rab_max", 4), jx(F7, "bonferroni_adjusted_p_rab", 2e-4, rel=5e-3),
          tx("None of 10,000 pseudo-panels reached the observed RAB GTPase or trafficking counts"),
          tx("approximately 2x10-4")]),

    spec("WP1-036", "CRISPR-tested universe sensitivity",
         "62.2x / 49.0x / 16.3x; P=5.4e-16 / 7.2e-12", A4, "postfreeze_a4_panel(*)", POST,
         "Results L65; Supp Fig. 6A",
         [tr(A4, {"label": "RAB_GTPase_family | panel | crispr_tested"},
             {"fold": 62.2367, "fisher_p": 5.4472e-16}),
          tr(A4, {"label": "RAB_associated_machinery | panel | crispr_tested"},
             {"fold": 48.9595}),
          tr(A4, {"label": "trafficking_universe | panel | crispr_tested"},
             {"fold": 16.2596, "fisher_p": 7.1804e-12}),
          tx("RAB GTPase 62.2-fold, P=5.4x10-16"),
          tx("trafficking 16.3-fold, P=7.2x10-12")]),

    spec("WP1-037", "RAB1A/RAB10 leave-two-out directional statistic",
         "52 genes; mean oriented z=3.78; P~1e-4", LO, "postfreeze_g_sensitivity(*)", POST,
         "Results L67",
         [tr(LO, {"scenario": "minus_RAB1A_RAB10"},
             {"panel_n": 52, "observed_mean_oriented_z": 3.7775, "empirical_p": 1e-4}),
          tx("Removing both genes left 52 genes with mean oriented z=3.78 and empirical P~1x10-4")]),

    spec("WP1-038", "CAMERA correlation-aware sensitivity",
         "54-gene P=1.87e-10; 52-gene P=6.34e-10", CAM, "postfreeze_g_camera(*)", POST,
         "Results L67; Methods L139",
         [jx(CAM, "results.0.camera_p", 1.8659e-10),
          jx(CAM, "results.3.panel_n", 52),
          jx(CAM, "results.3.camera_p", 6.3355e-10),
          tx("P=1.87x10-10"), tx("P=6.34x10-10")]),

    # ---------------- holdout ----------------
    spec("WP1-039", "Prespecified set tests passing", "8 of 15", HS,
         "stage7b_holdout_validation.py", ORIG, "Results L71; Abstract",
         [jx(S7B, "set_tests_pass_fdr05_n", 8), cx("holdout_set_test_passes", 8),
          cx("holdout_set_tests", 15),
          tx("Eight of 15 prespecified set tests passed endpoint-family FDR<0.05")]),

    spec("WP1-040", "Holdout-evaluable panel size", "56 -> 54 (AKAP2, HSPA2)",
         "09_holdout_validation/holdout_all_contrasts.tsv.gz",
         "stage7b_holdout_validation.py", ORIG, "Results L71; Methods L119, L131",
         [jx(S7B, "frozen_panel_n", 56), jx(S7B, "panel_measured_n", 54),
          cx("panel_genes", 56),
          tc(HC, 0, where={"human_symbol": "AKAP2"}),
          tc(HC, 0, where={"human_symbol": "HSPA2"}),
          tx("two panel genes (AKAP2 and HSPA2) fell below")]),

    spec("WP1-041", "Four-primary set test", "FDR=0.0603 (borderline, not positive)", HS,
         "stage7b_holdout_validation.py", ORIG, "Results L71; Fig. 4D",
         [tr(HS, {"set_name": "PRIMARY4_DIRECTIONAL"},
             {"set_n": 4, "empirical_fdr_within_family": 0.060294,
              "passes_frozen_endpoint": False}),
          tx("four-primary set was borderline (FDR=0.0603)")]),

    spec("WP1-042", "RAB1A holdout effect", "log2FC=+0.466, FDR=7.8e-6", HC,
         "stage7b_holdout_validation.py", POST, "Results L73; Fig. 4A",
         [tr(HC, {"human_symbol": "RAB1A", "contrast_id": "WT_vs_UI"},
             {"logFC": 0.466, "PValue": 7.335434e-07, "FDR": 7.800011e-06}),
          cx("rab1a_bh_rank", 1020),
          tx("log2FC of +0.466 (FDR=7.8x10-6)")]),

    spec("WP1-043", "RAB10 holdout effect", "log2FC=+0.531, FDR=1.5e-5", HC,
         "stage7b_holdout_validation.py", POST, "Results L73; Fig. 4B",
         [tr(HC, {"human_symbol": "RAB10", "contrast_id": "WT_vs_UI"},
             {"logFC": 0.5308, "FDR": 1.5133e-05}),
          tx("+0.531 (FDR=1.5x10-5)")]),

    spec("WP1-044", "Final mutually exclusive classes", "1 / 38 / 14 / 2 / 1", S8,
         "stage8_final_tier_matrix.py", POST, "Results L73",
         [jx(S8, "final_tier_counts.TIER1_EXTERNALLY_VALIDATED", 1),
          jx(S8, "final_tier_counts.TIER2_MULTIMODAL", 38),
          jx(S8, "final_tier_counts.TIER3_CONTEXT_DEPENDENT", 14),
          jx(S8, "final_tier_counts.TIER3_FITNESS_CAVEAT", 2),
          jx(S8, "final_tier_counts.REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION", 1),
          tx("38 Tier-2 candidates, 14 Tier-3 context-dependent candidates, "
             "two Tier-3 fitness-caveat candidates")]),

    spec("WP1-045", "Directional 100k matched null",
         "z=3.82 vs null 0.38+/-0.63, P~1e-5", F6B, "postfreeze_f6b_directional_100k.py",
         POST, "Results L79",
         [jx(F6B, "n_permutations", 100000), jx(F6B, "observed_directional_stat", 3.8183),
          jx(F6B, "matched_null_mean", 0.3767), jx(F6B, "matched_null_sd", 0.6258),
          jx(F6B, "empirical_p", 9.9999e-06),
          tx("signed z=3.82 versus matched-null mean 0.38+/-0.63")]),

    spec("WP1-046", "Global BH sensitivity", "8 of 15", BH, "postfreeze BH comparison script",
         POST, "Results L85",
         [cx("global_bh_passes", 8), cx("holdout_set_tests", 15),
          tx("with 8 of 15 tests passing")]),

    spec("WP1-047", "Module-decimation matched null",
         "median fold 6.1 (5.2-7.4); median Fisher P=5.9e-36; 200 replicates", CRS,
         "postfreeze_g_sensitivity(*)", POST, "Results L85; Fig. 5F",
         [jx(CRS, "B_decimation_fold_median", 6.1249),
          jx(CRS, "B_decimation_fold_p05", 5.1923),
          jx(CRS, "B_decimation_fold_p95", 7.4317),
          jx(CRS, "B_decimation_fisher_p_median", 5.9375e-36),
          jx(CRS, "B_decimation_B", 200),
          tx("median fold 6.1, 5th-95th percentile 5.2-7.4; median Fisher P=5.9x10-36")]),

    spec("WP1-048", "DepMap co-essentiality set", "20 Rab/ER-Golgi genes", DM,
         "postfreeze_depmap(*)", POST, "Results L85; Supp Fig. 4",
         [jx(DM, "n_machinery_genes", 20), jx(DM, "n_cell_lines", 1208),
          jx(DM, "empirical_p_pearson", 9.999e-05),
          tx("20 Rab/ER-Golgi genes")]),

    # ---------------- methods ----------------
    spec("WP1-049", "CRISPR tested universe", "20,561 genome-wide; 20,563 integrated", S2,
         "stage2_crispr_dependency.py", ORIG, "Methods L113",
         [jx(S2, "genomewide_genes", 20561),
          tr(A4, {"label": "RAB_GTPase_family | panel | crispr_tested"},
             {"n_universe": 20563}),
          tx("20,561 unique genes"), tx("20,563 genes")]),

    spec("WP1-050", "Ortholog mapping retention",
         "15,961/24,520 expression; 2,195/2,433 proteomics", f"{S5}; {S6B}",
         "stage3/stage5/stage6b", ORIG, "Methods L113",
         [jx(S5, "orthologs_mapped", 15961), jx(S5, "genes_tested", 24520),
          jx(S6B, "mouse_orthologs_mapped", 2195),
          jx(S6B, "mouse_symbols_submitted", 2433),
          tx("15,961 of 24,520"), tx("2,195 of 2,433")]),

    spec("WP1-051", "Holdout expression filter", "10,846 genes retained",
         "09_holdout_validation/holdout_all_contrasts.tsv.gz",
         "stage7b_holdout_validation.py", ORIG, "Methods L119, L131",
         [cx("holdout_genes", 10846), tx("10,846 genes retained")]),

    spec("WP1-052", "Holdout design", "15 groups x 3 replicates; MOI 10; 18 h",
         "08_evidence_atlas/FROZEN_HOLDOUT_SAMPLE_DESIGN.tsv",
         "stage6d_freeze_candidates.py", ORIG, "Methods L131",
         [cx("holdout_design_rows", 45),
          jx(D5B, "n_official_conditions", 15),
          tx("15 prespecified groups with three replicates per group"),
          tx("MOI 10 and 18 h")]),

    spec("WP1-053", "Permutation counts",
         "10,000 set tests; 100,000 directional; 5,000 multi-effector; 4,000 coupling",
         f"{S7B}; {F6B}; {D5B}; scripts/postfreeze_a1_coupling.py",
         "stage7b_holdout_validation.py; postfreeze_f6b_directional_100k.py; "
         "postfreeze_e3_matched_null.py; postfreeze_a1_coupling.py",
         "original pre-holdout freeze (endpoints) + original post-freeze (nulls)", "Methods L137",
         [jx(S7B, "permutations_per_set_test", 10000),
          jx(F6B, "n_permutations", 100000),
          jx(D5B, "n_matched_null_sets", 5000),
          # AST read of `n_iter = 4000`; the frozen script is parsed, never imported.
          cx("coupling_permutations", 4000),
          tx("10,000 gene-label permutations"),
          tx("Matched fold-enrichment coupling analyses used 4,000 label permutations")]),

    spec("WP1-054", "Gene-set definitions", "61 RAB; 77 RAB-associated; 17 families; 276 genes",
         TRAFFIC, "trafficking_universe.py", "frozen 2026-08-25 (one day before rules)",
         "Methods L155; ST21",
         [jx("13_postfreeze_interpretation/f7_selection_aware_enrichment.json",
             "universe_n", 26761),
          cx("rab_gtpase_size", 61), cx("rab_associated_size", 77),
          cx("trafficking_universe_size", 276), cx("trafficking_families", 17),
          tx("61 human RAB* GTPase genes"), tx("77 genes in total"), tx("17 families")]),

    spec("WP1-055", "Expected direction definition",
         "expected_wt_vs_uninfected = expression_consensus_direction", f"{PANEL}; {PARITY}",
         "stage6d_freeze_candidates.py", ORIG, "Methods L119",
         [cx("panel_predefined_direction", 55),
          cx("panel_not_predefined", "AKAP2"),
          cx("panel_genes", 56),
          tx("Expected directions were encoded by multiplying each panel row by +1")]),

    spec("WP1-056", "Implemented candidate-entry hierarchy and route rationale families",
         "four evidence_tier admission branches; two routes are rationale families",
         f"{ATLAS}; {PARITY}", "stage6c + stage6d", ORIG,
         "Methods L117; Results L73; Table 1 col. 2",
         [tx("Candidate inclusion was implemented by the frozen `evidence_tier` hierarchy rather than by two executable route filters"),
          tx("The two prespecified routes served as biological rationale families"),
          tx("one of four Tier-1/Tier-2 branches"),
          tx("Direct evidence comprised curated Class A/B PPI or high-confidence AP-MS evidence"),
          tx("lacked the translatome/mechanism support required by the strict final Tier-1 rule"),
          px(PARITY)]),

    # ---------------- post-review discovery-only rule sensitivity (WP4) ----------------
    spec("WP4-001", "Local rule-variant stability",
         "5 distinct panels; 33-86 genes; Jaccard 0.589-0.918; retention 13/40/0",
         f"{WP4_STATS}; {WP4_TABLE}", "render_wp4_rule_sensitivity.py",
         "post-review, algorithm-locked", "Results: rule sensitivity",
         [jx(WP4_STATS, "holdout_accessed", False),
          jx(WP4_STATS, "local_retention.n_local_observations", 5),
          jx(WP4_STATS, "local_retention.retained_in_all_local", 13),
          jx(WP4_STATS, "local_retention.retained_in_at_least_local_minus_one", 40),
          jx(WP4_STATS, "local_retention.local_retention_distribution.0", 0),
          jx(WP4_STATS, "local_retention.local_retention_distribution.1", 0),
          jx(WP4_STATS, "local_retention.local_retention_distribution.2", 0),
          tx("five distinct local one-change variants contained 33-86 genes"),
          tx("Jaccard similarities of 0.589-0.918"),
          tx("Thirteen frozen genes were retained by all five variants"),
          tx("40 were retained by at least four"),
          tx("none was retained by fewer than three")]),

    spec("WP4-002", "Structural domain and tier-only sensitivity",
         "domain panels 36 -> 14 -> 1; Tier 1 4 -> 1 at fixed n=56",
         WP4_TABLE, "render_wp4_rule_sensitivity.py",
         "post-review, algorithm-locked", "Results: rule sensitivity",
         [tj(WP4_TABLE, {"variant": "D3_panel_min_discovery_domains"}, "panel_n", 36),
          tj(WP4_TABLE, {"variant": "D4_panel_min_discovery_domains"}, "panel_n", 14),
          tj(WP4_TABLE, {"variant": "D5_panel_min_discovery_domains"}, "panel_n", 1),
          tj(WP4_TABLE, {"variant": "A3_T1_orthogonal_response_only"}, "panel_n", 56),
          tj(WP4_TABLE, {"variant": "A3_T1_orthogonal_response_only"}, "tier1_n", 1),
          tj(WP4_TABLE, {"variant": "A5_T1_require_response_and_mechanism"}, "panel_n", 56),
          tj(WP4_TABLE, {"variant": "A5_T1_require_response_and_mechanism"}, "tier1_n", 1),
          tx("reduced panel size from 36 to 14 to 1"),
          tx("changed Tier 1 membership from four genes to one without changing the 56-gene panel")]),

    spec("WP4-003", "Rule-variant enrichment is descriptive",
         "9 formal membership panels; all RAB/trafficking within-panel adjusted P < 0.05",
         f"{WP4_STATS}; {WP4_TABLE}", "render_wp4_rule_sensitivity.py",
         "post-review, algorithm-locked", "Results and Methods: rule sensitivity",
         [jx(WP4_STATS, "formal_membership_panels_n", 9),
          jx(WP4_STATS,
             "all_formal_membership_panels_rab_and_trafficking_bonferroni_p_lt_0_05",
             True),
          jx(WP4_STATS, "universe_n", 26761),
          jx(WP4_STATS, "curated_trafficking_gene_list_n", 276),
          jx(WP4_STATS, "gene_set_sizes_in_atlas.trafficking", 275),
          tx("26,761-gene atlas"),
          tx("curated trafficking list contains 276 genes, of which 275 are represented"),
          tx("descriptive sensitivity checks rather than independent discoveries")]),

    spec("WP4-004", "Locked WP4 interpretation",
         "local stability with structural and tier sensitivity; holdout not accessed",
         f"{WP4_STATS}; {WP4_TABLE}", "render_wp4_rule_sensitivity.py",
         "post-review, algorithm-locked", "Results; Methods; Supplementary legends",
         [jx(WP4_STATS, "status", "PASS"),
          jx(WP4_STATS, "gate_g2", "PASS"),
          jx(WP4_STATS, "holdout_accessed", False),
          jx(WP4_STATS, "a7_b1_counted_once", True),
          tx("The panel was therefore stable to several local rule modifications but sensitive to structural assumptions concerning cross-domain evidence requirements"),
          tx("Tier assignments were additionally sensitive to some rule changes that did not alter panel membership"),
          tx("The holdout was not accessed for this analysis")]),
]
