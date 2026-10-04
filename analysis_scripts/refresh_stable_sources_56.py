"""Refresh the canonical source TSVs under 14_submission_results/supplementary_tables
for the second-round editorial must-fix items (56-based discovery + RAB split).

Regenerates STable2 (enrichment), STable10 (selection null), STable18a/18b (annotation),
all from the current 13_postfreeze_interpretation outputs. The holdout-side tables are
untouched.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PF = ROOT / "13_postfreeze_interpretation"
DST = ROOT / "14_submission_results" / "supplementary_tables"
DST.mkdir(parents=True, exist_ok=True)

A4 = pd.read_csv(PF / "a4_panel_enrichment.tsv", sep="\t")
F7 = json.loads((PF / "f7_selection_aware_enrichment.json").read_text(encoding="utf-8"))
H2 = json.loads((PF / "h2_annotation_opportunity_null.json").read_text(encoding="utf-8"))

# STable2 -- enrichment (6 rows, split RAB)
A4.to_csv(DST / "STable2_panel_enrichment.tsv", sep="\t", index=False)

# STable10 -- selection null
rows = [
    ("universe", "26,761 unique human genes after collapsing 26,932 atlas rows by human symbol"),
    ("sampling", "%d-gene pseudo-panels sampled within propensity-score calipers; without "
                 "replacement within each panel" % F7["n_panel"]),
    ("matching_variables", "evidence-domain count; CRISPR tier; LCV detectability; PPI availability; "
                           "expression-study count; baseline UI expression; U937 Chronos gene effect; "
                           "network degree; cross-pathogen publication count; common-essential status"),
    ("propensity_caliper", "0.02 (widened to 0.05, then the full non-panel pool, if fewer than 20 "
                           "controls were available)"),
    ("n_pseudo_panels", "10,000"),
    ("empirical_p_convention", "(exceedances + 1)/(N + 1); no exceedances gives 1/10,001, "
                               "approximately 1.0e-4 per test"),
    ("multiple_testing", "Two prespecified enrichment tests; Bonferroni-adjusted empirical P "
                         "approximately 2.0e-4 per test"),
    ("observed", "RAB=%d/%d; trafficking=%d/%d" % (F7["obs_rab_in_panel"], F7["n_panel"],
                                                    F7["obs_trafficking_in_panel"], F7["n_panel"])),
    ("null_rab", "median %g; 95th percentile %g; maximum %g"
                 % (F7["null_rab_median"], F7["null_rab_p95"], F7["null_rab_max"])),
    ("null_trafficking", "median %g; 95th percentile %g; maximum %g"
                         % (F7["null_traf_median"], F7["null_traf_p95"], F7["null_traf_max"])),
    ("caveat", "Controls measured ascertainment features; unmeasured ascertainment bias remains possible"),
]
pd.DataFrame(rows, columns=["parameter", "value"]).to_csv(
    DST / "STable10_selection_null_details.tsv", sep="\t", index=False)

# STable18a -- annotation null
def _row(n, label):
    cal = "%d / %d / %d" % (n["pool_caliper_counts"]["0.02"], n["pool_caliper_counts"]["0.05"],
                            n["pool_caliper_counts"]["all"])
    return [label, n["n_covariates"], ", ".join(n["covariates"]), n["universe_n"], n["n_panel"],
            n["n_draws"], n["obs_rab_in_panel"], n["null_rab_median"], n["null_rab_p95"],
            n["null_rab_max"], "%.2e" % n["empirical_p_rab"], round(n["bonferroni_adjusted_p_rab"], 4),
            n["obs_trafficking_in_panel"], n["null_traf_median"], n["null_traf_p95"],
            n["null_traf_max"], "%.2e" % n["empirical_p_trafficking"],
            round(n["bonferroni_adjusted_p_trafficking"], 4),
            "%.3f" % n["max_abs_smd_before"], "%.3f" % n["max_abs_smd_after"], cal]

hdr = ["Matched null", "Covariates (n)", "Covariates", "Universe (n)", "Panel (n)",
       "Pseudo-panels (n)", "RAB observed in panel", "RAB null median", "RAB null 95th pct",
       "RAB null max", "Empirical P, RAB", "Bonferroni P, RAB", "Trafficking observed in panel",
       "Trafficking null median", "Trafficking null 95th pct", "Trafficking null max",
       "Empirical P, trafficking", "Bonferroni P, trafficking", "Max |SMD| before matching",
       "Max |SMD| after matching", "Controls at caliper 0.02 / 0.05 / full"]
pd.DataFrame([hdr,
              _row(H2["nulls"]["nullA_annotation_opportunity"], "Null A: annotation opportunity"),
              _row(H2["nulls"]["nullB_broad_membrane_context"], "Null B: Null A + broad membrane context")],
             ).to_csv(DST / "STable18a_annotation_opportunity_null.tsv", sep="\t",
                      index=False, header=False)

# STable18b -- panel annotation profile
p = H2["panel_annotation_profile"]
b_rows = [
    ["Annotation feature", "Statistic", "Panel (%d genes)" % H2["nulls"]["nullA_annotation_opportunity"]["n_panel"], "Universe"],
    ["GO biological-process terms (distinct)", "median", p["_go_bp_n"]["panel_median"], p["_go_bp_n"]["universe_median"]],
    ["GO cellular-component terms (distinct)", "median", p["_go_cc_n"]["panel_median"], p["_go_cc_n"]["universe_median"]],
    ["GO terms, all aspects (distinct)", "median", p["_go_terms_n"]["panel_median"], p["_go_terms_n"]["universe_median"]],
    ["GO:0016020 membrane descendant terms", "median", p["_n_membrane_terms"]["panel_median"], p["_n_membrane_terms"]["universe_median"]],
    ["GO:0043231 intracellular membrane-bounded organelle terms", "median", p["_n_imbo_terms"]["panel_median"], p["_n_imbo_terms"]["universe_median"]],
    ["Any GO:0016020 membrane annotation (fraction)", "fraction", round(p["_membrane"]["panel_frac"], 3), round(p["_membrane"]["universe_frac"], 3)],
    ["Any GO:0043231 IMBO annotation (fraction)", "fraction", round(p["_imbo"]["panel_frac"], 3), round(p["_imbo"]["universe_frac"], 3)],
]
pd.DataFrame(b_rows).to_csv(DST / "STable18b_panel_annotation_profile.tsv", sep="\t",
                            index=False, header=False)

print("refreshed: STable2, STable10, STable18a, STable18b")
