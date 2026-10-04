"""In-place patch of FINAL_SUBMISSION_2026-10-04/03_tables/Supplementary_Tables_1-21.xlsx
for the second-round editorial must-fix items.

Patches (no sheet is added or removed; the workbook keeps its 28 sheets):
  必须1 (56 vs 54)      -> ST2_Enrichment rebuilt from the regenerated 56-based
                           a4_panel_enrichment.tsv; ST18a/ST18b panel n 54->56.
  必须3 (RAB definition) -> ST21_GeneSetDefs rebuilt with the canonical RAB
                           GTPase (61) vs RAB-associated machinery (77) split.
  必须2 (abundance)      -> ST15d_Detectability gains Panel C (ordinal abundance
                           dose-response) and Panel D (within-detection internal
                           contrast + negative control + identifiability note).
  Contents descriptions  -> refreshed for the changed sheets.

The workbook is regenerated in place (a copy of the prior file is written first).
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import json
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "FINAL_SUBMISSION_2026-10-04"
XLSX = PKG / "03_tables" / "Supplementary_Tables_1-21.xlsx"
BACKUP = PKG / "03_tables" / "_backup_Supplementary_Tables_1-21_pre_56.xlsx"
PF = ROOT / "13_postfreeze_interpretation"

sys.path.insert(0, str(ROOT / "scripts"))
from trafficking_universe import (  # noqa: E402
    GENE_TO_FAMILY, RAB_GTPASE, RAB_ASSOCIATED, TRAFFICKING,
)

bold = Font(bold=True)
wrap = Alignment(wrap_text=True, vertical="top")

A4 = pd.read_csv(PF / "a4_panel_enrichment.tsv", sep="\t")
H2 = json.loads((PF / "h2_annotation_opportunity_null.json").read_text(encoding="utf-8"))
F7 = json.loads((PF / "f7_selection_aware_enrichment.json").read_text(encoding="utf-8"))
DR = pd.read_csv(PF / "p0b_abundance_dose_response.tsv", sep="\t")
DRS = json.loads((PF / "p0b_abundance_dose_response_summary.json").read_text(encoding="utf-8"))

if not BACKUP.exists():
    shutil.copy2(XLSX, BACKUP)
    print("backup written:", BACKUP.name)

wb = openpyxl.load_workbook(XLSX)


def clear_body(ws, first_row):
    """Delete every row from first_row to the end, leaving the title block intact."""
    if ws.max_row >= first_row:
        ws.delete_rows(first_row, ws.max_row - first_row + 1)


def put(ws, row, values, wrap_cols=()):
    for j, v in enumerate(values, start=1):
        c = ws.cell(row=row, column=j, value=v)
        if j in wrap_cols:
            c.alignment = wrap
    return row + 1


# ---------------------------------------------------------------- ST2
ws = wb["ST2_Enrichment"]
ws.cell(row=2, column=1, value=(
    "Canonical source: 13_postfreeze_interpretation\\a4_panel_enrichment.tsv "
    "(discovery enrichment computed on the frozen 56-gene candidate panel). "
    "RAB_GTPase_family = the 61 HGNC RAB* GTPases; RAB_associated_machinery additionally "
    "includes RAB GEF/effector and RAB GAP genes (77 total); trafficking_universe = the "
    "curated membrane-trafficking category (275 genes)."))
clear_body(ws, 4)
r = 4
r = put(ws, r, ["label", "n_query", "n_universe", "members_in_query", "members_in_background",
                "query_rate", "background_rate", "fold", "odds_ratio", "fisher_p"])
for _, row in A4.iterrows():
    r = put(ws, r, [row["label"], int(row["n_query"]), int(row["n_universe"]),
                    int(row["members_in_query"]), int(row["members_in_background"]),
                    float(row["query_rate"]), float(row["background_rate"]),
                    float(row["fold"]), float(row["odds_ratio"]), float(row["fisher_p"])])
print("ST2_Enrichment rows:", r - 4)

# ---------------------------------------------------------------- ST18a
ws = wb["ST18a_AnnotationNull"]
ws.cell(row=1, column=1, value=(
    "Supplementary Table 18a. Annotation-opportunity and broad-membrane matched nulls for "
    "the 56-gene frozen candidate panel. Null A matches GO annotation counts, annotation "
    "breadth, generic membrane and intracellular membrane-bounded organelle annotation, "
    "baseline expression, study count, network degree, publication count and CRISPR tier. "
    "Null B additionally matches the depth of membrane annotation. RAB membership, "
    "vesicle-mediated transport, ER-to-Golgi transport, Rab signal transduction and the "
    "membrane-trafficking category are deliberately NOT matched, since adjusting for the "
    "structure under test would be circular."))
clear_body(ws, 3)
header = ["Matched null", "Covariates (n)", "Covariates", "Universe (n)", "Panel (n)",
          "Pseudo-panels (n)", "RAB observed in panel", "RAB null median", "RAB null 95th pct",
          "RAB null max", "Empirical P, RAB", "Bonferroni P, RAB", "Trafficking observed in panel",
          "Trafficking null median", "Trafficking null 95th pct", "Trafficking null max",
          "Empirical P, trafficking", "Bonferroni P, trafficking", "Max |SMD| before matching",
          "Max |SMD| after matching", "Controls at caliper 0.02 / 0.05 / full"]
r = 3
r = put(ws, r, header)
_label_a = "Null A: annotation opportunity"
_label_b = "Null B: Null A + broad membrane context"
for key, label in (("nullA_annotation_opportunity", _label_a),
                   ("nullB_broad_membrane_context", _label_b)):
    n = H2["nulls"][key]
    cal = "%d / %d / %d" % (n["pool_caliper_counts"]["0.02"],
                            n["pool_caliper_counts"]["0.05"],
                            n["pool_caliper_counts"]["all"])
    r = put(ws, r, [
        label, n["n_covariates"], ", ".join(n["covariates"]), n["universe_n"], n["n_panel"],
        n["n_draws"], n["obs_rab_in_panel"], n["null_rab_median"], n["null_rab_p95"],
        n["null_rab_max"], "%.2e" % n["empirical_p_rab"], round(n["bonferroni_adjusted_p_rab"], 4),
        n["obs_trafficking_in_panel"], n["null_traf_median"], n["null_traf_p95"],
        n["null_traf_max"], "%.2e" % n["empirical_p_trafficking"],
        round(n["bonferroni_adjusted_p_trafficking"], 4),
        "%.3f" % n["max_abs_smd_before"], "%.3f" % n["max_abs_smd_after"], cal])
ws.cell(row=1, column=1).alignment = wrap
print("ST18a rows:", r - 3)

# ---------------------------------------------------------------- ST18b
ws = wb["ST18b_PanelAnnotation"]
prof = H2["panel_annotation_profile"]
ws.cell(row=1, column=1, value=(
    "Supplementary Table 18b. Annotation profile of the 56-gene frozen candidate panel versus "
    "the full enrichment universe, showing the annotation-opportunity imbalance that the "
    "matched nulls in 18a remove."))
clear_body(ws, 3)
r = 3
r = put(ws, r, ["Annotation feature", "Statistic", "Panel (56 genes)", "Universe"])
r = put(ws, r, ["GO biological-process terms (distinct)", "median",
                prof["_go_bp_n"]["panel_median"], prof["_go_bp_n"]["universe_median"]])
r = put(ws, r, ["GO cellular-component terms (distinct)", "median",
                prof["_go_cc_n"]["panel_median"], prof["_go_cc_n"]["universe_median"]])
r = put(ws, r, ["GO terms, all aspects (distinct)", "median",
                prof["_go_terms_n"]["panel_median"], prof["_go_terms_n"]["universe_median"]])
r = put(ws, r, ["GO:0016020 membrane descendant terms", "median",
                prof["_n_membrane_terms"]["panel_median"], prof["_n_membrane_terms"]["universe_median"]])
r = put(ws, r, ["GO:0043231 intracellular membrane-bounded organelle terms", "median",
                prof["_n_imbo_terms"]["panel_median"], prof["_n_imbo_terms"]["universe_median"]])
r = put(ws, r, ["Any GO:0016020 membrane annotation (fraction)", "fraction",
                round(prof["_membrane"]["panel_frac"], 3), round(prof["_membrane"]["universe_frac"], 3)])
r = put(ws, r, ["Any GO:0043231 IMBO annotation (fraction)", "fraction",
                round(prof["_imbo"]["panel_frac"], 3), round(prof["_imbo"]["universe_frac"], 3)])
ws.cell(row=1, column=1).alignment = wrap
print("ST18b rows:", r - 3)

# ---------------------------------------------------------------- ST21
ws = wb["ST21_GeneSetDefs"]
ws.cell(row=2, column=1, value=(
    "The canonical RAB GTPase family, the RAB-associated machinery, and the curated "
    "membrane-trafficking category used for the panel enrichment, with the Gene Ontology "
    "annotation source used for the functional-module partition. The canonical RAB GTPase "
    "family (the RAB* genes themselves) is the primary enrichment set; the RAB-associated "
    "machinery additionally includes RAB GEF/effector and RAB GAP genes."))
clear_body(ws, 5)
r = 5
r = put(ws, r, ["Panel A. Canonical RAB GTPase family (RAB* genes), n=%d" % len(RAB_GTPASE)])
r = put(ws, r, ["gene", "family"])
for g in sorted(RAB_GTPASE):
    r = put(ws, r, [g, GENE_TO_FAMILY.get(g, "?")])
r = put(ws, r, [None])
r = put(ws, r, ["Panel B. RAB-associated machinery additions (RAB GEF/effector + RAB GAP), n=%d"
                % (len(RAB_ASSOCIATED) - len(RAB_GTPASE))])
r = put(ws, r, ["gene", "family"])
for g in sorted(set(RAB_ASSOCIATED) - set(RAB_GTPASE)):
    r = put(ws, r, [g, GENE_TO_FAMILY.get(g, "?")])
r = put(ws, r, [None])
r = put(ws, r, ["Panel C. Curated membrane-trafficking category: %d families, %d genes (HGNC"
                % (len(TRAFFICKING), len(GENE_TO_FAMILY))])
r = put(ws, r, ["family", "n_genes", "genes"])
for fam in TRAFFICKING:
    genes = sorted(TRAFFICKING[fam])
    r = put(ws, r, [fam, len(genes), ", ".join(genes)], wrap_cols=(3,))
r = put(ws, r, [None])
r = put(ws, r, ["Panel D. Annotation source for the functional-module partition"])
r = put(ws, r, ["Item", "Value"])
r = put(ws, r, ["Gene Ontology Annotation file", "goa_human.gaf.gz (GAF version 2.2)"])
r = put(ws, r, ["Source database", "UniProt (Gene Ontology Annotation)"])
r = put(ws, r, ["Source association file date", "2026-04-30"])
r = put(ws, r, ["Retrieval / GOC date-generated", "2026-05-21"])
r = put(ws, r, ["Use", "Functional-module (GO biological process most-specific term) partition "
                       "for the correlation-robust analysis (Supplementary Table 15c)"])
ws.cell(row=5, column=1).font = bold
print("ST21 rows:", r - 5)

# ---------------------------------------------------------------- ST10
ws = wb["ST10_SelectionNull"]
vals = {
    "sampling": "%d-gene pseudo-panels sampled within propensity-score calipers; without "
                "replacement" % F7["n_panel"],
    "observed": "RAB=%d/%d; trafficking=%d/%d"
                % (F7["obs_rab_in_panel"], F7["n_panel"],
                   F7["obs_trafficking_in_panel"], F7["n_panel"]),
    "null_rab": "median %g; 95th percentile %g; maximum %g"
                % (F7["null_rab_median"], F7["null_rab_p95"], F7["null_rab_max"]),
    "null_trafficking": "median %g; 95th percentile %g; maximum %g"
                        % (F7["null_traf_median"], F7["null_traf_p95"], F7["null_traf_max"]),
}
for i, row in enumerate(ws.iter_rows(min_col=1, max_col=1), start=1):
    key = row[0].value
    if key in vals:
        ws.cell(row=i, column=2, value=vals[key])
print("ST10_SelectionNull patched:", sorted(vals))

# ---------------------------------------------------------------- ST15d Panel C/D
ws = wb["ST15d_Detectability"]
# idempotent: drop any previously appended Panel C/D block
for i, row in enumerate(ws.iter_rows(min_col=1, max_col=1), start=1):
    if isinstance(row[0].value, str) and row[0].value.startswith("Panel C. Ordinal abundance"):
        ws.delete_rows(i, ws.max_row - i + 1)
        break
last = ws.max_row
r = last + 2
r = put(ws, r, ["Panel C. Ordinal abundance dose-response of CRISPR-hit status "
                "(0 = not detected, 1-3 = LFQ tertiles among detections)"])
r = put(ws, r, ["Exposure level", "n", "CRISPR-clean rate"])
for _, row in DR.iterrows():
    r = put(ws, r, ["%d - %s" % (int(row["level"]), row["level_label"]),
                    int(row["n"]), round(float(row["crispr_clean_rate"]), 4)])
r = put(ws, r, [None])
r = put(ws, r, ["Dose-response trend OR per level", round(DRS["dose_response_trend_OR_per_level"], 3)])
r = put(ws, r, ["Dose-response trend 95% CI",
                "[%.3f, %.3f]" % tuple(DRS["dose_response_trend_CI"])])
r = put(ws, r, ["Dose-response trend P", "%.2e" % DRS["dose_response_trend_p"]])
r = put(ws, r, [None])
r = put(ws, r, ["Panel D. Within-detection abundance contrast, negative control and identifiability"])
r = put(ws, r, ["Item", "Value"])
hd = DRS["within_detected_lfq_halves_matched"]
r = put(ws, r, ["Within-detection LFQ high-vs-low matched fold", round(hd["fold"], 3)])
r = put(ws, r, ["Within-detection LFQ halves: Fisher P", "%.2e" % hd["fisher_p"]])
r = put(ws, r, ["Within-detection LFQ halves: paired-permutation P", "%.2e" % hd["perm_p"]])
r = put(ws, r, ["Within-detection LFQ halves: matched pairs (n)", hd["n_matched"]])
nc = DRS["negative_control_detection_on_response"]
r = put(ws, r, ["Negative control: detection vs holdout response fold", round(nc["fold"], 3)])
r = put(ws, r, ["Negative control: Fisher P", "%.3f" % nc["fisher_p"]])
r = put(ws, r, ["Fold reproducible from abundance IQR", round(DRS["fold_from_abundance_iqr"], 3)])
r = put(ws, r, ["Fold reproducible from abundance p5-p95 range",
                round(DRS["fold_from_abundance_full_range"], 3)])
r = put(ws, r, ["Identifiability", DRS["abundance_identifiability"]], wrap_cols=(2,))
ws.cell(row=r - 1, column=2).alignment = wrap
ws.cell(row=last + 2, column=1).font = bold
print("ST15d extended to row:", r - 1)

# ---------------------------------------------------------------- Contents
contents = wb["Contents"]
updates = {
    3: ("Frozen-panel (56-gene) RAB GTPase, RAB-associated machinery and trafficking enrichment",
        "13_postfreeze_interpretation\\a4_panel_enrichment.tsv"),
    19: ("Mouse BMM LCV detectability: within-detection gradient, high-confidence restriction, "
         "ordinal abundance dose-response and within-detection abundance contrast",
         "13_postfreeze_interpretation\\p0b_abundance_dose_response.tsv; p0b_abundance_dose_response_summary.json"),
    24: ("Annotation-opportunity and broad-membrane matched nulls for the 56-gene frozen candidate panel",
         "13_postfreeze_interpretation\\h2_annotation_opportunity_null.json"),
    25: ("Annotation profile of the 56-gene frozen candidate panel versus the universe",
         "13_postfreeze_interpretation\\h2_panel_annotation_profile.tsv"),
    28: ("Canonical RAB GTPase family, RAB-associated machinery and membrane-trafficking gene-set "
         "definitions and GO annotation source",
         "scripts\\trafficking_universe.py; 01_data_raw\\annotations\\goa_human.gaf.gz"),
}
for row, (desc, src) in updates.items():
    contents.cell(row=row, column=3, value=desc)
    contents.cell(row=row, column=4, value=src)

wb.save(XLSX)
print("wrote", XLSX.name)
print("sheets:", len(wb.sheetnames))
