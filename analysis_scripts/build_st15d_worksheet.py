"""Add Supplementary Table 15d (mouse BMM LCV detectability analysis) to the workbook.

Reads the P0-B outputs and appends a ST15d_Detectability worksheet plus a
matching Contents row, matching the formatting of ST15c.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "FINAL_SUBMISSION_2026-10-04"
OUT = ROOT / "13_postfreeze_interpretation"
XLSX = PKG / "03_tables" / "Supplementary_Tables_1-19.xlsx"

grad = pd.read_csv(OUT / "p0b_detectability_gradient.tsv", sep="\t")
import json
summary = json.loads((OUT / "p0b_detectability_adjusted_summary.json").read_text(encoding="utf-8"))
fm = summary["full_mouse_exposure"]
hc = summary["highconf_mouse_exposure"]

wb = openpyxl.load_workbook(XLSX)
if "ST15d_Detectability" in wb.sheetnames:
    del wb["ST15d_Detectability"]
ws = wb.create_sheet("ST15d_Detectability")

bold = Font(bold=True)
wrap = Alignment(wrap_text=True, vertical="top")

rows: list = []
rows.append(["Supplementary Table 15d. Detectability analysis of the mouse BMM LCV layer"])
rows.append(["The LCV-associated exposure is defined by mass-spectrometry detection, which rises with protein abundance. "
             "The mouse BMM LCV proteome (PMID 26829557) is quantitative (Average LFQ intensity and sequence coverage per protein), "
             "so the confound can be tested directly there. Panel A asks whether CRISPR-hit status rises with detectability inside the "
             "detected set (adjusted for the audited covariates). Panel B restricts the mouse exposure to high-confidence identifications "
             "(sequence coverage at or above the median among detections) and recomputes the confounder-matched enrichment against the same "
             "not-detected controls. The human U937 LCV proteomes (PMID 25369898) are presence-only and carry no abundance or coverage, so they "
             "cannot be detectability-adjusted."])
rows.append(["Canonical source: 13_postfreeze_interpretation\\p0b_detectability_gradient.tsv and p0b_detectability_adjusted_summary.json"])
rows.append([None])
rows.append(["Panel A. Within-detection logistic gradient of CRISPR-hit status on detectability (mouse-detected LCV proteins in the coupling universe)"])
rows.append(["term", "OR", "ci_lo", "ci_hi", "p"])
for _, r in grad.iterrows():
    rows.append([str(r["term"]), float(r["OR"]), float(r["ci_lo"]), float(r["ci_hi"]), float(r["p"])])
rows.append([None])
rows.append(["Panel B. Mouse BMM LCV versus fitness-filtered CRISPR matched enrichment (audited-covariate matching, paired permutation)"])
rows.append(["Item", "Value"])
rows.append(["Coupling universe (n)", summary["n_universe"]])
rows.append(["Mouse-detected in coupling universe", summary["n_mouse_detected_in_universe"]])
rows.append(["Median sequence coverage among detections (%)", summary["mouse_seq_cov_median_among_detected"]])
rows.append(["Full exposure - matched enrichment fold", round(fm["fold"], 4)])
rows.append(["Full exposure - Fisher P", fm["fisher_p"]])
rows.append(["Full exposure - paired-permutation P", fm["perm_p"]])
rows.append(["Full exposure - matched pairs (n)", fm["n_matched"]])
rows.append(["Full exposure - CRISPR-clean rate, mouse-detected query", round(fm["rate_query"], 4)])
rows.append(["Full exposure - CRISPR-clean rate, matched control", round(fm["rate_matched"], 4)])
rows.append(["High-confidence universe (n)", summary["n_universe_highconf"]])
rows.append(["High-confidence exposure (coverage >= median, n)", summary["n_query_highconf"]])
rows.append(["Restricted exposure - matched enrichment fold", round(hc["fold"], 4)])
rows.append(["Restricted exposure - Fisher P", hc["fisher_p"]])
rows.append(["Restricted exposure - paired-permutation P", hc["perm_p"]])
rows.append(["Restricted exposure - matched pairs (n)", hc["n_matched"]])
rows.append(["Restricted exposure - CRISPR-clean rate, mouse-detected query", round(hc["rate_query"], 4)])
rows.append(["Restricted exposure - CRISPR-clean rate, matched control", round(hc["rate_matched"], 4)])

for r_i, r in enumerate(rows, start=1):
    for c_i, v in enumerate(r, start=1):
        cell = ws.cell(row=r_i, column=c_i, value=v)
        if r_i in (1, 5, 9) and c_i == 1:
            cell.font = bold
        if r_i in (2,):
            cell.alignment = wrap

for col, w in zip("ABCDE", [54, 14, 14, 14, 16]):
    ws.column_dimensions[col].width = w
ws.column_dimensions["A"].width = 54

# Contents row: insert after the ST15c row (currently row 18)
contents = wb["Contents"]
st15c_row = None
for i, r in enumerate(contents.iter_rows(values_only=True), start=1):
    if r[1] == "ST15c_CorrRobust":
        st15c_row = i
        break
assert st15c_row is not None, "ST15c Contents row not found"
contents.insert_rows(st15c_row + 1)
contents.cell(row=st15c_row + 1, column=1, value="Supplementary Table 15")
contents.cell(row=st15c_row + 1, column=2, value="ST15d_Detectability")
contents.cell(row=st15c_row + 1, column=3, value=(
    "Mouse BMM LCV detectability analysis: within-detection CRISPR gradient on LFQ/coverage "
    "and high-confidence-detection-restricted matched enrichment"))
contents.cell(row=st15c_row + 1, column=4, value=(
    "13_postfreeze_interpretation\\p0b_detectability_gradient.tsv; p0b_detectability_adjusted_summary.json"))

wb.save(XLSX)
print("sheets:", wb.sheetnames)
print("contents rows after insert:")
for r in wb["Contents"].iter_rows(values_only=True):
    print(r[:3])
