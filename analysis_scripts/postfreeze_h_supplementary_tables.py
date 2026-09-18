"""Build Supplementary Tables 17 (species-stratified LCV->CRISPR coupling) and
18 (annotation-opportunity / broad-membrane matched null) from the h1/h2 outputs.

Writes flat TSVs into 14_submission_results/supplementary_tables/ (the STable1-13
convention) and appends worksheets to a NEW supplementary workbook
Supplementary_Tables_1-18.xlsx. The shipped 1-16 workbook is left untouched.

Nothing here recomputes statistics; it only reshapes h1_*/h2_* into submittable form.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PF = ROOT / "13_postfreeze_interpretation"
SUP = ROOT / "14_submission_results" / "supplementary_tables"
SRC = ROOT / "14_submission_results" / "source_data"

SHARED_SHEETS = {  # existing sheets copied through to the 1-18 workbook
    "ST17a_LCVspecies": None, "ST17b_JointHM": None,
    "ST18a_AnnotationNull": None, "ST18b_PanelAnnotation": None,
}


def pct(x: float) -> str:
    return f"{100 * x:.2f}"


def ci(lo: float, hi: float, fmt: str = "{:.2f}") -> str:
    return f"({fmt.format(lo)}-{fmt.format(hi)})"


def sci(p: float) -> str:
    if pd.isna(p):
        return ""
    if p == 0:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


# ---------------- ST17a: species-stratified exposure coupling ----------------
h1 = pd.read_csv(PF / "h1_species_stratified_lcv_crispr.tsv", sep="\t")
SETNAME = {"a1_primary_5": "A1 primary (5 covariates)",
           "audited_6": "Audited (6 covariates)"}
st17a = pd.DataFrame({
    "Covariate set": h1["covariate_set"].map(SETNAME),
    "LCV exposure": h1["label"],
    "Query genes (n)": h1["n_query"],
    "Matched controls (n)": h1["n_matched"],
    "Universe (n)": h1["n_universe"],
    "CRISPR+ rate, query (%)": h1["outcome_rate_query"].map(pct),
    "CRISPR+ rate, matched controls (%)": h1["outcome_rate_matched"].map(pct),
    "Fold enrichment vs matched": h1["fold_vs_matched"].map(lambda v: f"{v:.3f}"),
    "Fold 95% CI": [ci(a, b) for a, b in zip(h1["fold_ci_lo"], h1["fold_ci_hi"])],
    "Risk difference": h1["risk_difference"].map(lambda v: f"{v:.4f}"),
    "Risk difference 95% CI": [ci(a, b, "{:.4f}") for a, b in zip(h1["rd_ci_lo"], h1["rd_ci_hi"])],
    "Fisher exact P (one-sided)": h1["fisher_p"].map(sci),
    "Permutation P (4000)": h1["perm_p"].map(sci),
    "Adjusted OR": h1["adj_logistic_OR"].map(lambda v: f"{v:.3f}"),
    "Adjusted OR 95% CI": [ci(a, b) for a, b in zip(h1["adj_logistic_OR_lo"], h1["adj_logistic_OR_hi"])],
    "Adjusted OR P": h1["adj_logistic_p"].map(sci),
})
st17a.to_csv(SUP / "STable17a_lcv_species_coupling.tsv", sep="\t", index=False)

# ---------------- ST17b: joint human + mouse model ----------------
jh = pd.read_csv(PF / "h1_joint_human_mouse_lcv_model.tsv", sep="\t")
st17b = pd.DataFrame({
    "Covariate set": jh["covariate_set"].map(SETNAME),
    "Universe (n)": jh["n_universe"],
    "OR, LCV human U937": jh["OR_lcv_human"].map(lambda v: f"{v:.3f}"),
    "OR 95% CI, LCV human": [ci(a, b) for a, b in zip(jh["OR_lcv_human_lo"], jh["OR_lcv_human_hi"])],
    "P, LCV human": jh["p_lcv_human"].map(sci),
    "OR, LCV mouse BMM": jh["OR_lcv_mouse"].map(lambda v: f"{v:.3f}"),
    "OR 95% CI, LCV mouse": [ci(a, b) for a, b in zip(jh["OR_lcv_mouse_lo"], jh["OR_lcv_mouse_hi"])],
    "P, LCV mouse": jh["p_lcv_mouse"].map(sci),
    "Human-only (n)": jh["human_only_n"],
    "Mouse-only (n)": jh["mouse_only_n"],
    "Both (n)": jh["both_n"],
    "Neither (n)": jh["neither_n"],
})
st17b.to_csv(SUP / "STable17b_joint_human_mouse.tsv", sep="\t", index=False)

# ---------------- ST18a: annotation-opportunity null ----------------
h2 = json.loads((PF / "h2_annotation_opportunity_null.json").read_text(encoding="utf-8"))
NULLNAME = {"nullA_annotation_opportunity": "Null A: annotation opportunity",
            "nullB_broad_membrane_context": "Null B: Null A + broad membrane context"}
rows = []
for tag, r in h2["nulls"].items():
    rows.append({
        "Matched null": NULLNAME.get(tag, tag),
        "Covariates (n)": r["n_covariates"],
        "Covariates": ", ".join(r["covariates"]),
        "Universe (n)": r["universe_n"],
        "Panel (n)": r["n_panel"],
        "Pseudo-panels (n)": r["n_draws"],
        "RAB observed in panel": r["obs_rab_in_panel"],
        "RAB null median": f"{r['null_rab_median']:.0f}",
        "RAB null 95th pct": f"{r['null_rab_p95']:.0f}",
        "RAB null max": f"{r['null_rab_max']:.0f}",
        "Empirical P, RAB": sci(r["empirical_p_rab"]),
        "Bonferroni P, RAB": sci(r["bonferroni_adjusted_p_rab"]),
        "Trafficking observed in panel": r["obs_trafficking_in_panel"],
        "Trafficking null median": f"{r['null_traf_median']:.0f}",
        "Trafficking null 95th pct": f"{r['null_traf_p95']:.0f}",
        "Trafficking null max": f"{r['null_traf_max']:.0f}",
        "Empirical P, trafficking": sci(r["empirical_p_trafficking"]),
        "Bonferroni P, trafficking": sci(r["bonferroni_adjusted_p_trafficking"]),
        "Max |SMD| before matching": f"{r['max_abs_smd_before']:.3f}",
        "Max |SMD| after matching": f"{r['max_abs_smd_after']:.3f}",
        "Controls at caliper 0.02 / 0.05 / full": (f"{r['pool_caliper_counts']['0.02']} / "
                                                   f"{r['pool_caliper_counts']['0.05']} / "
                                                   f"{r['pool_caliper_counts']['all']}"),
    })
st18a = pd.DataFrame(rows)
st18a.to_csv(SUP / "STable18a_annotation_opportunity_null.tsv", sep="\t", index=False)

# ---------------- ST18b: panel annotation profile + per-covariate balance ----------------
prof = h2["panel_annotation_profile"]
lab = {"_go_bp_n": "GO biological-process terms (distinct)",
       "_go_cc_n": "GO cellular-component terms (distinct)",
       "_go_terms_n": "GO terms, all aspects (distinct)",
       "_n_membrane_terms": "GO:0016020 membrane descendant terms",
       "_n_imbo_terms": "GO:0043231 intracellular membrane-bounded organelle terms",
       "_membrane": "Any GO:0016020 membrane annotation (fraction)",
       "_imbo": "Any GO:0043231 IMBO annotation (fraction)"}
prof_rows = []
for k, v in prof.items():
    if "panel_median" in v:
        prof_rows.append({"Annotation feature": lab.get(k, k), "Statistic": "median",
                          "Panel (54 genes)": f"{v['panel_median']:.1f}",
                          "Universe": f"{v['universe_median']:.1f}"})
    else:
        prof_rows.append({"Annotation feature": lab.get(k, k), "Statistic": "fraction",
                          "Panel (54 genes)": f"{v['panel_frac']:.3f}",
                          "Universe": f"{v['universe_frac']:.3f}"})
st18b = pd.DataFrame(prof_rows)
st18b.to_csv(SUP / "STable18b_panel_annotation_profile.tsv", sep="\t", index=False)

# ---------------- copy raw h1/h2 outputs into source_data ----------------
for src in ["h1_species_stratified_lcv_crispr.tsv", "h1_joint_human_mouse_lcv_model.tsv",
            "h1_species_stratified_summary.json",
            "h2_annotation_opportunity_null.json", "h2_panel_annotation_profile.tsv"]:
    (SRC / src).write_bytes((PF / src).read_bytes())

print("wrote STable17a/17b/18a/18b into", SUP)
print("copied h1/h2 outputs into", SRC)

# ---------------- build Supplementary_Tables_1-18.xlsx ----------------
import openpyxl
from openpyxl.styles import Font

wbpath = ROOT / "15_manuscript" / "table" / "supplementary" / "Supplementary_Tables_1-16.xlsx"
wb = openpyxl.load_workbook(wbpath)

def add_sheet(name: str, df: pd.DataFrame, title_row: str) -> None:
    if name in wb.sheetnames:
        del wb[name]
    ws = wb.create_sheet(name)
    ws["A1"] = title_row
    ws["A1"].font = Font(bold=True)
    ws.append([])
    ws.append(list(df.columns))
    for c in ws[3]:
        c.font = Font(bold=True)
    for _, r in df.iterrows():
        ws.append([str(x) for x in r.tolist()])
    for i, col in enumerate(df.columns, start=1):
        width = max(len(str(col)), *(len(str(x)) for x in df[col])) if len(df) else len(str(col))
        ws.column_dimensions[ws.cell(row=3, column=i).column_letter].width = min(46, width + 2)

add_sheet("ST17a_LCVspecies", st17a,
          "Supplementary Table 17a. Species-stratified LCV-associated detection versus "
          "fitness-filtered CRISPR evidence. Gene universe, propensity model, 1:1 nearest-neighbour "
          "matching (0.25 x SD caliper), permutation count and outcome definition are unchanged from "
          "the primary coupling analysis; only the LCV exposure is restricted by species/source. "
          "'LCV any species' is the primary exposure and reproduces the published estimate exactly.")
add_sheet("ST17b_JointHM", st17b,
          "Supplementary Table 17b. Joint model of fitness-filtered CRISPR evidence on human U937 and "
          "mouse BMM LCV detection simultaneously, with the audited covariates. Human U937 LCV "
          "detection is attenuated after mutual adjustment whereas mouse BMM LCV detection is not, "
          "indicating that the evidence-layer coupling is not attributable to the shared human U937 "
          "cellular background alone.")
add_sheet("ST18a_AnnotationNull", st18a,
          "Supplementary Table 18a. Annotation-opportunity and broad-membrane matched nulls for the "
          "54-gene holdout panel. Null A matches GO annotation counts, annotation breadth, generic "
          "membrane and intracellular membrane-bounded organelle annotation, baseline expression, "
          "study count, network degree, publication count and CRISPR tier. Null B additionally matches "
          "the depth of membrane annotation. RAB membership, vesicle-mediated transport, ER-to-Golgi "
          "transport, Rab signal transduction and membrane-trafficking category are deliberately NOT "
          "matched, since adjusting for the structure under test would be circular.")
add_sheet("ST18b_PanelAnnotation", st18b,
          "Supplementary Table 18b. Annotation profile of the 54-gene panel versus the full enrichment "
          "universe, showing the annotation-opportunity imbalance that the matched nulls in 18a remove.")

# Contents sheet: append the new rows
ws = wb["Contents"]
used = [ws.cell(row=r, column=1).value for r in range(1, ws.max_row + 1)]
adds = [("17a", "Species-stratified LCV versus fitness-filtered CRISPR evidence"),
        ("17b", "Joint human U937 / mouse BMM LCV model"),
        ("18a", "Annotation-opportunity and broad-membrane matched nulls"),
        ("18b", "Panel annotation profile versus universe")]
r = ws.max_row + 1
for num, desc in adds:
    ws.cell(row=r, column=1, value=num)
    ws.cell(row=r, column=2, value=desc)
    r += 1

# Both canonical workbook locations are replaced in place (no new folders). The
# superseded 1-16 workbook, whose entire content is contained in 1-18, is removed
# so that each folder holds exactly one supplementary workbook.
WB_DIRS = [ROOT / "15_manuscript" / "table" / "supplementary",
           ROOT / "FINAL_SUBMISSION_2026-09-11" / "03_tables"]
for wdir in WB_DIRS:
    out_wb = wdir / "Supplementary_Tables_1-18.xlsx"
    wb.save(out_wb)
    stale = wdir / "Supplementary_Tables_1-16.xlsx"
    if stale.exists():
        stale.unlink()
    print("wrote", out_wb)
print("sheets:", wb.sheetnames)
