from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

# ============ 1. AUTHORITATIVE MUTANT IDENTITY TABLE ============
confirmed = {
    "A": ("lpw00071-lpw00221", 16, "cegC1", "CONFIRMED", "PRINT E-Table1 + 130b GBK FR687201.1"),
    "B": ("lpw04831-lpw04861", 4, "legA7", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
    "C": ("lpw11451-lpw11581", 14, "lem4;lem5", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
    "D": ("lpw16861-lpw16971", 12, "legL3", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
    "E": ("lpw19951-lpw20201", 26, "lem16;lem17;ralF;legC4;ltpG;lirA;lirB;lirE;lirF;pieE", "CONFIRMED", "PRINT E-Table1 + 130b GBK; Wibawa2022 confirms legC4=lpw20041"),
    "F": ("lpw23101-lpw23561", 45, "legK2;legAU13;sdeC;sde;sidJ;sdeB;sdeA;lem19;legS2", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
    "H": ("lpw26521-lpw26861", 35, "ceg30;legA14;legA15;sidM;sidD", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
    "I": ("lpw27241-lpw27551", 32, "ceg32;sdjA;sdeD;sdcA;sidC;lem26;lem27", "CONFIRMED", "PRINT E-Table1 + 130b GBK"),
}
status = {"A": "PASS", "B": "PASS", "C": "PASS", "D": "PASS", "E": "FAIL", "F": "FAIL", "H": "in_FGHI", "I": "FAIL"}

# SDRF official labels (from E-MTAB-16655.sdrf.txt)
sdrf_label = {"A": "Delta-A Lp130b infection", "B": "Delta-B Lp130b infection", "C": "Delta-C Lp130b infection",
              "D": "Delta-D Lp130b infection", "E": "Delta-E Lp130b infection", "F": "Delta-F Lp130b infection",
              "I": "Delta-I Lp130b infection", "FGHI": "Delta-FGHI Lp130b infection", "O": "Delta-O Lp130b infection",
              "P": "Delta-P Lp130b infection", "R": "Delta-R Lp130b infection", "S": "Delta-S Lp130b infection",
              "dotA": "Delta-dotA Lp 130b infection", "WT": "WT Lp 130b infection", "UI": "Uninfected"}

rows = []
for m in ["A", "B", "C", "D", "E", "F", "I"]:
    iv, n, eff, conf, src = confirmed[m]
    rows.append({"sdrf_official_label": sdrf_label[m], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "SINGLE-REGION DELETION",
                 "region_or_gene": iv, "n_genes_deleted": n, "deleted_genes_effectors": eff,
                 "holdout_status": status[m], "mapping_source": src, "confidence": conf,
                 "nested_note": "Also present within FGHI combined strain"})
# H (in library, not separate in RNA-seq)
iv, n, eff, conf, src = confirmed["H"]
rows.append({"sdrf_official_label": "N/A (not a separate RNA-seq condition)", "sample_n": 0, "in_rnaseq": "NO (within FGHI only)",
             "mutant_type": "SINGLE-REGION DELETION (library)", "region_or_gene": iv, "n_genes_deleted": n,
             "deleted_genes_effectors": eff, "holdout_status": "in_FGHI", "mapping_source": src, "confidence": conf,
             "nested_note": "Strain exists in library (ExtFig1/print E-Table1) but not assayed alone in E-MTAB-16655"})
# G
rows.append({"sdrf_official_label": "N/A (not a separate RNA-seq condition)", "sample_n": 0, "in_rnaseq": "NO (within FGHI only)",
             "mutant_type": "SINGLE-REGION DELETION (library)", "region_or_gene": "UNKNOWN", "n_genes_deleted": None,
             "deleted_genes_effectors": "UNKNOWN", "holdout_status": "in_FGHI", "mapping_source": "PRINT ExtFig1 (exists); coords PENDING",
             "confidence": "EXISTS_BUT_COORDS_PENDING", "nested_note": "Confirmed to exist (ExtFig1) but no coordinate/effector list available"})
# FGHI composite
rows.append({"sdrf_official_label": sdrf_label["FGHI"], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "COMBINED DELETION F+G+H+I",
             "region_or_gene": "F(lpw23101-23561)+G(?)+H(lpw26521-26861)+I(lpw27241-27551)", "n_genes_deleted": "45+?+35+32",
             "deleted_genes_effectors": "union of F/H/I; G pending", "holdout_status": "PASS",
             "mapping_source": "PRINT E-Table1 (F/H/I) + ExtFig1 (G exists)", "confidence": "COMPOSITE (F/H/I confirmed, G pending)",
             "nested_note": "F and I also assayed separately; G and H not separate conditions"})
# O/P/R/S
for m in ["O", "P", "R", "S"]:
    rows.append({"sdrf_official_label": sdrf_label[m], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "SINGLE-REGION DELETION",
                 "region_or_gene": "UNKNOWN", "n_genes_deleted": None, "deleted_genes_effectors": "UNKNOWN",
                 "holdout_status": status.get(m, "?"), "mapping_source": "SDRF only; content PENDING_PRIMARY_SOURCE",
                 "confidence": "MAPPING_PENDING_PRIMARY_SOURCE", "nested_note": "Present in E-MTAB-16655; deletion content not in any public source"})
# dotA / WT / UI
rows.append({"sdrf_official_label": sdrf_label["dotA"], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "T4SS-INACTIVE CONTROL",
             "region_or_gene": "lpw29401 (dotA gene)", "n_genes_deleted": 1, "deleted_genes_effectors": "dotA (Dot/Icm T4SS ATPase)",
             "holdout_status": "PASS (Dot/Icm subset)", "mapping_source": "PRINT E-Table1", "confidence": "CONFIRMED", "nested_note": ""})
rows.append({"sdrf_official_label": sdrf_label["WT"], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "WT", "region_or_gene": "full genome",
             "n_genes_deleted": 0, "deleted_genes_effectors": "none", "holdout_status": "baseline", "mapping_source": "SDRF", "confidence": "CONFIRMED", "nested_note": ""})
rows.append({"sdrf_official_label": sdrf_label["UI"], "sample_n": 3, "in_rnaseq": "YES", "mutant_type": "UNINFECTED", "region_or_gene": "n/a",
             "n_genes_deleted": 0, "deleted_genes_effectors": "none", "holdout_status": "baseline", "mapping_source": "SDRF", "confidence": "CONFIRMED", "nested_note": ""})

auth = pd.DataFrame(rows)
auth.to_csv(OUT / "E-MTAB-16655_mutant_identity_AUTHORITATIVE.tsv", sep="\t", index=False)

# ============ 2. PER-CANDIDATE EVIDENCE ANNOTATION (leading candidates) ============
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz", sep="\t", compression=None, low_memory=False)
conv = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
tier = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")
a2 = pd.read_csv(OUT / "a2_evidence_fragility_matrix.tsv", sep="\t")

stable = set(conv.loc[conv.stable_cross_cohort.astype(str).str.upper().eq("TRUE"), "human_symbol"])
ev["fdr05"] = ev["adj.P.Val"] < 0.05
n_geo = ev[ev.fdr05].groupby("human_symbol")["dataset"].nunique()
core_ge3 = set(n_geo[n_geo >= 3].index)
# human/mouse per candidate
hm = ev.groupby("human_symbol")["species"].apply(lambda s: (int((s == "HUMAN").sum()), int((s == "MOUSE").sum()))).to_dict()
# a2 leave-one-out collapse
a2_collapse = dict(zip(a2.candidate, a2.collapse_datasets.fillna("")))

leads = ["RAB1A", "RAB10", "ACTR3", "VAPA"] + tier.loc[tier.final_tier.eq("TIER2_MULTIMODAL"), "human_symbol"].tolist()
rows2 = []
for g in leads:
    r = tier[tier.human_symbol == g].iloc[0]
    nh, nm = hm.get(g, (0, 0))
    rows2.append({
        "gene": g, "final_tier": r.final_tier,
        "n_datasets_fdr05": int(n_geo.get(g, 0)),
        "n_human_studies": nh, "n_mouse_studies": nm,
        "in_meta_convergent_2601": g in stable,
        "in_ge3_replicated_490": g in core_ge3,
        "leave1_dataset_collapse": a2_collapse.get(g, ""),
        "human_only_stable": bool(a2.loc[a2.candidate == g, "human_only_stable"].iloc[0]) if (a2.candidate == g).any() else None,
    })
cand = pd.DataFrame(rows2)
cand.to_csv(OUT / "d5_candidate_evidence_annotation.tsv", sep="\t", index=False)

summary = {
    "meta_convergent_2601": len(stable),
    "ge3_replicated_core_490": len(core_ge3),
    "both_429": len(stable & core_ge3),
    "authoritative_mutant_table_n": int(len(auth)),
    "candidate_annotation_n": int(len(cand)),
}
(OUT / "e1_authoritative_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print("=== authoritative summary ===")
print(json.dumps(summary, indent=2, ensure_ascii=False))
print("\n=== candidate evidence annotation (leading) ===")
print(cand.head(12).to_string(index=False))
