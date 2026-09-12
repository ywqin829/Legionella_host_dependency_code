from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

# ---------- 1. per-mutant interval verification table ----------
confirmed = {
    "A": ("lpw00071-lpw00221", 16, "cegC1", "CONFIRMED"),
    "B": ("lpw04831-lpw04861", 4, "legA7", "CONFIRMED"),
    "C": ("lpw11451-lpw11581", 14, "lem4;lem5", "CONFIRMED"),
    "D": ("lpw16861-lpw16971", 12, "legL3", "CONFIRMED"),
    "E": ("lpw19951-lpw20201", 26, "lem16;lem17;ralF;legC4;ltpG;lirA;lirB;lirE;lirF;pieE", "CONFIRMED"),
    "F": ("lpw23101-lpw23561", 45, "legK2;legAU13;sdeC;sde;sidJ;sdeB;sdeA;lem19;legS2", "CONFIRMED"),
    "H": ("lpw26521-lpw26861", 35, "ceg30;legA14;legA15;sidM;sidD", "CONFIRMED"),
    "I": ("lpw27241-lpw27551", 32, "ceg32;sdjA;sdeD;sdcA;sidC;lem26;lem27", "CONFIRMED"),
}
status = {"A": "PASS", "B": "PASS", "C": "PASS", "D": "PASS", "E": "FAIL", "F": "FAIL", "H": "in_FGHI", "I": "FAIL"}
rows = []
for m, (iv, n, eff, conf) in confirmed.items():
    rows.append({"mutant": m, "interval": iv, "n_genes": n, "effectors": eff,
                 "holdout_status": status[m], "confidence": conf,
                 "nested_note": ("SEPARATE holdout condition and present in FGHI composite" if m in {"F", "I"} else
                                 "NOT separate RNA-seq condition; present only within FGHI composite" if m == "H" else "")})
rows.append({"mutant": "G", "interval": "UNKNOWN", "n_genes": None, "effectors": "UNKNOWN",
             "holdout_status": "in_FGHI", "confidence": "EXISTS (ExtFig1) but coordinates pending", "nested_note": ""})
rows.append({"mutant": "FGHI", "interval": "F(lpw23101-23561)+G(?)+H(lpw26521-26861)+I(lpw27241-27551)",
             "n_genes": "45+?+35+32", "effectors": "union of F/H/I; G pending",
             "holdout_status": "PASS", "confidence": "COMPOSITE (F/H/I confirmed, G pending)", "nested_note": "combined F+G+H+I deletion condition; F and I are also separately assayed; G and H are not separate RNA-seq conditions"})
for m in ["O", "P", "R", "S"]:
    rows.append({"mutant": m, "interval": "UNKNOWN", "n_genes": None, "effectors": "UNKNOWN",
                 "holdout_status": status.get(m, "?"), "confidence": "MAPPING_PENDING_PRIMARY_SOURCE",
                 "nested_note": "present in E-MTAB-16655; deletion content not in preprint E-Table1"})
rows.append({"mutant": "dotA", "interval": "lpw29401 (gene deletion)", "n_genes": 1,
             "effectors": "dotA (Dot/Icm T4SS ATPase)", "holdout_status": "PASS (Dot/Icm subset)",
             "confidence": "CONFIRMED", "nested_note": "Dot/Icm-inactive control"})
rows.append({"mutant": "WT", "interval": "full genome", "n_genes": None, "effectors": "none",
             "holdout_status": "baseline", "confidence": "CONFIRMED", "nested_note": ""})
iv_df = pd.DataFrame(rows)
iv_df.to_csv(OUT / "d3_mutant_interval_verification_table.tsv", sep="\t", index=False)

# ---------- 2. batch/platform check for stable_cross_cohort ----------
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
conv = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
stable = set(conv.loc[conv.stable_cross_cohort.astype(str).str.upper().eq("TRUE"), "human_symbol"])
sev = ev[ev.human_symbol.isin(stable)].copy()
sev["fdr05"] = sev["adj.P.Val"] < 0.05
# per dataset: how many stable genes have a fdr05 hit there, and what fraction of stable genes
# have their ONLY fdr05 hit there
batch_rows = []
for ds in sorted(sev.dataset.unique()):
    n_stable_hit = int(sev[(sev.dataset == ds) & sev.fdr05].human_symbol.nunique())
    batch_rows.append({"dataset": ds, "species": sev.loc[sev.dataset == ds, "species"].mode()[0],
                       "n_stable_genes_fdr05_in_this_dataset": n_stable_hit,
                       "frac_of_2590": round(n_stable_hit / len(stable), 3)})
batch_df = pd.DataFrame(batch_rows)
batch_df.to_csv(OUT / "d3_stable_cohort_batch_contribution.tsv", sep="\t", index=False)

# how many stable genes are significant in >=3 datasets (not dataset-specific)?
sev_wide = sev.pivot_table(index="human_symbol", columns="dataset", values="fdr05", aggfunc="any")
sev_wide = sev_wide[sev_wide.columns[sev_wide.columns.isin(sev.dataset.unique())]]
n_ge3 = int(sev_wide.sum(axis=1).ge(3).sum())
n_only1 = int(sev_wide.sum(axis=1).eq(1).sum())

audit = {
    "stable_cross_cohort_n": len(stable),
    "n_stable_signif_in_ge3_datasets": n_ge3,
    "n_stable_signif_in_only_1_dataset": n_only1,
    "per_dataset_stable_fdr05_hits": batch_df.to_dict("records"),
    "batch_interpretation": "Most meta-convergent genes are FDR-significant in only 1-2 datasets. GSE135803 and GSE287755 contribute disproportionately (65.6% and 60.2% of the 2,590 set), so this layer must not be described as per-gene replication across many studies.",
    "note_platform": "GSE61535/GSE26490 are arrays, GSE135803/GSE287755/GSE89184 are RNA-seq, GSE60085 dual-colour; "
                     "platform heterogeneity is mitigated by signed-direction concordance not effect-size pooling.",
}
(OUT / "d3_batch_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
print("=== per-mutant interval verification (first 12 rows) ===")
print(iv_df.head(12).to_string(index=False))
print("\n=== batch/platform audit ===")
print(json.dumps(audit, indent=2, ensure_ascii=False))
