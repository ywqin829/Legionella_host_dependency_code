from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import GENE_TO_FAMILY, ALL_GENES, RAB_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
holdout = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
lcv_src = pd.read_csv(ROOT / "08_evidence_atlas" / "lcv_spatial_source_records.tsv", sep="\t")


def B(s):
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


atlas["_ppi_AB"] = B(atlas["ppi_AB_flag"])
atlas["_apms"] = B(atlas["apms_high_confidence"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fit_ok"] = B(atlas["fitness_artifact_compatible"])
atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")

# per-gene LCV condition list
cond_map = lcv_src.groupby("human_symbol")["condition"].apply(lambda s: ",".join(sorted(set(s.astype(str))))).to_dict()

# per-gene holdout (WT vs UI)
hu = holdout[holdout.contrast_id == "WT_vs_UI"].set_index("human_symbol")
hu_dir = hu["observed_direction"].to_dict()
hu_fdr = hu["FDR"].to_dict()
hu_lfc = hu["logFC"].to_dict()

rows = []
for g in ALL_GENES:
    a = atlas[atlas.human_symbol == g]
    if len(a) == 0:
        continue
    a = a.iloc[0]
    mech = []
    for route in ["mechanism_dot_icm_convergent", "mechanism_rna_rpf_interaction_convergent",
                   "mechanism_single_effector_or_t2ss_convergent", "mechanism_translation_effector_module_convergent"]:
        if B(pd.Series([a[route]]))[0]:
            mech.append(route.replace("mechanism_", "").replace("_convergent", ""))
    rows.append({
        "gene": g,
        "family": GENE_TO_FAMILY[g],
        "in_frozen_panel": bool(g in set(hu.index)),
        "ppi_AB": bool(a["_ppi_AB"]),
        "apms_high": bool(a["_apms"]),
        "crispr_hit": bool(a["_hit"]),
        # fitness_artifact_compatible=True means hit+baseline-essential (caveat present)
        "crispr_clean": bool(a["_hit"] and (not a["_fit_ok"])),
        "fitness_caveat": bool(a["_fit_ok"] and a["_hit"]),
        "stable_cross_cohort": bool(a["_resp"]),
        "mechanism_routes": ";".join(mech),
        "lcv_spatial": bool(a["_lcv"]),
        "lcv_conditions": cond_map.get(g, ""),
        "final_tier": a["evidence_tier"] if pd.notna(a["evidence_tier"]) else "",
        "holdout_dir": hu_dir.get(g, ""),
        "holdout_FDR": round(float(hu_fdr[g]), 4) if g in hu_fdr else "",
        "holdout_logFC": round(float(hu_lfc[g]), 3) if g in hu_lfc else "",
    })

df = pd.DataFrame(rows)
df.to_csv(OUT / "a3_trafficking_gene_evidence_matrix.tsv", sep="\t", index=False)

# --- summaries ---
rab = df[df.family.isin(["RAB", "RAB_effector_GEF", "RAB_GAP"])].copy()
rab["n_evidence"] = rab[["ppi_AB", "apms_high", "crispr_hit", "stable_cross_cohort", "lcv_spatial"]].sum(axis=1)
rab = rab.sort_values(["n_evidence", "gene"], ascending=[False, True])

print("=== RAB-family evidence atlas (sorted by evidence count) ===")
print(rab[["gene", "family", "ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "fitness_caveat",
           "stable_cross_cohort", "lcv_spatial", "lcv_conditions", "final_tier", "holdout_dir", "holdout_FDR"]].to_string(index=False))

# per-family saturation
fam = df.groupby("family")[["ppi_AB", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]].mean().round(3)
fam = fam.reindex(sorted(fam.index))
fam.to_csv(OUT / "a3_family_evidence_saturation.tsv", sep="\t")

print("\n=== family saturation (fraction of family genes) ===")
print(fam.to_string())

stats = {
    "trafficking_genes_in_atlas": int(len(df)),
    "rab_genes_in_atlas": int(len(rab)),
    "rab_in_frozen_panel": int(rab.in_frozen_panel.sum()),
    "rab_ppi_AB": int(rab.ppi_AB.sum()),
    "rab_crispr_clean": int(rab.crispr_clean.sum()),
    "rab_lcv": int(rab.lcv_spatial.sum()),
    "rab_stable_cross": int(rab.stable_cross_cohort.sum()),
}
(OUT / "a3_stats.json").write_text(__import__("json").dumps(stats, indent=2), encoding="utf-8")
print("\n=== stats ===")
print(__import__("json").dumps(stats, indent=2))
