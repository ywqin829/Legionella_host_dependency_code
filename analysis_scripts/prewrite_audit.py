from __future__ import annotations

"""Pre-writing integrity audit: independently recompute headline numbers from RAW data
and compare against the values stated in reports/legends. PASS/FAIL output."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import RAB_GENES, ALL_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

results = []

def check(name, stated, computed, tol=0.02):
    ok = abs(stated - computed) <= max(tol, tol * abs(stated)) if isinstance(stated, (int, float)) else stated == computed
    results.append((name, stated, computed, "PASS" if ok else "FAIL"))

# ---- 1. evidence layers ----
conv = pd.read_csv(ROOT / "07_convergence" / "primary_infection_convergence.tsv", sep="\t")
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
n_meta = int((conv.stable_cross_cohort.astype(str).str.upper() == "TRUE").sum())
ev["fdr05"] = ev["adj.P.Val"] < 0.05
n_geo = ev[ev.fdr05].groupby("human_symbol")["dataset"].nunique()
n_core = int((n_geo >= 3).sum())
n_both = int((conv.loc[conv.stable_cross_cohort.astype(str).str.upper() == "TRUE", "human_symbol"].isin(n_geo[n_geo >= 3].index)).sum())
check("meta_convergent_2601", 2601, n_meta)
check("ge3_core_490", 490, n_core)
check("both_432", 432, n_both)

# ---- 2. four-primary z ----
h = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
w = h[(h.contrast_id == "WT_vs_UI") & (h.human_symbol.isin(["RAB1A", "RAB10", "ACTR3", "VAPA"]))].set_index("human_symbol")
zs = {g: float(np.sign(w.loc[g, "logFC"]) * stats.norm.ppf(np.clip(w.loc[g, "PValue"] / 2, 1e-300, 0.5 - 1e-300))) for g in w.index}
check("RAB1A_z_neg4.95", -4.95, round(zs["RAB1A"], 2), 0.01)
check("RAB10_z_neg4.81", -4.81, round(zs["RAB10"], 2), 0.01)
check("ACTR3_z_neg1.46", -1.46, round(zs["ACTR3"], 2), 0.01)
check("VAPA_z_neg0.88", -0.88, round(zs["VAPA"], 2), 0.01)

# ---- 3. A4 panel enrichment (independent Fisher) ----
panel = sorted(set(h.loc[h.contrast_id == "WT_vs_UI", "human_symbol"]))
atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
bg = set(atlas.human_symbol.dropna())
panel_set, rab_set, traf_set = set(panel), set(RAB_GENES), set(ALL_GENES)
for label, mem in [("RAB", rab_set), ("trafficking", traf_set)]:
    a = len(mem & panel_set); b = len(panel_set - mem); c = len((bg - panel_set) & mem); d = len((bg - panel_set) - mem)
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    fold = (a / len(panel_set)) / ((c + a) / len(bg))
    print(f"  A4 {label}: in_panel={a}, bg_rate={(c+a)/len(bg):.5f}, fold={fold:.1f}x, p={p:.2e}")

# ---- 4. B1 LCV-state rates (independent) ----
lcv_src = pd.read_csv(ROOT / "08_evidence_atlas" / "lcv_spatial_source_records.tsv", sep="\t")
cond = lcv_src.groupby("human_symbol")["condition"].apply(lambda s: set(s.astype(str))).to_dict()
atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_caveat"] = atlas["fitness_artifact_compatible"].astype(str).str.upper().isin(["TRUE", "1"])
atlas["_clean"] = atlas["_hit"] & (~atlas["_caveat"])
ankb_only = [g for g, c in cond.items() if any("ANKB" in x for x in c) and not any("WT_LCV" in x for x in c)]
wt_ankb = [g for g, c in cond.items() if any("WT_LCV" in x for x in c) and any("ANKB" in x for x in c)]
a_map = dict(zip(atlas.human_symbol, atlas._clean))
a_rate = sum(a_map.get(g, False) for g in ankb_only) / len(ankb_only)
w_rate = sum(a_map.get(g, False) for g in wt_ankb) / len(wt_ankb)
print(f"  B1 ankB_only clean={sum(a_map.get(g,False) for g in ankb_only)}/{len(ankb_only)}={a_rate:.3f}; WT+ankB={sum(a_map.get(g,False) for g in wt_ankb)}/{len(wt_ankb)}={w_rate:.3f}")

# ---- 5. matched-null claim ----
mn = json.load(open(OUT / "d2_matched_null_summary.json"))
check("matched_null_0_of_12", 0, mn["n_contrasts_still_enriched_after_matched_null_FDR"])

# ---- 6. RAB counts in panel / clean CRISPR ----
a3 = pd.read_csv(OUT / "a3_trafficking_gene_evidence_matrix.tsv", sep="\t")
rab_clean = int((a3.family.eq("RAB") & a3.crispr_clean.astype(str).str.upper().isin(["TRUE", "True", "1"])).sum())
rab_in_panel = int((a3.family.eq("RAB") & a3.in_frozen_panel.astype(str).str.upper().isin(["TRUE", "True", "1"])).sum())
check("RAB_clean_CRISPR_13", 13, rab_clean)
check("RAB_in_panel_10", 10, rab_in_panel)

# ---- output ----
df = pd.DataFrame(results, columns=["check", "stated", "computed", "result"])
print("\n==== AUDIT RESULT ====")
print(df.to_string(index=False))
n_fail = int((df.result == "FAIL").sum())
print(f"\n{n_fail} FAIL / {len(df)} checks")
df.to_csv(OUT / "prewrite_audit_numbers.tsv", sep="\t", index=False)
