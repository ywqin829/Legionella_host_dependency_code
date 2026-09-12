from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import GENE_TO_FAMILY, RAB_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
hold = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")

# measured frozen panel (54)
panel = sorted(set(hold.loc[hold.contrast_id == "WT_vs_UI", "human_symbol"]))
background = sorted(atlas.human_symbol.dropna().unique())

panel_set = set(panel)
bg_set = set(background)
rab_set = set(RAB_GENES)
traf_set = set(GENE_TO_FAMILY)

rab_in_panel = panel_set & rab_set
traf_in_panel = panel_set & traf_set


def fisher_enrichment(members: set, query: set, universe: set, label: str) -> dict:
    a = len(query & members)          # query genes that are members
    b = len(query - members)          # query genes not members
    c = len((universe - query) & members)  # background members
    d = len((universe - query) - members)  # background non-members
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "label": label,
        "n_query": len(query), "n_universe": len(universe),
        "members_in_query": a, "members_in_background": c + a,
        "query_rate": a / len(query) if query else np.nan,
        "background_rate": (c + a) / len(universe) if universe else np.nan,
        "fold": (a / len(query)) / ((c + a) / len(universe)) if query and universe and (c + a) > 0 else np.nan,
        "odds_ratio": odds, "fisher_p": float(p),
    }


# backgrounds: (a) all atlas genes, (b) CRISPR-tested genes
crispr_tested = set(atlas.loc[atlas.functional_tier.astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
     "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"]), "human_symbol"])

res = []
for bg_name, bg in [("atlas", set(background)), ("crispr_tested", crispr_tested)]:
    res.append(fisher_enrichment(rab_set, panel_set, bg, f"RAB_family | panel | {bg_name}"))
    res.append(fisher_enrichment(traf_set, panel_set, bg, f"trafficking_universe | panel | {bg_name}"))
# also: are panel trafficking genes higher in clean-CRISPR?
panel_traf = sorted(panel_set & traf_set)
panel_non_traf = sorted(panel_set - traf_set)
print(f"panel size: {len(panel)}; trafficking in panel: {len(panel_traf)} ({panel_traf})")

# ---- per-contrast: trafficking subset response vs rest ----
hold["signed_z"] = np.sign(hold["logFC"]) * stats.norm.ppf(np.clip(hold["PValue"] / 2, 1e-300, 0.5 - 1e-300))
multi = ["A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "E_vs_WT", "F_vs_WT",
         "FGHI_vs_WT", "I_vs_WT", "O_vs_WT", "P_vs_WT", "R_vs_WT", "S_vs_WT"]
passing = {"A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "FGHI_vs_WT", "S_vs_WT"}

contrast_rows = []
for c in multi:
    d = hold[hold.contrast_id == c].set_index("human_symbol").reindex(panel)
    traf_z = d.loc[panel_traf, "signed_z"].dropna()
    nont_z = d.loc[panel_non_traf, "signed_z"].dropna()
    if len(traf_z) >= 3 and len(nont_z) >= 3:
        p_w = stats.mannwhitneyu(traf_z.abs(), nont_z.abs(), alternative="greater").pvalue
    else:
        p_w = np.nan
    contrast_rows.append({
        "contrast": c, "passed": c in passing,
        "n_traf": int(traf_z.notna().sum()), "n_nontraf": int(nont_z.notna().sum()),
        "mean_abs_z_traf": round(float(traf_z.abs().mean()), 3),
        "mean_abs_z_nontraf": round(float(nont_z.abs().mean()), 3),
        "median_abs_z_traf": round(float(traf_z.abs().median()), 3),
        "median_abs_z_nontraf": round(float(nont_z.abs().median()), 3),
        "wilcoxon_p_traf_gt": float(p_w) if not np.isnan(p_w) else "",
    })

contrast_df = pd.DataFrame(contrast_rows)
contrast_df.to_csv(OUT / "a4_panel_trafficking_per_contrast.tsv", sep="\t", index=False)

# aggregate: passing vs non-passing
pass_rows = contrast_df[contrast_df.passed]
npass_rows = contrast_df[~contrast_df.passed]
agg = {
    "mean_mean_abs_z_traf_pass": float(pass_rows.mean_abs_z_traf.mean()),
    "mean_mean_abs_z_traf_nonpass": float(npass_rows.mean_abs_z_traf.mean()),
    "mean_mean_abs_z_nontraf_pass": float(pass_rows.mean_abs_z_nontraf.mean()),
    "mean_mean_abs_z_nontraf_nonpass": float(npass_rows.mean_abs_z_nontraf.mean()),
    "n_pass": int(len(pass_rows)), "n_nonpass": int(len(npass_rows)),
}

res_df = pd.DataFrame(res)
res_df.to_csv(OUT / "a4_panel_enrichment.tsv", sep="\t", index=False)

# combined stats json
stats_out = {
    "panel_n": len(panel),
    "rab_in_panel": len(rab_in_panel), "trafficking_in_panel": len(traf_in_panel),
    "enrichment": res, "contrast_summary": agg,
}
(OUT / "a4_stats.json").write_text(json.dumps(stats_out, indent=2), encoding="utf-8")

pd.set_option("display.width", 200)
print("\n=== Panel enrichment ===")
print(res_df[["label", "members_in_query", "query_rate", "background_rate", "fold", "odds_ratio", "fisher_p"]].to_string(index=False))
print("\n=== Per-contrast trafficking response ===")
print(contrast_df.to_string(index=False))
print("\n=== passing vs non-passing aggregate ===")
print(json.dumps(agg, indent=2))
