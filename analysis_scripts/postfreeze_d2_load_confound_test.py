from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

h = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz", sep="\t", low_memory=False)
panel = sorted(set(pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
                   .loc[lambda d: d.contrast_id == "WT_vs_UI", "human_symbol"]))

h["signed_z"] = np.sign(h["logFC"]) * stats.norm.ppf(np.clip(h["PValue"] / 2, 1e-300, 0.5 - 1e-300))

multi = ["A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "E_vs_WT", "F_vs_WT",
         "FGHI_vs_WT", "I_vs_WT", "O_vs_WT", "P_vs_WT", "R_vs_WT", "S_vs_WT"]
passing = {"A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "FGHI_vs_WT", "S_vs_WT"}

rows = []
for c in multi + ["WT_vs_UI", "WT_vs_dotA"]:
    d = h[h.contrast_id == c]
    d = d[d.FDR.notna() & d.logFC.notna()]
    n_deg = int((d.FDR < 0.05).sum())
    n_deg1 = int(((d.FDR < 0.05) & (d.logFC.abs() >= 1)).sum())
    global_mean_abs_z = float(d.signed_z.abs().mean())
    global_median_abs_z = float(d.signed_z.abs().median())
    global_mean_abs_lfc = float(d.logFC.abs().mean())
    # panel subset
    dp = d[d.human_symbol.isin(panel)]
    panel_mean_abs_z = float(dp.signed_z.abs().mean())
    panel_median_abs_z = float(dp.signed_z.abs().median())
    # random gene-set expectation ~ global mean |z|; panel enrichment = panel / global
    ratio = panel_mean_abs_z / global_mean_abs_z if global_mean_abs_z else np.nan
    # one-sample: is panel mean|z| > global mean|z|?
    p_panel_gt = stats.mannwhitneyu(dp.signed_z.abs(), d[~d.human_symbol.isin(panel)].signed_z.abs(), alternative="greater").pvalue
    rows.append({
        "contrast": c,
        "passes_frozen": c in passing,
        "n_genes": len(d), "n_deg_fdr05": n_deg, "n_deg_fdr05_abs1": n_deg1,
        "global_mean_abs_z": round(global_mean_abs_z, 4), "global_median_abs_z": round(global_median_abs_z, 4),
        "global_mean_abs_logFC": round(global_mean_abs_lfc, 4),
        "panel_mean_abs_z": round(panel_mean_abs_z, 4), "panel_median_abs_z": round(panel_median_abs_z, 4),
        "panel_over_global_ratio": round(ratio, 3),
        "mw_p_panel_gt_nonpanel": round(p_panel_gt, 5),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT / "d2_load_confound_covariate_test.tsv", sep="\t", index=False)

# ---- tests ----
multi_df = df[df.contrast.isin(multi)].copy()
multi_df["passed"] = multi_df.passes_frozen
# (1) Does panel response track global perturbation across contrasts?
r_panel_global, p_panel_global = stats.spearmanr(multi_df.panel_mean_abs_z, multi_df.global_mean_abs_z)
r_pass_global, p_pass_global = stats.spearmanr(multi_df.n_deg_fdr05, multi_df.global_mean_abs_z)
# (2) Do PASS groups have higher global perturbation than FAIL?
pass_g = multi_df[multi_df.passed].global_mean_abs_z
fail_g = multi_df[~multi_df.passed].global_mean_abs_z
p_global_pass_vs_fail = stats.mannwhitneyu(pass_g, fail_g, alternative="greater").pvalue
# (3) Is panel SPECIFICALLY enriched (panel>global) in each contrast?
n_enriched = int((multi_df.panel_over_global_ratio > 1.05).sum())
n_sig_specific = int((multi_df.mw_p_panel_gt_nonpanel < 0.05).sum())

summary = {
    "n_contrasts_multi": len(multi_df),
    "spearman_panel_vs_global": {"rho": round(r_panel_global, 3), "p": round(p_panel_global, 4)},
    "mw_global_load_pass_vs_fail_p": round(p_global_pass_vs_fail, 4),
    "mean_global_abs_z_PASS": round(float(pass_g.mean()), 4),
    "mean_global_abs_z_FAIL": round(float(fail_g.mean()), 4),
    "n_contrasts_panel_enriched_over_global": n_enriched,
    "n_contrasts_panel_signif_over_nonpanel": n_sig_specific,
    "interpretation": "If the frozen-panel response were purely a bacterial-load artifact, (1) panel response would track global perturbation AND (2) panel would NOT be specifically enriched vs random gene sets (which share the same load).",
}
(OUT / "d2_load_confound_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print(df.to_string(index=False))
print("\n=== summary ===")
print(json.dumps(summary, indent=2))
