from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
RNG = np.random.default_rng(42)

panel = sorted(set(pd.read_csv(ROOT / "09_holdout_validation" / "holdout_frozen_candidate_results.tsv", sep="\t")
                   .loc[lambda d: d.contrast_id == "WT_vs_UI", "human_symbol"]))
# frozen expected directions for the directional panel test
frozen = pd.read_csv(ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv", sep="\t")
exp_dir = dict(zip(frozen.human_symbol, frozen.expected_wt_vs_uninfected))

h = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz", sep="\t", low_memory=False)
# positive signed_z for significant UP genes (matches the frozen set-test convention)
h["signed_z"] = np.sign(h["logFC"]) * stats.norm.ppf(np.clip(h["PValue"] / 2, 1e-300, 0.5 - 1e-300)) * -1.0

# matching covariates from discovery
ev = pd.read_csv(ROOT / "07_convergence" / "primary_infection_gene_study_evidence.tsv.gz",
                 sep="\t", compression=None, low_memory=False)
evz = ev.assign(az=ev.signed_z.abs())
cov = evz.groupby("human_symbol")["az"].agg(["mean", "std", "count"]).rename(columns={"mean": "resp_mean", "std": "resp_sd", "count": "n_studies"})
cov["resp_mean"] = cov.resp_mean.fillna(0); cov["resp_sd"] = cov.resp_sd.fillna(cov.resp_sd.median())

def load_ui(mat, gsm, gmap, probe):
    m = pd.read_csv(ROOT / "05_expression" / "matrices" / mat, sep="\t", compression="gzip")
    idc = "feature_id" if "feature_id" in m.columns else "Geneid"
    cols = [c for c in gsm if c in m.columns]
    gm = pd.read_csv(ROOT / "05_expression" / "gene_maps" / gmap, sep="\t")
    gm = gm[gm.gene_symbol.notna()]
    if probe and "ambiguous_mapping" in gm.columns:
        gm = gm[~gm.ambiguous_mapping.astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in gm.columns and idc != "feature_id":
        gm = gm.rename(columns={"feature_id": idc})
    mg = m.merge(gm, on=idc, how="inner")
    sub = mg[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    return sub.rank(pct=True)[cols].mean(axis=1)
ui = pd.concat([
    load_ui("GSE61535_log2_expression.tsv.gz", ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    load_ui("GSE135803_fpkm.tsv.gz", ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    load_ui("GSE287755_raw_counts.tsv.gz", ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
], axis=1).mean(axis=1, skipna=True)
cov["base_expr"] = cov.index.map(ui).fillna(0)

uni = sorted(set(h.human_symbol) & set(cov.index))
cov = cov.loc[uni].copy()
panel_in = [g for g in panel if g in cov.index]
FEAT = ["resp_mean", "resp_sd", "base_expr", "n_studies"]
Z = (cov[FEAT] - cov[FEAT].mean()) / cov[FEAT].std().replace(0, 1)
flag = cov.index.isin(panel_in).astype(int)
clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
clf.fit(Z.values, flag); prop = clf.predict_proba(Z.values)[:, 1]
q_idx = np.where(flag == 1)[0]; ctl_idx = np.where(flag == 0)[0]
pools = {i: np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.05]) for i in q_idx}
for i, p in pools.items():
    if len(p) < 10: pools[i] = np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.15])
N = 100000
null_sets = []
for _ in range(N):
    used = set(); chosen = []
    for i in q_idx:
        avail = [j for j in pools[i] if j not in used]
        if not avail: break
        chosen.append(int(RNG.choice(avail))); used.add(chosen[-1])
    null_sets.append(chosen)
print(f"matched null sets: {N}; universe {len(uni)}; panel genes {len(panel_in)}")

d_ui = h[h.contrast_id == "WT_vs_UI"].set_index("human_symbol")["signed_z"]

# panel's frozen direction composition (for assigning directions to null genes)
panel_dirs = [exp_dir.get(g) for g in panel_in if exp_dir.get(g) in ("UP", "DOWN")]
panel_dir_pool = np.array(panel_dirs)

def panel_dir_stat(d_ui, genes, use_panel_dirs=True):
    zs = []
    if use_panel_dirs:
        # null genes: sample directions from the panel's frozen composition
        dirs = RNG.choice(panel_dir_pool, size=len(genes), replace=True)
        for g, d in zip(genes, dirs):
            if g in d_ui.index:
                zs.append((1 if d == "UP" else -1) * d_ui[g])
    else:
        for g in genes:
            if g in d_ui.index and g in exp_dir:
                d = exp_dir[g]
                if d in ("UP", "DOWN"):
                    zs.append((1 if d == "UP" else -1) * d_ui[g])
    return float(np.mean(zs)) if zs else np.nan

obs_dir = panel_dir_stat(d_ui, panel_in, use_panel_dirs=False)
null_dir = np.array([panel_dir_stat(d_ui, [cov.index[i] for i in s]) for s in null_sets])
null_dir = null_dir[~np.isnan(null_dir)]
p_dir = float((null_dir >= obs_dir).mean() + 1e-6)
print(f"WT_vs_UI directional: obs={obs_dir:.3f}, matched-null mean={null_dir.mean():.3f}+-{null_dir.std():.3f}, empirical P={p_dir:.4f}")

# multi-effector two-sided tests with matched null
multi = ["A_vs_WT", "B_vs_WT", "C_vs_WT", "D_vs_WT", "E_vs_WT", "F_vs_WT",
         "FGHI_vs_WT", "I_vs_WT", "O_vs_WT", "P_vs_WT", "R_vs_WT", "S_vs_WT"]
rows = []
for c in multi:
    dc = h[h.contrast_id == c].set_index("human_symbol")["signed_z"]
    obs = float(dc.loc[panel_in].abs().mean())
    null = []
    for s in null_sets:
        gz = dc.loc[[cov.index[i] for i in s if cov.index[i] in dc.index]].dropna()
        null.append(float(gz.abs().mean()))
    null = np.array(null)
    rows.append({"contrast": c, "panel_mean_abs_z": round(obs, 3),
                 "matched_null_mean": round(float(null.mean()), 3), "matched_null_sd": round(float(null.std()), 3),
                 "matched_empirical_p": round(float((null >= obs).mean() + 1e-6), 5)})
res = pd.DataFrame(rows)
res["matched_fdr"] = pd.Series(stats.false_discovery_control(pd.to_numeric(res.matched_empirical_p), method="bh"), index=res.index)
res["pass_after_matched"] = res.matched_fdr < 0.05
res.to_csv(OUT / "f6_holdout_matched_null.tsv", sep="\t", index=False)
summary = {"WT_vs_UI_directional_obs": round(obs_dir, 3), "WT_vs_UI_matched_null_p": round(p_dir, 5),
           "n_multieffector_pass_after_matched_fdr": int(res.pass_after_matched.sum())}
(OUT / "f6_holdout_matched_null_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(res.to_string(index=False))
print("\n" + json.dumps(summary, indent=2))
