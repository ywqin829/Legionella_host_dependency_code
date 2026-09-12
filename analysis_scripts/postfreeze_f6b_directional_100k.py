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
frozen = pd.read_csv(ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv", sep="\t")
exp_dir = dict(zip(frozen.human_symbol, frozen.expected_wt_vs_uninfected))

h = pd.read_csv(ROOT / "09_holdout_validation" / "holdout_all_contrasts.tsv.gz", sep="\t", low_memory=False)
# frozen set-test convention: signed_z = sign(logFC)*norm.isf(P/2); positive = concordant with UP
h["signed_z"] = np.sign(h["logFC"]) * stats.norm.isf(np.clip(h["PValue"] / 2, 1e-300, 0.5))
d_ui = h[h.contrast_id == "WT_vs_UI"].set_index("human_symbol")["signed_z"]

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

# direction pool from panel's frozen composition
panel_dirs = np.array([exp_dir.get(g) for g in panel_in if exp_dir.get(g) in ("UP", "DOWN")])

# vectorized matched sampling (with replacement from pools; large pools => valid approximation)
N = 100000
pools = []
for i in q_idx:
    pool = np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.05])
    if len(pool) < 10: pool = np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.15])
    pools.append(pool)
pool_lens = np.array([len(p) for p in pools])
# precompute z of all genes for fast lookup
z_arr = np.array([d_ui.get(g, np.nan) for g in cov.index])
dir_signs = np.array([1 if d == "UP" else -1 for d in panel_dirs])

# per-set: sample one gene from each panel gene's pool, then sample a direction
# vectorized: gene_samples shape (n_panel_genes, N)
gene_samples = np.empty((len(q_idx), N), dtype=int)
for k, pool in enumerate(pools):
    gene_samples[k] = pool[RNG.integers(0, len(pool), size=N)]
# for each sampled gene, get its z (0 if missing) and assign a sampled direction
samp_z = z_arr[gene_samples]
samp_z = np.nan_to_num(samp_z, nan=0.0)
dir_idx = RNG.integers(0, len(dir_signs), size=(len(q_idx), N))
samp_sign = dir_signs[dir_idx]
stats_null = (samp_sign * samp_z).mean(axis=0)

# observed directional statistic (frozen convention)
obs = float(np.mean([(1 if exp_dir.get(g) == "UP" else -1) * d_ui[g] for g in panel_in if g in d_ui.index]))
p_emp = float((stats_null >= obs).mean() + 1e-6)
print(f"WT_vs_UI directional: obs={obs:.3f}; matched-null mean={stats_null.mean():.3f}+-{stats_null.std():.3f}; empirical P={p_emp:.6f} (N={N})")
print(f"  -> exact lower bound on empirical P = 1/(N+1) = {1/(N+1):.2e}")

out = {"n_permutations": N, "observed_directional_stat": round(obs, 4),
       "matched_null_mean": round(float(stats_null.mean()), 4), "matched_null_sd": round(float(stats_null.std()), 4),
       "empirical_p": p_emp, "min_achievable_p": 1 / (N + 1),
       "n_obs_ge_null": int((stats_null >= obs).sum())}
(OUT / "f6b_directional_100k.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
np.save(OUT / "f6b_directional_null_100k.npy", stats_null)
print(json.dumps(out, indent=2))
