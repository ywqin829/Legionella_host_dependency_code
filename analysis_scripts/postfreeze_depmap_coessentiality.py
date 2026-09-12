from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
dep = ROOT / "01_data_raw" / "depmap" / "CRISPRGeneEffect.csv"
common_file = ROOT / "01_data_raw" / "depmap" / "CRISPRInferredCommonEssentials.csv"

MACHINERY = ["RAB1A", "RAB1B", "RAB10", "RABIF", "SEC22B", "SAR1A", "SEC23A", "SEC24C", "SEC31A",
             "COPA", "COPB1", "COPE", "TRAPPC3", "TRAPPC10", "STX5", "GOSR1", "GOLGA2", "VAPA", "RAB8A", "RAB8B"]
MIN_PERIODS = 50
N_PERM = 10000
N_PERM_SP = 1000
MATCH_POOL = 200
RNG = np.random.default_rng(42)

# ---- load full gene-effect matrix (gene x cell line) ----
t0 = time.time()
df = pd.read_csv(dep)
df = df.rename(columns={df.columns[0]: "cell_line"}).set_index("cell_line")
df = df.T  # gene x cell line
df.index.name = "gene"
# gene name cleanup: "SYMBOL (entrez)" -> "SYMBOL"
df.index = [str(x).split(" (")[0] if " (" in str(x) else str(x) for x in df.index]
print(f"loaded {df.shape[0]} genes x {df.shape[1]} cell lines in {time.time()-t0:.0f}s")

# keep genes with enough measurements
valid = df.notna().sum(axis=1)
gene_ok = valid[valid >= MIN_PERIODS]
print(f"genes with >= {MIN_PERIODS} measurements: {len(gene_ok)}")

# Matching variables for the permutation null. Matching avoids comparing the
# module with genes having very different baseline dependency, variance,
# measurement completeness, or common-essential status.
common = set()
if common_file.exists():
    ce = pd.read_csv(common_file)
    common = {
        str(x).split(" (")[0] if " (" in str(x) else str(x)
        for x in ce.iloc[:, 0].dropna()
    }
feat = pd.DataFrame(index=gene_ok.index)
feat["mean_effect"] = df.loc[gene_ok.index].mean(axis=1, skipna=True)
feat["sd_effect"] = df.loc[gene_ok.index].std(axis=1, skipna=True)
feat["missing_fraction"] = 1.0 - valid.loc[gene_ok.index] / df.shape[1]
feat["common_essential"] = feat.index.isin(common)
match_cols = ["mean_effect", "sd_effect", "missing_fraction"]
mu = feat[match_cols].mean()
sd = feat[match_cols].std().replace(0, 1)
feat_z = (feat[match_cols] - mu) / sd


def zscore_rows(M):
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M, axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (M - mu) / sd


def mean_pairwise_r(gene_vecs, method="pearson"):
    """gene_vecs: list of gene names. returns mean+median pairwise correlation."""
    G = df.loc[list(gene_vecs)]
    G = G.astype(float)
    if method == "pearson":
        Zn = zscore_rows(G.values)
        with np.errstate(invalid="ignore"):
            Zc = np.where(np.isnan(Zn), 0.0, Zn)
            C = (Zc @ Zc.T) / np.sum(~np.isnan(Zn), axis=1)[:, None]
            C = (C + C.T) / 2
        iu = np.triu_indices(len(G), k=1)
        vals = C[iu]
    else:  # spearman (slow; use on smaller sets)
        vals = []
        genes = list(gene_vecs)
        for a in range(len(genes)):
            for b in range(a + 1, len(genes)):
                x = G.loc[genes[a]].values; y = G.loc[genes[b]].values
                m = ~(np.isnan(x) | np.isnan(y))
                if m.sum() >= MIN_PERIODS:
                    vals.append(spearmanr(x[m], y[m])[0])
        vals = np.array(vals)
    vals = vals[~np.isnan(vals)]
    return float(vals.mean()), float(np.median(vals)), vals


mach = [g for g in MACHINERY if g in gene_ok.index]
if len(mach) < 5:
    raise SystemExit("too few machinery genes measured")
obs_mean, obs_med, _ = mean_pairwise_r(mach, "pearson")
obs_mean_sp, obs_med_sp, _ = mean_pairwise_r(mach, "spearman")

# ---- permutation null: covariate-matched 20-gene sets ----
# Each machinery gene receives a pool of nearest non-machinery genes with the
# same common-essential status. One unique match per machinery gene is drawn.
universe = [g for g in gene_ok.index if g not in set(mach)]
n = len(mach)
match_pools = {}
for g in mach:
    same_ce = feat.index[feat["common_essential"] == feat.loc[g, "common_essential"]]
    candidates = [x for x in same_ce if x not in set(mach)]
    delta = feat_z.loc[candidates, match_cols] - feat_z.loc[g, match_cols]
    dist = (delta ** 2).sum(axis=1).sort_values()
    match_pools[g] = dist.head(min(MATCH_POOL, len(dist))).index.to_numpy()
    if len(match_pools[g]) < n:
        raise SystemExit(f"too few matched controls for {g}")


def sample_matched_set():
    used = set()
    chosen = []
    for g in RNG.permutation(mach):
        available = [x for x in match_pools[g] if x not in used]
        if not available:
            raise RuntimeError(f"no unique matched control remains for {g}")
        x = RNG.choice(available)
        used.add(x)
        chosen.append(x)
    return chosen


perm_sets = [sample_matched_set() for _ in range(N_PERM)]
perm_means = np.empty(N_PERM)
t0 = time.time()
for k, samp in enumerate(perm_sets):
    m, _, _ = mean_pairwise_r(samp, "pearson")
    perm_means[k] = m
    if (k + 1) % 500 == 0:
        print(f"  perm {k+1}/{N_PERM} ({time.time()-t0:.0f}s)")
emp_p = float(((perm_means >= obs_mean).sum() + 1) / (N_PERM + 1))

# Spearman sensitivity analysis uses 1,000 identically matched sets.
perm_means_sp = np.empty(N_PERM_SP)
for k, samp in enumerate(perm_sets[:N_PERM_SP]):
    m, _, _ = mean_pairwise_r(samp, "spearman")
    perm_means_sp[k] = m
emp_p_sp = float(((perm_means_sp >= obs_mean_sp).sum() + 1) / (N_PERM_SP + 1))

print("\n=== Observed vs permutation null (matched 20-gene sets) ===")
print(f"Pearson: observed mean r = {obs_mean:.4f} (median {obs_med:.4f}); null mean = {perm_means.mean():.4f} +- {perm_means.std():.4f}; empirical P = {emp_p:.5f}")
print(f"Spearman: observed mean r = {obs_mean_sp:.4f} (median {obs_med_sp:.4f}); null mean = {perm_means_sp.mean():.4f}; empirical P = {emp_p_sp:.5f}")
print(f"null percentiles: 50%={np.percentile(perm_means,50):.4f}, 95%={np.percentile(perm_means,95):.4f}, 99%={np.percentile(perm_means,99):.4f}")

stats_out = {
    "n_cell_lines": int(df.shape[1]),
    "n_machinery_genes": len(mach), "machinery_genes": mach,
    "n_permutations_pearson": N_PERM, "n_permutations_spearman": N_PERM_SP,
    "min_periods": MIN_PERIODS, "match_pool_per_query_gene": MATCH_POOL,
    "matching_variables": ["mean_effect", "sd_effect", "missing_fraction", "common_essential_status"],
    "observed_mean_r_pearson": obs_mean, "observed_median_r_pearson": obs_med,
    "observed_mean_r_spearman": obs_mean_sp, "observed_median_r_spearman": obs_med_sp,
    "null_mean_r_pearson": float(perm_means.mean()), "null_sd_r_pearson": float(perm_means.std()),
    "null_p95_pearson": float(np.percentile(perm_means, 95)),
    "empirical_p_pearson": emp_p, "empirical_p_spearman": emp_p_sp,
    "interpretation": "Empirical permutation P from covariate-matched 20-gene sets with the plus-one correction; exploratory only",
    "caveat": "DepMap co-essentiality is exploratory background only; does NOT remove common-essential caveat or establish infection specificity",
}
(OUT / "depmap_coessentiality_stats.json").write_text(json.dumps(stats_out, indent=2), encoding="utf-8")
np.save(OUT / "depmap_null_distribution_pearson.npy", perm_means)
np.save(OUT / "depmap_null_distribution_spearman.npy", perm_means_sp)
print("\n" + json.dumps(stats_out, indent=2))
