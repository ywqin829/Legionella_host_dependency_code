from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

df = pd.read_csv(OUT / "a3_trafficking_gene_evidence_matrix.tsv", sep="\t")
fam = pd.read_csv(OUT / "a3_family_evidence_saturation.tsv", sep="\t", index_col=0)
enr = pd.read_csv(OUT / "a4_panel_enrichment.tsv", sep="\t")
ct = pd.read_csv(OUT / "a4_panel_trafficking_per_contrast.tsv", sep="\t")

tracks = ["ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]
track_lab = ["PPI(A/B)", "AP-MS", "CRISPR", "Clean\nCRISPR", "Cross-cohort\nresponse", "LCV"]
fam_tracks = ["ppi_AB", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]
fam_lab = ["PPI(A/B)", "CRISPR", "Clean\nCRISPR", "Cross-cohort\nresp.", "LCV"]

# final-size figure (~7.2 in wide, double column): fonts must stay >= 6pt
fig, axes = plt.subplots(2, 2, figsize=(7.4, 7.6), dpi=200, gridspec_kw={"hspace": 0.55, "wspace": 0.45})

# ---- Panel A: RAB genes with evidence (readable subset), final-size ----
rab = df[df.family.eq("RAB")].copy()
rab["n"] = rab[tracks].sum(axis=1)
rab_main = rab[(rab.n >= 3) | (rab.in_frozen_panel.astype(str).str.upper().isin(["TRUE", "1"]))].copy()
rab_main = rab_main.sort_values(["n", "gene"], ascending=[False, True])
ax = axes[0, 0]
M = rab_main[tracks].astype(int).values
ax.imshow(M, aspect="auto", cmap="RdPu", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(rab_main)))
ax.set_yticklabels(rab_main.gene, fontsize=7)
ax.set_xticks(range(len(tracks)))
ax.set_xticklabels(track_lab, fontsize=6)
for i, g in enumerate(rab_main.gene):
    if g in ("RAB1A", "RAB10", "RAB8A", "RAB8B", "RAB7A", "RAB1B", "RAB5C", "RAB35", "RAB13", "RAB22A", "RAB6A"):
        ax.get_yticklabels()[i].set_color("#C1121F"); ax.get_yticklabels()[i].set_fontweight("bold")
for xt, tl in zip(range(len(tracks)), ax.get_xticklabels()):
    tl.set_rotation(30); tl.set_ha("right")
ax.set_title("A. RAB family evidence atlas\n(genes with >=2 evidence domains or in frozen panel; full 61-gene map in Supp. Fig. 5)", fontsize=8)

# ---- Panel B: family saturation ----
ax = axes[0, 1]
fam_t = fam[fam_tracks].values
ax.imshow(fam_t, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(fam))); ax.set_yticklabels(fam.index, fontsize=6)
ax.set_xticks(range(len(fam_tracks))); ax.set_xticklabels(fam_lab, fontsize=6)
for i in range(fam_t.shape[0]):
    for j in range(fam_t.shape[1]):
        v = fam_t[i, j]
        if v > 0:
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.5, color="white" if v > 0.55 else "black")
for xt, tl in zip(range(len(fam_tracks)), ax.get_xticklabels()):
    tl.set_rotation(30); tl.set_ha("right")
ax.set_title("B. Trafficking-family evidence saturation\n(fraction of family genes positive)", fontsize=8)

# ---- Panel C: panel enrichment ----
ax = axes[1, 0]
sub = enr[enr["label"].str.contains("atlas$")]
bars = ax.bar(["RAB\nfamily", "Trafficking\nuniverse"], sub["fold"].values, color=["#C1121F", "#2E86AB"], edgecolor="black", lw=0.4)
for b, f, p in zip(bars, sub["fold"].values, sub["fisher_p"].values):
    star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
    ax.text(b.get_x() + b.get_width() / 2, f * 1.02, f"{f:.1f}x{star}", ha="center", fontsize=7)
ax.set_ylabel("Fold enrichment in 54-gene panel\nvs genome background", fontsize=7)
ax.set_title("C. Frozen panel is trafficking-centric\n(10 RAB / 12 trafficking of 54)", fontsize=8)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_ylim(0, 78)
ax.tick_params(axis="x", labelsize=6.5); ax.tick_params(axis="y", labelsize=6)

# ---- Panel D: per-contrast points ----
ax = axes[1, 1]
for _, r in ct.iterrows():
    c = "#2E86AB" if r.passed else "#8D99AE"
    ax.plot([0, 1], [r.mean_abs_z_traf, r.mean_abs_z_nontraf], "-", color=c, lw=0.7, alpha=0.5, zorder=1)
ax.scatter([0]*len(ct), ct.mean_abs_z_traf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=18, zorder=2, edgecolor="black", lw=0.3)
ax.scatter([1]*len(ct), ct.mean_abs_z_nontraf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=18, zorder=2, edgecolor="black", lw=0.3)
ax.plot([0, 1], [ct[ct.passed].mean_abs_z_traf.mean(), ct[ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#C1121F", lw=1.8, ms=5, label="Passing (6) mean")
ax.plot([0, 1], [ct[~ct.passed].mean_abs_z_traf.mean(), ct[~ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#222222", lw=1.8, ms=5, label="Non-passing (6) mean")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Trafficking\nsubset", "Non-traff.\nsubset"], fontsize=6.5)
ax.set_ylabel("Mean |signed z| of panel genes", fontsize=7)
ax.set_title("D. Response by trafficking membership\nper deletion contrast", fontsize=8)
ax.legend(fontsize=6, frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", labelsize=6)

fig.suptitle("The frozen host-dependency panel is enriched for RAB / membrane-trafficking machinery\nacross diverse effector-deletion backgrounds", fontsize=9.5, y=0.99)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Figure6_trafficking_atlas_PUBREADY.{ext}", bbox_inches="tight")
print("wrote", OUT / "Figure6_trafficking_atlas_PUBREADY.pdf")
print(f"Panel A shows {len(rab_main)} RAB genes")
