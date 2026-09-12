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

tracks = ["ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]
track_lab = ["Effector int.\n(A/B)", "AP-MS\nhigh-conf", "CRISPR\nhit", "Clean\nCRISPR", "Cross-cohort\nresponse", "LCV\nspatial"]
fam_tracks = ["ppi_AB", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]
fam_lab = ["Effector int.\n(A/B)", "CRISPR\nhit", "Clean\nCRISPR", "Cross-cohort\nresponse", "LCV\nspatial"]

fig = plt.figure(figsize=(13.5, 11.5), dpi=200)
gs = fig.add_gridspec(2, 3, height_ratios=[1.55, 0.8], width_ratios=[1.5, 1, 1], hspace=0.42, wspace=0.35)

# ---- Panel A: RAB-family genes (all 61, tall & readable) ----
rab = df[df.family.eq("RAB")].copy()
rab["n"] = rab[tracks].sum(axis=1)
rab = rab.sort_values(["n", "gene"], ascending=[True, False])
M = rab[tracks].astype(int).values
ax = fig.add_subplot(gs[0, 0])
im = ax.imshow(M, aspect="auto", cmap="RdPu", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(rab)))
ax.set_yticklabels(rab.gene, fontsize=6.5)
ax.set_xticks(range(len(tracks)))
ax.set_xticklabels(track_lab, fontsize=6.5)
for i, g in enumerate(rab.gene):
    if g in ("RAB1A", "RAB10", "RAB8A", "RAB8B", "RAB7A", "RAB33B", "RAB35", "RAB5C", "RAB1B", "RAB22A", "RAB13", "RAB6A"):
        ax.get_yticklabels()[i].set_color("#C1121F"); ax.get_yticklabels()[i].set_fontweight("bold")
for xt, tl in zip(range(len(tracks)), ax.get_xticklabels()):
    tl.set_rotation(35); tl.set_ha("right")
ax.set_title("A. RAB GTPase family evidence atlas\n(61 genes; red = key Legionella trafficking candidates)", fontsize=9)
cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_ticks([0, 1]); cbar.set_ticklabels(["absent", "present"])

# ---- Panel B: family saturation ----
ax = fig.add_subplot(gs[0, 1])
fam_t = fam[fam_tracks].values
im2 = ax.imshow(fam_t, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(fam)))
ax.set_yticklabels(fam.index, fontsize=7)
ax.set_xticks(range(len(fam_tracks)))
ax.set_xticklabels(fam_lab, fontsize=6)
for i in range(fam_t.shape[0]):
    for j in range(fam_t.shape[1]):
        v = fam_t[i, j]
        if v > 0:
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if v > 0.55 else "black")
for xt, tl in zip(range(len(fam_tracks)), ax.get_xticklabels()):
    tl.set_rotation(35); tl.set_ha("right")
ax.set_title("B. Trafficking-family evidence saturation\n(fraction of family genes positive)", fontsize=9)

# ---- Panel C: A4 panel enrichment (bars) ----
ax = fig.add_subplot(gs[0, 2])
sub = enr[enr["label"].str.contains("atlas$")]
labels = ["RAB GTPase family", "Full trafficking universe"]
folds = sub["fold"].values; ps = sub["fisher_p"].values
bars = ax.bar(labels, folds, color=["#C1121F", "#2E86AB"], edgecolor="black", lw=0.5)
for b, f, p in zip(bars, folds, ps):
    star = "" if p > 0.05 else ("*" if p > 0.01 else "***")
    ax.text(b.get_x() + b.get_width() / 2, f * 1.02, f"{f:.1f}x{star}", ha="center", fontsize=9)
ax.set_ylabel("Fold enrichment in 54-gene frozen panel\nvs genome background", fontsize=8)
ax.set_title("C. Frozen host-response panel is trafficking-centric\n(10 RAB / 12 trafficking genes of 54)", fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_ylim(0, 78)

# ---- Panel D: A4 per-contrast individual points ----
ct = pd.read_csv(OUT / "a4_panel_trafficking_per_contrast.tsv", sep="\t")
ax = fig.add_subplot(gs[1, :])
for _, r in ct.iterrows():
    c = "#2E86AB" if r.passed else "#8D99AE"
    ax.plot([0, 1], [r.mean_abs_z_traf, r.mean_abs_z_nontraf], "-", color=c, lw=0.8, alpha=0.55, zorder=1)
ax.scatter([0]*len(ct), ct.mean_abs_z_traf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=42, zorder=2, edgecolor="black", lw=0.3)
ax.scatter([1]*len(ct), ct.mean_abs_z_nontraf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=42, zorder=2, edgecolor="black", lw=0.3)
# group means
ax.plot([0, 1], [ct[ct.passed].mean_abs_z_traf.mean(), ct[ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#C1121F", lw=2.2, ms=7, label="Passing contrasts (A,B,C,D,FGHI,S) mean")
ax.plot([0, 1], [ct[~ct.passed].mean_abs_z_traf.mean(), ct[~ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#333333", lw=2.2, ms=7, label="Non-passing (E,F,I,O,P,R) mean")
ax.set_xticks([0, 1])
ax.set_xticklabels(["Trafficking subset\n(n=12 panel genes)", "Non-trafficking subset\n(n=42 panel genes)"], fontsize=8)
ax.set_ylabel("Mean |signed z| of frozen panel genes", fontsize=9)
ax.set_title("D. Per-effector-deletion-contrast response by trafficking membership\n(each line = one of 12 deletion contrasts; points = individual contrasts)", fontsize=9)
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.suptitle(
    "The frozen host-dependency panel is highly enriched for RAB / membrane-trafficking machinery,\n"
    "supporting a conserved host-trafficking core across diverse effector-deletion backgrounds",
    fontsize=11, y=0.995)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Figure6_trafficking_atlas_FIXED.{ext}", bbox_inches="tight")
fig.savefig(OUT / "a3_trafficking_atlas.pdf", bbox_inches="tight")
fig.savefig(OUT / "a3_trafficking_atlas.png", bbox_inches="tight")
print("wrote", OUT / "Figure6_trafficking_atlas_FIXED.pdf")
