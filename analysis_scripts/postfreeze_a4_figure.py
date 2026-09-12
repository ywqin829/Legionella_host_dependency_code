from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

enr = pd.read_csv(OUT / "a4_panel_enrichment.tsv", sep="\t")
ct = pd.read_csv(OUT / "a4_panel_trafficking_per_contrast.tsv", sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=200)

# ---- Panel A: enrichment ----
ax = axes[0]
sub = enr[enr["label"].str.contains("atlas$")]
labels = ["RAB GTPase family", "Full trafficking universe"]
folds = sub["fold"].values; ps = sub["fisher_p"].values
bars = ax.bar(labels, folds, color=["#C1121F", "#2E86AB"], edgecolor="black", lw=0.5)
for b, f, p in zip(bars, folds, ps):
    star = "" if p > 0.05 else ("*" if p > 0.01 else "***")
    ax.text(b.get_x() + b.get_width() / 2, f * 1.02, f"{f:.1f}x{star}", ha="center", fontsize=10)
ax.set_ylabel("Fold enrichment in 54-gene frozen panel vs genome background")
ax.set_title("A. Frozen host-response panel is trafficking-centric\n(10 RAB / 12 trafficking genes of 54)", fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_ylim(0, 78)

# ---- Panel B: per-contrast individual points ----
ax = axes[1]
for _, r in ct.iterrows():
    c = "#2E86AB" if r.passed else "#8D99AE"
    ax.plot([0, 1], [r.mean_abs_z_traf, r.mean_abs_z_nontraf], "-", color=c, lw=0.8, alpha=0.5, zorder=1)
ax.scatter([0]*len(ct), ct.mean_abs_z_traf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=40, zorder=2, edgecolor="black", lw=0.3)
ax.scatter([1]*len(ct), ct.mean_abs_z_nontraf, c=["#2E86AB" if p else "#8D99AE" for p in ct.passed], s=40, zorder=2, edgecolor="black", lw=0.3)
ax.plot([0, 1], [ct[ct.passed].mean_abs_z_traf.mean(), ct[ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#C1121F", lw=2.4, ms=8, label="Passing (6) mean")
ax.plot([0, 1], [ct[~ct.passed].mean_abs_z_traf.mean(), ct[~ct.passed].mean_abs_z_nontraf.mean()], "o-", color="#222222", lw=2.4, ms=8, label="Non-passing (6) mean")
ax.set_xticks([0, 1])
ax.set_xticklabels(["Trafficking subset\n(n=12)", "Non-trafficking subset\n(n=42)"], fontsize=8)
ax.set_ylabel("Mean |signed z| of panel genes")
ax.set_title("B. Per-contrast response by trafficking membership\n(12 individual deletion contrasts shown as lines/points)", fontsize=9)
ax.legend(fontsize=8, frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.suptitle("The frozen host-dependency panel is highly enriched for RAB / membrane-trafficking machinery", fontsize=10, y=1.02)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Figure6C_D_panel_enrichment_FIXED.{ext}", bbox_inches="tight")
fig.savefig(OUT / "a4_panel_enrichment.pdf", bbox_inches="tight")
fig.savefig(OUT / "a4_panel_enrichment.png", bbox_inches="tight")
print("wrote", OUT / "Figure6C_D_panel_enrichment_FIXED.pdf")
