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

df = pd.read_csv(OUT / "a3_trafficking_gene_evidence_matrix.tsv", sep="\t")
tracks = ["ppi_AB", "apms_high", "crispr_hit", "crispr_clean", "stable_cross_cohort", "lcv_spatial"]
track_lab = ["PPI(A/B)", "AP-MS", "CRISPR", "Clean\nCRISPR", "Cross-cohort\nresponse", "LCV"]

rab = df[df.family.eq("RAB")].copy()
rab["n"] = rab[tracks].sum(axis=1)
rab = rab.sort_values(["n", "gene"], ascending=[True, False])

fig, ax = plt.subplots(figsize=(7.4, 10.0), dpi=200)
M = rab[tracks].astype(int).values
ax.imshow(M, aspect="auto", cmap="RdPu", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(rab)))
ax.set_yticklabels(rab.gene, fontsize=6.5)
ax.set_xticks(range(len(tracks)))
ax.set_xticklabels(track_lab, fontsize=7)
for i, g in enumerate(rab.gene):
    if g in ("RAB1A", "RAB10", "RAB8A", "RAB8B", "RAB7A", "RAB33B", "RAB35", "RAB5C", "RAB1B", "RAB22A", "RAB13", "RAB6A"):
        ax.get_yticklabels()[i].set_color("#C1121F"); ax.get_yticklabels()[i].set_fontweight("bold")
for xt, tl in zip(range(len(tracks)), ax.get_xticklabels()):
    tl.set_rotation(30); tl.set_ha("right")
ax.set_title("Full RAB GTPase family evidence atlas (n=" + str(len(rab)) + ")\nred = key Legionella trafficking candidates", fontsize=9)
cbar = fig.colorbar(ax.images[0], ax=ax, shrink=0.5, pad=0.02)
cbar.set_ticks([0, 1]); cbar.set_ticklabels(["absent", "present"])
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"SupplementaryFigure5_RAB_full_atlas.{ext}", bbox_inches="tight")
print("wrote", OUT / "SupplementaryFigure5_RAB_full_atlas.pdf")
