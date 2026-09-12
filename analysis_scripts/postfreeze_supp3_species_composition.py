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

b2 = pd.read_csv(OUT / "b2_species_dependency_scores.tsv", sep="\t")
# Use the same frozen leading-candidate subset as the fragility analysis
# (RAB1A + 16 leading Tier-2 candidates), rather than plotting all 56 rows.
frag = pd.read_csv(OUT / "a2_evidence_fragility_matrix.tsv", sep="\t")
available = set(b2["human_symbol"])
order = [g for g in frag["candidate"] if g in available]
b2 = b2.set_index("human_symbol").loc[order].copy()
b2["lcv_h"] = b2["lcv_state"].str.startswith("HUMAN").astype(int)
b2["lcv_ankb"] = b2["lcv_state"].str.contains("ANKB").astype(int)
b2["lcv_m"] = b2["lcv_state"].str.contains("MOUSE").astype(int)

fig, axes = plt.subplots(1, 2, figsize=(7.2, 7.5), dpi=200, gridspec_kw={"width_ratios": [1.45, 1]})
fig.suptitle("Species composition of evidence for leading candidates\n(expression: human vs mouse study count; LCV: human-WT / human-ankB / mouse source presence)", fontsize=10, y=0.99)

# ---- Panel A: stacked expression studies ----
ax = axes[0]
y = np.arange(len(b2))
nh = b2["n_human"].values.astype(int)
nm = b2["n_mouse"].values.astype(int)
ax.barh(y, nm, color="#C9A227", edgecolor="black", lw=0.4, label="Mouse expression studies")
ax.barh(y, nh, left=nm, color="#1B4965", edgecolor="black", lw=0.4, label="Human expression studies")
ax.set_yticks(y)
ax.set_yticklabels(b2.index, fontsize=8)
for i, (h, m) in enumerate(zip(nh, nm)):
    ax.text(h + m + 0.05, i, f"{h}+{m}", va="center", fontsize=7)
ax.set_xlabel("n expression studies")
ax.set_title("A. Expression evidence by species\n(human+3 cross-cohort stable = "+str(int((b2.expression_stable_cross_cohort.sum()))) + " candidates)", fontsize=9)
ax.set_xlim(0, 7)
ax.legend(fontsize=8, frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ---- Panel B: LCV source type presence ----
ax = axes[1]
L = b2[["lcv_h", "lcv_ankb", "lcv_m"]].values
ax.imshow(L, aspect="auto", cmap="Oranges", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(b2)))
ax.set_yticklabels(b2.index, fontsize=8)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["Human WT LCV", "Human ankB LCV", "Mouse LCV"], fontsize=7)
# annotate PPI source count
for i in range(len(b2)):
    v = b2.loc[b2.index[i], "ppi_n_pmids"]
    v_int = 0 if pd.isna(v) else int(v)
    ax.text(3.2, i, f"PPI:{v_int}src", ha="left", va="center", fontsize=7,
            color="#C1121F" if v_int <= 1 else "#1B4965")
ax.set_xlim(-0.5, 5.2)
ax.set_title("B. LCV spatial source composition + PPI source count\n(red PPI text = single-source)", fontsize=9)
ax.set_xticks([0, 1, 2], ["Human WT LCV", "Human ankB LCV", "Mouse LCV"], fontsize=7)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"SupplementaryFigure3_species_composition.{ext}", bbox_inches="tight")
print("wrote", OUT / "SupplementaryFigure3_species_composition.pdf")
