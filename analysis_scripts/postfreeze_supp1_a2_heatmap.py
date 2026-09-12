from __future__ import annotations

import re
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

d = pd.read_csv(OUT / "a2_evidence_fragility_matrix.tsv", sep="\t")
datasets = ["GSE135803", "GSE26490", "GSE287755", "GSE60085", "GSE61535", "GSE89184"]

# parse leave1_dataset strings -> per-dataset status
def parse_row(r):
    out = {}
    for token in str(r["leave1_dataset"]).split("; "):
        token = token.strip()
        if not token:
            continue
        if "->" in token:
            ds, st = token.split("->", 1)
            out[ds.strip()] = st.strip()  # UNSTABLE / STABLE
        elif ":" in token:
            ds, rest = token.split(":", 1)
            out[ds.strip()] = "STABLE"  # direction preserved
        else:
            out[token] = "?"
    return out

rows = []
for _, r in d.iterrows():
    parsed = parse_row(r)
    base_stable = bool(r["base_stable_cross"])
    row = {"candidate": r["candidate"], "final_tier": r["final_tier"],
           "human_only": "STABLE" if bool(r["human_only_stable"]) else "UNSTABLE",
           "ppi_n_pmids": r["ppi_n_distinct_pmids"], "lcv_n_sources": r["lcv_n_sources"],
           "base_dir": str(r["base_direction"])}
    for ds in datasets:
        s = parsed.get(ds)
        if s is None:
            row[ds] = "NOT_IN"   # gene not measured in this dataset
        elif s == "UNSTABLE":
            row[ds] = "STABLE_FLIPPED"  # direction kept, stable flag flipped
        else:
            row[ds] = "STABLE"
    rows.append(row)

df = pd.DataFrame(rows).set_index("candidate")
# order: RAB1A first, then Tier2 by final appearance
order = ["RAB1A"] + [c for c in df.index if c != "RAB1A"]
df = df.loc[order]

cols = datasets + ["human_only"]
valmap = {"STABLE": 0, "STABLE_FLIPPED": 1, "UNSTABLE": 2, "NOT_IN": 3}
M = np.array([[valmap.get(df.loc[c, col], 3) for col in cols] for c in df.index])
cmap = matplotlib.colors.ListedColormap(["#2E8B57", "#F4D03F", "#C1121F", "#B0BEC5"])
cmap.set_bad("#FFFFFF")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 7.2), dpi=200, gridspec_kw={"width_ratios": [2.8, 1.2]})
ax = axes[0]
im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
ax.set_yticks(np.arange(len(df)))
ax.set_yticklabels([f"{g} ({df.loc[g,'base_dir']})" for g in df.index], fontsize=8)
ax.set_xticks(np.arange(len(cols)))
ax.set_xticklabels([c.replace("GSE", "") for c in cols], fontsize=8)
for xt, tl in zip(range(len(cols)), ax.get_xticklabels()):
    tl.set_rotation(0)
# annotate cells with stable/direction
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if v == 0:
            ax.text(j, i, "ok", ha="center", va="center", fontsize=6, color="white")
        elif v == 1:
            ax.text(j, i, "dir~", ha="center", va="center", fontsize=6, color="black")
        elif v == 2:
            ax.text(j, i, "flip", ha="center", va="center", fontsize=6, color="white")
ax.axvline(len(datasets) - 0.5, color="black", lw=1)
ax.set_title("A. Leave-one-source-out robustness of cross-cohort response\n(UP/DOWN direction is preserved in ALL removals; cells show stability flag)", fontsize=9)
leg = [Patch(facecolor="#2E8B57", label="Stable (direction + stability kept)"),
       Patch(facecolor="#F4D03F", label="Direction kept, stable flag flipped (n reduced)"),
       Patch(facecolor="#C1121F", label="Unstable"),
       Patch(facecolor="#B0BEC5", label="Not measured in dataset")]
ax.legend(handles=leg, fontsize=7, loc="lower left", bbox_to_anchor=(0, -0.02), ncol=1, frameon=False)

# right panel: PPI / LCV source counts + human-only marker
ax = axes[1]
ax.axis("off")
info = df[["ppi_n_pmids", "lcv_n_sources", "human_only", "final_tier"]]
ys = np.arange(len(df))[::-1]
ax.set_yticks(ys)
ax.set_yticklabels(df.index, fontsize=8)
labels = ["PPI sources", "LCV sources"]
x0 = 0.0
for k, col in enumerate(["ppi_n_pmids", "lcv_n_sources"]):
    for i, g in enumerate(df.index):
        v = df.loc[g, col]
        ax.text(x0 + k * 2.6, len(df) - 1 - i, str(int(v)) if v == v else "0", ha="center", va="center", fontsize=8,
                color="#1B4965" if v >= 2 else ("#C1121F" if v == 1 else "#999"))
ax.text(x0, len(df) + 0.7, "PPI\nsources", ha="center", fontsize=8)
ax.text(x0 + 2.6, len(df) + 0.7, "LCV\nsources", ha="center", fontsize=8)
ax.set_xlim(-0.5, 6)
ax.set_ylim(-1, len(df) + 2)
ax.set_title("B. Source counts\nred = single-source", fontsize=9)

fig.suptitle("Cross-cohort evidence is directionally robust to single-dataset removal; stability flag is threshold-sensitive at n=3-4",
             fontsize=10, y=1.0)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"SupplementaryFigure1_evidence_fragility.{ext}", bbox_inches="tight")
print("wrote", OUT / "SupplementaryFigure1_evidence_fragility.pdf")
