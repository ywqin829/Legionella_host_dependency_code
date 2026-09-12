from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import ALL_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

dep = ROOT / "01_data_raw" / "depmap" / "CRISPRGeneEffect.csv"
hdr = pd.read_csv(dep, nrows=0)
sym_map = {c.split(" (")[0]: c for c in hdr.columns[1:] if " (" in c}
machinery = ["RAB1A", "RAB1B", "RAB10", "RABIF", "SEC22B", "SAR1A", "SEC23A", "SEC24C", "SEC31A",
             "COPA", "COPB1", "COPE", "TRAPPC3", "TRAPPC10", "STX5", "GOSR1", "GOLGA2", "VAPA", "RAB8A", "RAB8B"]
mach = [g for g in machinery if g in sym_map]
df = pd.read_csv(dep, usecols=["Unnamed: 0"] + [sym_map[g] for g in mach])
df = df.rename(columns={"Unnamed: 0": "cl"}).set_index("cl"); df.columns = mach
C = df.corr(method="pearson", min_periods=50)
# cluster by simple order (families grouped)
order = [g for g in machinery if g in mach]
M = C.loc[order, order].values

fig, ax = plt.subplots(figsize=(7.2, 6.8), dpi=200)
im = ax.imshow(M, vmin=-0.15, vmax=0.35, cmap="RdBu_r", interpolation="nearest")
ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
ax.set_xticklabels(order, rotation=55, ha="right", fontsize=7)
ax.set_yticklabels(order, fontsize=7)
ax.set_title("DepMap pan-cancer co-essentiality of Rab / ER-Golgi trafficking machinery\n(mean r=0.045 vs matched null 0.0016; Pearson empirical P=0.00010; exploratory, does not remove common-essential caveat)", fontsize=9)
cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("Pearson r (Chronos gene effect)")
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"SupplementaryFigure4_depmap_coessentiality.{ext}", bbox_inches="tight")
print("wrote", OUT / "SupplementaryFigure4_depmap_coessentiality.pdf")
