from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

# ASCII-only text to avoid any Windows/GBK rendering pitfalls
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

df = pd.read_csv(OUT / "a1_evidence_coupling_matched.tsv", sep="\t")
joint = pd.read_csv(OUT / "a1_joint_ppi_lcv_crispr.tsv", sep="\t")

df["log2fold"] = np.log2(df["fold_vs_matched"].replace([np.inf, -np.inf], np.nan))
df["star"] = df["perm_p"].apply(lambda p: "" if p > 0.1 else ("*" if p > 0.01 else ("**" if p > 0.001 else "***")))

short = {
    "PPI(A/B) -> CRISPR hit": "PPI(A/B) -> CRISPR hit",
    "PPI(A/B) -> clean CRISPR": "PPI(A/B) -> clean CRISPR",
    "PPI(any) -> clean CRISPR": "PPI(any) -> clean CRISPR",
    "LCV -> CRISPR hit": "LCV -> CRISPR hit",
    "LCV -> clean CRISPR": "LCV -> clean CRISPR",
    "PPI(A/B) -> LCV": "PPI(A/B) -> LCV",
    "CRISPR hit -> response": "CRISPR hit -> response",
    "clean CRISPR -> response": "clean CRISPR -> response",
    "response -> clean CRISPR": "response -> clean CRISPR",
    "LCV -> response": "LCV -> response",
    "PPI(A/B) -> response": "PPI(A/B) -> response",
}
df["plot_label"] = df["label"].map(short).fillna(df["label"])
order = list(df["plot_label"])  # keep TSV order
df["plot_label"] = pd.Categorical(df["plot_label"], categories=order, ordered=True)
df = df.sort_values("plot_label")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 6.0), dpi=200, gridspec_kw={"width_ratios": [1.45, 1.0]})

# ---- Panel A: matched fold with significance ----
ax = axes[0]
y = np.arange(len(df))
folds = df["fold_vs_matched"].replace(np.inf, np.nan).fillna(0)
colors = ["#2E86AB" if f >= 1.25 else "#C1121F" if f <= 0.8 else "#8D99AE" for f in folds]
ax.barh(y, df["log2fold"].fillna(0), color=colors, edgecolor="none")
ax.set_yticks(y)
ax.set_yticklabels(df["plot_label"], fontsize=8.5)
for i, (_, r) in enumerate(df.iterrows()):
    if r["log2fold"] is not None and np.isfinite(r["log2fold"]):
        ax.text(r["log2fold"] + 0.03, i, r["star"], va="center", fontsize=11)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("log2 fold-enrichment vs confounder-matched background")
ax.set_title("A. Multimodal evidence coupling (matched)", fontsize=10)
ax.set_xlim(-1.2, 3.3)
ax.axvspan(-0.1, 0.1, color="grey", alpha=0.15)

# ---- Panel B: joint increment (y-axis must reach 0.45 > 0.389) ----
ax = axes[1]
j = joint.iloc[0]
labels = ["Universe", "PPI(A/B)", "LCV", "PPI(A/B)+LCV"]
rates = [j["P(clean|universe)"], j["P(clean|PPI)"], j["P(clean|LCV)"], j["P(clean|PPI&LCV)"]]
bars = ax.bar(labels, rates, color=["#8D99AE", "#2E86AB", "#2E86AB", "#1B4965"], edgecolor="black", lw=0.5)
for b, v in zip(bars, rates):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.008, f"{v*100:.1f}%", ha="center", fontsize=9)
ax.set_ylabel("P(clean CRISPR dependency | evidence)")
ax.set_title("B. Joint PPI+LCV -> clean CRISPR\n(n = " + str(int(j["n_ppi_and_lcv"])) + " PPI AND LCV genes)", fontsize=10)
ax.set_ylim(0, 0.45)  # 0.389 must be visible
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig.suptitle(
    "Physical contact, spatial recruitment, functional dependency and transcriptional response\n"
    "are partially overlapping but NON-EQUIVALENT layers of Legionella-host biology",
    fontsize=10, y=1.02)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Figure5_evidence_coupling_FIXED.{ext}", bbox_inches="tight")
# also overwrite the old name for back-compat
fig.savefig(OUT / "a1_evidence_coupling.pdf", bbox_inches="tight")
fig.savefig(OUT / "a1_evidence_coupling.png", bbox_inches="tight")
print("wrote", OUT / "Figure5_evidence_coupling_FIXED.pdf")
