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

# region -> effector list (confirmed A-I from preprint E-Table1 x FR687201.1 GBK)
regions = {
    "A": ["cegC1"],
    "B": ["legA7"],
    "C": ["lem4", "lem5"],
    "D": ["legL3"],
    "E": ["lem16", "lem17", "ralF", "legC4", "ltpG", "lirA", "lirB", "lirE", "lirF", "pieE"],
    "F": ["legK2", "legAU13", "sdeC", "sde", "sidJ", "sdeB", "sdeA", "lem19", "legS2"],
    "H": ["ceg30", "legA14", "legA15", "sidM", "sidD"],
    "I": ["ceg32", "sdjA", "sdeD", "sdcA", "sidC", "lem26", "lem27"],
}
status = {"A": "PASS", "B": "PASS", "C": "PASS", "D": "PASS", "E": "FAIL", "F": "FAIL",
          "H": "in_FGHI", "I": "FAIL"}

# representative effector columns, grouped roughly by module
cols = ["cegC1", "legA7", "lem4", "lem5", "legL3",
        "ralF", "legC4", "ltpG", "lirA", "lirB", "lirE", "lirF", "pieE", "lem16", "lem17",
        "legK2", "legAU13", "sdeA", "sdeB", "sdeC", "sde", "sidJ", "lem19", "legS2",
        "ceg30", "legA14", "legA15", "sidM", "sidD",
        "ceg32", "sdjA", "sdeD", "sidC", "sdcA", "lem26", "lem27"]

order = ["A", "B", "C", "D", "E", "F", "H", "I"]
M = np.zeros((len(order), len(cols)))
for r_i, mut in enumerate(order):
    for eff in regions[mut]:
        if eff in cols:
            M[r_i, cols.index(eff)] = 1

fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=200)
im = ax.imshow(M, aspect="auto", cmap="Greens", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(np.arange(len(order)))
ax.set_yticklabels([f"{m} ({status[m]})" for m in order], fontsize=8)
ax.set_xticks(np.arange(len(cols)))
ax.set_xticklabels(cols, fontsize=6.5, rotation=55, ha="right")
# module separators
for sep in [5, 15, 24, 29]:
    ax.axvline(sep - 0.5, color="#666666", lw=0.6)
for i, m in enumerate(order):
    s = status[m]
    ax.get_yticklabels()[i].set_color("#2E8B57" if s == "PASS" else ("#C1121F" if s == "FAIL" else "#666666"))
    ax.get_yticklabels()[i].set_fontweight("bold" if s in ("PASS", "FAIL") else "normal")
ax.set_title("Deletion regions vs representative Dot/Icm effectors (A-F, H, I mapped from preprint E-Table1 + 130b annotation; holdout-letter equivalence inferred)\n"
             "green = effector deleted in that region; text color = frozen holdout set-test status (PASS/FAIL); region 'F' deletes the SidE family AND SidJ metaeffector", fontsize=9)
ax.text(-0.5, len(order), "Grey bands separate effector modules:\nRAB/isolated | ARF1/RNA-degrad. | SidE/SidJ | Rab1-AMP | Rab10/SidC", fontsize=6.5, va="bottom")
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"SupplementaryFigure2_effector_deletion_context.{ext}", bbox_inches="tight")
print("wrote", OUT / "SupplementaryFigure2_effector_deletion_context.pdf")
