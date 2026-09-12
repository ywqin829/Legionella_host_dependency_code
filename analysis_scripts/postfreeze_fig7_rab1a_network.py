from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
tier = pd.read_csv(ROOT / "10_final_tiers" / "final_candidate_tier_matrix.tsv", sep="\t")


def B(s):
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


atlas["_ppi"] = B(atlas["ppi_AB_flag"])
atlas["_apms"] = B(atlas["apms_high_confidence"])
atlas["_crispr"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])

# FINAL tier (only RAB1A is Tier 1)
final_map = dict(zip(tier.human_symbol, tier.final_tier))
atlas["final_tier"] = atlas["human_symbol"].map(final_map).fillna("background")
atlas = atlas.drop_duplicates("human_symbol").set_index("human_symbol")


def tier_color(t):
    t = str(t)
    if t.startswith("TIER1"): return "#C1121F"
    if t.startswith("TIER2"): return "#E07A1F"
    if t.startswith("TIER3"): return "#F4D03F"
    if t == "REJECT": return "#333333"
    return "#B0BEC5"


# ---- Panel A: direct effector -> RAB1A (literature-confirmed edges only) ----
effectors = [
    ("DrrA/SidM", "GEF + AMPylator", "GEF PMID:16906144; AMPylation PMID:20651120"),
    ("SidD", "de-AMPylase", "PMID:21680813;21734656"),
    ("LepB", "Rab1 GAP (inactivator)", "PMID:17952054"),
    ("AnkX", "phosphocholination", "PMID:21822290"),
    ("Lem3", "de-phosphocholinase", "PMID:22158903"),
    ("SetA", "Rab1 glycosylase", "PMID:30323948"),
]

fig, axes = plt.subplots(2, 1, figsize=(7.2, 9.0), dpi=200, gridspec_kw={"height_ratios": [1.0, 1.35]})
fig.suptitle("RAB1A: the sole Tier-1 host dependency combines direct multi-effector manipulation\nwith a host trafficking program carrying concordant genetic/spatial/response evidence",
             fontsize=10.5, y=0.99)

# Panel A
ax = axes[0]; ax.axis("off"); ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.35, 1.05)
ne = len(effectors)
for i, (nm, role, ev) in enumerate(effectors):
    x = 0.02; y = 0.95 - 0.17 * i
    ax.scatter(x, y, s=1500, color="#D35400", edgecolor="black", linewidth=0.6, zorder=3)
    ax.text(x, y, nm, ha="center", va="center", fontsize=8.5, color="white", zorder=4, fontweight="bold")
    ax.plot([x + 0.045, 0.72], [y, 0.32], color="#D35400", lw=1.3, alpha=0.85, zorder=1)
    # role + evidence label on the edge
    ax.text(0.40, (y + 0.32) / 2 + 0.02, role + "\n" + ev, ha="center", va="bottom", fontsize=6.8, color="#7B3F00", zorder=2)
ax.scatter([0.78], [0.32], s=2600, color="#C1121F", edgecolor="black", linewidth=1.5, zorder=5)
ax.text(0.78, 0.32, "RAB1A\nTIER 1", ha="center", va="center", fontsize=9, fontweight="bold", color="white", zorder=6)
ax.set_title("A. Direct effector manipulation of RAB1A\n(literature-confirmed; edges = biochemistry/PPI with evidence source)", fontsize=8.5, pad=6)

# Panel B: host trafficking evidence matrix (5-slot bars) with FINAL tier
ax = axes[1]
genes = ["RAB1A", "RAB10", "RABIF", "SEC22B", "VAPA", "ARF1", "RAB8A", "RAB8B", "RAB7A",
         "RAB1B", "RAB5C", "RAB35", "RAB13", "RAB22A", "SAR1A", "COPA", "STX5", "TRAPPC3"]
tracks = ["_ppi", "_apms", "_crispr", "_lcv", "_resp"]
track_lab = ["PPI(A/B)", "AP-MS", "CRISPR", "LCV", "Response"]
rows = []
for g in genes:
    if g not in atlas.index:
        continue
    a = atlas.loc[g]
    values = [int(a[t]) for t in tracks]
    rows.append({"gene": g, "tier": a["final_tier"],
                 **dict(zip(track_lab, values)), "n": sum(values)})
host = pd.DataFrame(rows)
tier_rank = {"TIER1_EXTERNALLY_VALIDATED": 0, "TIER2_MULTIMODAL": 1, "TIER3_CONTEXT_DEPENDENT": 2}
host["tier_rank"] = host["tier"].map(tier_rank).fillna(9)
host = host.sort_values(["tier_rank", "n"], ascending=[True, False])

ys = np.arange(len(host))[::-1]
for i, (_, r) in enumerate(host.iterrows()):
    y = ys[i]
    ax.add_patch(Rectangle((-0.02, y - 0.34), 1.04, 0.68, facecolor=tier_color(r["tier"]), edgecolor="black", lw=0.5, zorder=1))
    ax.text(-0.10, y, r["gene"], ha="right", va="center", fontsize=8.5,
            color="#C1121F" if r["gene"] == "RAB1A" else "black",
            fontweight="bold" if r["gene"] == "RAB1A" else "normal", zorder=3)
    for j, col in enumerate(track_lab):
        v = int(r[col])
        ax.add_patch(Rectangle((j, y - 0.30), 0.88, 0.60, facecolor="#2E8B57" if v else "#EEEEEE",
                               edgecolor="#BBBBBB", lw=0.4, zorder=2))
        if v:
            ax.text(j + 0.44, y, "", ha="center", va="center", fontsize=6)
ax.set_xticks(np.arange(len(track_lab)))
ax.set_xticklabels(track_lab, fontsize=8)
ax.set_ylim(-1, len(host))
ax.set_xlim(-0.35, len(tracks) + 0.1)
ax.set_title("B. Host trafficking program: evidence domains with FINAL tier\n(green = evidence present; left bar = final Tier: red=Tier1, orange=Tier2, yellow=Tier3)", fontsize=8.5, pad=6)
for xt, tl in zip(range(len(track_lab)), ax.get_xticklabels()):
    tl.set_ha("center")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
legend = [Patch(facecolor="#C1121F", label="Tier 1 (RAB1A)"),
          Patch(facecolor="#E07A1F", label="Tier 2 multimodal"),
          Patch(facecolor="#F4D03F", label="Tier 3 context-dependent"),
          Patch(facecolor="#2E8B57", label="evidence present")]
ax.legend(handles=legend, fontsize=7, loc="lower left", bbox_to_anchor=(0, -0.18), ncol=2, frameon=False)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Figure7_RAB1A_evidence_architecture_FIXED.{ext}", bbox_inches="tight")
host.to_csv(OUT / "fig7_host_trafficking_evidence_finaltier.tsv", sep="\t", index=False)
pd.DataFrame(effectors, columns=["effector", "role", "evidence_source"]).to_csv(
    OUT / "fig7_effector_rab1a_edges.tsv", sep="\t", index=False)
print("wrote", OUT / "Figure7_RAB1A_evidence_architecture_FIXED.pdf")
print(host[["gene", "tier", "n"]].to_string(index=False))
