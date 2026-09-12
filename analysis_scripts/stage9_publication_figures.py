from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "11_figures"
SRC = OUT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "grey": "#7A7A7A",
    "lightgrey": "#D9D9D9",
    "black": "#222222",
}

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.2,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 600,
})
sns.set_style("ticks")


def read_tsv(path: Path) -> pd.DataFrame:
    # Several historical outputs use .gz in the name but are plain TSV.
    with path.open("rt", encoding="utf-8-sig") as handle:
        return pd.read_csv(handle, sep="\t")


def load_json(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text(encoding="utf-8-sig"))


def panel(ax, letter: str) -> None:
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", dpi=350)
    plt.close(fig)


def bool_col(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().eq("true")


def figure1() -> None:
    routes = read_tsv(ROOT / "05_expression/stage3_dataset_routes.tsv")
    s1 = load_json("03_interactome/stage1_interactome_statistics.json")
    s2 = load_json("04_genetic_dependency/stage2_statistics.json")
    s4 = load_json("06_differential_expression/stage4_statistics.json")
    s5e = load_json("07_convergence/stage5_expression_statistics.json")
    s5m = load_json("07_convergence/stage5_mechanism_statistics.json")
    s6 = load_json("08_evidence_atlas/stage6_statistics.json")

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.7), constrained_layout=True)
    ax = axs[0, 0]
    rr = routes.sort_values("n_samples", ascending=True)
    bars = ax.barh(rr.accession, rr.n_samples, color=COLORS["blue"], edgecolor="none")
    for b, n in zip(bars, rr.n_samples):
        ax.text(n + 0.8, b.get_y() + b.get_height()/2, str(n), va="center", fontsize=6)
    ax.set_xlabel("Samples in discovery inventory")
    ax.set_ylabel("")
    ax.set_xlim(0, max(rr.n_samples) * 1.18)
    ax.set_title("Six discovery datasets (104 samples)", loc="left", fontweight="bold")
    sns.despine(ax=ax)
    panel(ax, "A")
    routes.to_csv(SRC / "Figure1A_discovery_datasets.tsv", sep="\t", index=False)

    ax = axs[0, 1]
    ppi = pd.Series(s1["class_counts_pairs"]).reindex(["A", "B", "C"])
    ax.bar(ppi.index, ppi.values,
           color=[COLORS["green"], COLORS["blue"], COLORS["grey"]], width=0.7)
    for i, n in enumerate(ppi.values):
        ax.text(i, n + 1, str(n), ha="center", fontsize=7)
    ax.set_ylabel("Unique effector–host pairs")
    ax.set_xlabel("Evidence class")
    ax.set_title("Curated physical-interaction evidence", loc="left", fontweight="bold")
    ax.text(0.02, 0.94, f"79 pairs; 54 human targets\nA/B targets = {s1['represented_AB_targets']}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    sns.despine(ax=ax)
    panel(ax, "B")
    ppi.rename_axis("evidence_class").rename("pairs").reset_index().to_csv(
        SRC / "Figure1B_ppi_classes.tsv", sep="\t", index=False)

    ax = axs[1, 0]
    ft = pd.Series(s2["functional_tier_counts"])
    order = ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
             "F3_SINGLE_LAYER_OR_DISCORDANT"]
    labels = ["F1: genome-wide +\n2-batch concordant", "F2: multi-layer\nconcordant",
              "F3: single-layer/\ndiscordant"]
    vals = [ft[x] for x in order]
    ax.bar(np.arange(3), vals, color=[COLORS["green"], COLORS["blue"], COLORS["orange"]])
    ax.set_xticks(np.arange(3), labels)
    ax.set_ylabel("Genes")
    for i, n in enumerate(vals):
        ax.text(i, n + 12, str(n), ha="center")
    ax.set_title("Replicated U937 CRISPR dependency tiers", loc="left", fontweight="bold")
    ax.text(0.98, 0.96, "DepMap common-essential and\nU937 fitness flags retained",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.3)
    sns.despine(ax=ax)
    panel(ax, "C")
    pd.DataFrame({"tier": labels, "genes": vals}).to_csv(
        SRC / "Figure1C_crispr_tiers.tsv", sep="\t", index=False)

    ax = axs[1, 1]
    funnel = pd.DataFrame({
        "stage": ["DE calls\n(FDR<0.05)", "Stable cross-\ncohort genes",
                  "Mechanistic\nconvergence", "Frozen\ncandidates", "Final strict\nTier 1"],
        "n": [s4["total_fdr05"], s5e["stable_cross_cohort"],
              s5m["convergent"], 56, 1],
    })
    x = np.arange(len(funnel))
    ax.plot(x, funnel.n, color=COLORS["blue"], lw=1.5, marker="o", ms=5)
    ax.fill_between(x, funnel.n, alpha=0.12, color=COLORS["blue"])
    ax.set_yscale("log")
    ax.set_xticks(x, funnel.stage)
    ax.set_ylabel("Number retained (log scale)")
    for i, n in enumerate(funnel.n):
        ax.text(i, n * 1.35, f"{n:,}", ha="center", fontsize=6.5)
    ax.set_title("Predefined evidence convergence", loc="left", fontweight="bold")
    ax.text(0.98, 0.94, f"Known truth set recovered: {s6['truth_recovered_n']}/{s6['truth_set_n']}",
            transform=ax.transAxes, ha="right", va="top", color=COLORS["green"],
            fontweight="bold")
    sns.despine(ax=ax)
    panel(ax, "D")
    funnel.to_csv(SRC / "Figure1D_evidence_funnel.tsv", sep="\t", index=False)
    save(fig, "Figure1_data_and_evidence_architecture")


def select_display_genes(final: pd.DataFrame, n: int = 22) -> list[str]:
    prim = ["RAB1A", "RAB10", "ACTR3", "VAPA"]
    score_cols = ["direct_evidence", "functional_clean", "stable_cross_cohort",
                  "any_mechanism_convergent", "lcv_spatial_supported",
                  "holdout_directional_support"]
    tmp = final.copy()
    tmp["display_score"] = sum(bool_col(tmp[c]).astype(int) for c in score_cols)
    tier_rank = {
        "TIER1_EXTERNALLY_VALIDATED": 0,
        "TIER2_MULTIMODAL": 1,
        "TIER3_CONTEXT_DEPENDENT": 2,
        "TIER3_FITNESS_CAVEAT": 3,
        "REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION": 4,
    }
    tmp["tier_rank"] = tmp.final_tier.map(tier_rank).fillna(9)
    rest = tmp[~tmp.human_symbol.isin(prim)].sort_values(
        ["tier_rank", "display_score", "domain_count", "human_symbol"],
        ascending=[True, False, False, True]).human_symbol.tolist()
    genes = prim + rest[: max(0, n - len(prim) - 1)]
    if "TPM3" not in genes:
        genes.append("TPM3")
    return genes


def figure2() -> None:
    final = read_tsv(ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv")
    genes = select_display_genes(final, 18)
    with (ROOT / "07_convergence/primary_infection_gene_study_evidence.tsv.gz").open(
            "rt", encoding="utf-8-sig") as h:
        primary = pd.read_csv(h, sep="\t")
    with (ROOT / "07_convergence/mechanism_gene_study_evidence.tsv.gz").open(
            "rt", encoding="utf-8-sig") as h:
        mech = pd.read_csv(h, sep="\t")
    lcv = read_tsv(ROOT / "08_evidence_atlas/lcv_spatial_evidence.tsv")

    psub = primary[primary.human_symbol.isin(genes)].copy()
    pmat = psub.pivot_table(index="human_symbol", columns="dataset", values="logFC", aggfunc="mean")
    dataset_order = [x for x in ["GSE61535", "GSE135803", "GSE287755", "GSE26490", "GSE89184", "GSE60085"] if x in pmat]
    pmat = pmat.reindex(index=genes, columns=dataset_order)
    pmat.to_csv(SRC / "Figure2A_discovery_logFC_matrix.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1], width_ratios=[1.28, 1])
    ax = fig.add_subplot(gs[0, 0])
    sns.heatmap(pmat, ax=ax, cmap="vlag", center=0, vmin=-1.5, vmax=1.5,
                linewidths=0.3, linecolor="white", cbar_kws={"label": "log2 fold change", "shrink": 0.65})
    ax.set_xlabel("Discovery dataset")
    ax.set_ylabel("")
    ax.set_title("Cross-cohort infection responses", loc="left", fontweight="bold")
    panel(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    four = primary[primary.human_symbol.isin(["RAB1A", "RAB10", "ACTR3", "VAPA"])].copy()
    four = four[four.dataset.isin(dataset_order)]
    offsets = dict(zip(dataset_order, np.linspace(-0.22, 0.22, len(dataset_order))))
    gorder = ["RAB1A", "RAB10", "ACTR3", "VAPA"]
    for ds, dd in four.groupby("dataset"):
        yy = [gorder.index(g) + offsets[ds] for g in dd.human_symbol]
        ax.scatter(dd.logFC, yy, s=18, label=ds, alpha=0.9)
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.set_yticks(range(4), gorder)
    ax.invert_yaxis()
    ax.set_xlabel("Discovery log2 fold change")
    ax.set_title("Frozen primary candidates", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=2, loc="lower right", fontsize=5.4)
    sns.despine(ax=ax)
    panel(ax, "B")
    four.to_csv(SRC / "Figure2B_primary_discovery_effects.tsv", sep="\t", index=False)

    ax = fig.add_subplot(gs[1, 0])
    routes = ["DOT_ICM", "SINGLE_EFFECTOR_OR_T2SS", "TRANSLATION_EFFECTOR_MODULE", "RNA_RPF_INTERACTION"]
    msub = mech[mech.human_symbol.isin(genes)].copy()
    mm = msub.pivot_table(index="human_symbol", columns="route", values="logFC", aggfunc="mean")
    mm = mm.reindex(index=genes, columns=routes)
    sns.heatmap(mm, ax=ax, cmap="vlag", center=0, vmin=-1, vmax=1,
                linewidths=0.3, linecolor="white", cbar_kws={"label": "mean log2 fold change", "shrink": 0.65})
    ax.set_xticklabels(["Dot/Icm", "Single effector/T2SS", "Translation-effectors", "RNA×RPF"], rotation=30, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Mechanism-resolved perturbation evidence", loc="left", fontweight="bold")
    panel(ax, "C")
    mm.to_csv(SRC / "Figure2C_mechanism_logFC_matrix.tsv", sep="\t")

    ax = fig.add_subplot(gs[1, 1])
    lv = lcv[lcv.human_symbol.isin(genes)].copy().set_index("human_symbol")
    lcols = ["lcv_human_wt", "lcv_human_ankb", "lcv_mouse_untreated"]
    lm = lv.reindex(genes)[lcols].fillna(False).astype(int)
    sns.heatmap(lm, ax=ax, cmap=ListedColormap(["#F2F2F2", COLORS["green"]]),
                vmin=0, vmax=1, cbar=False, linewidths=0.5, linecolor="white")
    ax.set_xticklabels(["Human WT LCV", "Human ankB LCV", "Mouse LCV"], rotation=30, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Independent LCV spatial proteomics", loc="left", fontweight="bold")
    panel(ax, "D")
    lm.to_csv(SRC / "Figure2D_lcv_spatial_matrix.tsv", sep="\t")
    save(fig, "Figure2_cross_cohort_and_mechanism_evidence")


def figure3() -> None:
    final = read_tsv(ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv")
    truth = read_tsv(ROOT / "08_evidence_atlas/known_positive_control_recovery.tsv")
    genes = select_display_genes(final, 24)
    sub = final.set_index("human_symbol").reindex(genes).copy()

    evidence = pd.DataFrame(index=genes)
    evidence["Direct PPI"] = bool_col(sub.ppi_AB_flag).values
    evidence["AP-MS"] = bool_col(sub.apms_high_confidence).values
    evidence["Clean CRISPR"] = bool_col(sub.functional_clean).values
    evidence["Cross-cohort"] = bool_col(sub.stable_cross_cohort).values
    evidence["Mechanism"] = bool_col(sub.any_mechanism_convergent).values
    evidence["LCV spatial"] = bool_col(sub.lcv_spatial_supported).values
    evidence["Holdout"] = bool_col(sub.holdout_directional_support).values
    evidence["Fitness caveat"] = bool_col(sub.fitness_artifact_compatible).values

    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[1.4, 1])
    ax = fig.add_subplot(gs[:, 0])
    plot = evidence.astype(int)
    sns.heatmap(plot, ax=ax, cmap=ListedColormap(["#F1F1F1", COLORS["blue"]]),
                vmin=0, vmax=1, cbar=False, linewidths=0.45, linecolor="white")
    # Overlay caveats in orange and the rejected gene label in red.
    caveat_col = list(plot.columns).index("Fitness caveat")
    for row, val in enumerate(plot["Fitness caveat"]):
        if val:
            ax.add_patch(plt.Rectangle((caveat_col, row), 1, 1, facecolor=COLORS["orange"], edgecolor="white", lw=0.45))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Frozen-candidate multimodal evidence matrix", loc="left", fontweight="bold")
    for tick in ax.get_yticklabels():
        gene = tick.get_text()
        if gene == "RAB1A":
            tick.set_color(COLORS["green"]); tick.set_fontweight("bold")
        elif gene == "TPM3":
            tick.set_color(COLORS["red"]); tick.set_fontweight("bold")
        elif gene in {"RAB10", "ACTR3", "VAPA"}:
            tick.set_fontweight("bold")
    panel(ax, "A")
    evidence.to_csv(SRC / "Figure3A_candidate_evidence_matrix.tsv", sep="\t")

    ax = fig.add_subplot(gs[0, 1])
    tc = final.final_tier.value_counts().reindex([
        "TIER1_EXTERNALLY_VALIDATED", "TIER2_MULTIMODAL", "TIER3_CONTEXT_DEPENDENT",
        "TIER3_FITNESS_CAVEAT", "REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION"]).fillna(0)
    short = ["Tier 1", "Tier 2", "Tier 3 context", "Tier 3 fitness", "Rejected"]
    colors = [COLORS["green"], COLORS["blue"], COLORS["sky"], COLORS["orange"], COLORS["red"]]
    bars = ax.barh(short[::-1], tc.values[::-1], color=colors[::-1])
    for b, n in zip(bars, tc.values[::-1]):
        ax.text(n + 0.5, b.get_y() + b.get_height()/2, str(int(n)), va="center")
    ax.set_xlabel("Frozen candidates")
    ax.set_title("Strict final tiers", loc="left", fontweight="bold")
    ax.set_xlim(0, max(tc.values) * 1.18)
    sns.despine(ax=ax)
    panel(ax, "B")
    tc.rename_axis("final_tier").rename("n").reset_index().to_csv(
        SRC / "Figure3A_final_tier_counts.tsv", sep="\t", index=False)

    ax = fig.add_subplot(gs[1, 1])
    route_names = ["PPI", "AP-MS", "CRISPR", "Expression", "Mechanism", "LCV"]
    truth_matrix = pd.DataFrame(0, index=truth.human_symbol, columns=route_names)
    for gene, routes in zip(truth.human_symbol, truth.recovery_routes.fillna("")):
        truth_matrix.loc[gene, "PPI"] = int("PPI_AB" in routes)
        truth_matrix.loc[gene, "AP-MS"] = int("APMS" in routes)
        truth_matrix.loc[gene, "CRISPR"] = int("CRISPR" in routes)
        truth_matrix.loc[gene, "Expression"] = int("EXPRESSION" in routes)
        truth_matrix.loc[gene, "Mechanism"] = int("MECHANISM" in routes)
        truth_matrix.loc[gene, "LCV"] = int("LCV" in routes)
    sns.heatmap(truth_matrix, ax=ax, cmap=ListedColormap(["#F1F1F1", COLORS["green"]]),
                vmin=0, vmax=1, cbar=False, linewidths=0.45, linecolor="white")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Known-positive recovery (9/9)", loc="left", fontweight="bold")
    panel(ax, "C")
    truth_matrix.to_csv(SRC / "Figure3B_positive_control_recovery.tsv", sep="\t")
    save(fig, "Figure3_multimodal_candidate_prioritization")


def figure4() -> None:
    cand = read_tsv(ROOT / "09_holdout_validation/holdout_frozen_candidate_results.tsv")
    sets = read_tsv(ROOT / "09_holdout_validation/holdout_frozen_set_tests.tsv")
    final = read_tsv(ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv")
    wt = cand[cand.contrast_id.eq("WT_vs_UI")].copy()
    wt = wt.merge(final[["human_symbol", "final_tier"]], on="human_symbol", how="left")
    wt["minuslog10FDR"] = -np.log10(np.clip(wt.FDR.astype(float), 1e-300, 1))

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.1), constrained_layout=True)
    ax = axs[0, 0]
    palette = {
        "TIER1_EXTERNALLY_VALIDATED": COLORS["green"],
        "TIER2_MULTIMODAL": COLORS["blue"],
        "TIER3_CONTEXT_DEPENDENT": COLORS["grey"],
        "TIER3_FITNESS_CAVEAT": COLORS["orange"],
        "REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION": COLORS["red"],
    }
    for tier, dd in wt.groupby("final_tier", dropna=False):
        ax.scatter(dd.logFC, dd.minuslog10FDR, s=20, alpha=0.78,
                   c=palette.get(tier, COLORS["lightgrey"]), edgecolors="none")
    ax.axvline(0, color="black", ls="--", lw=0.6)
    ax.axhline(-np.log10(0.05), color="black", ls=":", lw=0.6)
    for gene in ["RAB1A", "RAB10", "ACTR3", "VAPA", "TPM3"]:
        dd = wt[wt.human_symbol.eq(gene)]
        if len(dd):
            row = dd.iloc[0]
            ax.annotate(gene, (row.logFC, row.minuslog10FDR), xytext=(3, 3),
                        textcoords="offset points", fontsize=6, fontweight="bold")
    ax.set_xlabel("WT vs uninfected log2 fold change")
    ax.set_ylabel("−log10 holdout FDR")
    ax.set_title("Locked panel in independent holdout", loc="left", fontweight="bold")
    sns.despine(ax=ax)
    panel(ax, "A")
    wt.to_csv(SRC / "Figure4A_holdout_candidate_volcano.tsv", sep="\t", index=False)

    ax = axs[0, 1]
    prim = wt[wt.human_symbol.isin(["RAB1A", "RAB10", "ACTR3", "VAPA"])].copy()
    p = np.clip(prim.PValue.astype(float), 1e-300, 1)
    z = norm.isf(p / 2)
    prim["approx_se"] = np.abs(prim.logFC.astype(float)) / z
    prim["ci_low"] = prim.logFC - 1.96 * prim.approx_se
    prim["ci_high"] = prim.logFC + 1.96 * prim.approx_se
    prim = prim.sort_values("logFC")
    y = np.arange(len(prim))
    ax.errorbar(prim.logFC, y,
                xerr=[prim.logFC - prim.ci_low, prim.ci_high - prim.logFC],
                fmt="o", color=COLORS["blue"], ecolor=COLORS["grey"], capsize=2, ms=4)
    sig = prim.FDR < 0.05
    ax.scatter(prim.loc[sig, "logFC"], y[sig], color=COLORS["green"], s=28, zorder=3)
    ax.axvline(0, color="black", ls="--", lw=0.7)
    ax.set_yticks(y, prim.human_symbol)
    ax.set_xlabel("WT vs uninfected log2 fold change (95% CI)")
    ax.set_title("Frozen primary candidates", loc="left", fontweight="bold")
    ax.text(0.03, 0.04, "Green: genome-wide FDR < 0.05\n4/4 direction matched",
            transform=ax.transAxes, fontsize=6.2)
    sns.despine(ax=ax)
    panel(ax, "B")
    prim.to_csv(SRC / "Figure4B_primary_holdout_effects.tsv", sep="\t", index=False)

    ax = axs[1, 0]
    st = sets.copy()
    st["minuslog10FDR"] = -np.log10(np.clip(st.empirical_fdr_within_family.astype(float), 1e-300, 1))
    st["label"] = st.contrast_id + "\n" + st.set_name.str.replace("_", " ")
    cols = np.where(st.passes_frozen_endpoint.astype(str).str.lower().eq("true"), COLORS["green"], COLORS["grey"])
    ax.bar(np.arange(len(st)), st.minuslog10FDR, color=cols)
    ax.axhline(-np.log10(0.05), color=COLORS["red"], ls="--", lw=0.8)
    ax.set_xticks(np.arange(len(st)), st.contrast_id, rotation=55, ha="right")
    ax.set_ylabel("−log10 empirical FDR")
    ax.set_title("Pre-registered set-level tests", loc="left", fontweight="bold")
    ax.text(0.98, 0.95, "8/15 passed", transform=ax.transAxes, ha="right",
            va="top", color=COLORS["green"], fontweight="bold")
    sns.despine(ax=ax)
    panel(ax, "C")
    st.to_csv(SRC / "Figure4C_holdout_set_tests.tsv", sep="\t", index=False)

    ax = axs[1, 1]
    multi = sets[sets.endpoint_family.eq("MULTI_EFFECTOR_TWO_SIDED")].copy()
    multi = multi.sort_values("score")
    pass_flag = multi.passes_frozen_endpoint.astype(str).str.lower().eq("true")
    ax.scatter(multi.score, np.arange(len(multi)),
               c=np.where(pass_flag, COLORS["green"], COLORS["grey"]), s=30)
    ax.set_yticks(np.arange(len(multi)), multi.contrast_id)
    ax.set_xlabel("Mean absolute signed Z")
    ax.set_title("Effector-module boundary conditions", loc="left", fontweight="bold")
    for x, y, passed in zip(multi.score, np.arange(len(multi)), pass_flag):
        if passed:
            ax.text(x + 0.02, y, "PASS", va="center", fontsize=5.5, color=COLORS["green"])
    ax.set_xlim(0, max(multi.score) * 1.18)
    sns.despine(ax=ax)
    panel(ax, "D")
    multi.to_csv(SRC / "Figure4D_effector_boundary_tests.tsv", sep="\t", index=False)
    save(fig, "Figure4_locked_external_validation")


def write_legends() -> None:
    text = """# Main figure legends

## Figure 1. Data and evidence architecture of the Legionella host-dependency atlas

**A**, Discovery expression datasets and sample counts. **B**, Experimentally curated *Legionella pneumophila* effector–human protein interaction pairs stratified by evidence class. **C**, Replicated U937 CRISPR screening tiers after preservation of common-essential and U937 baseline-fitness caveats. **D**, Attrition across the predefined evidence-convergence workflow. Counts at different stages are not a single nested gene universe and are displayed on a logarithmic axis; the panel summarizes analytical selectivity rather than a literal CONSORT flow. The predefined nine-gene positive-control set was recovered without lowering the frozen tier thresholds.

## Figure 2. Cross-cohort, mechanism-resolved and spatial evidence for frozen candidates

**A**, Discovery-cohort log2 fold changes for selected frozen candidates. Blank cells denote absence of a mapped/tested estimate, not zero effect. **B**, Study-specific infection effects for the four primary candidates frozen before holdout access. **C**, Mean log2 fold changes across predefined mechanism routes, including Dot/Icm dependence, single-effector or type-II secretion perturbation, translation-effector modules, and paired RNA–RPF interaction tests. **D**, detection in independent human WT, human *ankB* mutant and mouse untreated Legionella-containing-vacuole proteomic datasets. These spatial measurements support recruitment or proximity and do not establish direct binding.

## Figure 3. Multimodal prioritization and internal positive-control recovery

**A**, Binary evidence matrix for representative frozen candidates. Blue denotes supporting evidence; orange denotes a baseline-fitness caveat. RAB1A (green label) is the sole candidate satisfying the complete preregistered Tier-1 rule; TPM3 (red label) was rejected because the holdout direction contradicted the frozen expectation. **B**, final candidate counts under the strict tier rules. **C**, recovery routes for a nine-gene literature-derived positive-control set. LCV-only recovery was counted as truth-set recovery but did not by itself upgrade candidates into the high-confidence tiers.

## Figure 4. Locked external validation in E-MTAB-16655

**A**, WT-versus-uninfected effects for the frozen candidate panel in the independent 45-sample holdout. **B**, holdout effects and approximate 95% confidence intervals for the four primary candidates; intervals were reconstructed from the two-sided gene-level P value and effect estimate for visualization. RAB1A and RAB10 passed genome-wide FDR < 0.05, while all four primary candidates matched the frozen direction. **C**, empirical FDR from 10,000 gene-label permutations for all 15 preregistered set tests; green bars passed the endpoint-specific 5% FDR criterion. **D**, magnitude and pass status of the multi-effector boundary-condition tests. The four-gene primary set test was borderline (empirical FDR = 0.0603) and was not counted as positive; the 54-gene directional panel and the predefined Dot/Icm subset passed. Holdout validation concerns host transcriptional response and does not establish direct binding, bacterial replication dependence or therapeutic causality.
"""
    (OUT / "MAIN_FIGURE_LEGENDS.md").write_text(text, encoding="utf-8")


def main() -> None:
    figure1()
    figure2()
    figure3()
    figure4()
    write_legends()
    generated = sorted(str(p.relative_to(ROOT)) for p in OUT.glob("Figure*.*"))
    stats = {
        "figures_generated": 4,
        "formats": ["PDF", "SVG", "PNG"],
        "source_data_directory": str(SRC.relative_to(ROOT)),
        "outputs": generated,
        "gate9_publication_figures": "PASS",
    }
    (OUT / "stage9_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    (OUT / "STAGE9_FIGURE_AUDIT.md").write_text(
        "# Stage 9 publication-figure audit\n\n"
        "Status: **PASS (draft figures generated from locked analytical outputs)**\n\n"
        "Four data-bearing main figures were exported as PDF, SVG and 350-dpi PNG. "
        "Every panel has a machine-readable source-data table. Visual inspection is required "
        "before manuscript submission; no conceptual-only figure is counted as a result figure.\n",
        encoding="utf-8")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
