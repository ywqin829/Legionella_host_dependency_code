"""P0-B: is the mouse BMM LCV -> CRISPR coupling explained by MS detectability?

Motivation
----------
The LCV-associated proteome is an MS detection: a host protein is called
"LCV-associated" only if it is enriched enough to be identified. Detection
probability rises with protein abundance, and abundant proteins can differ in
CRISPR essentiality, so a detection-based exposure is in principle confounded by
detectability. The mouse BMM LCV proteome (PMID 26829557) is quantitative
(Average LFQ intensity and sequence coverage per protein), so the confound can be
tested there directly. The human U937 LCV proteomes (PMID 25369898) are
presence-only and carry no abundance or coverage, so they cannot be adjusted and
this limitation is disclosed rather than modelled.

This script does NOT alter candidate selection, tiers, frozen endpoints or the
published a1/h1 coupling values. It writes new outputs only.

Analyses
--------
(1) Detectability gradient within the detected proteome: among mouse-detected LCV
    proteins in the coupling universe, test whether CRISPR-hit status is
    associated with log10 LFQ intensity or sequence coverage, adjusted for the
    audited covariates. A null result shows the mouse LCV signal is not an
    MS-abundance gradient inside the detected set.

(2) High-confidence-detection-restricted coupling: recompute the mouse LCV ->
    fitness-filtered CRISPR matched enrichment with the exposure restricted to
    unambiguously identified mouse proteins (sequence coverage >= median among
    detections), against the same not-detected controls, using the a1 matched
    design. Persistence shows the coupling is not carried by marginal detections.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
SRC = ROOT / "08_evidence_atlas" / "lcv_spatial_source_records.tsv"
MOUSE_XLSX = ROOT / "01_data_raw" / "lcv_proteomics" / "PLOS_S2_LCV_untreated.xlsx"
RNG = np.random.default_rng(42)
N_PERM = 4000


def B(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


# ---------------- Atlas + covariates (identical to postfreeze_a1_coupling.py) ----
atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz",
                    sep="\t", low_memory=False)
deg = (pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv",
                   sep="\t", low_memory=False)
       [["gene", "cross_pathogen_degree", "cross_pathogen_records", "cross_pathogen_pmids"]]
       .drop_duplicates())

atlas["_ppi_AB"] = B(atlas["ppi_AB_flag"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_lcv_mouse"] = B(atlas["lcv_mouse_untreated"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_fitness_ok"] = B(atlas["fitness_artifact_compatible"])
atlas["_crispr_tested"] = atlas["functional_tier"].astype(str).isin(
    ["F1_GW_AND_BOTH_BATCH_CONCORDANT", "F2_MULTILAYER_CONCORDANT",
     "F3_SINGLE_LAYER_OR_DISCORDANT", "NO_FUNCTIONAL_HIT"])
atlas["_crispr_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_crispr_clean"] = atlas["_crispr_hit"] & (~atlas["_fitness_ok"])
atlas["_common_ess"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas["_nstudies"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)


def load_ui_percentile(matrix_name, ui_gsms, map_name, use_probe_map):
    mat = pd.read_csv(ROOT / "05_expression" / "matrices" / matrix_name, sep="\t", compression="gzip")
    cols = [c for c in ui_gsms if c in mat.columns]
    if not cols:
        return None
    id_col = "feature_id" if "feature_id" in mat.columns else "Geneid"
    gm = pd.read_csv(ROOT / "05_expression" / "gene_maps" / map_name, sep="\t")
    gm = gm[gm["gene_symbol"].notna()]
    if use_probe_map and "ambiguous_mapping" in gm.columns:
        gm = gm[~gm["ambiguous_mapping"].astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in gm.columns and id_col != "feature_id":
        gm = gm.rename(columns={"feature_id": id_col})
    m = mat.merge(gm, on=id_col, how="inner")
    sub = m[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    return sub.rank(pct=True)[cols].mean(axis=1)


ui_percentiles = {
    "GSE61535": load_ui_percentile("GSE61535_log2_expression.tsv.gz",
                                   ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    "GSE135803": load_ui_percentile("GSE135803_fpkm.tsv.gz",
                                    ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    "GSE287755": load_ui_percentile("GSE287755_raw_counts.tsv.gz",
                                    ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
}
ui_stack = pd.concat(ui_percentiles.values(), axis=1, keys=list(ui_percentiles.keys()))
atlas["_ui_expr"] = atlas["human_symbol"].map(ui_stack.mean(axis=1, skipna=True))
atlas = atlas.merge(deg, left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)
atlas["_ui_expr"] = atlas["_ui_expr"].fillna(0.0)

# ---------------- Mouse protein detectability (LFQ + sequence coverage) --------
src = pd.read_csv(SRC, sep="\t", low_memory=False)
mouse = src[src["host_species"] == "MOUSE"].copy()
cov_x = pd.read_excel(MOUSE_XLSX, sheet_name="identified proteins (host)")
cov_x = cov_x[["Uniprot identifier", "Sequence coverage [%]"]].rename(columns={"Uniprot identifier": "accession"})
cov_x["accession"] = cov_x["accession"].astype(str).str.strip()
mouse["accession"] = mouse["accession"].astype(str).str.strip()
mouse = mouse.merge(cov_x, on="accession", how="left")
mouse["_lfq"] = pd.to_numeric(mouse["quantitative_value"], errors="coerce")
mouse["_cov"] = pd.to_numeric(mouse["Sequence coverage [%]"], errors="coerce")
prot = mouse.groupby("human_symbol", as_index=False).agg(_lfq=("_lfq", "max"), _cov=("_cov", "max"))
prot["_log_lfq"] = np.log10(prot["_lfq"].clip(lower=1.0))

atlas = atlas.groupby("human_symbol", as_index=False, sort=True).agg({
    "_ppi_AB": "max", "_lcv": "max", "_lcv_mouse": "max", "_resp": "max", "_fitness_ok": "max",
    "_crispr_tested": "max", "_crispr_hit": "max", "_crispr_clean": "max",
    "_common_ess": "max", "_chronos": "median", "_nstudies": "max",
    "_ui_expr": "median", "_deg": "max", "_pmids": "max",
})
atlas = atlas.merge(prot, on="human_symbol", how="left")
CONF = ["_deg", "_pmids", "_chronos", "_common_ess", "_ui_expr"]


def universe():
    d = atlas.copy()
    d = d[d["_crispr_tested"]]
    for c in CONF:
        d = d[d[c].notna()]
    return d.reset_index(drop=True)


def matched_fold(d, exposure_col, outcome_col="_crispr_clean"):
    """a1 confounder-matched enrichment for a given exposure column."""
    d = d.reset_index(drop=True)
    X = np.column_stack([
        np.log1p(d["_deg"].clip(lower=0).replace(0, 1e-9)),
        np.log1p(d["_pmids"].clip(lower=0).replace(0, 1e-9)),
        d["_chronos"].astype(float),
        d["_common_ess"].astype(float),
        d["_ui_expr"].astype(float),
    ])
    y = d[exposure_col].astype(int).values
    if y.sum() < 5:
        return None
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(X, y)
    prop = clf.predict_proba(X)[:, 1]
    q_idx = np.where(y == 1)[0]
    ctl_idx = np.where(y == 0)[0]
    if len(ctl_idx) == 0:
        return None
    caliper = 0.25 * float(np.std(prop))
    used = np.zeros(len(ctl_idx), dtype=bool)
    positions = np.arange(len(ctl_idx))
    pairs = []  # (query_idx, control_idx), 1:1 matched
    for i in q_idx[np.argsort(prop[q_idx])]:
        avail = positions[~used]
        if avail.size == 0:
            break
        diff = np.abs(prop[ctl_idx[avail]] - prop[i])
        b = int(np.argmin(diff))
        if diff[b] <= caliper:
            used[avail[b]] = True
            pairs.append((i, ctl_idx[avail[b]]))
    if len(pairs) < 5:
        return None
    q_pairs = np.array([p[0] for p in pairs])
    m_pairs = np.array([p[1] for p in pairs])
    o = d[outcome_col].astype(int).values
    o_q, o_m = int(o[q_pairs].sum()), int(o[m_pairs].sum())
    rate_q, rate_m = o_q / len(q_pairs), o_m / len(m_pairs)
    _, p_fisher = stats.fisher_exact([[o_q, len(q_pairs) - o_q], [o_m, len(m_pairs) - o_m]], alternative="greater")
    # matched-pair permutation: labels flipped within each 1:1 pair
    obs = rate_q - rate_m
    perm = []
    q_vals, m_vals = o[q_pairs], o[m_pairs]
    for _ in range(N_PERM):
        flip = RNG.random(len(q_pairs)) < 0.5
        a = np.where(flip, m_vals, q_vals)
        b = np.where(flip, q_vals, m_vals)
        perm.append(a.mean() - b.mean())
    perm_p = float((1 + np.sum(np.array(perm) >= obs)) / (1 + len(perm)))
    return {"fold": rate_q / rate_m if rate_m > 0 else np.nan, "fisher_p": p_fisher,
            "perm_p": perm_p, "n_query": len(q_pairs), "n_matched": len(m_pairs),
            "rate_query": rate_q, "rate_matched": rate_m}


def gradient_within_detected(d):
    """Logistic: CRISPR outcome ~ detectability, among mouse-detected proteins."""
    sub = d[d["_lcv_mouse"] & d["_cov"].notna() & d["_log_lfq"].notna()].copy()
    y = sub["_crispr_clean"].astype(int).values
    cols = ["_log_lfq", "_cov"] + CONF
    X = sm.add_constant(sub[cols].astype(float).values)
    names = ["const", "log10_LFQ", "seq_coverage_pct"] + CONF
    res = sm.Logit(y, X).fit(disp=0, maxiter=200)
    rows = []
    for i, nm in enumerate(names):
        rows.append({"analysis": "detectability gradient within mouse-detected proteome", "term": nm,
                     "OR": float(np.exp(res.params[i])),
                     "ci_lo": float(np.exp(res.params[i] - 1.96 * res.bse[i])),
                     "ci_hi": float(np.exp(res.params[i] + 1.96 * res.bse[i])),
                     "p": float(res.pvalues[i])})
    return sub, rows


# ---------------- Run ----------------
d = universe()
n_detected_in_universe = int(d["_lcv_mouse"].sum())
cov_med = float(np.nanmedian(d.loc[d["_lcv_mouse"], "_cov"]))

full_mouse = matched_fold(d, "_lcv_mouse")

# high-confidence exposure: detections with coverage >= median; low-cov detections dropped
d_hc = d.copy()
low_conf = d_hc["_lcv_mouse"] & (d_hc["_cov"] < cov_med)
d_hc["_lcv_mouse_hc"] = d_hc["_lcv_mouse"] & (d_hc["_cov"] >= cov_med)
d_hc = d_hc[~low_conf].reset_index(drop=True)
hc_mouse = matched_fold(d_hc, "_lcv_mouse_hc")

grad_sub, grad_rows = gradient_within_detected(d)

# ---------------- Abundance dose-response / abundance-adjustment (must-fix) ----
# A strict abundance-adjusted detection effect is not identifiable: the not-detected
# proteins have no LFQ value, so abundance cannot be held fixed across the detection
# boundary. The strongest feasible substitutes are (i) a graded dose-response over an
# ordinal exposure that starts at "not detected" and steps up through LFQ tertiles among
# detections, (ii) an internal abundance contrast among detections alone, and (iii) a
# comparison of how much of the detection-level fold a within-detection abundance range
# can reproduce.
det = d[d["_lcv_mouse"] & d["_log_lfq"].notna()].copy()
lfq_cuts = det["_log_lfq"].quantile([1 / 3, 2 / 3]).tolist()
lfq_med = float(det["_log_lfq"].median())
lfq_iqr = float(det["_log_lfq"].quantile(0.75) - det["_log_lfq"].quantile(0.25))
lfq_rng = float(det["_log_lfq"].quantile(0.95) - det["_log_lfq"].quantile(0.05))

d_dr = d.copy()
lvl = pd.Series(0, index=d_dr.index)
lvl.loc[det.index] = np.digitize(det["_log_lfq"], lfq_cuts) + 1  # 1,2,3
d_dr["_abund_level"] = lvl
level_labels = {0: "Not detected", 1: "LFQ tertile 1 (low)", 2: "LFQ tertile 2 (mid)", 3: "LFQ tertile 3 (high)"}
dr_rows = []
for L in [0, 1, 2, 3]:
    sub = d_dr[d_dr["_abund_level"] == L]
    rate = float(sub["_crispr_clean"].mean()) if len(sub) else np.nan
    dr_rows.append({"level": L, "level_label": level_labels[L], "n": int(len(sub)),
                    "crispr_clean_rate": rate})
df_dr = pd.DataFrame(dr_rows)

# Ordinal trend: logistic CRISPR-clean ~ abundance level + audited covariates.
Xtr = sm.add_constant(np.column_stack([
    d_dr["_abund_level"].astype(float).values,
    np.log1p(d_dr["_deg"].clip(lower=0).replace(0, 1e-9)),
    np.log1p(d_dr["_pmids"].clip(lower=0).replace(0, 1e-9)),
    d_dr["_chronos"].astype(float),
    d_dr["_common_ess"].astype(float),
    d_dr["_ui_expr"].astype(float),
]))
ytr = d_dr["_crispr_clean"].astype(int).values
trend = sm.Logit(ytr, Xtr).fit(disp=0, maxiter=200)
trend_or = float(np.exp(trend.params[1]))
trend_ci = (float(np.exp(trend.params[1] - 1.96 * trend.bse[1])),
            float(np.exp(trend.params[1] + 1.96 * trend.bse[1])))
trend_p = float(trend.pvalues[1])

# Internal abundance contrast among detections: high-LFQ half vs low-LFQ half, matched.
d_det = det.copy()
d_det["_lfq_high"] = (d_det["_log_lfq"] >= lfq_med).astype(int)
lfq_split = matched_fold(d_det, "_lfq_high")

# Negative control: does mouse LCV detection also predict an unrelated outcome
# (cross-cohort transcriptional response) after the same matching? If detection were a
# generic detectability artefact it would inflate any correlated outcome.
neg_ctrl = matched_fold(d, "_lcv_mouse", outcome_col="_resp")

# How much of the detection-level fold a within-detection abundance range reproduces,
# extrapolating the fitted per-unit log10-LFQ OR across the observed detected range.
grad_or_lfq = float(pd.DataFrame(grad_rows).set_index("term").loc["log10_LFQ", "OR"])
fold_abund_iqr = grad_or_lfq ** lfq_iqr
fold_abund_range = grad_or_lfq ** lfq_rng

df_grad = pd.DataFrame(grad_rows)
df_grad.to_csv(OUT / "p0b_detectability_gradient.tsv", sep="\t", index=False)

# Abundance dose-response table + summary
df_dr.to_csv(OUT / "p0b_abundance_dose_response.tsv", sep="\t", index=False)

summary = {
    "n_universe": int(len(d)),
    "n_mouse_detected_in_universe": n_detected_in_universe,
    "n_human_layer": "presence-only (no LFQ/coverage) - cannot be detectability-adjusted",
    "mouse_seq_cov_median_among_detected": cov_med,
    "full_mouse_exposure": full_mouse,
    "highconf_mouse_exposure": hc_mouse,
    "n_universe_highconf": int(len(d_hc)),
    "n_query_highconf": int(d_hc["_lcv_mouse_hc"].sum()),
    "gradient_log_lfq_OR": float(df_grad.loc[df_grad.term == "log10_LFQ", "OR"].iloc[0]),
    "gradient_log_lfq_p": float(df_grad.loc[df_grad.term == "log10_LFQ", "p"].iloc[0]),
    "gradient_seq_cov_OR": float(df_grad.loc[df_grad.term == "seq_coverage_pct", "OR"].iloc[0]),
    "gradient_seq_cov_p": float(df_grad.loc[df_grad.term == "seq_coverage_pct", "p"].iloc[0]),
    "n_detected_for_gradient": int(len(grad_sub)),
    # ---- abundance dose-response / adjustability ----
    "lfq_tertile_cuts_log10": [float(x) for x in lfq_cuts],
    "dose_response_levels": df_dr.to_dict("records"),
    "dose_response_trend_OR_per_level": trend_or,
    "dose_response_trend_CI": list(trend_ci),
    "dose_response_trend_p": trend_p,
    "within_detected_lfq_halves_matched": lfq_split,
    "negative_control_detection_on_response": neg_ctrl,
    "log_lfq_iqr_detected": lfq_iqr,
    "log_lfq_p5_p95_range_detected": lfq_rng,
    "fold_from_abundance_iqr": fold_abund_iqr,
    "fold_from_abundance_full_range": fold_abund_range,
    "abundance_identifiability": (
        "A strict abundance-held-fixed detection contrast is not identifiable: the "
        "not-detected stratum has no LFQ value, so abundance cannot be conditioned across "
        "the detection boundary and no whole-cell reference proteome is available. The "
        "dose-response and internal-contrast analyses bound how much of the detection-level "
        "fold a within-detection abundance range can reproduce."),
}
(OUT / "p0b_detectability_adjusted_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
(OUT / "p0b_abundance_dose_response_summary.json").write_text(
    json.dumps({k: summary[k] for k in summary if k not in ("full_mouse_exposure", "highconf_mouse_exposure")},
               indent=2, default=str), encoding="utf-8")

print(json.dumps(summary, indent=2, default=str))
print(df_grad.to_string(index=False))
print(df_dr.to_string(index=False))
