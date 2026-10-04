"""Annotation-opportunity / broad-membrane matched null for the holdout panel (post-freeze #3).

Purpose
-------
The frozen 56-gene candidate panel is enriched for RAB GTPases and membrane-trafficking
genes (10 canonical RAB GTPases / 12 trafficking observed; selection-aware null medians
1 / 4; f7_selection_aware_enrichment.json). A trivial alternative explanation is
*annotation opportunity*: genes that are heavily annotated in GO are more likely to be
recorded in both the LCV proteome layer and the CRISPR hit layer, and GO annotation
density correlates with how "membrane-ish" a gene looks. This script asks whether the
RAB-GTPase / trafficking enrichment survives matching on that opportunity.

Two nested nulls
----------------
Null A : annotation-opportunity null. Pseudo-panels of 56 unique genes matched on
         GO annotation counts (BP, CC, total distinct terms), generic membrane
         (GO:0016020 + descendants) yes/no, generic intracellular membrane-bounded
         organelle (GO:0043231 + descendants) yes/no, baseline UI expression, expression
         study count, cross-pathogen network degree, cross-pathogen publication count and
         CRISPR functional tier.
Null B : Null A + broad membrane context, i.e. the *depth/breadth* of the membrane
         annotation itself (number of distinct membrane and IMBO descendant terms, and the
         membrane fraction of the gene's annotation profile). This is the "the gene simply
         looks membranous in every annotation system" adjustment.

Deliberately NOT matched (would be circular)
--------------------------------------------
RAB membership, vesicle-mediated transport, ER-to-Golgi transport, Rab protein signal
transduction, or any membrane-trafficking category. Matching on the structure under test
would make the enrichment untestable by construction.

Method
------
Identical machinery to postfreeze_f7_selection_aware_enrichment.py: standardised
covariates -> LogisticRegression propensity for panel membership -> per-query control
pool at caliper 0.02, relaxed to 0.05, then the full control set -> 10,000 unique-gene
pseudo-panel draws -> empirical P((null count) >= observed) = (exceedances + 1)/(N + 1).
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trafficking_universe import RAB_GTPASE, ALL_GENES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
GAF = ROOT / "01_data_raw" / "annotations" / "goa_human.gaf.gz"
OBO = ROOT / "01_data_raw" / "annotations" / "go-basic.obo"
RNG = np.random.default_rng(42)
N_DRAW = 10000

MEMBRANE_ROOT = "GO:0016020"           # membrane
IMBO_ROOT = "GO:0043231"               # intracellular membrane-bounded organelle


# ---------------- GO ontology: parents, then descendants of the two roots ----------------
def load_obo(path: Path) -> dict[str, list[str]]:
    """Return {term_id: [parent_id, ...]} from is_a and part_of relations."""
    parents: dict[str, list[str]] = {}
    cur: str | None = None
    obsolete: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
                continue
            if line.startswith("["):
                cur = None
                continue
            if not line:
                continue
            if line.startswith("id: GO:"):
                cur = line[4:].strip()
                parents.setdefault(cur, [])
            elif cur is None:
                continue
            elif line.startswith("is_a: GO:"):
                parents[cur].append(line[6:].split("!")[0].strip())
            elif line.startswith("relationship: part_of GO:"):
                parents[cur].append(line[len("relationship: part_of "):].split("!")[0].strip())
            elif line.startswith("is_obsolete: true"):
                obsolete.add(cur)
    for t in obsolete:
        parents.pop(t, None)
    return parents


def descendants(parents: dict[str, list[str]], root: str) -> set[str]:
    children: dict[str, list[str]] = {}
    for child, pars in parents.items():
        for p in pars:
            children.setdefault(p, []).append(child)
    out: set[str] = set()
    stack = [root]
    while stack:
        t = stack.pop()
        for c in children.get(t, []):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


# ---------------- GAF: per-gene distinct GO terms (NOT-qualified rows dropped) ----------------
def load_gaf(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return (gene -> distinct GO terms by aspect P/C, gene -> all distinct GO terms).

    Only human (taxon:9606) protein-coding associations whose qualifier does not contain
    NOT are used. Distinct terms are counted so that reference density cannot inflate the
    annotation-opportunity measure.
    """
    by_aspect: dict[str, set[str]] = {"P": set(), "C": set()}
    all_terms: dict[str, set[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            qual = f[3]
            if qual and "NOT" in qual.split("|"):
                continue
            gene, go_id, aspect = f[2], f[4], f[8]
            if not gene or not go_id.startswith("GO:"):
                continue
            all_terms.setdefault(gene, set()).add(go_id)
            if aspect in by_aspect:
                by_aspect[aspect].add(gene + "\x00" + go_id)
    per_gene_aspect = {"P": {}, "C": {}}
    for asp, keys in by_aspect.items():
        for k in keys:
            g, go = k.split("\x00")
            per_gene_aspect[asp].setdefault(g, set()).add(go)
    return per_gene_aspect, all_terms


print("parsing go-basic.obo ...")
parents = load_obo(OBO)
membrane_desc = descendants(parents, MEMBRANE_ROOT)
imbo_desc = descendants(parents, IMBO_ROOT)
print(f"  GO:0016020 descendants: {len(membrane_desc)}")
print(f"  GO:0043231 descendants: {len(imbo_desc)}")

print("parsing goa_human.gaf.gz ...")
gaf_aspect, gaf_all = load_gaf(GAF)
print(f"  genes with GO annotation: {len(gaf_all)}")


# ---------------- Atlas universe (identical to f7) ----------------
atlas = pd.read_csv(ROOT / "08_evidence_atlas" / "host_dependency_evidence_atlas.tsv.gz", sep="\t", low_memory=False)
deg = pd.read_csv(ROOT / "03_interactome" / "cross_pathogen_degree_target_excluded.tsv", sep="\t", low_memory=False)
# Discovery/candidate null uses the truly frozen 56-gene panel.
panel = sorted(set(pd.read_csv(ROOT / "08_evidence_atlas" / "FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv", sep="\t")
                   .loc[:, "human_symbol"].dropna().astype(str)))


def B(s):
    return s.astype(str).str.upper().isin(["TRUE", "1", "YES"])


atlas["_ppi"] = B(atlas["ppi_AB_flag"])
atlas["_apms"] = B(atlas["apms_high_confidence"])
atlas["_lcv"] = B(atlas["lcv_spatial_supported"])
atlas["_resp"] = B(atlas["stable_cross_cohort"])
atlas["_hit"] = atlas["functional_tier"].astype(str).str.startswith("F")
atlas["_nstud"] = pd.to_numeric(atlas["expression_n_studies"], errors="coerce").fillna(0)
atlas["_common"] = B(atlas["common_essential_flag"])
atlas["_chronos"] = pd.to_numeric(atlas["u937_chronos_gene_effect"], errors="coerce")
atlas = atlas.merge(deg[["gene", "cross_pathogen_degree", "cross_pathogen_pmids"]],
                    left_on="human_symbol", right_on="gene", how="left")
atlas["_deg"] = pd.to_numeric(atlas["cross_pathogen_degree"], errors="coerce").fillna(0)
atlas["_pmids"] = pd.to_numeric(atlas["cross_pathogen_pmids"], errors="coerce").fillna(0)


def load_ui(mat, gsm, gmap, probe):
    m = pd.read_csv(ROOT / "05_expression" / "matrices" / mat, sep="\t", compression="gzip")
    idc = "feature_id" if "feature_id" in m.columns else "Geneid"
    cols = [c for c in gsm if c in m.columns]
    gm = pd.read_csv(ROOT / "05_expression" / "gene_maps" / gmap, sep="\t")
    gm = gm[gm.gene_symbol.notna()]
    if probe and "ambiguous_mapping" in gm.columns:
        gm = gm[~gm.ambiguous_mapping.astype(str).str.upper().isin(["TRUE", "1"])]
    if "feature_id" in gm.columns and idc != "feature_id":
        gm = gm.rename(columns={"feature_id": idc})
    mg = m.merge(gm, on=idc, how="inner")
    sub = mg[["gene_symbol"] + cols].groupby("gene_symbol").mean()
    return sub.rank(pct=True)[cols].mean(axis=1)


ui = pd.concat([
    load_ui("GSE61535_log2_expression.tsv.gz", ["GSM1507677", "GSM1507678", "GSM1507679"], "GSE61535_gene_map.tsv", True),
    load_ui("GSE135803_fpkm.tsv.gz", ["GSM4030445", "GSM4030446", "GSM4030447"], "GSE135803_gene_map.tsv", False),
    load_ui("GSE287755_raw_counts.tsv.gz", ["GSM8750874", "GSM8750875", "GSM8750876"], "GSE287755_gene_map.tsv", False),
], axis=1).mean(axis=1, skipna=True)
atlas["_ui"] = atlas.human_symbol.map(ui).fillna(0)

atlas["_crispr_tier"] = atlas["functional_tier"].astype(str).map(
    {"F1_GW_AND_BOTH_BATCH_CONCORDANT": 3, "F2_MULTILAYER_CONCORDANT": 2,
     "F3_SINGLE_LAYER_OR_DISCORDANT": 1, "NO_FUNCTIONAL_HIT": 0, "NO_FUNCTIONAL_RECORD": 0}).fillna(0)

# ---- GO annotation-opportunity features ----
sym = atlas["human_symbol"].astype(str)
atlas["_go_bp_n"] = sym.map(lambda g: len(gaf_aspect["P"].get(g, ()))).astype(float)
atlas["_go_cc_n"] = sym.map(lambda g: len(gaf_aspect["C"].get(g, ()))).astype(float)
atlas["_go_terms_n"] = sym.map(lambda g: len(gaf_all.get(g, ()))).astype(float)
atlas["_n_membrane_terms"] = sym.map(lambda g: len(gaf_all.get(g, set()) & membrane_desc)).astype(float)
atlas["_n_imbo_terms"] = sym.map(lambda g: len(gaf_all.get(g, set()) & imbo_desc)).astype(float)
atlas["_membrane"] = (atlas["_n_membrane_terms"] > 0).astype(float)
atlas["_imbo"] = (atlas["_n_imbo_terms"] > 0).astype(float)
atlas["_membrane_frac"] = np.where(atlas["_go_terms_n"] > 0,
                                   atlas["_n_membrane_terms"] / atlas["_go_terms_n"], 0.0)

collapse = {
    "_ppi": "max", "_apms": "max", "_hit": "max", "_lcv": "max", "_resp": "max",
    "_crispr_tier": "max", "_nstud": "max", "_ui": "median", "_deg": "max",
    "_pmids": "max", "_common": "max", "_chronos": "median",
    "_go_bp_n": "max", "_go_cc_n": "max", "_go_terms_n": "max",
    "_n_membrane_terms": "max", "_n_imbo_terms": "max",
    "_membrane": "max", "_imbo": "max", "_membrane_frac": "max",
}
d = atlas.groupby("human_symbol", as_index=False, sort=True).agg(collapse)
d["_chronos"] = d["_chronos"].fillna(d["_chronos"].median())
# genes absent from the GAF are genuinely annotation-poor, not missing: keep them at 0.
for c in ["_go_bp_n", "_go_cc_n", "_go_terms_n", "_n_membrane_terms", "_n_imbo_terms", "_membrane_frac"]:
    d[c] = d[c].fillna(0.0)
for c in ["_membrane", "_imbo"]:
    d[c] = d[c].fillna(0.0)

BASE = ["_go_bp_n", "_go_cc_n", "_go_terms_n", "_membrane", "_imbo",
        "_ui", "_nstud", "_deg", "_pmids", "_crispr_tier"]
NULLS = {
    "nullA_annotation_opportunity": BASE,
    "nullB_broad_membrane_context": BASE + ["_n_membrane_terms", "_n_imbo_terms", "_membrane_frac"],
}

d["_in_panel"] = d.human_symbol.isin(panel).astype(int)
if not d.human_symbol.is_unique:
    raise RuntimeError("Annotation-opportunity null universe must contain unique human symbols")
if int(d["_in_panel"].sum()) != len(panel):
    missing = sorted(set(panel) - set(d.loc[d["_in_panel"] == 1, "human_symbol"]))
    raise RuntimeError(f"Panel genes missing from null universe: {missing}")

rab_set = set(RAB_GTPASE)
traf_set = set(ALL_GENES)
obs_rab = len(set(panel) & rab_set)
obs_traf = len(set(panel) & traf_set)


def run_null(feat: list[str], tag: str) -> dict:
    Z = d[feat].astype(float)
    Z = (Z - Z.mean()) / Z.std().replace(0, 1)
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
    clf.fit(Z.values, d["_in_panel"].values)
    prop = clf.predict_proba(Z.values)[:, 1]
    q_idx = np.where(d["_in_panel"].values == 1)[0]
    ctl_idx = np.where(d["_in_panel"].values == 0)[0]
    pools = {int(i): np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.02]) for i in q_idx}
    calib = {"0.02": 0, "0.05": 0, "all": 0}
    for i, p in pools.items():
        if len(p) < 20:
            pools[i] = np.array([j for j in ctl_idx if abs(prop[j] - prop[i]) <= 0.05])
            calib["0.05"] += 1
        if len(p) < 20:
            pools[i] = ctl_idx
            calib["all"] += 1
        else:
            calib["0.02"] += 1
    pool_sizes = np.array([len(pools[int(i)]) for i in q_idx], dtype=float)

    # balance check: standardised mean difference of the panel vs universe, and the mean
    # SMD of the drawn pseudo-panels vs the panel (matched "after" balance).
    Zf = d[feat].astype(float).values
    pan_mean = Zf[q_idx].mean(axis=0)
    pan_sd = Zf[q_idx].std(axis=0, ddof=1)
    uni_mean = Zf[ctl_idx].mean(axis=0)
    uni_sd = Zf[ctl_idx].std(axis=0, ddof=1)
    sd_pool = np.sqrt((pan_sd ** 2 + uni_sd ** 2) / 2)
    sd_pool[sd_pool == 0] = 1.0
    smd_before = (pan_mean - uni_mean) / sd_pool

    counts_rab = np.empty(N_DRAW)
    counts_traf = np.empty(N_DRAW)
    matched_mean = np.zeros((N_DRAW, len(feat)))
    for k in range(N_DRAW):
        used = set()
        chosen = []
        for i in q_idx:
            avail = [j for j in pools[int(i)] if j not in used]
            if not avail:
                break
            c = int(RNG.choice(avail))
            chosen.append(c)
            used.add(c)
        genes = d.iloc[chosen].human_symbol
        counts_rab[k] = len(genes[genes.isin(rab_set)])
        counts_traf[k] = len(genes[genes.isin(traf_set)])
        matched_mean[k] = Zf[np.array(chosen)].mean(axis=0)
        if (k + 1) % 2500 == 0:
            print(f"    {tag}: draw {k + 1}/{N_DRAW}")

    smd_after = (pan_mean - matched_mean.mean(axis=0)) / sd_pool
    p_rab = float(((counts_rab >= obs_rab).sum() + 1) / (N_DRAW + 1))
    p_traf = float(((counts_traf >= obs_traf).sum() + 1) / (N_DRAW + 1))
    print(f"  {tag}: median RAB={np.median(counts_rab):.1f} p95={np.percentile(counts_rab, 95):.0f} "
          f"max={counts_rab.max():.0f} | median traf={np.median(counts_traf):.1f} "
          f"p95={np.percentile(counts_traf, 95):.0f} max={counts_traf.max():.0f}")
    print(f"  {tag}: P(#RAB>={obs_rab})={p_rab:.5f}  P(#traf>={obs_traf})={p_traf:.5f}  "
          f"max|SMD_before|={np.abs(smd_before).max():.3f}  max|SMD_after|={np.abs(smd_after).max():.3f}")

    return {
        "null": tag, "covariates": feat, "n_covariates": len(feat), "n_draws": N_DRAW,
        "universe_n": int(len(d)), "n_panel": len(panel),
        "pool_caliper_counts": calib,
        "pool_size_median": float(np.median(pool_sizes)), "pool_size_min": float(pool_sizes.min()),
        "obs_rab_in_panel": obs_rab, "obs_trafficking_in_panel": obs_traf,
        "null_rab_median": float(np.median(counts_rab)),
        "null_rab_p95": float(np.percentile(counts_rab, 95)),
        "null_rab_max": float(counts_rab.max()),
        "null_traf_median": float(np.median(counts_traf)),
        "null_traf_p95": float(np.percentile(counts_traf, 95)),
        "null_traf_max": float(counts_traf.max()),
        "empirical_p_rab": p_rab, "empirical_p_trafficking": p_traf,
        "empirical_p_convention": "(exceedances + 1)/(N + 1)",
        "bonferroni_adjusted_p_rab": min(1.0, 2 * p_rab),
        "bonferroni_adjusted_p_trafficking": min(1.0, 2 * p_traf),
        "max_abs_smd_before": float(np.abs(smd_before).max()),
        "max_abs_smd_after": float(np.abs(smd_after).max()),
        "smd_before": {f: float(v) for f, v in zip(feat, smd_before)},
        "smd_after": {f: float(v) for f, v in zip(feat, smd_after)},
        "_counts_rab": counts_rab, "_counts_traf": counts_traf,
    }


results = {}
for tag, feat in NULLS.items():
    print(f"\n[{tag}]")
    results[tag] = run_null(feat, tag)

# ---- observed annotation profile of the panel vs universe (honest context) ----
panel_rows = d[d["_in_panel"] == 1]
univ_rows = d[d["_in_panel"] == 0]
prof = {}
for c in ["_go_bp_n", "_go_cc_n", "_go_terms_n", "_n_membrane_terms", "_n_imbo_terms"]:
    prof[c] = {"panel_median": float(panel_rows[c].median()),
               "universe_median": float(univ_rows[c].median())}
for c in ["_membrane", "_imbo"]:
    prof[c] = {"panel_frac": float(panel_rows[c].mean()),
               "universe_frac": float(univ_rows[c].mean())}

for tag in results:
    results[tag].pop("_counts_rab", None)
    results[tag].pop("_counts_traf", None)

summary = {
    "analysis": "Annotation-opportunity / broad-membrane matched null for the 56-gene frozen candidate panel",
    "panel_source": "08_evidence_atlas/FROZEN_CANDIDATE_PANEL_PRE_HOLDOUT.tsv (56-gene frozen panel)",
    "annotation_source": {
        "gaf": "goa_human.gaf.gz (GOC, date-generated 2026-05-21)",
        "obo": "go-basic.obo (releases/2026-07-26)",
        "membrane_root": MEMBRANE_ROOT, "membrane_descendants": len(membrane_desc),
        "imbo_root": IMBO_ROOT, "imbo_descendants": len(imbo_desc),
        "not_qualified_rows_dropped": True,
        "counts_are_distinct_terms": True,
    },
    "deliberately_not_matched": [
        "RAB membership", "vesicle-mediated transport", "ER-to-Golgi transport",
        "Rab protein signal transduction", "any membrane-trafficking category",
    ],
    "matching_machinery": ("identical to postfreeze_f7_selection_aware_enrichment.py: "
                           "standardised covariates -> logistic propensity -> per-query control "
                           "pool (caliper 0.02 -> 0.05 -> all) -> 10000 unique-gene pseudo-panels"),
    "panel_annotation_profile": prof,
    "nulls": results,
    "interpretation": ("Null A removes annotation opportunity; Null B additionally removes the "
                       "depth/breadth of generic membrane annotation. If the RAB/trafficking "
                       "enrichment survives both, it is not an annotation-opportunity artefact. "
                       "Unmeasured ascertainment and the curated nature of the trafficking "
                       "universe remain possible limitations."),
}
(OUT / "h2_annotation_opportunity_null.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

pd.DataFrame({c: prof[c] for c in prof}).to_csv(
    OUT / "h2_panel_annotation_profile.tsv", sep="\t")

print("\n" + json.dumps(summary, indent=2, ensure_ascii=False))
print("\nwrote", OUT / "h2_annotation_opportunity_null.json")
