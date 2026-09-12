from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
INFILE = ROOT / "02_metadata/intact_legionella_human_records.tsv"
OUT = ROOT / "03_interactome"
CP = ROOT / "checkpoints"
MAP = Path(r"${PA_HOST_PPI_ROOT}\data\processed\human_uniprot_gene_map.tsv")
CROSS = Path(r"${PA_HOST_PPI_ROOT}\data\processed\cross_species_raw.tsv")
TARGET_TAXIDS = {446, 272624, 297245, 297246}
SEED = 42


def as_bool(value: object) -> bool:
    """Parse MITAB boolean fields without treating the string 'false' as True."""
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def bool_mask(series: pd.Series) -> pd.Series:
    return series.map(as_bool).astype(bool)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def evidence_class(row: pd.Series) -> str:
    typ = str(row["type"]).lower()
    method = str(row["method"]).lower()
    if as_bool(row["negative"]):
        return "NEGATIVE"
    reaction = any(x in typ for x in ["phosphorylation", "glycosylation", "enzymatic reaction", "ampylation", "deamidation"])
    direct = "direct interaction" in typ
    strong_method = any(x in method for x in ["x-ray crystallography", "isothermal titration", "surface plasmon resonance",
                                               "competition binding", "protein kinase assay", "glycosylase assay",
                                               "gdp/gtp exchange assay", "enzymatic study"])
    if reaction or direct or strong_method:
        return "A"
    if "physical association" in typ:
        return "B"
    return "C"


def join_unique(s: pd.Series) -> str:
    vals = sorted({str(x) for x in s if str(x) and str(x) != "nan"})
    return " | ".join(vals)


def best_class(classes: pd.Series) -> str:
    order = {"A": 0, "B": 1, "C": 2, "NEGATIVE": 3}
    return min((str(x) for x in classes), key=lambda x: order.get(x, 9))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True); CP.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INFILE, sep="\t")
    raw["evidence_class"] = raw.apply(evidence_class, axis=1)
    gene_map = pd.read_csv(MAP, sep="\t").drop_duplicates("uniprot")
    raw = raw.merge(gene_map, left_on="human_accession", right_on="uniprot", how="left")
    raw["gene"] = raw["gene"].fillna(raw["human_accession"])
    raw.to_csv(OUT / "legionella_human_records_classified.tsv", sep="\t", index=False)

    positive_mask = ~bool_mask(raw["negative"])
    pair = raw[positive_mask].groupby(["bacterial_accession", "human_accession", "gene"], as_index=False).agg(
        evidence_class=("evidence_class", best_class), n_records=("human_accession", "size"), n_methods=("method", "nunique"),
        n_pmids=("pmid", lambda x: len({p for v in x for p in str(v).split("|") if p})), methods=("method", join_unique),
        interaction_types=("type", join_unique), sources=("source", join_unique), pmids=("pmid", join_unique))
    pair.to_csv(OUT / "legionella_human_pairs_evidence.tsv", sep="\t", index=False)

    host = pair.groupby(["human_accession", "gene"], as_index=False).agg(
        best_class=("evidence_class", best_class), bacterial_partners=("bacterial_accession", "nunique"),
        interaction_records=("n_records", "sum"), methods=("methods", join_unique), pmids=("pmids", join_unique), sources=("sources", join_unique))
    host["n_pmids"] = host.pmids.map(lambda x: len({p for v in str(x).split(" | ") for p in v.split("|") if p}))
    host.to_csv(OUT / "legionella_host_target_summary.tsv", sep="\t", index=False)

    source_rows = []
    positive = raw[positive_mask].copy()
    for source in sorted(positive.source.unique()):
        kept = positive[positive.source != source]
        source_rows.append({"removed_source": source, "remaining_records": len(kept),
                            "remaining_pairs": len(kept.drop_duplicates(["bacterial_accession", "human_accession"])),
                            "remaining_host_targets": kept.human_accession.nunique(),
                            "remaining_AB_targets": kept.loc[kept.evidence_class.isin(["A", "B"]), "human_accession"].nunique()})
    pd.DataFrame(source_rows).to_csv(OUT / "source_leave_one_out.tsv", sep="\t", index=False)

    cross = pd.read_csv(CROSS, sep="\t", low_memory=False)
    cross = cross[~bool_mask(cross["negative"])].copy()
    cross = cross[~cross.taxid_a.astype(str).isin({str(x) for x in TARGET_TAXIDS})]
    cross = cross[~cross.species_a.astype(str).str.contains("legionella pneumophila", case=False, na=False)]
    degree = cross.groupby("uniprot_b").agg(cross_pathogen_degree=("taxid_a", "nunique"),
                                             cross_pathogen_records=("taxid_a", "size"),
                                             cross_pathogen_pmids=("pmid", "nunique")).reset_index().rename(columns={"uniprot_b": "human_accession"})
    degree = degree.merge(gene_map, left_on="human_accession", right_on="uniprot", how="left")
    target_set = set(host.loc[host.best_class.isin(["A", "B"]), "human_accession"])
    degree["legionella_target"] = degree.human_accession.isin(target_set)
    degree["gene"] = degree.gene.fillna(degree.human_accession)
    degree.to_csv(OUT / "cross_pathogen_degree_target_excluded.tsv", sep="\t", index=False)

    represented = degree[degree.legionella_target]
    other = degree[~degree.legionella_target]
    threshold = float(degree.cross_pathogen_degree.quantile(0.80))
    a = int((represented.cross_pathogen_degree >= threshold).sum()); b = int((represented.cross_pathogen_degree < threshold).sum())
    c = int((other.cross_pathogen_degree >= threshold).sum()); d = int((other.cross_pathogen_degree < threshold).sum())
    fish = fisher_exact([[a, b], [c, d]])
    mw = mannwhitneyu(represented.cross_pathogen_degree, other.cross_pathogen_degree, alternative="greater")

    degree["record_bin"] = pd.qcut(degree.cross_pathogen_records.rank(method="first"), q=10, labels=False, duplicates="drop")
    rng = np.random.default_rng(SEED); observed = float(represented.cross_pathogen_degree.mean()); null = []
    for _ in range(10000):
        chosen = []
        for bin_id, n in degree[degree.legionella_target].record_bin.value_counts().items():
            pool = degree[(~degree.legionella_target) & (degree.record_bin == bin_id)].cross_pathogen_degree.to_numpy()
            if len(pool):
                chosen.extend(rng.choice(pool, size=int(n), replace=len(pool) < int(n)).tolist())
        null.append(float(np.mean(chosen)) if chosen else np.nan)
    null_arr = np.asarray([x for x in null if np.isfinite(x)])
    empirical_p = float((1 + np.sum(null_arr >= observed)) / (1 + len(null_arr)))
    pd.DataFrame({"permuted_mean_degree": null_arr}).to_csv(OUT / "network_intensity_matched_null.tsv", sep="\t", index=False)

    stats = {"input_sha256": digest(INFILE), "positive_records": int(len(positive)), "unique_pairs": int(len(pair)),
             "unique_hosts": int(len(host)), "class_counts_pairs": pair.evidence_class.value_counts().to_dict(),
             "class_counts_hosts": host.best_class.value_counts().to_dict(), "represented_AB_targets": int(len(represented)),
             "top20_degree_threshold": threshold, "fisher_table": [[a, b], [c, d]], "fisher_odds_ratio": float(fish.statistic),
             "fisher_p": float(fish.pvalue), "mann_whitney_u": float(mw.statistic), "mann_whitney_p": float(mw.pvalue),
             "observed_mean_degree": observed, "other_mean_degree": float(other.cross_pathogen_degree.mean()),
             "matched_permutation_mean": float(np.mean(null_arr)), "matched_permutation_p": empirical_p}
    (OUT / "stage1_interactome_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = ["# Stage 1 Interactome and Bias Audit", "", f"- Positive records: {stats['positive_records']}",
              f"- Unique pairs: {stats['unique_pairs']}", f"- Unique host targets: {stats['unique_hosts']}",
              f"- Pair evidence classes: {stats['class_counts_pairs']}", f"- Host evidence classes: {stats['class_counts_hosts']}",
              f"- Target-excluded Fisher OR: {stats['fisher_odds_ratio']:.3f} (p={stats['fisher_p']:.4g})",
              f"- Mean degree, represented A/B targets: {observed:.3f}", f"- Mean degree, other hosts: {stats['other_mean_degree']:.3f}",
              f"- Record-intensity-matched empirical p: {empirical_p:.4g}", "",
              "Interpretation is limited to source- and ascertainment-audited prioritization; degree is not causal evidence."]
    (OUT / "STAGE1_INTERACTOME_AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage1.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
