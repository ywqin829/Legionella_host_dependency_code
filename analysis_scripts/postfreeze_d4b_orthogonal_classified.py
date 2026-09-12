from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

# Truth-set / discovery source PMIDs used within the pipeline (so these overlap a discovery domain)
TRUTH_SOURCE = {"15117975", "24244168", "28835546", "28842705", "31540829", "PMC5364298"}
PPI_SOURCE = {"19942850", "30323948", "22228731"}
RELATION_SOURCE = {"PMC8605245"}

cols = ["candidate", "study", "pmid_pmcid", "experimental_model", "evidence_type",
        "used_in_discovery_evidence", "used_in_truth_set", "classification",
        "evidence_domains_overlapped", "used_for_tier_assignment"]
rows = [
    ["RAB1A", "Rab1 recruitment / dominant-negative replication phenotype", "15117975", "macrophage", "function/replication",
     "YES", "YES", "MECHANISTIC_CONTEXT", "truth-set/known positive control", "NO"],
    ["RAB1A", "DrrA/SidM binding + AMPylation + GEF", "19942850", "biochemistry", "direct interaction",
     "YES", "NO", "SOURCE_OVERLAPPING", "PPI(A/B) track", "NO"],
    ["RAB1A", "SidM/SidD/AnkX/Lem3 reversible Rab1 modification", "PMC8605245", "review/biochemistry", "PTM",
     "YES", "NO", "SOURCE_OVERLAPPING", "effector-relationship table (PPI context)", "NO"],
    ["SEC22B", "Sec22b LCV localisation + replicative-organelle requirement", "15117975", "macrophage", "function",
     "YES", "YES", "MECHANISTIC_CONTEXT", "truth-set/known positive control", "NO"],
    ["ARF1", "RalF-dependent Arf1 recruitment/activation at LCV", "24244168", "macrophage/biochemistry", "direct interaction",
     "YES", "YES", "MECHANISTIC_CONTEXT", "truth-set/known positive control", "NO"],
    ["RAB10", "Genome-wide U937 CRISPR screen + infection phenotyping", "31540829", "U937 CRISPR/infection", "genetic screen",
     "YES", "YES", "SOURCE_OVERLAPPING", "CRISPR track (primary source) + truth set", "NO"],
    ["RABIF", "RAB10 stabilisation and LCV replication (same screen)", "31540829", "U937 CRISPR/infection", "genetic + AP-MS",
     "YES", "YES", "SOURCE_OVERLAPPING", "CRISPR/AP-MS track", "NO"],
    ["RTN4", "LCV proteomics localisation + vacuole-maturation perturbation", "28835546", "LCV proteomics/function", "localisation + perturbation",
     "YES", "YES", "SOURCE_OVERLAPPING", "LCV spatial track + truth set", "NO"],
    ["ATL3", "LCV proteomics localisation + replication-support perturbation", "28835546", "LCV proteomics/function", "localisation + perturbation",
     "YES", "YES", "SOURCE_OVERLAPPING", "LCV spatial track + truth set", "NO"],
    ["LAMTOR1", "WipB co-purification; lysosomal nutrient-sensing machinery", "28842705", "AP-MS/localisation", "interaction + localisation",
     "YES", "YES", "SOURCE_OVERLAPPING", "AP-MS/PPI track + truth set", "NO"],
    ["TBK1", "TBK1 knockdown reduces Legionella-induced IL-6 (human macrophage)", "PMC5364298", "human macrophage KD", "function (cytokine)",
     "YES", "YES", "MECHANISTIC_CONTEXT", "response track context + truth set", "NO"],
    ["VAPA", "Dictyostelium Vap localizes to the LCV-ER contact site and affects vacuole remodeling and intracellular replication", "36588479", "Dictyostelium infection and host-gene perturbation", "ortholog localization and perturbation", "NO", "NO", "ORTHOLOGOUS_EXPERIMENTAL_CONTEXT", "post-freeze ortholog evidence; not human VAPA validation", "NO"],
    ["ACBD3", "No independent Legionella-specific experimental evidence identified in the focused post-freeze search", "", "not identified", "none identified", "NO", "NO", "NO_INDEPENDENT_EVIDENCE_IDENTIFIED", "none", "NO"],
]
df = pd.DataFrame(rows, columns=cols)
df.to_csv(OUT / "d4_orthogonal_literature_evidence_CLASSIFIED.tsv", sep="\t", index=False)
print(df[["candidate", "classification", "used_in_discovery_evidence", "used_for_tier_assignment"]].to_string(index=False))
print("\nNOTE: entries sharing a source with discovery/truth-set evidence are MECHANISTIC_CONTEXT or SOURCE_OVERLAPPING, "
      "not independent orthogonal evidence. Independent evidence is reported only when a distinct PMID-verified experimental lineage was identified; no pending classifications remain.")
