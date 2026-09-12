from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"

# Post-freeze literature triangulation. These studies were not used to choose
# candidates or assign tiers. Metadata and evidence statements were checked
# against PubMed/PMC records; unsupported positives are not inferred.
cols = [
    "candidate", "independent_study", "pmid_pmcid", "doi",
    "experimental_model", "evidence_type", "genetic_perturbation",
    "protein_interaction", "lcv_localization", "intracellular_replication",
    "microscopy", "rescue", "overlaps_our_data", "used_for_tier_assignment",
    "verification_note",
]

rows = [
    ["RAB1A", "Kagan et al. — Legionella subvert Rab1 and Sec22b to create a replicative organelle",
     "15117975; PMC2211909", "10.1084/jem.20031706", "mammalian-cell infection",
     "localization + genetic inhibition + intracellular growth", "yes (Rab1 functional inhibition)",
     "not tested as an effector-Rab1 binding assay", "yes", "yes", "yes", "not reported",
     "LCV/function", "NO", "PubMed abstract directly supports rapid LCV recruitment, ER-vesicle recruitment and intracellular-growth phenotypes."],
    ["SEC22B", "Kagan et al. — Legionella subvert Rab1 and Sec22b to create a replicative organelle",
     "15117975; PMC2211909", "10.1084/jem.20031706", "mammalian-cell infection",
     "localization + functional perturbation", "yes (Sec22b functional perturbation)",
     "not a direct bacterial-effector binding assay", "yes", "yes", "yes", "not reported",
     "LCV/function", "NO", "PubMed abstract supports LCV delivery and a role in replicative-organelle biogenesis."],
    ["RAB1A", "Suh et al. — structural insights into SidM/DrrA nucleotide exchange and GDI displacement",
     "19942850; PMC2824451", "10.1038/emboj.2009.347", "biochemistry/structure",
     "direct biochemical/structural mechanism", "no", "yes (SidM/DrrA-Rab1 mechanism)",
     "not applicable", "not applicable", "no", "not applicable", "PPI-A", "NO",
     "Supports GEF/GDI-displacement activity; AMPylation is not assigned to this PMID."],
    ["RAB1A", "Wang et al. — SetA glucosyltransferase regulation of Rab1",
     "30323948; PMC6175885", "10.1038/s41421-018-0055-9", "biochemistry/cell biology",
     "direct post-translational modification", "no", "yes (SetA-Rab1 mechanism)",
     "not applicable", "not applicable", "yes", "not reported", "PPI-A/mechanism", "NO",
     "Directly supports SetA-mediated Rab1 glucosylation."],
    ["ARF1", "Folly-Klan et al. — membrane sensor controls RalF localization and ArfGEF activity",
     "24244168; PMC3828167", "10.1371/journal.ppat.1003747", "cell biology/biochemistry",
     "effector localization + ArfGEF mechanism", "yes (RalF variants)", "yes (RalF-Arf pathway)",
     "yes", "not primary endpoint", "yes", "not reported", "PPI/LCV", "NO",
     "Supports RalF localization and ArfGEF activity, not a clean-CRISPR dependency."],
    ["RAB10", "Jeng et al. — genome-wide CRISPR screen of Legionella pathogenesis",
     "31540829; PMC6800164", "10.1016/j.chom.2019.08.017", "U937 CRISPR/infection",
     "genetic dependency + infection phenotyping", "yes (CRISPR)", "yes (SidC/SdcA-related evidence)",
     "yes", "yes", "yes", "yes", "CRISPR/LCV/AP-MS", "NO",
     "Supports RAB10/RABIF requirement for optimal replication, Rab10 LCV recruitment and SidC/SdcA-associated ubiquitination."],
    ["RABIF", "Jeng et al. — genome-wide CRISPR screen of Legionella pathogenesis",
     "31540829; PMC6800164", "10.1016/j.chom.2019.08.017", "U937 CRISPR/infection",
     "genetic dependency + infection phenotyping", "yes (CRISPR)", "yes (Rab10 stabilization pathway)",
     "not separately established in abstract", "yes", "reported in study", "yes",
     "CRISPR/AP-MS", "NO", "Independent anchor; not used to set the candidate tier."],
    ["RAB8A", "Hoffmann et al. — functional analysis of Rab GTPases in purified LCVs",
     "24373249", "10.1111/cmi.12256", "RAW264.7 macrophage LCV proteomics/RNAi",
     "LCV localization + RNAi replication phenotype", "yes (RNAi)", "no direct effector binding",
     "yes", "yes", "yes", "not reported", "LCV/function", "NO",
     "PubMed abstract reports that Rab8a depletion impairs replication and reduces SidC-positive LCVs."],
    ["RAB10", "Hoffmann et al. — functional analysis of Rab GTPases in purified LCVs",
     "24373249", "10.1111/cmi.12256", "RAW264.7 macrophage LCV proteomics/RNAi",
     "LCV localization + RNAi replication phenotype", "yes (RNAi)", "no direct effector binding",
     "yes", "yes", "yes", "not reported", "LCV/function", "NO",
     "PubMed abstract reports LCV localization and a replication-promoting RNAi phenotype."],
    ["RTN4/ATL", "Steiner et al. — ER remodeling by atlastin promotes vacuolar growth",
     "28835546; PMC5623866", "10.15252/embr.201743903", "cell infection/ER remodeling",
     "localization + functional perturbation", "yes", "no direct effector binding",
     "yes", "yes", "yes", "not reported", "LCV/function", "NO",
     "Supports an atlastin-dependent ER-remodeling module; isoform claims must follow the paper's experiments."],
    ["LAMTOR1", "Prevost et al. — WipB targets lysosomal nutrient-sensing machinery",
     "28842705; PMC5572681", "10.1038/s41598-017-10249-6", "biochemistry/cell biology",
     "effector interaction/localization", "no host knockout claimed", "yes (WipB-associated pathway)",
     "reported for WipB-associated compartment", "not primary endpoint", "yes", "not reported",
     "AP-MS/localization", "NO", "Does not remove common-essential or infection-specificity caveats."],
    ["VAPA", "No candidate-specific Legionella functional validation located in the targeted audit", "", "",
     "not located", "no verified independent anchor", "no", "no", "no", "no", "no", "no",
     "computational-only", "NO", "Absence from this targeted audit is not proof of biological absence."],
    ["ACTR3", "No candidate-specific Legionella functional validation located in the targeted audit", "", "",
     "not located", "no verified independent anchor", "no", "no", "no", "no", "no", "no",
     "computational-only", "NO", "The frozen four-gene signal must not imply equal experimental support."],
    ["ACBD3", "No candidate-specific Legionella functional validation located in the targeted audit", "", "",
     "not located", "no verified independent anchor", "no", "no", "no", "no", "no", "no",
     "computational-only", "NO", "Mouse-only LCV support remains a species-extrapolation caveat."],
]

df = pd.DataFrame(rows, columns=cols)
assert (df["used_for_tier_assignment"] == "NO").all()
df.to_csv(OUT / "d4_orthogonal_literature_evidence_v2.tsv", sep="\t", index=False)
print(df[["candidate", "pmid_pmcid", "evidence_type", "used_for_tier_assignment"]].to_string(index=False))
