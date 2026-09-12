from __future__ import annotations

from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "13_postfreeze_interpretation"

# 8-item formal Biological Confound Audit (structure: Evidence / Audit result / Residual / Consequences)
rows = [
    ["Intracellular replication / bacterial burden",
     "Originating investigators (2024/2025 ViIN) report all genomic deletion mutants replicate efficiently in PMA-THP-1 with no gross defect; Wibawa 2022 warns growth/burden can confound; d2 shows panel response correlates with global perturbation (rho=0.685, P=0.014)",
     "PARTIALLY RESOLVED / RESIDUAL QUANTITATIVE UNCERTAINTY",
     "No 18-h per-strain quantitative CFU; Wibawa2022 mouse-lung CFU is a different estimand and cannot be applied to THP-1/18h",
     "Genomic-deletion-context responses are supported; subtle burden differences cannot be formally excluded; claim requires conservative wording"],
    ["Global host-response amplitude",
     "Panel response tracks global transcriptome perturbation across 12 contrasts (Spearman rho=0.685, P=0.014); PASS groups have slightly higher global |z| (1.029 vs 0.936, MW P=0.021)",
     "CORRELATED / NOT FULLY INDEPENDENT",
     "Causal direction between panel and global amplitude unresolved; matched-expression null not yet run",
     "Qualify effector-specific claims; describe as response-amplitude sensitivity rather than load-proven confound"],
    ["Cell context (THP-1 vs discovery models)",
     "Holdout uses PMA-differentiated THP-1 (cell line); discovery uses MDM/other macrophage models and mouse; cross-context by design",
     "ORTHOGONAL CONTEXT (strengthens external validity)",
     "THP-1 is not primary macrophage; differentiation differences remain",
     "Call it cross-context external validation, not physiological validation"],
    ["Species composition (human vs mouse)",
     "b2_species_dependency_scores.tsv annotates all 56 candidates; RAB1A LCV support is human-only core; ACBD3/ATL2 LCV support is mouse-only; NAIP/NLRC4 axis species-boundary flagged",
     "EXPLICIT / SYSTEMATICALLY ANNOTATED",
     "Mouse-only-supported candidates cannot be extrapolated to human infection",
     "Human-only sensitivity reported; species extrapolation boundary stated"],
    ["Single-cohort dominance",
     "GSE135803 provides FDR hits for 65.6% and GSE287755 for 60.2% of the 2,601 meta-convergent genes; only 432 of 2,601 are significant in >=3 datasets",
     "CONFIRMED CONCENTRATION",
     "Meta-convergence layer is dataset-collapsed, not per-gene multi-study replication",
     "Rename layer to META-CONVERGENT RESPONSE SET (2,601) vs >=3-DATASET REPLICATED CORE (490); do not claim per-gene replication"],
    ["Mutant genotype mapping",
     "SDRF official labels (Delta-A/B/...) confirmed for 15 groups; A-F,H,I gene content mapped from preprint E-Table1 + 130b annotation; O/P/R/S and G content unresolved",
     "PARTIALLY RESOLVED",
     "O/P/R/S and G gene content not in any public source",
     "Prohibit gene-specific inference for unmapped regions; keep genomic-deletion-context response terminology"],
    ["Effector redundancy / antagonism",
     "effector_functional_relationships.tsv + region mapping show Delta-F deletes SidE family AND SidJ; Delta-H deletes SidM+SidD; Wibawa2022 lists redundancy/metaeffector antagonism as official challenge",
     "CONTEXT-DEPENDENT / PLAUSIBLE HYPOTHESIS",
     "Single-mutant negatives may be redundancy-masked; no formal enrichment test (n=7 mapped regions)",
     "No hard negative claims; redundancy-masked = inconclusive, not functional absence"],
    ["Cellular heterogeneity (bulk limitation)",
     "Bulk RNA-seq and pooled CRISPR reflect population means; permissive vs non-permissive subpopulations cannot be resolved",
     "ACKNOWLEDGED LIMITATION",
     "No matched single-cell data",
     "State as explicit limitation and future direction; no subpopulation inference"],
]
df = pd.DataFrame(rows, columns=["Confound", "Evidence", "Audit result", "Residual uncertainty", "Consequences for claims"])
df.to_csv(OUT / "e2_biological_confound_audit_table.tsv", sep="\t", index=False)
print(df.to_string(index=False))
print("\nwrote", OUT / "e2_biological_confound_audit_table.tsv")
