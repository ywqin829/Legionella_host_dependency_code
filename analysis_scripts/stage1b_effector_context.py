from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
UNIPROT = ROOT / "01_data_raw" / "metadata" / "uniprot_legionella_philadelphia1.tsv"
META = ROOT / "02_metadata" / "dataset_metadata_audit.tsv"
OUT_META = ROOT / "02_metadata"
OUT_INT = ROOT / "03_interactome"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sequence_hash(sequence: object) -> str:
    seq = str(sequence).strip().upper()
    return hashlib.sha256(seq.encode("ascii")).hexdigest() if seq and seq != "NAN" else ""


def first_locus(gene_names: object) -> str:
    matches = re.findall(r"\b(?:lpg|lpp|lpl|lpc)\d+\b", str(gene_names), flags=re.I)
    return matches[0].lower() if matches else ""


def primary_gene(gene_names: object) -> str:
    tokens = str(gene_names).strip().split()
    return tokens[0] if tokens and tokens[0].lower() != "nan" else ""


RELATIONS = [
    # effector_a, effector_b, relationship, pathway, direction, evidence, source
    ("AnkX", "Lem3", "antagonistic_host_modification", "Rab1/Rab35 phosphocholination", "add_vs_remove", "direct_biochemistry", "PMC8605245"),
    ("DrrA/SidM", "SidD", "antagonistic_host_modification", "Rab1 AMPylation", "add_vs_remove", "direct_biochemistry", "PMC8605245"),
    ("DrrA/SidM", "LepB", "shared_pathway_non_equivalent", "Rab1 activation cycle", "activate_vs_inactivate", "direct_biochemistry", "PMC8605245"),
    ("Lgt1", "Lgt2", "paralog_redundancy", "eEF1A-dependent translation inhibition", "same_direction", "family_biochemistry", "10.1128/JB.01798-07"),
    ("Lgt1", "Lgt3", "paralog_redundancy", "eEF1A-dependent translation inhibition", "same_direction", "family_biochemistry", "10.1128/JB.01798-07"),
    ("Lgt2", "Lgt3", "paralog_redundancy", "eEF1A-dependent translation inhibition", "same_direction", "family_biochemistry", "10.1128/JB.01798-07"),
    ("SidI", "MesI/Lpg2505", "direct_metaeffector_regulation", "translation inhibition", "inhibit_effector", "direct_binding_and_function", "PMC7171249"),
    ("SdeA", "SidE", "paralog_redundancy", "phosphoribosyl ubiquitination", "same_direction", "family_biochemistry", "PMC6941232"),
    ("SdeA", "SidJ", "direct_metaeffector_regulation", "phosphoribosyl ubiquitination", "inhibit_effector", "direct_biochemistry", "PMC8605245"),
    ("SdeA", "DupA/Lpg2154", "antagonistic_host_modification", "phosphoribosyl ubiquitination", "add_vs_remove", "direct_biochemistry", "PMC6941232"),
    ("SdeA", "DupB/Lpg2509", "antagonistic_host_modification", "phosphoribosyl ubiquitination", "add_vs_remove", "direct_biochemistry", "PMC6941232"),
    ("MavC/Lpg2147", "MvcA/Lpg2148", "antagonistic_host_modification", "UBE2N ubiquitination", "add_vs_remove", "direct_biochemistry", "PMC7024838"),
    ("MavC/Lpg2147", "Lpg2149", "direct_metaeffector_regulation", "UBE2N ubiquitination", "inhibit_effector", "direct_binding_structure", "PDB:7BXF"),
    ("MvcA/Lpg2148", "Lpg2149", "direct_metaeffector_regulation", "UBE2N ubiquitination", "inhibit_effector", "direct_binding_structure", "PDB:7BXF"),
    ("LphD", "RomA", "cooperative_paraeffector", "histone H3K14 modification", "sequential_cooperation", "functional_biochemistry", "PMC11651115"),
    ("LegC2", "LegC3", "cooperative_paraeffector", "SNARE-like membrane fusion", "same_complex", "physical_and_functional", "PMC11651115"),
    ("LegC2", "LegC7", "cooperative_paraeffector", "SNARE-like membrane fusion", "same_complex", "physical_and_functional", "PMC11651115"),
    ("LegC3", "LegC7", "cooperative_paraeffector", "SNARE-like membrane fusion", "same_complex", "physical_and_functional", "PMC11651115"),
]

ALIASES = {
    "ANKX": "lpg0695", "LEM3": "lpg0696", "DRRA/SIDM": "lpg2464",
    "SIDD": "lpg2465", "LEPB": "lpg2490", "LGT1": "lpg1368",
    "LGT2": "lpg2862", "LGT3": "lpg1488", "SIDI": "lpg2504",
    "MESI/LPG2505": "lpg2505", "SDEA": "lpg2157", "SIDE": "lpg0234",
    "SIDJ": "lpg2155", "DUPA/LPG2154": "lpg2154", "DUPB/LPG2509": "lpg2509",
    "MAVC/LPG2147": "lpg2147", "MVCA/LPG2148": "lpg2148", "LPG2149": "lpg2149",
    "LEGC2": "lpg1884", "LEGC3": "lpg1701", "LEGC7": "lpg2298",
}


def main() -> int:
    OUT_META.mkdir(parents=True, exist_ok=True)
    OUT_INT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    if not UNIPROT.exists():
        raise FileNotFoundError(UNIPROT)

    up = pd.read_csv(UNIPROT, sep="\t", low_memory=False)
    up = up.rename(columns={"Entry": "uniprot_accession", "Gene Names": "gene_names", "Protein names": "protein_name", "Sequence": "sequence", "Length": "length"})
    up["locus_tag"] = up.gene_names.map(first_locus)
    up["primary_gene"] = up.gene_names.map(primary_gene)
    up["protein_sequence_sha256"] = up.sequence.map(sequence_hash)
    up["strain"] = "Philadelphia-1"
    up["taxid"] = 272624
    up["assembly_accession_version"] = "AE017354.1"
    up["mapping_method"] = "UniProt organism_id:272624 stream; locus parsed from Gene Names"
    up["ambiguity_status"] = up.locus_tag.map(lambda x: "RESOLVED" if x else "NO_LOCUS_TAG_IN_UNIPROT")

    relation_rows = []
    relevant_loci = set()
    for a, b, rel, pathway, direction, evidence, source in RELATIONS:
        la = ALIASES.get(a.upper(), "")
        lb = ALIASES.get(b.upper(), "")
        relevant_loci.update(x for x in (la, lb) if x)
        relation_rows.append({
            "effector_a": a,
            "effector_b": b,
            "locus_tag_a_philadelphia1": la,
            "locus_tag_b_philadelphia1": lb,
            "relationship_type": rel,
            "shared_host_target_or_pathway": pathway,
            "direction": direction,
            "evidence_system": evidence,
            "source_identifier": source,
            "evidence_strength": "CURATED_PUBLISHED",
            "candidate_rule": "null single-mutant result is inconclusive when redundancy/antagonism is relevant; relation itself is not host-dependency evidence",
        })
    rel_df = pd.DataFrame(relation_rows)
    rel_df.to_csv(OUT_INT / "effector_functional_relationships.tsv", sep="\t", index=False)

    core = up[up.locus_tag.isin(relevant_loci)].copy()
    core["common_effector_name"] = core.locus_tag.map({v: k for k, v in ALIASES.items()})
    core["orthogroup"] = ""
    core["identity"] = 1.0
    core["coverage"] = 1.0
    crosswalk_cols = [
        "strain", "taxid", "assembly_accession_version", "common_effector_name", "locus_tag",
        "uniprot_accession", "gene_names", "protein_name", "length", "protein_sequence_sha256",
        "orthogroup", "mapping_method", "identity", "coverage", "ambiguity_status",
    ]
    crosswalk = core[crosswalk_cols].sort_values(["locus_tag", "uniprot_accession"])

    unresolved = []
    for locus in sorted(relevant_loci - set(crosswalk.locus_tag)):
        unresolved.append({
            "strain": "Philadelphia-1", "taxid": 272624, "assembly_accession_version": "AE017354.1",
            "common_effector_name": {v: k for k, v in ALIASES.items()}.get(locus, ""), "locus_tag": locus,
            "uniprot_accession": "", "gene_names": "", "protein_name": "", "length": "",
            "protein_sequence_sha256": "", "orthogroup": "", "mapping_method": "curated locus; UniProt name unresolved",
            "identity": "", "coverage": "", "ambiguity_status": "NOT_FOUND_IN_CURRENT_UNIPROT_EXPORT",
        })
    if unresolved:
        crosswalk = pd.concat([crosswalk, pd.DataFrame(unresolved)], ignore_index=True)

    derivative_rows = []
    for strain, note in [
        ("Lp01", "Philadelphia-1 laboratory derivative; known lineage-specific mutations/deletions require strain-level caveat"),
        ("Lp02", "Philadelphia-1 rpsL hsdR thyA derivative; do not assume complete identity to reference"),
        ("JR32", "Philadelphia-1-derived laboratory lineage; do not assume complete identity to Lp01/Lp02"),
        ("AA100/130b", "non-Philadelphia lineage; independent assembly and orthology mapping required"),
    ]:
        derivative_rows.append({
            "strain": strain, "taxid": "", "assembly_accession_version": "", "common_effector_name": "ALL_CORE_EFFECTORS",
            "locus_tag": "", "uniprot_accession": "", "gene_names": "", "protein_name": "", "length": "",
            "protein_sequence_sha256": "", "orthogroup": "", "mapping_method": note, "identity": "", "coverage": "",
            "ambiguity_status": "REFERENCE_RELATION_REQUIRES_SEQUENCE_RESOLUTION",
        })
    crosswalk = pd.concat([crosswalk, pd.DataFrame(derivative_rows)], ignore_index=True)
    crosswalk.to_csv(OUT_META / "legionella_strain_locus_tag_crosswalk.tsv", sep="\t", index=False)

    dataset_strains = pd.DataFrame([
        {"dataset": "GSE61535", "reported_strain": "AA100", "reference_relation": "130b/AA100 non-Philadelphia lineage", "status": "SEPARATE_LINEAGE"},
        {"dataset": "GSE135803", "reported_strain": "not resolved in GEO series metadata", "reference_relation": "publication methods required", "status": "PENDING_PUBLICATION_AUDIT"},
        {"dataset": "GSE287755", "reported_strain": "not resolved in GEO series metadata", "reference_relation": "publication methods required", "status": "PENDING_PUBLICATION_AUDIT"},
        {"dataset": "GSE26490", "reported_strain": "Philadelphia-derived experimental strain expected", "reference_relation": "confirm publication strain table", "status": "PENDING_PUBLICATION_AUDIT"},
        {"dataset": "GSE89184", "reported_strain": "flaA/dotA/7-effector mutant backgrounds", "reference_relation": "confirm publication strain table", "status": "PENDING_PUBLICATION_AUDIT"},
        {"dataset": "GSE60085", "reported_strain": "L. pneumophila; exact laboratory lineage pending", "reference_relation": "confirm publication strain table", "status": "PENDING_PUBLICATION_AUDIT"},
        {"dataset": "PRJEB108343", "reported_strain": "LOCKED", "reference_relation": "do not inspect outcome/candidate results before freeze", "status": "HOLDOUT_LOCKED"},
    ])
    dataset_strains.to_csv(OUT_META / "dataset_strain_resolution.tsv", sep="\t", index=False)

    resolved_core = int(crosswalk.uniprot_accession.astype(str).ne("").sum())
    unresolved_core = int((crosswalk.ambiguity_status == "NOT_FOUND_IN_CURRENT_UNIPROT_EXPORT").sum())
    stats = {
        "uniprot_source_sha256": sha256(UNIPROT),
        "uniprot_records": int(len(up)),
        "curated_effector_relationships": int(len(rel_df)),
        "relationship_type_counts": rel_df.relationship_type.value_counts().to_dict(),
        "core_philadelphia_mappings_with_uniprot": resolved_core,
        "core_philadelphia_loci_unresolved_in_export": unresolved_core,
        "gate1b": "BORDERLINE",
        "gate1b_reason": "Philadelphia-1 core mapping is available; non-Philadelphia and laboratory-derivative sequence-level crosswalk remains unresolved.",
        "holdout_locked": True,
    }
    (OUT_META / "stage1b_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = [
        "# Stage 1B Strain and Effector Context Audit", "",
        f"- UniProt Philadelphia-1 records: {stats['uniprot_records']}",
        f"- Curated effector relationships: {stats['curated_effector_relationships']}",
        f"- Core mappings with UniProt accession: {resolved_core}",
        f"- Core curated loci unresolved in current export: {unresolved_core}",
        f"- Gate 1B: {stats['gate1b']}", "",
        "The borderline status is deliberate: Philadelphia-1 locus mapping is usable, but AA100/130b and laboratory-derivative equivalence must not be silently assumed.",
        "Effector redundancy/antagonism changes interpretation of null single-mutant contrasts to inconclusive; it does not convert them into positive evidence.",
        "PRJEB108343 remains locked.",
    ]
    (OUT_META / "STAGE1B_STRAIN_EFFECTOR_AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage1b.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
