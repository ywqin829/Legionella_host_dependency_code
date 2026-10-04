"""Add Supplementary Table 20 (dataset metadata/independence) and 21 (gene-set
definitions) to the workbook, and rename it to Supplementary_Tables_1-21.xlsx.

ST20 reuses 02_metadata/dataset_metadata_audit.tsv and dataset_independence_matrix.tsv.
ST21 records the RAB/trafficking definitions from scripts/trafficking_universe.py
and the GO annotation source version.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "FINAL_SUBMISSION_2026-10-04"
TAB = PKG / "03_tables"
SRC = TAB / "Supplementary_Tables_1-19.xlsx"
DST = TAB / "Supplementary_Tables_1-21.xlsx"

sys.path.insert(0, str(ROOT / "scripts"))
from trafficking_universe import GENE_TO_FAMILY, RAB_GTPASE, RAB_ASSOCIATED, TRAFFICKING  # noqa: E402

META = pd.read_csv(ROOT / "02_metadata" / "dataset_metadata_audit.tsv", sep="\t")
IND = pd.read_csv(ROOT / "02_metadata" / "dataset_independence_matrix.tsv", sep="\t")

bold = Font(bold=True)
wrap = Alignment(wrap_text=True, vertical="top")

wb = openpyxl.load_workbook(SRC)

# ---------------- ST20 ----------------
if "ST20_Datasets" in wb.sheetnames:
    del wb["ST20_Datasets"]
ws = wb.create_sheet("ST20_Datasets")
rows = []
rows.append(["Supplementary Table 20. Discovery dataset metadata and independence"])
rows.append(["Metadata for the six discovery expression datasets and the locked E-MTAB-16655 holdout, "
             "and the pairwise dataset-independence matrix. Independence was assessed by shared samples and "
             "shared PubMed records; all discovery pairs were sample-disjoint in this compilation."])
rows.append(["Canonical source: 02_metadata\\dataset_metadata_audit.tsv and dataset_independence_matrix.tsv"])
rows.append([None])
rows.append(["Panel A. Dataset metadata"])
rows.append(["accession", "role", "model_planned", "contrast_family", "heldout", "sample_count",
             "organisms", "pubmed_ids", "title"])
for _, r in META.iterrows():
    rows.append([r.get("accession"), r.get("role"), r.get("model_planned"), r.get("contrast_family"),
                 bool(r.get("heldout")) if pd.notna(r.get("heldout")) else None,
                 int(r["sample_count"]) if pd.notna(r.get("sample_count")) else None,
                 r.get("organisms"), (r.get("pubmed_ids") if pd.notna(r.get("pubmed_ids")) else ""),
                 r.get("title")])
rows.append([None])
rows.append(["Panel B. Pairwise dataset-independence matrix"])
rows.append(["dataset_a", "dataset_b", "independent_preliminary", "shared_samples", "shared_pmids"])
for _, r in IND.iterrows():
    rows.append([r.get("dataset_a"), r.get("dataset_b"), bool(r.get("independent_preliminary")),
                 r.get("shared_samples"), (r.get("shared_pmids") if pd.notna(r.get("shared_pmids")) else "")])
for i, r in enumerate(rows, start=1):
    for j, v in enumerate(r, start=1):
        c = ws.cell(row=i, column=j, value=v)
        if i == 1 and j == 1:
            c.font = bold
        if i in (2,):
            c.alignment = wrap
for col, w in zip("ABCDEFGHI", [14, 22, 30, 26, 8, 12, 16, 20, 60]):
    ws.column_dimensions[col].width = w

# ---------------- ST21 ----------------
if "ST21_GeneSetDefs" in wb.sheetnames:
    del wb["ST21_GeneSetDefs"]
ws = wb.create_sheet("ST21_GeneSetDefs")
fam_of = {f: sorted(g) for f, g in TRAFFICKING.items()}
rows = []
rows.append(["Supplementary Table 21. Gene-set and annotation definitions"])
rows.append(["The canonical RAB GTPase family, the RAB-associated machinery, and the curated "
             "membrane-trafficking category used for the panel enrichment, with the Gene Ontology "
             "annotation source used for the functional-module partition. The canonical RAB GTPase "
             "family (the RAB* genes themselves) is the primary enrichment set; the RAB-associated "
             "machinery additionally includes RAB GEF/effector and RAB GAP genes."])
rows.append(["Canonical source: scripts\\trafficking_universe.py; 01_data_raw\\annotations\\goa_human.gaf.gz"])
rows.append([None])
rows.append(["Panel A. Canonical RAB GTPase family (RAB* genes), n=%d" % len(RAB_GTPASE)])
rows.append(["gene", "family"])
for g in sorted(RAB_GTPASE):
    fam = GENE_TO_FAMILY.get(g, "?")
    rows.append([g, fam])
rows.append([None])
rows.append(["Panel B. RAB-associated machinery additions (RAB GEF/effector + RAB GAP), n=%d"
             % (len(RAB_ASSOCIATED) - len(RAB_GTPASE))])
rows.append(["gene", "family"])
for g in sorted(set(RAB_ASSOCIATED) - set(RAB_GTPASE)):
    fam = GENE_TO_FAMILY.get(g, "?")
    rows.append([g, fam])
rows.append([None])
rows.append(["Panel C. Curated membrane-trafficking category: %d families, %d genes (HGNC nomenclature)"
             % (len(TRAFFICKING), len(GENE_TO_FAMILY))])
rows.append(["family", "n_genes", "genes"])
for fam in TRAFFICKING:
    genes = fam_of[fam]
    rows.append([fam, len(genes), ", ".join(genes)])
rows.append([None])
rows.append(["Panel D. Annotation source for the functional-module partition"])
rows.append(["Item", "Value"])
rows.append(["Gene Ontology Annotation file", "goa_human.gaf.gz (GAF version 2.2)"])
rows.append(["Source database", "UniProt (Gene Ontology Annotation)"])
rows.append(["Source association file date", "2026-04-30"])
rows.append(["Retrieval / GOC date-generated", "2026-05-21"])
rows.append(["Use", "Functional-module (GO biological process most-specific term) partition for the correlation-robust analysis (Supplementary Table 15c)"])
for i, r in enumerate(rows, start=1):
    for j, v in enumerate(r, start=1):
        c = ws.cell(row=i, column=j, value=v)
        if i == 1 and j == 1:
            c.font = bold
        if i == 2:
            c.alignment = wrap
ws.column_dimensions["A"].width = 30
ws.column_dimensions["B"].width = 14
ws.column_dimensions["C"].width = 110
ws.column_dimensions["C"].alignment = wrap

# ---------------- Contents ----------------
contents = wb["Contents"]
# remove any prior ST20/ST21 rows
for row in list(contents.iter_rows(min_row=2)):
    for c in row:
        if isinstance(c.value, str) and (c.value == "ST20_Datasets" or c.value == "ST21_GeneSetDefs"):
            contents.delete_rows(c.row)
            break
last = contents.max_row + 1
contents.cell(row=last, column=1, value="Supplementary Table 20")
contents.cell(row=last, column=2, value="ST20_Datasets")
contents.cell(row=last, column=3, value="Discovery dataset metadata and pairwise independence matrix")
contents.cell(row=last, column=4, value="02_metadata\\dataset_metadata_audit.tsv; dataset_independence_matrix.tsv")
contents.cell(row=last + 1, column=1, value="Supplementary Table 21")
contents.cell(row=last + 1, column=2, value="ST21_GeneSetDefs")
contents.cell(row=last + 1, column=3, value="RAB family and membrane-trafficking gene-set definitions and GO annotation source")
contents.cell(row=last + 1, column=4, value="scripts\\trafficking_universe.py; 01_data_raw\\annotations\\goa_human.gaf.gz")

wb.save(DST)
print("wrote", DST.name)
print("sheets:", wb.sheetnames)
print("contents last rows:")
for r in wb["Contents"].iter_rows(min_row=wb["Contents"].max_row - 3, values_only=True):
    print(r[:3])
