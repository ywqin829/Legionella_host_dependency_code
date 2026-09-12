from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "05_expression" / "stage3_sample_design_inventory.tsv"
OUT = ROOT / "02_metadata"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return h


def integer_at_end(text: str) -> str:
    match = re.search(r"(?:rep|replicate)[ _-]*(\d+)\s*$", text, re.I)
    return match.group(1) if match else ""


def classify(row: pd.Series) -> dict[str, object]:
    ds, title, chars = row.dataset, row.title, row.characteristics
    low = f"{title} {chars}".lower()
    result: dict[str, object] = {
        "assay_modality": "ARRAY" if ds in {"GSE61535", "GSE26490", "GSE60085"} else "RNA",
        "design_unit": "BIOLOGICAL_SAMPLE",
        "infection_group": "UNRESOLVED",
        "host_genotype": "WT",
        "drug_treatment": "NONE",
        "biological_replicate": integer_at_end(title),
        "technical_replicate": "",
        "primary_eligible": True,
        "boundary_note": "",
    }
    if ds == "GSE61535":
        result["infection_group"] = "ANKB_MUTANT" if "ankb mutant" in low else "T2SS_MUTANT" if "t2ss mutant" in low else "UNINFECTED" if "uninfected" in low else "WT"
    elif ds == "GSE135803":
        result["infection_group"] = "UNINFECTED" if "control" in low else "ANKH_MUTANT" if "ankh" in low else "WT"
        result["biological_replicate"] = title.split("-")[-1].split("_")[0]
    elif ds == "GSE287755":
        result["infection_group"] = "UNINFECTED" if "uninfected" in low else "DOTA_MUTANT" if "deltadota" in low else "DELTA4_TRANSLATION_EFFECTOR" if "delta4" in low else "WT"
        result["biological_replicate"] = title.rsplit(" ", 1)[-1]
    elif ds == "GSE26490":
        result["infection_group"] = "UNINFECTED" if "uninfected" in low else "DELTA5_TRANSLATION_EFFECTOR" if "delta5" in low or "δ5" in low else "WT"
        bio = re.search(r"biological replicate (\d+)", low)
        tech = re.search(r"tech(?:ical|nical) replicate (\d+|[ab])", low)
        result["biological_replicate"] = bio.group(1) if bio else integer_at_end(title)
        result["technical_replicate"] = tech.group(1) if tech else ""
        result["boundary_note"] = "Collapse technical replicates before biological inference"
    elif ds == "GSE89184":
        result["assay_modality"] = "RNA" if "_rna" in title.lower() else "RPF"
        result["infection_group"] = "UNINFECTED" if "uninf" in low else "DOTA_FLAA_MUTANT" if "dotfla" in low else "SEVEN_EFFECTOR_FLAA_MUTANT" if "7fla" in low else "FLAA_MUTANT"
        result["host_genotype"] = "MYD88_KO" if "myd88" in low else "WT"
        result["drug_treatment"] = "HARRINGTONINE" if "harr" in low else "NONE"
        result["primary_eligible"] = result["drug_treatment"] == "NONE"
        result["boundary_note"] = "Mouse NAIP5/NLRC4 boundary; flaA background is the principal infection comparator"
    elif ds == "GSE60085":
        result["design_unit"] = "DUAL_COLOR_ARRAY_COMPARISON"
        result["infection_group"] = "ENCODED_IN_CHANNEL_PAIR"
        result["host_genotype"] = "IFNAR_KO" if "ifnar ko" in chars.lower() else "IFNGR_KO" if "ifngr ko" in chars.lower() else "IFN_DKO" if "ifn dko" in chars.lower() else "WT"
        result["boundary_note"] = "Dye-swap paired comparison; orient log-ratios before inference"
    return result


CONTRASTS = [
    ("A1", "GSE61535", "WT", "UNINFECTED", "infection", "PRIMARY"),
    ("A2", "GSE135803", "WT", "UNINFECTED", "infection", "PRIMARY"),
    ("A3", "GSE287755", "WT", "UNINFECTED", "infection", "PRIMARY"),
    ("A4", "GSE26490", "WT", "UNINFECTED", "infection", "PRIMARY"),
    ("A5", "GSE89184", "FLAA_MUTANT", "UNINFECTED", "infection_mouse_flaA", "PRIMARY_WITH_SPECIES_BOUNDARY"),
    ("A6", "GSE60085", "WT_LP", "WT_PBS", "infection_dual_color", "PRIMARY"),
    ("B1", "GSE287755", "WT", "DOTA_MUTANT", "Dot_Icm", "PRIMARY"),
    ("B2", "GSE89184", "FLAA_MUTANT", "DOTA_FLAA_MUTANT", "Dot_Icm_mouse_flaA", "PRIMARY_WITH_SPECIES_BOUNDARY"),
    ("C1", "GSE287755", "WT", "DELTA4_TRANSLATION_EFFECTOR", "translation_effector_module", "PRIMARY"),
    ("C2", "GSE26490", "WT", "DELTA5_TRANSLATION_EFFECTOR", "translation_effector_module", "PRIMARY"),
    ("C3", "GSE89184", "FLAA_MUTANT", "SEVEN_EFFECTOR_FLAA_MUTANT", "translation_effector_module", "PRIMARY_WITH_SPECIES_BOUNDARY"),
    ("D1", "GSE61535", "WT", "ANKB_MUTANT", "single_effector_AnkB", "SECONDARY_REDUNDANCY_AWARE"),
    ("D2", "GSE61535", "WT", "T2SS_MUTANT", "T2SS", "SECONDARY"),
    ("D3", "GSE135803", "WT", "ANKH_MUTANT", "single_effector_AnkH", "SECONDARY_REDUNDANCY_AWARE"),
    ("E1", "GSE89184", "RPF", "RNA", "translation_efficiency_interaction", "ORTHOGONAL_GLOBAL_SHUTOFF_AWARE"),
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(SOURCE, sep="\t", dtype=str).fillna("")
    annotations = pd.DataFrame([classify(row) for _, row in source.iterrows()])
    design = pd.concat([source.reset_index(drop=True), annotations], axis=1)
    design_path = OUT / "stage3_frozen_sample_design.tsv"
    design.to_csv(design_path, sep="\t", index=False)
    contrasts = pd.DataFrame(CONTRASTS, columns=["contrast_id", "dataset", "numerator", "denominator", "question", "evidence_role"])
    contrast_path = OUT / "stage3_frozen_contrasts.tsv"
    contrasts.to_csv(contrast_path, sep="\t", index=False)
    stats = {
        "input_sha256": sha256(SOURCE),
        "n_samples": int(len(design)),
        "n_datasets": int(design.dataset.nunique()),
        "modality_counts": design.assay_modality.value_counts().to_dict(),
        "n_contrasts": int(len(contrasts)),
        "primary_contrasts": int(contrasts.evidence_role.str.startswith("PRIMARY").sum()),
        "holdout_accessed": False,
        "design_sha256": sha256(design_path),
        "contrast_sha256": sha256(contrast_path),
        "gate3_design": "PASS" if len(design) == 104 and len(contrasts) >= 10 else "FAIL",
    }
    (OUT / "stage3_design_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = [
        "# Stage 3 Frozen Discovery Design", "",
        f"- Samples: {stats['n_samples']}", f"- Datasets: {stats['n_datasets']}",
        f"- Contrasts: {stats['n_contrasts']} ({stats['primary_contrasts']} primary)",
        f"- Modalities: {stats['modality_counts']}", f"- Gate: {stats['gate3_design']}", "",
        "Technical replicates are not biological replicates. GSE60085 is a dye-swap comparison design. Harringtonine samples are excluded from primary steady-state TE inference. Mouse flaA contrasts carry a pre-registered NAIP5/NLRC4 boundary flag.",
        "The external holdout was not accessed.",
    ]
    (OUT / "STAGE3_FROZEN_DESIGN.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage3_design.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0 if stats["gate3_design"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
