from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_data_raw" / "metadata"
OUT = ROOT / "05_expression"
CP = ROOT / "checkpoints"
ACCESSIONS = ["GSE61535", "GSE135803", "GSE287755", "GSE26490", "GSE89184", "GSE60085"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strip_value(line: str) -> str:
    return line.split("=", 1)[1].strip() if "=" in line else ""


def classify(url: str) -> str:
    value = url.lower()
    if any(x in value for x in ["count", "fpkm", "tpm", "matrix", "normalized", "expression", "qexpr", "tabular"]):
        return "PROCESSED_MATRIX_CANDIDATE"
    if value.endswith((".cel.gz", ".gpr.gz", ".txt.gz", ".csv.gz", ".tsv.gz")):
        return "SAMPLE_OR_ARRAY_FILE"
    if value.endswith((".tar", ".tar.gz", ".tgz")):
        return "ARCHIVE_REQUIRES_INSPECTION"
    if value.endswith((".bam", ".bam.gz", ".fastq.gz", ".fq.gz")):
        return "LARGE_RAW_OR_ALIGNED"
    return "OTHER"


def parse_soft(accession: str, path: Path) -> tuple[dict, list[dict], list[dict]]:
    series: dict[str, list[str]] = {}
    samples: list[dict] = []
    supplementary: list[dict] = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SERIES"):
                current = None
            elif line.startswith("^SAMPLE"):
                if current:
                    samples.append(current)
                current = {"dataset": accession, "gsm": strip_value(line), "title": "", "organism": "", "platform": "", "characteristics": []}
            elif line.startswith("^PLATFORM"):
                if current:
                    samples.append(current)
                    current = None
            elif line.startswith("!Series_"):
                key = line.split("=", 1)[0].replace("!Series_", "").strip()
                value = strip_value(line)
                series.setdefault(key, []).append(value)
                if key.startswith("supplementary_file") and value and value.upper() != "NONE":
                    supplementary.append({"dataset": accession, "scope": "SERIES", "gsm": "", "url": value, "classification": classify(value)})
            elif current is not None and line.startswith("!Sample_"):
                key = line.split("=", 1)[0].replace("!Sample_", "").strip()
                value = strip_value(line)
                if key == "title":
                    current["title"] = value
                elif key == "organism_ch1":
                    current["organism"] = value
                elif key == "platform_id":
                    current["platform"] = value
                elif key.startswith("characteristics_ch1"):
                    current["characteristics"].append(value)
                elif key.startswith("supplementary_file") and value and value.upper() != "NONE":
                    supplementary.append({"dataset": accession, "scope": "SAMPLE", "gsm": current["gsm"], "url": value, "classification": classify(value)})
        if current:
            samples.append(current)
    for sample in samples:
        sample["characteristics"] = " | ".join(sample["characteristics"])
    compact = {k: " | ".join(v) for k, v in series.items()}
    return compact, samples, supplementary


def infer_route(series: dict, supp: list[dict]) -> tuple[str, str]:
    exp = series.get("type", "").lower()
    classes = {x["classification"] for x in supp}
    if "ribosome profiling" in (series.get("summary", "") + series.get("overall_design", "")).lower():
        return "PAIRED_RNA_RPF", "Inspect processed per-sample counts first; raw Ribo-seq only if QC/interaction model cannot be supported"
    if "high throughput sequencing" in exp:
        if "PROCESSED_MATRIX_CANDIDATE" in classes:
            return "RNASEQ_PROCESSED_COUNTS_FIRST", "Validate integer/raw-count status; use DESeq2/edgeR only for counts"
        return "RNASEQ_RAW_REQUIRED", "Use ENA/SRA FASTQ with fixed reference and counting pipeline"
    if "array" in exp:
        return "MICROARRAY_RAW_OR_SERIES_MATRIX", "Prefer CEL/GPR re-normalization; otherwise audit normalized matrix and use limma"
    return "MANUAL_INSPECTION", "Assay type unresolved"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    series_rows, sample_rows, supp_rows, failures = [], [], [], []
    hashes = {}
    for accession in ACCESSIONS:
        path = RAW / f"{accession}_family.soft.gz"
        if not path.exists():
            failures.append({"dataset": accession, "reason": "SOFT_MISSING"})
            continue
        try:
            hashes[path.name] = sha256(path)
            series, samples, supp = parse_soft(accession, path)
            route, note = infer_route(series, supp)
            series_rows.append({
                "accession": accession,
                "title": series.get("title", ""),
                "organism": series.get("taxid_ch1", ""),
                "experiment_type": series.get("type", ""),
                "platforms": series.get("platform_id", ""),
                "n_samples": len(samples),
                "n_supplementary": len(supp),
                "analysis_route": route,
                "route_note": note,
                "pubmed_id": series.get("pubmed_id", ""),
            })
            sample_rows.extend(samples)
            supp_rows.extend(supp)
        except Exception as exc:
            failures.append({"dataset": accession, "reason": f"{type(exc).__name__}: {exc}"})

    series_df = pd.DataFrame(series_rows)
    sample_df = pd.DataFrame(sample_rows)
    supp_df = pd.DataFrame(supp_rows)
    failures_df = pd.DataFrame(failures, columns=["dataset", "reason"])
    series_df.to_csv(OUT / "stage3_dataset_routes.tsv", sep="\t", index=False)
    sample_df.to_csv(OUT / "stage3_sample_design_inventory.tsv", sep="\t", index=False)
    supp_df.to_csv(OUT / "stage3_supplementary_file_inventory.tsv", sep="\t", index=False)
    failures_df.to_csv(OUT / "stage3_inventory_failures.tsv", sep="\t", index=False)

    stats = {
        "soft_sha256": hashes,
        "datasets_expected": len(ACCESSIONS),
        "datasets_parsed": int(len(series_df)),
        "samples_parsed": int(len(sample_df)),
        "supplementary_links": int(len(supp_df)),
        "route_counts": series_df.analysis_route.value_counts().to_dict() if len(series_df) else {},
        "failures": failures,
        "holdout_accessed": False,
        "gate3_inventory": "PASS" if len(series_df) >= 4 else "FAIL",
    }
    (OUT / "stage3_inventory_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = [
        "# Stage 3 Discovery Data Inventory", "",
        f"- Datasets parsed: {stats['datasets_parsed']}/{stats['datasets_expected']}",
        f"- Samples parsed: {stats['samples_parsed']}",
        f"- Supplementary links: {stats['supplementary_links']}",
        f"- Analysis routes: {stats['route_counts']}",
        f"- Gate: {stats['gate3_inventory']}", "",
        "Processed matrices are preferred only after confirming whether values are raw counts, normalized abundance, or array intensities.",
        "PRJEB108343 was not opened or parsed by this stage.",
    ]
    (OUT / "STAGE3_DATA_INVENTORY.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage3_inventory.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0 if stats["gate3_inventory"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
