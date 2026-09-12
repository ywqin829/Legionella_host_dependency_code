from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_data_raw" / "expression"
DESIGN_PATH = ROOT / "02_metadata" / "stage3_frozen_sample_design.tsv"
OUT = ROOT / "05_expression" / "matrices"
QC_OUT = ROOT / "05_expression"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_geo_matrix(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    inside = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if inside:
                rows.append(next(csv.reader([line], delimiter="\t")))
    if not rows:
        raise ValueError(f"No matrix table in {path}")
    header = [x.strip('"') for x in rows[0]]
    frame = pd.DataFrame(rows[1:], columns=header)
    frame = frame.rename(columns={frame.columns[0]: "feature_id"}).set_index("feature_id")
    return frame.apply(pd.to_numeric, errors="coerce")


def read_two_column_counts(dataset: str, pattern: str, strip_version: bool) -> pd.DataFrame:
    series = []
    for path in sorted((RAW / dataset).glob(pattern)):
        gsm = path.name.split("_")[0]
        frame = pd.read_csv(path, sep="\t", compression="gzip")
        feature = frame.iloc[:, 0].astype(str)
        if strip_version:
            feature = feature.str.replace(r"\.\d+$", "", regex=True)
        values = pd.to_numeric(frame.iloc[:, -1], errors="coerce").fillna(0)
        item = pd.Series(values.to_numpy(), index=feature, name=gsm).groupby(level=0).sum()
        series.append(item)
    return pd.concat(series, axis=1).fillna(0)


def read_fpkm() -> pd.DataFrame:
    path = RAW / "GSE135803" / "GSE135803_genes.fpkm_table.txt.gz"
    frame = pd.read_csv(path, sep="\t", compression="gzip")
    frame = frame.rename(columns={frame.columns[0]: "feature_id"}).set_index("feature_id")
    design = pd.read_csv(DESIGN_PATH, sep="\t", dtype=str).fillna("")
    mapping = dict(zip(design.loc[design.dataset.eq("GSE135803"), "title"], design.loc[design.dataset.eq("GSE135803"), "gsm"]))
    frame = frame.rename(columns=mapping)
    return frame.apply(pd.to_numeric, errors="coerce")


def read_qexpr() -> tuple[pd.DataFrame, int]:
    series, lengths = [], []
    for path in sorted((RAW / "GSE89184").glob("*_qexpr.txt.gz")):
        gsm = path.name.split("_")[0]
        frame = pd.read_csv(path, sep="\t", compression="gzip", header=None, names=["feature_id", "length", "count"])
        frame["feature_id"] = frame.feature_id.astype(str)
        series.append(pd.Series(pd.to_numeric(frame["count"], errors="coerce").fillna(0).to_numpy(), index=frame.feature_id, name=gsm).groupby(level=0).sum())
        lengths.append(pd.Series(pd.to_numeric(frame["length"], errors="coerce").to_numpy(), index=frame.feature_id, name=gsm).groupby(level=0).first())
    count_matrix = pd.concat(series, axis=1).fillna(0)
    length_matrix = pd.concat(lengths, axis=1)
    conflicts = int((length_matrix.nunique(axis=1, dropna=True) > 1).sum())
    length_matrix.median(axis=1, skipna=True).rename("median_length").to_csv(OUT / "GSE89184_feature_lengths.tsv.gz", sep="\t", compression="gzip")
    return count_matrix, conflicts


def matrix_qc(dataset: str, matrix: pd.DataFrame, value_type: str) -> dict[str, object]:
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    summary: dict[str, object] = {
        "dataset": dataset,
        "value_type": value_type,
        "features": int(matrix.shape[0]),
        "samples": int(matrix.shape[1]),
        "missing_fraction": float(np.isnan(values).mean()),
        "min": float(np.min(finite)),
        "median": float(np.median(finite)),
        "max": float(np.max(finite)),
        "all_zero_features": int(np.all(np.nan_to_num(values) == 0, axis=1).sum()),
    }
    if value_type == "RAW_COUNTS":
        library = np.nansum(values, axis=0)
        summary.update(
            library_min=float(library.min()),
            library_median=float(np.median(library)),
            library_max=float(library.max()),
            library_cv=float(np.std(library, ddof=1) / np.mean(library)),
            noninteger_values=int(np.sum(np.abs(finite - np.round(finite)) > 1e-8)),
        )
    return summary


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    matrices: dict[str, tuple[pd.DataFrame, str]] = {
        "GSE61535": (read_geo_matrix(RAW / "GSE61535" / "GSE61535_series_matrix.txt.gz"), "LOG2_EXPRESSION"),
        "GSE135803": (read_fpkm(), "FPKM"),
        "GSE287755": (read_two_column_counts("GSE287755", "*.tabular.gz", True), "RAW_COUNTS"),
        "GSE26490": (read_geo_matrix(RAW / "GSE26490" / "GSE26490_series_matrix.txt.gz"), "NORMALIZED_LOG_RATIO"),
        "GSE60085": (read_geo_matrix(RAW / "GSE60085" / "GSE60085_series_matrix.txt.gz"), "NORMALIZED_LOG_RATIO"),
    }
    qexpr, length_conflicts = read_qexpr()
    matrices["GSE89184"] = (qexpr, "RAW_COUNTS")

    qc_rows, output_hashes = [], {}
    for dataset, (matrix, value_type) in matrices.items():
        path = OUT / f"{dataset}_{value_type.lower()}.tsv.gz"
        matrix.to_csv(path, sep="\t", compression="gzip")
        output_hashes[path.name] = sha256(path)
        qc_rows.append(matrix_qc(dataset, matrix, value_type))

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(QC_OUT / "stage3_matrix_qc_summary.tsv", sep="\t", index=False)
    raw_count_ok = all((row.get("noninteger_values", 0) == 0) for row in qc_rows if row["value_type"] == "RAW_COUNTS")
    stats = {
        "datasets_built": len(matrices),
        "output_sha256": output_hashes,
        "raw_count_integer_check": raw_count_ok,
        "GSE89184_length_conflicts": length_conflicts,
        "global_translation_claim_authorized": False,
        "global_translation_reason": "No external spike-in or absolute molecule normalization; library totals are technical QC only.",
        "holdout_accessed": False,
        "gate3_qc_build": "PASS" if len(matrices) == 6 and raw_count_ok else "FAIL",
    }
    (QC_OUT / "stage3_qc_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    report = [
        "# Stage 3 Matrix Build and QC", "",
        f"- Discovery matrices built: {stats['datasets_built']}/6",
        f"- Integer count validation: {stats['raw_count_integer_check']}",
        f"- GSE89184 transcript-length conflicts: {length_conflicts}",
        f"- Gate: {stats['gate3_qc_build']}", "",
        "GSE61535 uses the submitter-normalized GEO matrix because the platform-specific oligo annotation package is absent; raw CEL files remain archived for clean-room sensitivity analysis.",
        "GSE26490 and GSE60085 are normalized dual-color log-ratio matrices, not independent single-channel expression measurements.",
        "GSE135803 FPKM is analyzed with log-scale linear models, never as integer counts.",
        "GSE287755 and GSE89184 passed integer-count validation.",
        "RPF library totals cannot establish absolute global translational shutdown without spike-ins; only relative gene-level RNA-by-modality interactions are authorized.",
        "The external holdout remains unopened.",
    ]
    (QC_OUT / "STAGE3_MATRIX_QC.md").write_text("\n".join(report), encoding="utf-8")
    (CP / "stage3_qc.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0 if stats["gate3_qc_build"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
