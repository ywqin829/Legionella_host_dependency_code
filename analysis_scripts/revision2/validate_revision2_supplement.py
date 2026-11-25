"""Validate the assembled revision-2 Supplementary Tables 22-23 workbook."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "revision2_2026-11-23"
WORKBOOK = OUT_DIR / "Supplementary_Tables_22-23_revision2_working.xlsx"
MANIFEST = OUT_DIR / "REVISION2_SUPPLEMENT_22-23_MANIFEST.json"
REPORT = OUT_DIR / "REVISION2_SUPPLEMENT_22-23_VALIDATION.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_shape(path: Path) -> tuple[int, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader) if path.stat().st_size < 5_000_000 else reader
        if isinstance(rows, list):
            return len(rows), max((len(row) for row in rows), default=0)
        row_count = 0
        column_count = 0
        for row in rows:
            row_count += 1
            column_count = max(column_count, len(row))
        return row_count, column_count


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    workbook = load_workbook(WORKBOOK, read_only=True, data_only=True)
    sheet_shapes: dict[str, dict[str, object]] = {}
    for sheet_name in workbook.sheetnames:
        row_count = 0
        column_count = 0
        header = None
        for row in workbook[sheet_name].iter_rows(values_only=True):
            row_count += 1
            column_count = max(column_count, len(row))
            if header is None:
                header = list(row)
        sheet_shapes[sheet_name] = {
            "rows": row_count,
            "columns": column_count,
            "header": header,
        }

    shape_checks = {}
    for sheet_name, meta in manifest["sources"].items():
        source_path = Path(meta["path"])
        expected_rows, expected_columns = source_shape(source_path)
        observed = sheet_shapes[sheet_name]
        shape_checks[sheet_name] = {
            "expected_rows": expected_rows,
            "observed_rows": observed["rows"],
            "expected_columns": expected_columns,
            "observed_columns": observed["columns"],
            "pass": expected_rows == observed["rows"]
            and expected_columns == observed["columns"],
        }

    original = (
        ROOT
        / "FINAL_SUBMISSION_2026-10-04"
        / "03_tables"
        / "Supplementary_Tables_1-21.xlsx"
    )
    report = {
        "status": "PASS"
        if sha256(WORKBOOK) == manifest["output_sha256"]
        and all(item["pass"] for item in shape_checks.values())
        and all(value is False for value in manifest["guardrails"].values())
        else "FAIL",
        "workbook": str(WORKBOOK),
        "workbook_sha256": sha256(WORKBOOK),
        "manifest_sha256_match": sha256(WORKBOOK) == manifest["output_sha256"],
        "original_tables_1_21": str(original),
        "original_tables_1_21_sha256": sha256(original),
        "guardrails": manifest["guardrails"],
        "sheet_shapes": sheet_shapes,
        "source_shape_checks": shape_checks,
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
