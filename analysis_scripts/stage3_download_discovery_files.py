from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "05_expression" / "stage3_supplementary_file_inventory.tsv"
RAW = ROOT / "01_data_raw" / "expression"
OUT = ROOT / "05_expression"
CP = ROOT / "checkpoints"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def selected(row: dict[str, str]) -> bool:
    url = row["url"]
    name = Path(urlparse(url).path).name
    if row["classification"] == "ARCHIVE_REQUIRES_INSPECTION":
        return False
    if row["dataset"] == "GSE89184" and "_chrM_qexpr" in name:
        return False
    return row["classification"] in {
        "PROCESSED_MATRIX_CANDIDATE",
        "SAMPLE_OR_ARRAY_FILE",
    } or name == "GSE26490_Figure3.xls.gz"


def download(url: str, target: Path, retries: int = 2) -> tuple[str, int]:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        return "EXISTING", target.stat().st_size
    https_url = url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
    part = target.with_suffix(target.suffix + ".part")
    last_error = ""
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(https_url, headers={"User-Agent": "Legionella-atlas/2.0"})
            with urllib.request.urlopen(request, timeout=120) as response, part.open("wb") as output:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    output.write(block)
            part.replace(target)
            return "DOWNLOADED", target.stat().st_size
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if part.exists():
                part.unlink()
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(last_error)


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    CP.mkdir(parents=True, exist_ok=True)
    with INVENTORY.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if selected(row)]

    manifest: list[dict[str, object]] = []
    for index, row in enumerate(rows, start=1):
        name = Path(urlparse(row["url"]).path).name
        target = RAW / row["dataset"] / name
        record: dict[str, object] = {
            "dataset": row["dataset"], "scope": row["scope"], "gsm": row["gsm"],
            "url": row["url"], "local_path": str(target),
        }
        try:
            status, size = download(row["url"], target)
            record.update(status=status, bytes=size, sha256=sha256(target), error="")
        except Exception as exc:
            record.update(status="FAILED", bytes=0, sha256="", error=str(exc))
        manifest.append(record)
        print(f"[{index}/{len(rows)}] {row['dataset']} {name}: {record['status']}", flush=True)

    fields = ["dataset", "scope", "gsm", "url", "local_path", "status", "bytes", "sha256", "error"]
    with (OUT / "stage3_download_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)

    failures = [x for x in manifest if x["status"] == "FAILED"]
    stats = {
        "selected_files": len(rows),
        "available_files": len(rows) - len(failures),
        "failed_files": len(failures),
        "total_bytes": sum(int(x["bytes"]) for x in manifest),
        "holdout_accessed": False,
        "gate3_download": "PASS" if not failures else "COMPLETE_WITH_WARNINGS",
    }
    (OUT / "stage3_download_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    (CP / "stage3_download.done").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
