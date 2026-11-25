"""Stage the revision-2 release payload for GitHub and Zenodo.

Builds, under ``revision2_2026-11-23/RELEASE_STAGING/``:

* ``github/analysis_scripts/revision2/`` — the revision-2 analysis and figure
  scripts, ready to drop into ``GITHUB_CODE_PUBLIC_2026-09-12``;
* ``data/`` and ``Legionella_host_dependency_revision2_data.zip`` — revision-2
  outputs, figure source data and Supplementary Tables 22-23;
* ``MANIFEST_SHA256.txt`` — checksums of every staged file.

Internal work-in-progress documents are deliberately excluded.
"""

from __future__ import annotations

import hashlib
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revision2_2026-11-23"
STAGE = OUT / "RELEASE_STAGING"
SCRIPTS = ROOT / "scripts" / "revision2"

EXCLUDE_DOCS = {
    # internal work-in-progress documents
    "REVISION2_HANDOFF_WORKLOG.md",
    "REVISION2_ADDITIONAL_AUDIT_AND_EXECUTION_PLAN.md",
    "REVISION2_PRESUBMISSION_THREE_ITEM_CHECK.md",
    "REVISION2_REVIEWER_RESPONSE_STATUS_AND_DECISION_BRIEF.md",
    "REVISION2_MANUSCRIPT_CHANGE_MAP.md",
    # submission-facing correspondence and manuscript text are distributed through
    # the journal, not the public data archive; reviewer correspondence is not public
    "REVISION2_COVER_LETTER.md",
    "REVISION2_RESPONSE_TO_REVIEWERS.md",
    "REVISION2_MANUSCRIPT_WORKING_DRAFT.md",
    "REVISION2_SUPPLEMENTARY_LEGENDS.md",
    "REVISION2_G4_AUDIT.md",
    "REVISION2_G4_AUDIT.tsv",
    "REVISION2_G4_STATS.json",
    "REVISION2_UPLOAD_CHECKLIST.md",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def copy_code() -> list[Path]:
    target = STAGE / "github" / "analysis_scripts" / "revision2"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(SCRIPTS, target, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    return sorted(p for p in target.rglob("*") if p.is_file())


def copy_data() -> list[Path]:
    data = STAGE / "data"
    if data.exists():
        shutil.rmtree(data)
    data.mkdir(parents=True)
    copied: list[Path] = []

    for path in sorted(OUT.iterdir()):
        if not path.is_file():
            continue
        if path.name in EXCLUDE_DOCS or path.name.startswith("_render"):
            continue
        if path.name == "Supplementary_Tables_22-23_revision2_working.xlsx":
            dest = data / "Supplementary_Tables_22-23.xlsx"
            shutil.copy2(path, dest)
            copied.append(dest)
            continue
        if path.suffix.lower() in {".md", ".tsv", ".gz", ".json"}:
            dest = data / path.name
            shutil.copy2(path, dest)
            copied.append(dest)

    for old_fig in sorted((OUT / "figures").glob("*_source_data.tsv")):
        dest = data / "figure_source_data" / old_fig.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(old_fig, dest)
        copied.append(dest)

    print_dir = OUT / "figures_v2_print"
    for print_fig in sorted(print_dir.iterdir()):
        if print_fig.is_file() and print_fig.suffix.lower() in {".pdf", ".png", ".md"}:
            dest = data / "figures_print" / print_fig.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(print_fig, dest)
            copied.append(dest)

    return copied


def write_manifest(files: list[Path]) -> Path:
    manifest = STAGE / "MANIFEST_SHA256.txt"
    lines = [f"{sha256(p)}  {p.relative_to(STAGE).as_posix()}" for p in files]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def zip_data(files: list[Path]) -> Path:
    archive = STAGE / "Legionella_host_dependency_revision2_data.zip"
    if archive.exists():
        archive.unlink()
    inner_manifest = STAGE / "MANIFEST_DATA_SHA256.txt"
    inner_manifest.write_text(
        "\n".join(f"{sha256(p)}  {p.relative_to(STAGE / 'data').as_posix()}"
                  for p in files) + "\n",
        encoding="utf-8",
    )
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            zf.write(path, (Path("Legionella_host_dependency_revision2_data") /
                            path.relative_to(STAGE / "data")).as_posix())
        zf.write(inner_manifest,
                 "Legionella_host_dependency_revision2_data/MANIFEST_SHA256.txt")
    inner_manifest.unlink()
    return archive


def main() -> None:
    STAGE.mkdir(parents=True, exist_ok=True)
    code_files = copy_code()
    data_files = copy_data()
    manifest = write_manifest(code_files + data_files)
    archive = zip_data(data_files)
    print(f"code files: {len(code_files)}")
    print(f"data files: {len(data_files)}")
    print(f"manifest: {manifest}")
    print(f"archive: {archive} ({archive.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
