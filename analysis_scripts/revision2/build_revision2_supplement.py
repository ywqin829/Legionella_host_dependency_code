"""Build the revision-2 Supplementary Tables 22-23 workbook.

The workbook is assembled only from archived WP4-WP6/PD1 outputs.  It does not
recalculate panels or access raw holdout data.  Failed opportunity-matching
audits are included as gate records only; no invalid endpoint estimate is
exported.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment, Font, PatternFill


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "revision2_2026-11-23"
OUT_XLSX = OUT_DIR / "Supplementary_Tables_22-23_revision2_working.xlsx"
OUT_MANIFEST = OUT_DIR / "REVISION2_SUPPLEMENT_22-23_MANIFEST.json"

SOURCES = {
    "ST22_RuleVariants": "REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv",
    "ST22_GeneRetention": "REVISION2_VARIANT_GENE_RETENTION.tsv",
    "ST22_GateAudit": "REVISION2_VARIANT_EFFECT_AUDIT.tsv",
    "ST23_PanelDefinitions": "REVISION2_WP5_COMPARATOR_SUMMARY.tsv",
    "ST23_PanelMembership": "REVISION2_WP5_COMPARATOR_PANELS.tsv.gz",
    "ST23_StrategySummary": "REVISION2_WP6_STRATEGY_SUMMARY.tsv",
    "ST23_MainContrasts": "REVISION2_WP6_MAIN_CONTRASTS.tsv",
    "ST23_EqualN": "REVISION2_WP6_EQUAL_N_SENSITIVITY.tsv",
    "ST23_DrawResults": "REVISION2_WP6_PANEL_RESULTS.tsv.gz",
    "ST23_PD1A_Summary": "REVISION2_PD1A_STRATEGY_SUMMARY.tsv",
    "ST23_PD1A_Contrasts": "REVISION2_PD1A_CONTRASTS.tsv",
    "ST23_PD1B_Gate": "REVISION2_PD1B_AUDIT.tsv",
    "ST23_PD1C_Gate": "REVISION2_PD1C_STAGE1_AUDIT.tsv",
    "ST23_G3Audit": "REVISION2_WP6_G3_AUDIT.tsv",
}

DESCRIPTIONS = {
    "ST22_RuleVariants": "Discovery-only rule perturbations, panel overlap, tier counts, and descriptive RAB/trafficking enrichment.",
    "ST22_GeneRetention": "Per-gene retention across five distinct local membership-changing panels and the structural domain series.",
    "ST22_GateAudit": "Machine-audited WP4/G2 checks; the transcriptional holdout was not accessed.",
    "ST23_PanelDefinitions": "Frozen definitions and panel-level properties for the seven comparator strategies, including expression-eligible sensitivities.",
    "ST23_PanelMembership": "Complete 56-gene membership of every frozen comparator panel and tie-resolution draw (4,003 panels).",
    "ST23_StrategySummary": "Same-holdout strategy summaries; stochastic strategies are summarized across 1,000 frozen tie draws.",
    "ST23_MainContrasts": "Frozen-minus-comparator contrasts and sensitivity intervals; calibration P values are not direct superiority tests.",
    "ST23_EqualN": "No-backfill equal-n=50 sensitivity analysis.",
    "ST23_DrawResults": "One row per evaluated frozen panel/draw, including evaluable n and direction-oriented statistic.",
    "ST23_PD1A_Summary": "Post-result endpoint-alignment diagnostic summaries; explanatory and not a validation benchmark.",
    "ST23_PD1A_Contrasts": "Post-result endpoint-alignment diagnostic contrasts.",
    "ST23_PD1B_Gate": "Opportunity-matching audit gate record. Balance failed; no adjusted endpoint result is valid or included.",
    "ST23_PD1C_Gate": "Alternative opportunity audit Stage-1 gate record. Optimizer tolerance failed before outcomes; Stage 2 was not run.",
    "ST23_G3Audit": "Machine-audited WP6/G3 one-time same-holdout evaluation checks.",
}

HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
HEADER_FONT = Font(color="FFFFFF", bold=True)
NOTE_FILL = PatternFill("solid", fgColor="D9EAF7")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def count_rows(path: Path) -> int:
    with open_text(path) as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def coerce(value: str):
    value = value.strip()
    if value == "" or value.lower() in {"nan", "na", "none"}:
        return None
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    try:
        if value.isdigit() or (value.startswith("-") and value[1:].isdigit()):
            return int(value)
        return float(value)
    except ValueError:
        return value


def styled_row(ws, values, *, header=False, note=False):
    row = []
    for value in values:
        cell = WriteOnlyCell(ws, value=value)
        if header:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        elif note:
            cell.fill = NOTE_FILL
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        row.append(cell)
    ws.append(row)


def add_notes_sheet(wb: Workbook, title: str, rows: list[tuple[str, str]]):
    ws = wb.create_sheet(title)
    ws.freeze_panes = "A2"
    styled_row(ws, ["item", "detail"], header=True)
    for item, detail in rows:
        styled_row(ws, [item, detail], note=True)


def add_tsv_sheet(wb: Workbook, sheet: str, path: Path) -> int:
    ws = wb.create_sheet(sheet)
    ws.freeze_panes = "A2"
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        styled_row(ws, header, header=True)
        n = 0
        for fields in reader:
            ws.append([coerce(value) for value in fields])
            n += 1
    ws.auto_filter.ref = f"A1:{chr(64 + min(len(header), 26))}{n + 1}" if len(header) <= 26 else None
    return n


def main() -> None:
    paths = {sheet: OUT_DIR / name for sheet, name in SOURCES.items()}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing source files: " + "; ".join(missing))

    source_meta = {
        sheet: {
            "path": str(path),
            "rows": count_rows(path),
            "sha256": sha256(path),
            "description": DESCRIPTIONS[sheet],
        }
        for sheet, path in paths.items()
    }

    wb = Workbook(write_only=True)
    contents = wb.create_sheet("Contents")
    styled_row(contents, ["table", "sheet", "rows", "purpose", "source_file", "source_sha256"], header=True)
    for sheet, meta in source_meta.items():
        table = "Supplementary Table 22" if sheet.startswith("ST22") else "Supplementary Table 23"
        styled_row(
            contents,
            [table, sheet, meta["rows"], meta["description"], Path(meta["path"]).name, meta["sha256"]],
        )

    add_notes_sheet(
        wb,
        "ST22_ReadMe",
        [
            ("Scope", "Discovery-only candidate-rule sensitivity analysis. The transcriptional holdout was not accessed."),
            ("Duplicate handling", "A7 and B1 generate the same panel and are counted once in gene-retention summaries."),
            ("Inference boundary", "Repeated family-enrichment tests reuse the discovery evidence and are descriptive sensitivity checks, not independent discoveries."),
            ("Background", "The atlas contains 26,761 genes. The curated trafficking list contains 276 genes, of which 275 occur in the atlas."),
        ],
    )

    for sheet in ["ST22_RuleVariants", "ST22_GeneRetention", "ST22_GateAudit"]:
        add_tsv_sheet(wb, sheet, paths[sheet])

    add_notes_sheet(
        wb,
        "ST23_ReadMe",
        [
            ("Scope", "Post-review comparator analysis using the same PMA-THP-1 transcriptional holdout already examined in the original study."),
            ("Interpretation", "These are secondary comparisons of directional transcriptional corroboration, not independent validation of functional dependency or proof of general ranking superiority."),
            ("Panel construction", "All strategies use a common 26,761-gene discovery universe and select 56 genes. Stochastic ties were resolved in 1,000 frozen draws before WP6 holdout evaluation."),
            ("Observed ordering", "The frozen panel exceeded the CRISPR-only and LCV-only comparators but was below the unweighted discovery-domain-count comparator on this endpoint."),
            ("Coverage", "The principal strategies had 51, 54, 55 or 56 holdout-evaluable genes; the equal-n=50 no-backfill analysis retained the ordering."),
            ("Endpoint alignment", "PD1A was performed only after the comparator result was known and is explanatory. It suggests endpoint alignment as a plausible contributor but does not establish causation."),
            ("Opportunity bias", "PD1B failed balance and PD1C failed its Stage-1 optimizer tolerance before outcome evaluation. No opportunity-adjusted endpoint claim is valid; annotation/study-intensity bias remains unresolved."),
            ("Prohibited interpretation", "Do not treat permutation calibration P values as direct tests that one ranking strategy is superior to another."),
        ],
    )

    for sheet in [
        "ST23_PanelDefinitions",
        "ST23_PanelMembership",
        "ST23_StrategySummary",
        "ST23_MainContrasts",
        "ST23_EqualN",
        "ST23_DrawResults",
        "ST23_PD1A_Summary",
        "ST23_PD1A_Contrasts",
        "ST23_PD1B_Gate",
        "ST23_PD1C_Gate",
        "ST23_G3Audit",
    ]:
        add_tsv_sheet(wb, sheet, paths[sheet])

    wb.save(OUT_XLSX)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "builder": str(Path(__file__).resolve()),
        "output": str(OUT_XLSX),
        "output_sha256": sha256(OUT_XLSX),
        "sources": source_meta,
        "guardrails": {
            "recalculated_panels": False,
            "raw_holdout_accessed": False,
            "pd1b_adjusted_endpoint_exported": False,
            "pd1c_stage2_exported": False,
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(OUT_XLSX), "sha256": manifest["output_sha256"], "manifest": str(OUT_MANIFEST)}, indent=2))


if __name__ == "__main__":
    main()
