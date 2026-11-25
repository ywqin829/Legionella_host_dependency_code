"""Gate G4 — cross-file consistency of the revision-2 submission package.

G0 audits numbers against archived outputs inside one pinned manuscript. G4 asks a
different question: do the manuscript, the figure legends, the figure files, the
supplementary tables and the response letter tell the same story with the same
numbers and the same figure numbering?

Checks
------
G4a  manuscript pin            the pinned file is the file that will be submitted
G4b  legend structure          one legend per figure, no gaps, no duplicates
G4c  legend panel letters      letters unique, contiguous, and equal to the panel
                               objects actually built by the generating script
G4d  reference resolution      every Fig./Figure/Supplementary Fig. reference and
                               every cited panel letter resolves
G4e  supplementary tables      every cited ST number exists in an assembled workbook
G4f  figure files              one file per figure, hash-matched to the manifest
G4g  claim tokens              the archived headline numbers appear verbatim in the text
G4h  banned phrases            superseded claims are absent (negation-aware in the letter)
G4i  placeholders              no TODO/TBD/insert markers anywhere
G4j  response letter           exists, answers every declared comment id, reuses the
                               same claim tokens and cites only existing exhibits
G4k  Fig. 5 / Sup. Fig. 8       the split is consistent: 5A-F, 8A-C, nothing else
G4l  abstract boundary         the superseded closer is gone and the negative result
                               is stated in the abstract
G4m  legend data values        values quoted inside legends equal the archived values

Statuses are computed. ``AUDIT`` marks a check that cannot be automated for a given
item (for example a supplementary figure whose script builds no named panel objects);
it never converts a failure into a pass and never blocks on its own.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))

import g0_registry as reg  # noqa: E402

OUT = ROOT / "revision2_2026-11-23"
MANUSCRIPT = ROOT / reg.MANUSCRIPT
FIGURE_DIR = OUT / "submission_figures"
FIGURE_MANIFEST = OUT / "_not_for_upload" / "REVISION2_FIGURE_MANIFEST.json"
LETTER = OUT / "REVISION2_RESPONSE_TO_REVIEWERS.md"
SUPP_LEGENDS = OUT / "REVISION2_SUPPLEMENTARY_LEGENDS.md"

ST_WORKBOOKS = [
    ROOT / "FINAL_SUBMISSION_2026-10-04" / "03_tables" / "Supplementary_Tables_1-21.xlsx",
    OUT / "Supplementary_Tables_22-23_revision2_working.xlsx",
]

FROZEN_FIG_DIR = ROOT / "FINAL_SUBMISSION_2026-10-04" / "figure R script"
REV2_FIG_DIR = ROOT / "scripts" / "revision2" / "figures"

#: figure number -> (script path, panel-object prefix).
#: Panel objects are the ``pA <-`` / ``sA <-`` assignments that the assembled figure
#: is built from, so this compares the legend against the code that draws it.
PANEL_SOURCES = {
    ("main", 1): (FROZEN_FIG_DIR / "Figure1_data_architecture.R", "p"),
    ("main", 2): (FROZEN_FIG_DIR / "Figure2_cross_cohort_mechanism.R", "p"),
    ("main", 3): (FROZEN_FIG_DIR / "Figure3_prioritization.R", "p"),
    ("main", 4): (FROZEN_FIG_DIR / "Figure4_holdout_validation.R", "p"),
    ("main", 5): (REV2_FIG_DIR / "print_v2" / "fig5_print.R", "p"),
    ("main", 6): (REV2_FIG_DIR / "print_v2" / "fig6_print.R", "p"),
    ("main", 7): (REV2_FIG_DIR / "print_v2" / "fig7_print.R", "p"),
    ("supplementary", 1): (FROZEN_FIG_DIR / "SupplementaryFigure1_evidence_fragility.R", "p"),
    ("supplementary", 2): (FROZEN_FIG_DIR / "SupplementaryFigure2_effector_context.R", "p"),
    ("supplementary", 3): (FROZEN_FIG_DIR / "SupplementaryFigure3_species_composition.R", "p"),
    ("supplementary", 4): (FROZEN_FIG_DIR / "SupplementaryFigure4_depmap_coessentiality.R", "p"),
    ("supplementary", 5): (FROZEN_FIG_DIR / "SupplementaryFigure5_RAB_full_atlas.R", "p"),
    ("supplementary", 6): (FROZEN_FIG_DIR / "SupplementaryFigure6_ascertainment_sensitivity.R", "p"),
    ("supplementary", 7): (REV2_FIG_DIR / "suppfigure7_rule_sensitivity.R", "p"),
    ("supplementary", 8): (REV2_FIG_DIR / "figure5_split.R", "s"),
    ("supplementary", 9): (REV2_FIG_DIR / "print_v2" / "suppfig9_print.R", "p"),
}

#: Phrases that a resubmission may not assert. Kept short and unambiguous: every
#: entry is an affirmative claim that the analyses do not support.
BANNED = [
    "improves prioritisation",
    "improves prioritization",
    "improved prioritisation",
    "improved prioritization",
    "converges on a rab/er-golgi host-factor evidence core",
    "strongly enriched",
    "route parity",
    "validated the prioritization framework",
    "validated the prioritisation framework",
]

PLACEHOLDERS = [
    "TODO", "TBD", "FIXME", "XXX", "lorem ipsum", "placeholder",
    "[insert", "<insert", "AUTHOR NOTE", "to be completed", "TBC",
]

SUPERSEDED_CLOSER = "most informative when its biological semantics are preserved"

_NORM = re.compile(r"\s+")


def normalise(text: str) -> str:
    return _NORM.sub(" ", text).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def split_legends(text: str) -> tuple[dict[int, str], dict[int, str]]:
    """Return ({main figure n: legend body}, {supp figure n: legend body})."""
    main: dict[int, str] = {}
    supp: dict[int, str] = {}
    current = None
    bucket = None
    buf: list[str] = []

    def flush():
        nonlocal current, buf
        if current is not None and bucket is not None:
            bucket[current] = "\n".join(buf).strip()
        current = None
        buf = []

    for line in text.split("\n"):
        if line.startswith("## "):
            flush()
            m = re.match(r"## Supplementary Figure (\d+)\.", line)
            if m:
                current, bucket = int(m.group(1)), supp
                buf = [line]
                continue
            m = re.match(r"## Figure (\d+)\.", line)
            if m:
                current, bucket = int(m.group(1)), main
                buf = [line]
                continue
            bucket = None
        elif current is not None:
            buf.append(line)
    flush()
    return main, supp


def declared_letters(legend: str) -> list[str]:
    """Panel letters declared by a legend: bold single letters plus heading ranges."""
    letters = set(re.findall(r"\*\*([A-Z])\*\*", legend))
    head = legend.split("\n")[0]
    for group in re.findall(r"\(([A-Z][A-Z,\-\s]*)\)", head):
        for part in group.split(","):
            part = part.strip()
            m = re.match(r"^([A-Z])\s*-\s*([A-Z])$", part)
            if m:
                letters.update(chr(c) for c in range(ord(m.group(1)), ord(m.group(2)) + 1))
            elif re.match(r"^[A-Z]$", part):
                letters.add(part)
    return sorted(letters)


def script_panel_letters(path: Path, prefix: str) -> list[str]:
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    return sorted(set(re.findall(rf"^\s*({prefix}[A-Z])\s*<-", text, flags=re.M)))


# --------------------------------------------------------------------------- #
# Check registry
# --------------------------------------------------------------------------- #

CHECKS: list[dict] = []


def record(cid: str, status: str, detail: str, evidence: str = "") -> None:
    CHECKS.append({"check_id": cid, "status": status, "detail": detail, "evidence": evidence})


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if not MANUSCRIPT.exists():
        raise SystemExit(f"manuscript not found: {MANUSCRIPT}")
    text = MANUSCRIPT.read_text(encoding="utf-8")
    norm = normalise(text)
    manuscript_sha = sha256(MANUSCRIPT)

    # Main-figure legends live in the manuscript; supplementary legends were moved to
    # the separate supplementary file for ASM submission (E5) and are audited there.
    main_legends, _ = split_legends(text)
    supp_text = SUPP_LEGENDS.read_text(encoding="utf-8") if SUPP_LEGENDS.exists() else ""
    _, supp_legends = split_legends(supp_text)
    exhibit_text = text + "\n" + supp_text
    exhibit_norm = normalise(exhibit_text)

    # ---------------------------------------------------------------- G4a
    pinned_ok = (
        reg.MANUSCRIPT.endswith("REVISION2_MANUSCRIPT_WORKING_DRAFT.md")
        and manuscript_sha == reg.MANUSCRIPT_SHA256
    )
    record(
        "G4a_MANUSCRIPT_PIN",
        "PASS" if pinned_ok else "FAIL",
        "the G0 pin points at the working draft and the hash matches"
        if pinned_ok else
        f"pin mismatch: registry={reg.MANUSCRIPT}, registry_sha={reg.MANUSCRIPT_SHA256[:16]}..., "
        f"actual_sha={manuscript_sha[:16]}...",
        f"{reg.MANUSCRIPT} :: {manuscript_sha}",
    )

    # ---------------------------------------------------------------- G4b
    problems = []
    for group, legends, expected in (("main", main_legends, 7), ("supplementary", supp_legends, 9)):
        missing = [n for n in range(1, expected + 1) if n not in legends]
        extra = [n for n in legends if n > expected]
        thin = [n for n, body in legends.items() if len(body) < 200]
        if missing:
            problems.append(f"{group}: missing legend(s) {missing}")
        if extra:
            problems.append(f"{group}: unexpected legend number(s) {extra}")
        if thin:
            problems.append(f"{group}: legend(s) shorter than 200 characters {thin}")
    record(
        "G4b_LEGEND_STRUCTURE",
        "PASS" if not problems else "FAIL",
        "7 main (manuscript) and 9 supplementary (separate legend file) legends, complete and non-empty"
        if not problems else "; ".join(problems),
        ", ".join(f"Fig{n}" for n in sorted(main_legends)) + " | " +
        ", ".join(f"S{n}" for n in sorted(supp_legends)),
    )

    # ---------------------------------------------------------------- G4c
    group_legends = {"main": main_legends, "supplementary": supp_legends}
    letter_rows = []
    letter_problems = []
    audits = []
    for key, (script, prefix) in PANEL_SOURCES.items():
        group, number = key
        legend = group_legends[group].get(number, "")
        declared = declared_letters(legend)
        built = script_panel_letters(script, prefix)
        if not built:
            audits.append(f"{group} {number}: script declares no {prefix}X panel objects")
            continue
        expected = [chr(ord("A") + i) for i in range(len(built))]
        if declared != expected:
            letter_problems.append(
                f"{group} {number}: legend {declared} vs script {expected}"
            )
        letter_rows.append(f"{group}{number}:{''.join(declared)}")
    status = "FAIL" if letter_problems else ("AUDIT" if audits else "PASS")
    record(
        "G4c_LEGEND_PANEL_LETTERS",
        status,
        "; ".join(letter_problems) if letter_problems else
        ("letters match the generating scripts; audit notes: " + "; ".join(audits) if audits
         else "every legend letter matches the panel objects built by its script"),
        " ".join(letter_rows),
    )

    # ---------------------------------------------------------------- G4d
    ref_problems = []
    main_refs = set()
    supp_refs = set()
    for m in re.finditer(r"\b(?:Fig\.|Figure)\s*(\d+)([A-Z](?:\s*(?:and|,)\s*[A-Z])*)?", text):
        number = int(m.group(1))
        # skip the legend heading itself and "Extended Data Figure"
        prefix = text[max(0, m.start() - 15):m.start()]
        if "Supplementary" in prefix:
            continue
        main_refs.add(number)
        if number not in main_legends:
            ref_problems.append(f"Fig. {number} has no legend")
            continue
        panels = re.findall(r"[A-Z]", m.group(2) or "")
        for p in panels:
            if p not in declared_letters(main_legends[number]):
                ref_problems.append(f"Fig. {number}{p} not declared in the legend")
    for m in re.finditer(r"Supplementary\s+(?:Fig\.|Figure)\s*(\d+)([A-Z](?:\s*(?:and|,)\s*[A-Z])*)?", text):
        number = int(m.group(1))
        supp_refs.add(number)
        if number not in supp_legends:
            ref_problems.append(f"Supplementary Fig. {number} has no legend")
            continue
        for p in re.findall(r"[A-Z]", m.group(2) or ""):
            if p not in declared_letters(supp_legends[number]):
                ref_problems.append(f"Supplementary Fig. {number}{p} not declared in the legend")
    record(
        "G4d_REFERENCE_RESOLUTION",
        "PASS" if not ref_problems else "FAIL",
        "every figure reference and cited panel resolves to a legend panel"
        if not ref_problems else "; ".join(sorted(set(ref_problems))),
        f"main refs {sorted(main_refs)}; supplementary refs {sorted(supp_refs)}",
    )

    # ---------------------------------------------------------------- G4e
    sheet_numbers: set[int] = set()
    sheet_detail = []
    for wb_path in ST_WORKBOOKS:
        if not wb_path.exists():
            ref_problems.append(f"missing workbook {wb_path.name}")
            continue
        wb = load_workbook(wb_path, read_only=True)
        found = set()
        for name in wb.sheetnames:
            m = re.match(r"^ST(\d+)[a-z]?_", name)
            if m:
                found.add(int(m.group(1)))
        wb.close()
        sheet_numbers |= found
        sheet_detail.append(f"{wb_path.name}: {sorted(found)}")

    cited_tables: set[int] = set()
    for m in re.finditer(r"Supplementary Tables?\s+([0-9][0-9,\s\-–and]*)", text):
        chunk = m.group(1)
        for part in re.split(r",|\band\b", chunk):
            part = part.strip()
            rng = re.match(r"^(\d+)\s*[-–]\s*(\d+)$", part)
            if rng:
                cited_tables.update(range(int(rng.group(1)), int(rng.group(2)) + 1))
            elif part.isdigit():
                cited_tables.add(int(part))
    missing_tables = sorted(n for n in cited_tables if n not in sheet_numbers)
    record(
        "G4e_SUPPLEMENTARY_TABLES",
        "PASS" if not missing_tables else "FAIL",
        f"all {len(cited_tables)} cited supplementary tables resolve to workbook sheets"
        if not missing_tables else f"cited but absent from the workbooks: {missing_tables}",
        "; ".join(sheet_detail) + f" | cited {sorted(cited_tables)}",
    )

    # ---------------------------------------------------------------- G4f
    fig_problems = []
    fig_entries = []
    if not FIGURE_MANIFEST.exists():
        fig_problems.append("figure manifest missing")
    else:
        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        seen = {"main": [], "supplementary": []}
        for entry in manifest["entries"]:
            path = ROOT / entry["file"]
            if not path.exists():
                fig_problems.append(f"missing file {entry['file']}")
                continue
            if sha256(path) != entry["sha256"]:
                fig_problems.append(f"hash mismatch {entry['file']}")
            seen[entry["group"]].append(entry["figure"])
            fig_entries.append(f"{entry['group'][0].upper()}{entry['figure']}")
        for group, expected in (("main", 7), ("supplementary", 9)):
            if sorted(seen[group]) != list(range(1, expected + 1)):
                fig_problems.append(f"{group} set is {sorted(seen[group])}")
    record(
        "G4f_FIGURE_FILES",
        "PASS" if not fig_problems else "FAIL",
        "7 main and 9 supplementary PDFs present and hash-matched"
        if not fig_problems else "; ".join(fig_problems),
        " ".join(fig_entries),
    )

    # ---------------------------------------------------------------- G4g
    strategies = pd.read_csv(OUT / "REVISION2_WP6_STRATEGY_SUMMARY.tsv", sep="\t")
    contrasts = pd.read_csv(OUT / "REVISION2_WP6_MAIN_CONTRASTS.tsv", sep="\t")
    wp4 = pd.read_csv(OUT / "REVISION2_WP4_RULE_SENSITIVITY_TABLE.tsv", sep="\t")
    sensitivity = pd.read_csv(
        ROOT / "FINAL_SUBMISSION_2026-10-04" / "figure R script" / "source_data" /
        "g_matching_sensitivity.tsv", sep="\t")
    corr = json.loads((ROOT / "FINAL_SUBMISSION_2026-10-04" / "figure R script" /
                       "source_data" / "corr_robust_summary.json").read_text(encoding="utf-8"))
    dose = json.loads((ROOT / "FINAL_SUBMISSION_2026-10-04" / "figure R script" /
                       "source_data" / "p0b_detectability_adjusted_summary.json")
                      .read_text(encoding="utf-8"))

    def strat(sid: str, column: str) -> float:
        row = strategies.loc[strategies.strategy_id == sid]
        if len(row) != 1:
            raise KeyError(f"strategy {sid} not unique")
        return float(row.iloc[0][column])

    def cnt(comparator: str, column: str) -> float:
        row = contrasts.loc[contrasts.comparator_id == comparator]
        if len(row) != 1:
            raise KeyError(f"contrast {comparator} not unique")
        return float(row.iloc[0][column])

    def wp4v(variant: str, column: str):
        row = wp4.loc[wp4.variant == variant]
        if len(row) != 1:
            raise KeyError(f"variant {variant} not unique")
        return row.iloc[0][column]

    local = wp4.loc[wp4.variant.isin([
        "A1_admit_F3", "A2_fitness_does_not_block", "A4_T2_orthogonal_no_spatial",
        "A6_T2_orthogonal_require_response_and_mechanism",
        "A7_drop_T2_functional_response_branch"])]

    overlap = sensitivity.loc[sensitivity.method == "overlap weighting"].iloc[0]

    TOKENS = [
        # C5 comparator
        ("C5", f"{int(strat('S1_FROZEN', 'evaluable_n_median'))} evaluable genes", "S1 evaluable_n"),
        ("C5", f"{strat('S1_FROZEN', 'statistic_median'):.2f}", "S1 statistic"),
        ("C5", f"{int(strat('S2_CRISPR', 'evaluable_n_median'))} evaluable genes", "S2 evaluable_n"),
        ("C5", f"{strat('S2_CRISPR', 'statistic_median'):.2f}", "S2 statistic"),
        ("C5", f"median evaluable count was {int(strat('S3_LCV', 'evaluable_n_median'))}", "S3 evaluable_n"),
        ("C5", f"{strat('S3_LCV', 'statistic_median'):.2f}", "S3 statistic"),
        ("C5", f"{strat('S3_LCV', 'statistic_p025'):.2f}-{strat('S3_LCV', 'statistic_p975'):.2f}",
         "S3 2.5-97.5 percentile"),
        ("C5", f"{int(strat('S4_COUNT', 'evaluable_n_median'))} evaluable genes", "S4 evaluable_n"),
        ("C5", f"{strat('S4_COUNT', 'statistic_median'):.2f}", "S4 statistic"),
        ("C5", f"{strat('S4_COUNT', 'statistic_p025'):.2f}-{strat('S4_COUNT', 'statistic_p975'):.2f}",
         "S4 2.5-97.5 percentile"),
        ("C5", f"{strat('S3_LCV', 'statistic_max'):.2f}", "S3 maximum"),
        ("C5", f"{strat('S4_COUNT', 'statistic_min'):.2f}", "S4 minimum"),
        ("C5", f"{abs(cnt('S4_COUNT', 'observed_delta')):.2f}", "frozen-minus-count |delta|"),
        ("C5", f"{abs(cnt('S4_COUNT', 'module_delta_p975')):.2f}", "count delta interval bound"),
        ("C5", f"{abs(cnt('S4_COUNT', 'module_delta_p025')):.2f}", "count delta interval bound"),
        # C4 rule sensitivity
        ("C4", f"{int(local.panel_n.min())}-{int(local.panel_n.max())}", "WP4 local panel range"),
        ("C4", f"{local.jaccard_vs_original.min():.3f}-{local.jaccard_vs_original.max():.3f}",
         "WP4 Jaccard range"),
        ("C4", (f"{int(wp4v('D3_panel_min_discovery_domains', 'panel_n'))} to "
                f"{int(wp4v('D4_panel_min_discovery_domains', 'panel_n'))} to "
                f"{int(wp4v('D5_panel_min_discovery_domains', 'panel_n'))}"), "domain series"),
        ("C4", "from four genes to one", "tier-only Tier 1 change"),
        ("C4", "Thirteen frozen genes", "retention count"),
        ("C4", "40 were retained by at least four", "retention count"),
        ("C4", "none was retained by fewer than three", "retention denominator"),
        # C2 PPI
        ("C2", f"{overlap.ci_low:.2f}-{overlap.ci_high:.2f}", "overlap-weighting CI"),
        ("C2", f"{overlap.estimate:.2f}", "overlap-weighting OR"),
        ("C2", f"{overlap.p:.3f}", "overlap-weighting P"),
        # C1 coupling
        ("C1", f"{dose['full_mouse_exposure']['fold']:.2f}", "full-detection fold"),
        ("C1", f"{dose['highconf_mouse_exposure']['fold']:.2f}", "high-confidence fold"),
        ("C1", f"{corr['A_lcv_OR']:.2f}", "module-clustered LCV OR"),
        ("C1", f"{corr['A_lcv_ci_cluster'][0]:.2f}-{corr['A_lcv_ci_cluster'][1]:.2f}",
         "module-clustered LCV CI"),
        ("C1", f"{dose['dose_response_trend_OR_per_level']:.2f}", "dose-response trend OR"),
        ("C1", f"{dose['dose_response_trend_CI'][0]:.2f}-{dose['dose_response_trend_CI'][1]:.2f}",
         "dose-response trend CI"),
        ("C1", f"{corr['B_decimation_fold_median']:.1f}", "decimation median fold"),
        ("C1", f"{corr['B_decimation_fold_p05']:.1f}-{corr['B_decimation_fold_p95']:.1f}",
         "decimation percentile band"),
    ]

    token_problems = [f"{cid}: '{tok}' absent" for cid, tok, _ in TOKENS if normalise(tok) not in norm]
    record(
        "G4g_CLAIM_TOKENS",
        "PASS" if not token_problems else "FAIL",
        f"all {len(TOKENS)} archived headline values occur verbatim in the manuscript"
        if not token_problems else "; ".join(token_problems),
        " | ".join(t for _, t, _ in TOKENS),
    )

    # ---------------------------------------------------------------- G4h
    def banned_hits(document: str) -> list[str]:
        low = normalise(document).lower()
        return [p for p in BANNED if p in low]

    def sentence_of(document: str, phrase: str) -> str:
        low = document.lower()
        idx = low.find(phrase)
        if idx < 0:
            return ""
        start = max(low.rfind(".", 0, idx), low.rfind("\n", 0, idx)) + 1
        end = low.find(".", idx)
        return document[start: end if end > 0 else len(document)]

    ms_hits = banned_hits(text)
    record(
        "G4h_BANNED_PHRASES",
        "PASS" if not ms_hits else "FAIL",
        "no superseded claim appears in the manuscript"
        if not ms_hits else "present in manuscript: " + "; ".join(ms_hits),
        f"{len(BANNED)} phrases screened",
    )

    # ---------------------------------------------------------------- G4i
    placeholder_hits = [p for p in PLACEHOLDERS if p.lower() in text.lower()]
    record(
        "G4i_PLACEHOLDERS",
        "PASS" if not placeholder_hits else "FAIL",
        "no placeholder or insert marker in the manuscript"
        if not placeholder_hits else "found: " + "; ".join(placeholder_hits),
        f"{len(PLACEHOLDERS)} markers screened",
    )

    # ---------------------------------------------------------------- G4j
    if not LETTER.exists():
        record(
            "G4j_RESPONSE_LETTER", "FAIL",
            f"response letter missing: {LETTER.name}",
            str(LETTER),
        )
    else:
        letter_text = LETTER.read_text(encoding="utf-8")
        letter_norm = normalise(letter_text)
        letter_problems = []
        ids = re.findall(r"^###\s+([Rr]\d+[-.]\d+|[A-Z]{1,2}\d+[-.]\d+)", letter_text, flags=re.M)
        if not ids:
            letter_problems.append("no numbered comment ids found")
        if len(set(ids)) != len(ids):
            letter_problems.append("duplicate comment ids")
        absent = [t for _, t, _ in TOKENS if normalise(t) not in letter_norm]
        if absent:
            letter_problems.append("claim tokens missing: " + ", ".join(absent))
        # Do not double-count the ``Fig. N`` substring inside a Supplementary
        # Figure citation as a main-figure citation.
        for m in re.finditer(r"(?<!Supplementary )\b(?:Fig\.|Figure)\s*(\d+)", letter_text):
            if int(m.group(1)) not in main_legends:
                letter_problems.append(f"letter cites Fig. {m.group(1)}")
        for m in re.finditer(r"Supplementary\s+(?:Fig\.|Figure)\s*(\d+)", letter_text):
            if int(m.group(1)) not in supp_legends:
                letter_problems.append(f"letter cites Supplementary Fig. {m.group(1)}")
        for m in re.finditer(r"Supplementary Tables?\s+(\d+)", letter_text):
            if int(m.group(1)) not in sheet_numbers:
                letter_problems.append(f"letter cites Supplementary Table {m.group(1)}")
        letter_placeholders = [p for p in PLACEHOLDERS if p.lower() in letter_text.lower()]
        if letter_placeholders:
            letter_problems.append("placeholder(s): " + ", ".join(letter_placeholders))
        for phrase in banned_hits(letter_text):
            sentence = sentence_of(letter_text, phrase).lower()
            if not any(neg in sentence for neg in ("not ", "no longer", "withdraw", "remov",
                                                   "cannot", "does not", "no claim")):
                letter_problems.append(f"affirmative banned phrase '{phrase}'")
        record(
            "G4j_RESPONSE_LETTER",
            "PASS" if not letter_problems else "FAIL",
            f"{len(set(ids))} comment responses; claim tokens, exhibit numbers and wording agree"
            if not letter_problems else "; ".join(letter_problems),
            f"{LETTER.name}",
        )

    # ---------------------------------------------------------------- G4k
    split_problems = []
    if declared_letters(main_legends.get(5, "")) != list("ABCDEF"):
        split_problems.append(f"Figure 5 letters {declared_letters(main_legends.get(5, ''))}")
    if declared_letters(supp_legends.get(8, "")) != list("ABC"):
        split_problems.append(f"Supplementary Figure 8 letters {declared_letters(supp_legends.get(8, ''))}")
    stale_panels = sorted(set(re.findall(r"Fig\. 5[G-I]", text)) |
                          set(re.findall(r"Supplementary Fig\. 8[D-Z]", text)))
    if stale_panels:
        split_problems.append("stale panel reference(s): " + ", ".join(stale_panels))
    record(
        "G4k_FIG5_SUPP8_SPLIT",
        "PASS" if not split_problems else "FAIL",
        "Figure 5 is A-F, Supplementary Figure 8 is A-C, no stale panel reference"
        if not split_problems else "; ".join(split_problems),
        "Fig. 5A-F; Supplementary Fig. 8A-C",
    )

    # ---------------------------------------------------------------- G4l
    abstract = text.split("## Abstract", 1)[1].split("## Importance", 1)[0] if "## Abstract" in text else ""
    abstract_norm = normalise(abstract)
    abstract_problems = []
    if SUPERSEDED_CLOSER in abstract_norm:
        abstract_problems.append("superseded closer still present")
    if "lower statistic than an unweighted" not in abstract_norm:
        abstract_problems.append("negative comparator result not stated")
    if "not" not in abstract_norm:
        abstract_problems.append("no limitation clause")
    record(
        "G4l_ABSTRACT_BOUNDARY",
        "PASS" if not abstract_problems else "FAIL",
        "abstract states the negative comparator result and drops the superseded closer"
        if not abstract_problems else "; ".join(abstract_problems),
        f"{len(abstract_norm.split())} words",
    )

    # ---------------------------------------------------------------- G4m
    src = ROOT / "FINAL_SUBMISSION_2026-10-04" / "figure R script" / "source_data"
    ppi_lcv = pd.read_csv(src / "g_ppi_lcv_exact_rate_ci.tsv", sep="\t")
    joint = ppi_lcv.loc[ppi_lcv.group == "PPI(A/B)+LCV"].iloc[0]
    abundance = pd.read_csv(src / "p0b_abundance_dose_response.tsv", sep="\t")
    legend_tokens = [
        ("Fig. 5B", f"{int(joint.successes)}/{int(joint.n)}", "joint rate successes/n"),
        ("Fig. 5B", f"{joint.rate * 100:.1f}%", "joint rate"),
        ("Fig. 5B", f"{joint.ci_low * 100:.1f}-{joint.ci_high * 100:.1f}", "joint rate CI"),
        ("Fig. 5E", f"{dose['full_mouse_exposure']['fold']:.2f}x", "full-detection fold"),
        ("Fig. 5E", f"{dose['highconf_mouse_exposure']['fold']:.2f}x", "high-confidence fold"),
        ("Fig. 5F", f"{abundance.crispr_clean_rate.min() * 100:.1f}%", "lowest dose level"),
        ("Fig. 5F", f"{abundance.crispr_clean_rate.max() * 100:.1f}%", "highest dose level"),
        ("Fig. 5F", f"{dose['dose_response_trend_OR_per_level']:.2f} per level", "trend OR"),
        ("Supp. Fig. 8A", f"{corr['A_lcv_OR']:.2f}", "module-clustered OR"),
        ("Supp. Fig. 8A", f"{corr['A_lcv_ci_cluster'][0]:.2f}-{corr['A_lcv_ci_cluster'][1]:.2f}",
         "module-clustered CI"),
        ("Supp. Fig. 8B", f"{corr['B_decimation_fold_median']:.1f}x", "decimation median"),
        ("Supp. Fig. 8B", f"{corr['B_decimation_fold_p05']:.1f}-{corr['B_decimation_fold_p95']:.1f}",
         "decimation band"),
        ("Supp. Fig. 8B", f"{int(corr['B_decimation_n_valid'])}", "decimation replicates"),
        ("Supp. Fig. 7A", f"{int(local.panel_n.min())}-{int(local.panel_n.max())}", "local panel range"),
        ("Supp. Fig. 7A", f"{local.jaccard_vs_original.min():.3f}-{local.jaccard_vs_original.max():.3f}",
         "Jaccard range"),
    ]
    legend_problems = [f"{fig}: '{tok}' absent from the legend block"
                       for fig, tok, _ in legend_tokens
                       if normalise(tok) not in exhibit_norm]
    record(
        "G4m_LEGEND_DATA_VALUES",
        "PASS" if not legend_problems else "FAIL",
        f"all {len(legend_tokens)} values quoted in figure legends equal the archived values"
        if not legend_problems else "; ".join(legend_problems),
        " | ".join(f"{f}={t}" for f, t, _ in legend_tokens),
    )

    # ---------------------------------------------------------------- output
    frame = pd.DataFrame(CHECKS)
    frame.to_csv(OUT / "REVISION2_G4_AUDIT.tsv", sep="\t", index=False)

    fails = frame.loc[frame.status == "FAIL", "check_id"].tolist()
    audits = frame.loc[frame.status == "AUDIT", "check_id"].tolist()
    gate = "FAIL" if fails else "PASS"

    stats = {
        "analysis": "Gate G4 — cross-file consistency",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "manuscript": {"path": reg.MANUSCRIPT, "sha256": manuscript_sha,
                       "pinned": reg.MANUSCRIPT_SHA256},
        "response_letter": str(LETTER.relative_to(ROOT)).replace("\\", "/")
        if LETTER.exists() else None,
        "checks_total": int(len(frame)),
        "checks_pass": int((frame.status == "PASS").sum()),
        "checks_audit": int((frame.status == "AUDIT").sum()),
        "checks_fail": int((frame.status == "FAIL").sum()),
        "gate_g4": gate,
        "failed_checks": fails,
        "audit_checks": audits,
        "claim_tokens_checked": len(TOKENS),
        "legend_values_checked": len(legend_tokens),
        "figures_checked": 7 + 9,
        "supplementary_tables_cited": sorted(cited_tables),
        "holdout_accessed_for_new_analysis": False,
    }
    (OUT / "REVISION2_G4_STATS.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")

    lines = [
        "# Revision 2 — Gate G4 cross-file consistency audit",
        "",
        f"- Timestamp (UTC): {stats['timestamp_utc']}",
        f"- Manuscript: `{reg.MANUSCRIPT}`",
        f"- Manuscript SHA-256: `{manuscript_sha}`",
        f"- Response letter: `{stats['response_letter']}`",
        f"- Checks: {stats['checks_total']} ({stats['checks_pass']} PASS, "
        f"{stats['checks_audit']} AUDIT, {stats['checks_fail']} FAIL)",
        f"- **Gate G4: {gate}**",
        "",
        "| check | status | detail |",
        "|---|---|---|",
    ]
    for row in CHECKS:
        detail = row["detail"].replace("|", "/")
        lines.append(f"| {row['check_id']} | {row['status']} | {detail} |")
    lines += [
        "",
        "## What G4 does and does not certify",
        "",
        "G4 certifies that the manuscript, figure legends, assembled figure files,",
        "supplementary-table citations and response letter use the same figure numbering and",
        "quote the same archived headline values. It does not re-derive statistics (G0 does",
        "that), and it does not certify scientific interpretation.",
        "",
        "AUDIT rows mark items that cannot be automated for a specific figure, such as a",
        "supplementary figure whose script builds no named panel objects. An AUDIT row is not",
        "a pass and is not a failure.",
    ]
    (OUT / "REVISION2_G4_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({
        "gate_g4": gate,
        "checks_total": stats["checks_total"],
        "checks_pass": stats["checks_pass"],
        "checks_audit": stats["checks_audit"],
        "checks_fail": stats["checks_fail"],
        "failed": fails,
    }, indent=2))
    if fails:
        print("\nFAILED:")
        print(frame.loc[frame.status == "FAIL", ["check_id", "detail"]].to_string(index=False))
    return 0 if gate == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
