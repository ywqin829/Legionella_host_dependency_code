"""Assemble the revision-2 submission figure set.

Only copies existing rendered files. Unchanged main figures (1-4, 6, 7) and unchanged
supplementary figures (1-6) are copied from the frozen 2026-10-04 submission; the
re-drawn Figure 5 and the new Supplementary Figures 7 and 8 come from the revision-2
figure workspace. Nothing is re-rendered here, so no panel content can drift.

Outputs (``revision2_2026-11-23/submission_figures/``):
    main/Figure1..Figure7 (PDF)
    supplementary/SupplementaryFigure1..8 (PDF)
    REVISION2_FIGURE_MANIFEST.json
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "FINAL_SUBMISSION_2026-10-04" / "02_figures_final"
REV2_FIG = ROOT / "revision2_2026-11-23" / "figures"
OUT = ROOT / "revision2_2026-11-23" / "submission_figures"

PLAN = {
    "main": [
        (1, FROZEN / "main" / "Figure1_data_architecture.pdf",
         "Figure1_data_architecture.pdf", "frozen 2026-10-04", "unchanged"),
        (2, FROZEN / "main" / "Figure2_cross_cohort_mechanism.pdf",
         "Figure2_cross_cohort_mechanism.pdf", "frozen 2026-10-04", "unchanged"),
        (3, FROZEN / "main" / "Figure3_prioritization.pdf",
         "Figure3_prioritization.pdf", "frozen 2026-10-04", "unchanged"),
        (4, FROZEN / "main" / "Figure4_holdout_validation.pdf",
         "Figure4_holdout_validation.pdf", "frozen 2026-10-04", "unchanged"),
        (5, REV2_FIG / "Figure5_evidence_coupling_main.pdf",
         "Figure5_evidence_coupling.pdf",
         "scripts/revision2/figures/figure5_split.R", "re-drawn: nine panels split into 5A-F + Sup Fig 8A-C"),
        (6, FROZEN / "main" / "Figure6_trafficking_atlas.pdf",
         "Figure6_trafficking_atlas.pdf", "frozen 2026-10-04", "unchanged"),
        (7, FROZEN / "main" / "Figure7_RAB1A_architecture.pdf",
         "Figure7_RAB1A_architecture.pdf", "frozen 2026-10-04", "unchanged"),
    ],
    "supplementary": [
        (1, FROZEN / "supplementary" / "SupplementaryFigure1_evidence_fragility.pdf",
         "SupplementaryFigure1_evidence_fragility.pdf", "frozen 2026-10-04", "unchanged"),
        (2, FROZEN / "supplementary" / "SupplementaryFigure2_effector_context.pdf",
         "SupplementaryFigure2_effector_context.pdf", "frozen 2026-10-04", "unchanged"),
        (3, FROZEN / "supplementary" / "SupplementaryFigure3_species_composition.pdf",
         "SupplementaryFigure3_species_composition.pdf", "frozen 2026-10-04", "unchanged"),
        (4, FROZEN / "supplementary" / "SupplementaryFigure4_depmap_coessentiality.pdf",
         "SupplementaryFigure4_depmap_coessentiality.pdf", "frozen 2026-10-04", "unchanged"),
        (5, FROZEN / "supplementary" / "SupplementaryFigure5_RAB_full_atlas.pdf",
         "SupplementaryFigure5_RAB_full_atlas.pdf", "frozen 2026-10-04", "unchanged"),
        (6, FROZEN / "supplementary" / "SupplementaryFigure6_ascertainment_sensitivity.pdf",
         "SupplementaryFigure6_ascertainment_sensitivity.pdf", "frozen 2026-10-04", "unchanged"),
        (7, REV2_FIG / "SupplementaryFigure7_rule_sensitivity.pdf",
         "SupplementaryFigure7_rule_sensitivity.pdf",
         "scripts/revision2/figures/suppfigure7_rule_sensitivity.R",
         "re-drawn to the A-D legend from archived WP4 tables"),
        (8, REV2_FIG / "SupplementaryFigure8_coupling_diagnostics.pdf",
         "SupplementaryFigure8_coupling_diagnostics.pdf",
         "scripts/revision2/figures/figure5_split.R",
         "new: correlation-robust diagnostics moved out of Figure 5"),
    ],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    missing = [str(src) for _, src, *_ in PLAN["main"] + PLAN["supplementary"] if not src.exists()]
    if missing:
        raise SystemExit("missing rendered figure(s):\n" + "\n".join(missing))

    entries = []
    for group in ("main", "supplementary"):
        target_dir = OUT / group
        target_dir.mkdir(parents=True, exist_ok=True)
        for number, src, name, origin, note in PLAN[group]:
            dst = target_dir / name
            shutil.copyfile(src, dst)
            entries.append({
                "group": group,
                "figure": number,
                "file": str(dst.relative_to(ROOT)).replace("\\", "/"),
                "source_file": str(src.relative_to(ROOT)).replace("\\", "/"),
                "origin": origin,
                "note": note,
                "bytes": dst.stat().st_size,
                "sha256": sha256(dst),
            })

    manifest = {
        "analysis": "WP8 revision-2 submission figure set",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_output_directory": str(OUT),
        "rendered_here": False,
        "main_figures": sum(1 for e in entries if e["group"] == "main"),
        "supplementary_figures": sum(1 for e in entries if e["group"] == "supplementary"),
        "entries": entries,
    }
    (OUT / "REVISION2_FIGURE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "main": manifest["main_figures"],
        "supplementary": manifest["supplementary_figures"],
        "directory": str(OUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
