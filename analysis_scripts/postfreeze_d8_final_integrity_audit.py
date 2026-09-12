from __future__ import annotations

import hashlib
import json
import py_compile
from pathlib import Path

import fitz
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
POST = ROOT / "13_postfreeze_interpretation"
REPORT = ROOT / "12_reports"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check(name: str, passed: bool, detail: str) -> dict:
    return {"check": name, "status": "PASS" if passed else "FAIL", "detail": detail}


def main() -> int:
    checks: list[dict] = []

    s5 = json.loads((ROOT / "07_convergence/stage5_expression_statistics.json").read_text(encoding="utf-8-sig"))
    checks.append(check("dataset-collapsed meta-convergent count", s5["stable_cross_cohort"] == 2601,
                        f"observed={s5['stable_cross_cohort']} expected=2601"))

    evidence = pd.read_csv(ROOT / "07_convergence/primary_infection_gene_study_evidence.tsv.gz", sep="\t", compression=None, low_memory=False)
    dup = int(evidence.duplicated(["human_symbol", "dataset"]).sum())
    checks.append(check("one human gene x dataset unit", dup == 0, f"duplicate_units={dup}"))

    legacy = ROOT / "07_convergence/legacy_pre_dataset_collapse/primary_infection_convergence.tsv"
    checks.append(check("pre-correction result archived", legacy.exists() and legacy.stat().st_size > 0, str(legacy)))

    tier = pd.read_csv(ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv", sep="\t")
    tier1 = tier.loc[tier.final_tier.eq("TIER1_EXTERNALLY_VALIDATED"), "human_symbol"].tolist()
    checks.append(check("sole Tier-1 identity", tier1 == ["RAB1A"], f"tier1={tier1}"))
    checks.append(check("frozen candidate count", len(tier) == 56, f"n={len(tier)}"))

    d7 = json.loads((POST / "d7_dataset_collapse_audit_summary.json").read_text(encoding="utf-8"))
    checks.append(check("frozen candidates unchanged by correction",
                        d7["frozen_candidate_expression_flag_changed_n"] == 0,
                        f"changed={d7['frozen_candidate_expression_flag_changed_n']}"))
    checks.append(check("correction impact disclosed", d7["stable_status_changed_n"] == 17,
                        f"gene_status_changes={d7['stable_status_changed_n']}"))

    cond = pd.read_csv(POST / "d5_holdout_official_condition_audit.tsv", sep="\t")
    counts_ok = len(cond) == 15 and set(cond.n_biological_samples.astype(int)) == {3}
    checks.append(check("official holdout design", counts_ok,
                        f"conditions={len(cond)} sample_counts={sorted(cond.n_biological_samples.unique())}"))
    expected = {"UI", "WT", "dotA", "A", "B", "C", "D", "E", "F", "FGHI", "I", "O", "P", "R", "S"}
    checks.append(check("F/I/FGHI official conditions", set(cond.official_condition) == expected,
                        ";".join(cond.official_condition.astype(str))))

    null = pd.read_csv(POST / "d5_global_response_amplitude_matched_null.tsv", sep="\t")
    n_fdr = int((null.matched_empirical_fdr < 0.05).sum())
    checks.append(check("matched-null claim boundary", n_fdr == 0, f"FDR-significant contrasts={n_fdr}/12"))

    lit = pd.read_csv(POST / "d4_orthogonal_literature_evidence_v2.tsv", sep="\t", keep_default_na=False)
    checks.append(check("literature not used for tier assignment",
                        set(lit.used_for_tier_assignment) == {"NO"},
                        f"rows={len(lit)}"))
    bad_lit = lit.independent_study.str.contains("Machner & Isberg", regex=False).any() or (
        lit.loc[lit.pmid_pmcid.str.contains("19942850", regex=False), "verification_note"]
        .str.contains("Supports GEF/GDI-displacement").sum() != 1
    )
    checks.append(check("corrected PMID attribution", not bad_lit,
                        "15117975=Kagan et al.; 19942850=GEF/GDI displacement, not AMPylation"))

    edge = pd.read_csv(POST / "fig7_effector_rab1a_edges.tsv", sep="\t")
    source_text = " ".join(edge.evidence_source.astype(str))
    required_pmids = ["16906144", "20651120", "21680813", "21734656", "17952054", "21822290", "22158903", "30323948"]
    checks.append(check("Figure 7 direct-source mapping", all(x in source_text for x in required_pmids), source_text))

    pdf_names = [
        "Figure5_evidence_coupling_FIXED.pdf",
        "Figure6_trafficking_atlas_PUBREADY.pdf",
        "Figure7_RAB1A_evidence_architecture_FIXED.pdf",
        "SupplementaryFigure1_evidence_fragility.pdf",
        "SupplementaryFigure2_effector_deletion_context.pdf",
        "SupplementaryFigure3_species_composition.pdf",
        "SupplementaryFigure4_depmap_coessentiality.pdf",
        "SupplementaryFigure5_RAB_full_atlas.pdf",
    ]
    pdf_rows = []
    for name in pdf_names:
        path = POST / name
        ok = path.exists() and path.stat().st_size > 5000
        pages = 0
        width = height = 0.0
        text_chars = 0
        if ok:
            with fitz.open(path) as doc:
                pages = doc.page_count
                page = doc[0]
                width, height = page.rect.width, page.rect.height
                text_chars = len(page.get_text())
            ok = pages >= 1 and width > 0 and height > 0 and text_chars > 20
        pdf_rows.append({"file": name, "bytes": path.stat().st_size if path.exists() else 0,
                         "pages": pages, "width_pt": width, "height_pt": height,
                         "text_chars": text_chars, "sha256": sha256(path) if path.exists() else "",
                         "status": "PASS" if ok else "FAIL"})
    pdf_df = pd.DataFrame(pdf_rows)
    pdf_df.to_csv(POST / "d8_submission_figure_integrity.tsv", sep="\t", index=False)
    checks.append(check("submission PDF integrity", (pdf_df.status == "PASS").all(),
                        f"passed={(pdf_df.status == 'PASS').sum()}/{len(pdf_df)}"))

    compile_failures = []
    for path in sorted((ROOT / "scripts").glob("postfreeze_*.py")):
        try:
            py_compile.compile(str(path), doraise=True)
        except Exception as exc:  # pragma: no cover - audit path
            compile_failures.append(f"{path.name}: {exc}")
    checks.append(check("post-freeze Python syntax", not compile_failures,
                        "none" if not compile_failures else " | ".join(compile_failures)))

    controller = (ROOT / "run_pipeline_v2.ps1").read_text(encoding="utf-8-sig")
    checks.append(check("Stage-10 dependency includes Stage-5 statistics",
                        '07_convergence\\stage5_expression_statistics.json' in controller,
                        "report checkpoint invalidates when Stage 5 changes"))

    table = pd.DataFrame(checks)
    table.to_csv(POST / "d8_final_integrity_checks.tsv", sep="\t", index=False)
    failures = table.loc[table.status.eq("FAIL"), "check"].tolist()
    result = {
        "checks_total": len(table),
        "checks_passed": int((table.status == "PASS").sum()),
        "checks_failed": len(failures),
        "failures": failures,
        "fatal_error_status": "NO_DETECTED_FATAL_ERROR" if not failures else "FAILED_INTEGRITY_AUDIT",
        "submission_status": "MANUSCRIPT_DRAFTING_AUTHORIZED_FINAL_FREEZE_PENDING_EXTERNAL_G1_G2",
    }
    (POST / "d8_final_integrity_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(table.to_string(index=False))
    print(json.dumps(result, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
