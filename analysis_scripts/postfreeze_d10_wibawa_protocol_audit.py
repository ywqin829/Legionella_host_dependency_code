from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import fitz
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "Using Genomic Deletion Mutants to Investigate Effector-Triggered Immunity During Legionella pneumophila Infection.pdf"
OUT = ROOT / "13_postfreeze_interpretation"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def main() -> int:
    if not PDF.exists():
        raise FileNotFoundError(PDF)
    doc = fitz.open(PDF)
    pages = [page.get_text() for page in doc]
    terms = {
        "growth_rate_confound": ["growth rate", "growth rates"],
        "bacterial_number_inflammation": ["bacterial number", "bacterial numbers"],
        "bacterial_load_alternative": ["bacterial load"],
        "delta_e": ["ΔE", "Delta E", "delta-E"],
        "delta_legc4": ["ΔlegC4", "Delta legC4", "delta-legC4"],
        "legc4_locus": ["LPW_20041", "lpw_20041"],
        "redundancy": ["functional redundancy", "redundancy"],
        "metaeffector": ["metaeffector", "metaeffectors"],
        "primer": ["primer", "primers"],
    }
    rows = []
    for label, patterns in terms.items():
        for page_no, page_text in enumerate(pages, start=1):
            for pattern in patterns:
                for match in re.finditer(re.escape(pattern), page_text, flags=re.I):
                    start = max(0, match.start() - 350)
                    end = min(len(page_text), match.end() + 650)
                    rows.append({"concept": label, "matched_text": pattern, "page": page_no,
                                 "excerpt": clean(page_text[start:end])})
    hits = pd.DataFrame(rows).drop_duplicates(["concept", "page", "excerpt"])
    hits.to_csv(OUT / "d10_wibawa_protocol_keyword_evidence.tsv", sep="\t", index=False)

    full_text = "\n\n".join(f"[[PAGE {i + 1}]]\n{text}" for i, text in enumerate(pages))
    table_like = [line.strip() for line in full_text.splitlines()
                  if re.search(r"Table|Δ[A-I]|Delta-[A-I]|LPW_", line, flags=re.I)]
    (OUT / "d10_wibawa_protocol_table_locus_lines.txt").write_text("\n".join(table_like), encoding="utf-8")

    counts = hits.groupby("concept").size().to_dict() if not hits.empty else {}
    result = {
        "pdf": str(PDF.relative_to(ROOT)),
        "sha256": digest(PDF),
        "pages": doc.page_count,
        "title": doc.metadata.get("title", ""),
        "author": doc.metadata.get("author", ""),
        "keyword_hit_counts": {k: int(v) for k, v in counts.items()},
        "complete_A_to_I_strain_table_detected": False,
        "interpretation": {
            "supports": [
                "mutant growth/fitness is a recognized confound in genomic-deletion infection experiments",
                "genomic-deletion phenotypes cannot be assigned automatically to a single effector",
                "Delta-E versus Delta-legC4 is an experimental exemplar of burden-aware interpretation",
                "redundancy and metaeffector antagonism are biologically plausible explanations",
            ],
            "does_not_support": [
                "quantitative burden adjustment for all E-MTAB-16655 mutants at 18 h in PMA-THP-1",
                "a complete Delta-A-to-I or O/P/R/S effector-content crosswalk",
                "closing G1 or G2 as fully resolved",
            ],
        },
        "recommended_status": "POSTFREEZE_INTERPRETATION_UPDATE_ONLY_CORE_PIPELINE_REMAINS_FROZEN",
    }
    (OUT / "d10_wibawa_protocol_audit_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print(hits[["concept", "page"]].value_counts().sort_index().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
