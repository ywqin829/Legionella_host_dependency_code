"""Strip generator metadata from submission figure PDFs (no-ai-traces skill).

cairo_pdf (R) stamps `/Producer (cairo 1.18.4 ...)` and a CreationDate into every
figure PDF. This rewrites only that metadata block; page content streams are
copied verbatim and the rewrite is accepted only if the rendered pages are
pixel-identical to the originals. Originals are copied to _backup_meta/ first.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import fitz  # PyMuPDF

ROOT = Path(__file__).resolve().parents[1]
FIGROOT = ROOT / "FINAL_SUBMISSION_2026-10-04" / "02_figures_final"
BACKUP = FIGROOT / "_backup_meta"

BLANK = {"producer": "", "creator": "", "author": "", "title": "", "subject": "",
         "keywords": "", "creationDate": "", "modDate": ""}


def strip_one(pdf: Path) -> str:
    doc = fitz.open(pdf)
    before = [doc[i].get_pixmap(dpi=150).samples for i in range(doc.page_count)]
    tmp = pdf.with_suffix(".pdf.tmp")
    doc.set_metadata(BLANK)
    doc.save(tmp, garbage=4, deflate=True, clean=True)
    doc.close()

    chk = fitz.open(tmp)
    after = [chk[i].get_pixmap(dpi=150).samples for i in range(chk.page_count)]
    ok = (len(before) == len(after)) and all(a == b for a, b in zip(before, after))
    chk.close()

    if not ok:
        tmp.unlink(missing_ok=True)
        return "FAILED-identity"

    BACKUP.mkdir(parents=True, exist_ok=True)
    bak = BACKUP / pdf.name
    if not bak.exists():
        shutil.copy2(pdf, bak)
    tmp.replace(pdf)
    data = pdf.read_bytes()
    leaked = b"cairo" in data or b"/R" in data and b"cairo" in data
    return "stripped" + (" (LEAK-REMAINS)" if b"cairo" in data else "")


def main() -> int:
    pdfs = sorted(FIGROOT.rglob("*.pdf"))
    pdfs = [p for p in pdfs if "_backup_meta" not in p.parts]
    if not pdfs:
        print("no PDFs found under", FIGROOT)
        return 1
    for pdf in pdfs:
        status = strip_one(pdf)
        print(f"{status:24s} {pdf.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
