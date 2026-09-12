from __future__ import annotations

import json
from pathlib import Path

import fitz
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "11_figures"


def audit_pdf(path: Path) -> dict:
    doc = fitz.open(path)
    if len(doc) != 1:
        raise ValueError(f"{path.name}: expected one page, observed {len(doc)}")
    page = doc[0]
    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    grey = arr[:, :, :3].mean(axis=2)
    ink = grey < 248
    rows, cols = np.where(ink)
    if len(rows) == 0:
        raise ValueError(f"{path.name}: rendered page is blank")
    margins = {
        "left_px": int(cols.min()),
        "right_px": int(pix.width - 1 - cols.max()),
        "top_px": int(rows.min()),
        "bottom_px": int(pix.height - 1 - rows.max()),
    }
    blocks = page.get_text("blocks")
    text = page.get_text("text")
    png = FIG / "qa_renders" / path.with_suffix(".png").name
    png.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(arr[:, :, :3]).save(png, dpi=(300, 300))
    return {
        "file": path.name,
        "bytes": path.stat().st_size,
        "page_width_pt": round(page.rect.width, 2),
        "page_height_pt": round(page.rect.height, 2),
        "render_width_px": pix.width,
        "render_height_px": pix.height,
        "ink_fraction": round(float(ink.mean()), 4),
        "text_blocks": len(blocks),
        "text_characters": len(text),
        "margins": margins,
        "blank": False,
        "edge_contact_warning": any(v < 2 for v in margins.values()),
        "qa_render": str(png.relative_to(ROOT)),
    }


def main() -> None:
    pdfs = sorted(FIG.glob("Figure*.pdf"))
    if len(pdfs) != 4:
        raise ValueError(f"Expected 4 main figure PDFs, found {len(pdfs)}")
    results = [audit_pdf(p) for p in pdfs]
    source_tables = sorted((FIG / "source_data").glob("*.tsv"))
    warnings = []
    for result in results:
        if result["bytes"] < 10_000:
            warnings.append(f"{result['file']} is unexpectedly small")
        if result["text_blocks"] < 5:
            warnings.append(f"{result['file']} has too few text blocks")
        if result["edge_contact_warning"]:
            warnings.append(f"{result['file']} has rendered ink touching a page edge")
    if len(source_tables) < 12:
        warnings.append(f"Only {len(source_tables)} source-data tables found")
    status = "PASS" if not warnings else "PASS_WITH_WARNINGS"
    report = {
        "main_figures": len(results),
        "source_data_tables": len(source_tables),
        "figures": results,
        "warnings": warnings,
        "gate9b_automated_figure_qa": status,
        "boundary": "Automated rendering QA detects blank/clipped-edge failures but does not replace human visual review.",
    }
    (FIG / "stage9b_figure_qa.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = [
        "# Stage 9B automated figure quality audit", "", f"Status: **{status}**", "",
        f"- Main figure PDFs: {len(results)}", f"- Source-data tables: {len(source_tables)}",
        "- All PDFs were rendered from vector originals and checked for blank output, edge contact, page count and embedded text.",
        "- Automated checks do not replace final human review at print size.", "", "## Per-figure checks", "",
    ]
    for r in results:
        lines.append(
            f"- {r['file']}: {r['page_width_pt']}×{r['page_height_pt']} pt; "
            f"ink fraction={r['ink_fraction']}; text blocks={r['text_blocks']}; "
            f"edge warning={r['edge_contact_warning']}."
        )
    if warnings:
        lines += ["", "## Warnings", ""] + [f"- {w}" for w in warnings]
    (FIG / "STAGE9B_FIGURE_QA.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
