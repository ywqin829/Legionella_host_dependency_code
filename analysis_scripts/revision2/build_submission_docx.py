from __future__ import annotations

import difflib
import re
import subprocess
import tempfile
from pathlib import Path

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_LINE_SPACING, WD_PARAGRAPH_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[2]
PANDOC = Path(r"C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools\pandoc.exe")
BASELINE = ROOT / "FINAL_SUBMISSION_2026-10-04" / "01_manuscript" / "Legionella_host_dependency_manuscript_revised.md"
REVISED = ROOT / "revision2_2026-11-23" / "REVISION2_MANUSCRIPT_WORKING_DRAFT.md"
RESPONSE = ROOT / "revision2_2026-11-23" / "REVISION2_RESPONSE_TO_REVIEWERS.md"
SUPP_LEGENDS = ROOT / "revision2_2026-11-23" / "REVISION2_SUPPLEMENTARY_LEGENDS.md"
OUTDIR = ROOT / "revision2_2026-11-23" / "submission_files"


def set_east_asia_font(style, name: str) -> None:
    style.font.name = name
    style._element.rPr.rFonts.set(qn("w:eastAsia"), name)


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_PARAGRAPH_ALIGNMENT.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.append(begin)
    run._r.append(instruction)
    run._r.append(end)


def make_reference_docx(path: Path) -> None:
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)

    sect_pr = section._sectPr
    line_numbers = sect_pr.find(qn("w:lnNumType"))
    if line_numbers is None:
        line_numbers = OxmlElement("w:lnNumType")
        sect_pr.append(line_numbers)
    line_numbers.set(qn("w:countBy"), "1")
    line_numbers.set(qn("w:restart"), "continuous")
    line_numbers.set(qn("w:distance"), "360")

    normal = doc.styles["Normal"]
    set_east_asia_font(normal, "Times New Roman")
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    normal.paragraph_format.space_after = Pt(0)

    for name in ["Title", "Subtitle", "Heading 1", "Heading 2", "Heading 3", "Heading 4"]:
        style = doc.styles[name]
        set_east_asia_font(style, "Times New Roman")
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.font.size = Pt(12 if name != "Title" else 14)
        style.font.bold = True

    addition = doc.styles.add_style("Addition", WD_STYLE_TYPE.CHARACTER)
    set_east_asia_font(addition, "Times New Roman")
    addition.font.color.rgb = RGBColor(0, 0, 255)
    addition.font.underline = True

    deletion = doc.styles.add_style("Deletion", WD_STYLE_TYPE.CHARACTER)
    set_east_asia_font(deletion, "Times New Roman")
    deletion.font.color.rgb = RGBColor(192, 0, 0)
    deletion.font.strike = True

    add_glossary_box_style(doc)

    add_page_number(section.footer.paragraphs[0])
    doc.save(path)


def add_glossary_box_style(doc) -> None:
    """Box 1 (glossary) paragraph style; the border/shading are applied later by
    wrapping the block in a single-cell table so it reads as one box."""
    style = doc.styles.add_style("GlossaryBox", WD_STYLE_TYPE.PARAGRAPH)
    set_east_asia_font(style, "Times New Roman")
    style.font.size = Pt(10.5)
    style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    style.paragraph_format.space_before = Pt(0)
    style.paragraph_format.space_after = Pt(4)


def _set_cell_borders_and_shading(table) -> None:
    tbl_pr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right"):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), "8")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), "404040")
        borders.append(element)
    cell_margins = OxmlElement("w:tblCellMar")
    for edge, width in (("top", "100"), ("left", "120"), ("bottom", "100"), ("right", "120")):
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:w"), width)
        element.set(qn("w:type"), "dxa")
        cell_margins.append(element)
    anchor = tbl_pr.find(qn("w:tblLook"))
    if anchor is not None:
        anchor.addprevious(borders)
        anchor.addprevious(cell_margins)
    else:
        tbl_pr.append(borders)
        tbl_pr.append(cell_margins)

    shading = OxmlElement("w:shd")
    shading.set(qn("w:val"), "clear")
    shading.set(qn("w:color"), "auto")
    shading.set(qn("w:fill"), "F5F5F5")
    table.cell(0, 0)._tc.get_or_add_tcPr().append(shading)


def wrap_glossary_box(path: Path) -> None:
    """Wrap the contiguous GlossaryBox paragraphs in a single bordered cell."""
    doc = Document(str(path))
    paras = doc.paragraphs
    indices = [i for i, p in enumerate(paras) if p.style.name == "GlossaryBox"]
    if not indices:
        return
    first_i, last_i = indices[0], indices[-1]
    block = paras[first_i:last_i + 1]
    body = doc.element.body
    position = list(body).index(block[0]._p)

    table = doc.add_table(rows=1, cols=1)
    cell = table.cell(0, 0)
    while cell.paragraphs and not cell.paragraphs[0].text:
        cell._tc.remove(cell.paragraphs[0]._p)
    for para in block:
        cell._tc.append(para._p)
    _set_cell_borders_and_shading(table)

    body.remove(table._tbl)
    body.insert(position, table._tbl)
    doc.save(str(path))


def protect_span_text(text: str) -> str:
    return text.replace("[", r"\[").replace("]", r"\]")


def mark_inline(text: str, style: str) -> str:
    if not text.strip():
        return text
    return f'[{protect_span_text(text)}]{{custom-style="{style}"}}'


def mark_line(line: str, style: str) -> str:
    if not line.strip():
        return line
    if line.lstrip().startswith(":::"):
        return line
    if re.fullmatch(r"\s*\|?(?:\s*:?-+:?\s*\|)+\s*", line):
        return line
    if line.lstrip().startswith("|") and line.rstrip().endswith("|"):
        cells = line.split("|")
        for i in range(1, len(cells) - 1):
            cells[i] = mark_inline(cells[i], style)
        return "|".join(cells)
    heading = re.match(r"^(#{1,6}\s+)(.*)$", line)
    if heading:
        return heading.group(1) + mark_inline(heading.group(2), style)
    item = re.match(r"^(\s*(?:[-+*]|\d+[.)])\s+)(.*)$", line)
    if item:
        return item.group(1) + mark_inline(item.group(2), style)
    return mark_inline(line, style)


def make_marked_markdown(old_text: str, new_text: str) -> str:
    old_lines = old_text.splitlines()
    new_lines = new_text.splitlines()
    matcher = difflib.SequenceMatcher(a=old_lines, b=new_lines, autojunk=False)
    output = [
        "# Marked-up manuscript",
        "",
        '[Additions are shown in blue underline; deletions are shown in red strikethrough.]{custom-style="Addition"}',
        "",
    ]
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            output.extend(old_lines[i1:i2])
        elif tag == "delete":
            output.extend(mark_line(line, "Deletion") for line in old_lines[i1:i2])
        elif tag == "insert":
            output.extend(mark_line(line, "Addition") for line in new_lines[j1:j2])
        else:
            output.extend(mark_line(line, "Deletion") for line in old_lines[i1:i2])
            if i2 > i1 and j2 > j1:
                output.append("")
            output.extend(mark_line(line, "Addition") for line in new_lines[j1:j2])
    return "\n".join(output) + "\n"


def convert(source: Path, target: Path, reference: Path) -> None:
    command = [
        str(PANDOC),
        str(source),
        "--from=markdown+smart",
        "--to=docx",
        f"--reference-doc={reference}",
        "--wrap=none",
        f"--output={target}",
    ]
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    if not PANDOC.exists():
        raise FileNotFoundError(PANDOC)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="revision2_docx_") as tmp:
        tmpdir = Path(tmp)
        reference = tmpdir / "reference.docx"
        marked_md = tmpdir / "marked.md"
        make_reference_docx(reference)
        marked_md.write_text(
            make_marked_markdown(
                BASELINE.read_text(encoding="utf-8"),
                REVISED.read_text(encoding="utf-8"),
            ),
            encoding="utf-8",
        )
        convert(REVISED, OUTDIR / "Legionella_host_dependency_clean.docx", reference)
        convert(marked_md, OUTDIR / "Legionella_host_dependency_marked_up.docx", reference)
        convert(RESPONSE, OUTDIR / "Response_to_Reviewers.docx", reference)
        convert(SUPP_LEGENDS, OUTDIR / "Supplementary_Legends.docx", reference)
        wrap_glossary_box(OUTDIR / "Legionella_host_dependency_clean.docx")
        wrap_glossary_box(OUTDIR / "Legionella_host_dependency_marked_up.docx")


if __name__ == "__main__":
    main()
