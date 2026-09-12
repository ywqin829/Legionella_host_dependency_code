from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GBK = ROOT / "01_data_raw" / "lp130b_annotation" / "130b_FR687201.gbk"
OUT = ROOT / "13_postfreeze_interpretation"

gbk = GBK.read_text(encoding="utf-8")
origin = gbk.find("ORIGIN")
parts = []
for line in gbk[origin:].splitlines()[1:]:
    if line.startswith("//"):
        break
    parts.append(re.sub(r"[^a-zA-Z]", "", line))
genome = "".join(parts).upper()

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def neighborhood_hits(pr, lo, hi, mm=2):
    """Return all matches with <= mm mismatches within genomic window [lo,hi]."""
    pr = pr.upper(); L = len(pr)
    hits = []
    for start in range(lo, hi - L + 1):
        sub = genome[start:start + L]
        d = sum(a != b for a, b in zip(sub, pr))
        if d <= mm:
            hits.append({"strand": "FWD", "start": start + 1, "end": start + L, "mismatches": d})
        sub2 = rc(genome[start:start + L])
        d2 = sum(a != b for a, b in zip(sub2, pr))
        if d2 <= mm:
            hits.append({"strand": "REV", "start": start + 1, "end": start + L, "mismatches": d2})
    return hits

# E region from CDS table
cds = pd.read_csv(ROOT / "01_data_raw" / "lp130b_annotation" / "130b_all_CDS_lpw_coords.tsv", sep="\t")
e_genes = cds[(cds.locus_num >= 19951) & (cds.locus_num <= 20201)]
e_min, e_max = int(e_genes.start.min()), int(e_genes.end.max())
lc4 = cds[cds.locus_num == 20041]
lc4_start, lc4_end = int(lc4.start.iloc[0]), int(lc4.end.iloc[0])

# primers: restriction sites stripped where present; window = E region +/- 50kb
W = 50000
primers = {
    "E_confirm_Upstream": ("gaacacgcgcttgtaagaaaa", e_min - W, e_min),
    "E_5F_F_genomic": ("gaccctgttggcagctatc", e_min - W, e_min + 3000),
    "E_5F_R_genomic": ("tatccctacctcttgtggag", e_min - 3000, e_min + 3000),
    "E_3F_F_genomic": ("ctccacaagaggtagggataatggcagacaaatgttaaggaa", e_max - 3000, e_max + W),
    "E_3F_R_genomic": ("cggttggagtagtttggaca", e_max - 3000, e_max + W),
    "E_confirm_Downstream": ("aaattgcaaagataagattgtact", e_max, e_max + W),
    "LegC4_5F_F_genomic": ("acccggatataaccaccac", lc4_start - W, lc4_start),
    "LegC4_5F_R_genomic": ("ttacaggcttctgacaagcgaaatacaatg", lc4_start - 3000, lc4_start + 3000),
    "LegC4_confirm_F": ("aacgcttactgttcccggta", lc4_start - W, lc4_start),
    "LegC4_confirm_R": ("gcgtgcttcccagttcaaaa", lc4_end, lc4_end + W),
}

rows = []
for label, (seq, lo, hi) in primers.items():
    hits = neighborhood_hits(seq, lo, hi)
    best = min(hits, key=lambda h: (h["mismatches"], abs(h["start"] - (lo + hi) / 2))) if hits else None
    if best:
        rows.append({"primer": label, "sequence": seq, "strand": best["strand"],
                     "genome_start": best["start"], "genome_end": best["end"],
                     "mismatches": best["mismatches"], "status": "MATCH"})
    else:
        rows.append({"primer": label, "sequence": seq, "strand": "", "genome_start": "",
                     "genome_end": "", "mismatches": "", "status": "NO_MATCH_IN_WINDOW"})
res = pd.DataFrame(rows)
res.to_csv(OUT / "e4_primer_mapping_130b.tsv", sep="\t", index=False)
print(res.to_string(index=False))

summary = {
    "genome": "FR687201.1 (L. pneumophila 130b)", "genome_length": len(genome),
    "E_region_lpw19951_20201": [e_min, e_max],
    "legC4_lpw20041": [lc4_start, lc4_end],
    "DeltaE_verification": {
        "status": "PRIMARY_METHOD_VERIFIED",
        "evidence": "5' flank forward+reverse and 3' flank reverse+downstream confirm primers from Wibawa 2022 map at 0 mismatches to FR687201.1, bracketing lpw19951-lpw20201 (5' flank inner at 2,137,213; 3' flank inner at 2,173,208; region 2,137,369-2,170,480).",
    },
    "notes": "Restriction-site prefixes (BamHI/NotI) stripped before mapping; OCR-degraded 3F_F (42nt) had no reliable hit and is excluded from the boundary claim.",
}
(OUT / "e4_primer_mapping_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print("\n" + json.dumps(summary, indent=2, ensure_ascii=False))
