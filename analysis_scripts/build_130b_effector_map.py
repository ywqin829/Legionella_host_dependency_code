from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANNO = ROOT / "01_data_raw" / "lp130b_annotation"
OUT = ROOT / "13_postfreeze_interpretation"

gbk = (ANNO / "130b_FR687201.gbk").read_text(encoding="utf-8")

# ---- robust GBK feature + qualifier parser ----
FEAT_RE = re.compile(r"^\s{5}[A-Za-z_]+\s+(.*)$")
QUAL_RE = re.compile(r"^\s{21}/(\w+)=(.+)$")

cds = []
cur = None
for line in gbk.splitlines():
    fm = FEAT_RE.match(line)
    if fm:
        ftype = line.split()[0].strip()
        if ftype in ("CDS", "tRNA", "rRNA", "gene"):
            if cur and (cur.get("locus") or cur.get("gene")):
                cds.append(cur)
            cur = {"type": ftype, "loc": fm.group(1).strip(),
                   "locus": "", "gene": "", "product": "", "note": ""}
        else:
            if cur and (cur.get("locus") or cur.get("gene")):
                cds.append(cur)
            cur = None
        continue
    qm = QUAL_RE.match(line)
    if qm and cur is not None:
        k, v = qm.group(1), qm.group(2)
        v = v.strip().strip('"')
        if k == "locus_tag":
            cur["locus"] = v
        elif k == "gene":
            cur["gene"] = v
        elif k == "product":
            cur["product"] = v
        elif k == "note":
            cur["note"] = v
if cur and (cur.get("locus") or cur.get("gene")):
    cds.append(cur)

# ---- parse location to numeric range ----
def loc_range(loc: str):
    # strip join(123..456, ...) and complement()
    pos = loc.replace("complement(", "").replace("join(", "").replace(")", "")
    nums = re.findall(r"\d+", pos)
    if not nums:
        return None, None
    return int(nums[0]), int(nums[-1])

rows = []
for c in cds:
    if c["type"] != "CDS":
        continue
    s, e = loc_range(c["loc"])
    if s is None:
        continue
    m = re.search(r"(\d+)", c["locus"])
    num = int(m.group(1)) if m else None
    rows.append({
        "locus_tag": c["locus"], "locus_num": num,
        "start": s, "end": e,
        "gene": c["gene"], "product": c["product"], "note": c["note"],
    })
df = pd.DataFrame(rows)
df = df.drop_duplicates(subset=["locus_tag", "start"])

# ---- monotonicity check of locus numbers along genome ----
df = df.sort_values("start")
df["num_prev"] = df["locus_num"].shift(1)
df["mono_increase"] = (df["locus_num"] > df["num_prev"]).astype(int)
mono_frac = df["mono_increase"].mean()
print(f"CDS parsed: {len(df)}; locus-number monotonic-increase fraction along genome: {mono_frac:.3f}")
# count wraps (genome is circular)
wraps = (df["mono_increase"] == 0).sum()
print(f"wraps/decreases: {wraps}")

# ---- deletion regions from preprint Extended Data Table 1 ----
regions = {
    "A": (71, 221),
    "B": (4831, 4861),
    "C": (11451, 11581),
    "D": (16861, 16971),
    "E": (19951, 20201),
    "F": (23101, 23561),
    "H": (26521, 26861),
    "I": (27241, 27551),
}

out_rows = []
for mut, (lo, hi) in regions.items():
    genes = df[(df["locus_num"] >= lo) & (df["locus_num"] <= hi)].copy()
    genes = genes.sort_values("start")
    contig = genes["start"].max() - genes["start"].min()
    out_rows.append({
        "mutant": mut, "region_lpw": f"lpw{lo:05d}-lpw{hi:05d}",
        "n_genes": len(genes),
        "genome_span_nt": int(contig) if len(genes) else "",
        "gene_list": ";".join(genes["gene"].fillna("").astype(str)),
        "effector_candidates": ";".join(genes.loc[genes["gene"].astype(str).str.match(r"(?i)(leg|ceg|sid|sde|mav|rav|lem|lep|ank|vip|lid|set|lir|lyt|rlf|cag|piro|sidC|sidD)"), "gene"]),
        "products_short": ";".join(genes["product"].str[:40]),
    })
    print(f"\n=== mutant {mut} (lpw{lo:05d}-lpw{hi:05d}): {len(genes)} genes, span {contig} nt ===")
    for _, g in genes.iterrows():
        nm = g["gene"] if g["gene"] else f"[{g['locus_tag']}]"
        print(f"  {g['start']:>7,}-{g['end']:>7,}  {g['locus_tag']:<12} {nm:<10} {g['product'][:55]}")

map_df = pd.DataFrame(out_rows)
map_df.to_csv(OUT / "E-MTAB-16655_mutant_effector_mapping_A-I.tsv", sep="\t", index=False)
df.to_csv(ANNO / "130b_all_CDS_lpw_coords.tsv", sep="\t", index=False)
print("\nwrote", OUT / "E-MTAB-16655_mutant_effector_mapping_A-I.tsv")
