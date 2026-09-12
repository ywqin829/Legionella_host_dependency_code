from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_postfreeze_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

# ---- effectors per confirmed region (from preprint E-Table1 x FR687201.1 GBK) ----
regions = {
    "A": {"lpw": "lpw00071-lpw00221", "status": "PASS", "effectors": ["cegC1"], "n": 16},
    "B": {"lpw": "lpw04831-lpw04861", "status": "PASS", "effectors": ["legA7"], "n": 4},
    "C": {"lpw": "lpw11451-lpw11581", "status": "PASS", "effectors": ["lem4", "lem5"], "n": 14},
    "D": {"lpw": "lpw16861-lpw16971", "status": "PASS", "effectors": ["legL3"], "n": 12},
    "E": {"lpw": "lpw19951-lpw20201", "status": "FAIL", "effectors": ["lem16", "lem17", "ralF", "legC4", "ltpG", "lirA", "lirB", "lirE", "lirF", "pieE"], "n": 26},
    "F": {"lpw": "lpw23101-lpw23561", "status": "FAIL", "effectors": ["legK2", "legAU13", "sdeC", "sde", "sidJ", "sdeB", "sdeA", "lem19", "legS2"], "n": 45},
    "H": {"lpw": "lpw26521-lpw26861", "status": "in_FGHI", "effectors": ["ceg30", "legA14", "legA15", "sidM", "sidD"], "n": 35},
    "I": {"lpw": "lpw27241-lpw27551", "status": "FAIL", "effectors": ["ceg32", "sdjA", "sdeD", "sdcA", "sidC", "lem26", "lem27"], "n": 32},
}

# ---- cross-reference: effector relationship table (Philadelphia-1) -> 130b regions ----
# known redundancy / antagonism / metaeffector modules by member effector names
REL_MODULES = [
    {
        "module": "Rab1 AMPylation add/remove",
        "members": ["sidM", "sidD"],
        "relation": "antagonistic_host_modification",
        "regions": ["H"],
    },
    {
        "module": "Rab1 activation cycle (SidM <-> LepB)",
        "members": ["sidM", "lepB"],
        "relation": "shared_pathway_non_equivalent",
        "regions": ["H"],
    },
    {
        "module": "SidE PR-ubiquitination family (paralogs)",
        "members": ["sdeA", "sdeB", "sdeC", "sde"],
        "relation": "paralog_redundancy",
        "regions": ["F"],
    },
    {
        "module": "SidE <-> SidJ metaeffector (PR-Ub writer vs inactivator)",
        "members": ["sdeA", "sdeB", "sdeC", "sde", "sidJ"],
        "relation": "direct_metaeffector_regulation",
        "regions": ["F"],
    },
    {
        "module": "SidC/SdcA Rab10-PI4P tethers (paralogs)",
        "members": ["sidC", "sdcA"],
        "relation": "paralog_redundancy",
        "regions": ["I"],
    },
    {
        "module": "Lgt translation-inhibition paralogs",
        "members": ["lgt1", "lgt2", "lgt3"],
        "relation": "paralog_redundancy",
        "regions": [],
    },
    {
        "module": "MavC/MvcA UBE2N add/remove",
        "members": ["mavC", "mvcA"],
        "relation": "antagonistic_host_modification",
        "regions": [],
    },
    {
        "module": "AnkX/Lem3 phosphocholination add/remove",
        "members": ["ankX", "lem3"],
        "relation": "antagonistic_host_modification",
        "regions": [],
    },
]

# annotate each region with which relationship modules are FULLY contained
rows = []
for mut, info in regions.items():
    eff_low = {e.lower() for e in info["effectors"]}
    contained = []
    partial = []
    for mod in REL_MODULES:
        if not mod["regions"]:
            continue
        members = {m.lower() for m in mod["members"]}
        if members <= eff_low:
            contained.append(mod["module"])
        elif members & eff_low:
            partial.append(mod["module"])
    # redundancy-masked flag: region deletes BOTH sides of an add/remove or a full paralog family
    has_full_antag_module = any("add/remove" in c or "metaeffector" in c or "family" in c for c in contained)
    rows.append({
        "mutant": mut,
        "region_lpw": info["lpw"],
        "n_genes": info["n"],
        "holdout_status": info["status"],
        "effectors_deleted": ";".join(info["effectors"]),
        "relationship_modules_fully_deleted": ";".join(contained) if contained else "none",
        "relationship_modules_partially_deleted": ";".join(partial) if partial else "none",
        "redundancy_or_antagonism_module": "YES" if has_full_antag_module else "no",
    })

df = pd.DataFrame(rows)
df.to_csv(OUT / "c1_region_redundancy_x_holdout.tsv", sep="\t", index=False)

# ---- descriptive 2x2 (KNOWN regions only; O/P/R/S and FGHI composite excluded) ----
known = df[df.holdout_status.isin(["PASS", "FAIL"])]
tbl = pd.crosstab(known.holdout_status, known.redundancy_or_antagonism_module)
print("=== redundancy module x holdout status (confirmed single regions only) ===")
print(tbl.to_string())
print()
print(df[["mutant", "holdout_status", "effectors_deleted", "relationship_modules_fully_deleted", "redundancy_or_antagonism_module"]].to_string(index=False))

summary = {
    "regions_with_mapping": list(regions),
    "regions_with_holdout_single_status": known.mutant.tolist(),
    "redundancy_crosstab": tbl.to_dict(),
    "note": "A-F,H,I mapping inferred from preprint E-Table1 x FR687201.1 GBK; O/P/R/S unmapped; FGHI composite (F+H+I+G) excluded from 2x2; descriptive only, not a formal enrichment test",
}
(OUT / "c1_redundancy_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print("\n" + json.dumps(summary, indent=2, ensure_ascii=False))
