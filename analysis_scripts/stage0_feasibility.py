from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config/project_config.json").read_text(encoding="utf-8"))
RAW = ROOT / "01_data_raw/metadata"
OUT = ROOT / "02_metadata"
CP = ROOT / "checkpoints"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "LegionellaHostDependencyAtlas/1.0"})


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get(url: str, path: Path | None = None, retries: int = 2) -> bytes:
    if path and path.exists() and path.stat().st_size:
        return path.read_bytes()
    errors = []
    for n in range(retries + 1):
        try:
            r = SESSION.get(url, timeout=90)
            r.raise_for_status()
            if path:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(r.content)
            return r.content
        except Exception as e:  # noqa: BLE001
            errors.append(f"{type(e).__name__}: {e}")
            if n < retries:
                time.sleep(2**n)
    raise RuntimeError(" | ".join(errors))


def add(obj: dict[str, Any], key: str, value: str) -> None:
    if key not in obj:
        obj[key] = value
    elif isinstance(obj[key], list):
        obj[key].append(value)
    else:
        obj[key] = [obj[key], value]


def geo_soft(acc: str) -> dict[str, Any]:
    prefix = acc[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/soft/{acc}_family.soft.gz"
    path = RAW / f"{acc}_family.soft.gz"
    text = gzip.decompress(get(url, path)).decode("utf-8", errors="replace")
    series: dict[str, Any] = {}
    samples: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    kind = ""
    for line in text.splitlines():
        if line.startswith("^SERIES"):
            current, kind = series, "Series"
            series["accession"] = line.split("=", 1)[1].strip()
        elif line.startswith("^SAMPLE"):
            current, kind = {"accession": line.split("=", 1)[1].strip()}, "Sample"
            samples.append(current)
        elif line.startswith("^"):
            current, kind = None, ""
        elif current is not None and line.startswith("!") and "=" in line:
            key, value = line[1:].split("=", 1)
            key = key.strip()
            if key.startswith(kind + "_"):
                key = key[len(kind) + 1 :]
            add(current, key.lower(), value.strip())
    return {"series": series, "samples": samples, "path": str(path), "sha256": sha(path)}


def values(records: list[dict[str, Any]], key: str) -> list[str]:
    out = []
    for r in records:
        v = r.get(key, [])
        out.extend(map(str, v if isinstance(v, list) else [v]))
    return [x for x in out if x]


def sra_info(acc: str) -> tuple[pd.DataFrame, Path]:
    path = RAW / f"{acc}_SRA_RunInfo.csv"
    raw = get(f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={acc}", path)
    try:
        return pd.read_csv(io.BytesIO(raw)), path
    except Exception:  # noqa: BLE001
        return pd.DataFrame(), path


def summarize_geo(spec: dict[str, Any]) -> tuple[dict[str, Any], pd.DataFrame]:
    parsed = geo_soft(spec["accession"])
    s, samples = parsed["series"], parsed["samples"]
    runs, run_path = sra_info(spec["accession"])
    organisms = sorted(set(values(samples, "organism_ch1")))
    sources = sorted(set(values(samples, "source_name_ch1")))
    blob_parts = list(map(str, [s.get("title", ""), s.get("summary", ""), s.get("overall_design", "")]))
    blob_parts += values(samples, "title") + sources + values(samples, "characteristics_ch1")
    blob = " ".join(blob_parts).lower()
    pmids = s.get("pubmed_id", [])
    pmids = pmids if isinstance(pmids, list) else [pmids]
    relations = s.get("relation", [])
    relations = relations if isinstance(relations, list) else [relations]
    size_mb = float(pd.to_numeric(runs.get("size_MB", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
    row = {
        "accession": spec["accession"], "role": spec["role"], "model_planned": spec["model"],
        "contrast_family": spec["contrast_family"], "heldout": spec["heldout"], "title": s.get("title", ""),
        "status": s.get("status", ""), "sample_count": len(samples), "sra_runs": len(runs),
        "estimated_raw_mb": round(size_mb, 2), "organisms": " | ".join(organisms), "sources": " | ".join(sources),
        "human_model": any("homo sapiens" in x.lower() for x in organisms),
        "mouse_model": any("mus musculus" in x.lower() for x in organisms),
        "mutant_or_perturbation": any(x in blob for x in ["mutant", "delta", "deletion", "knockout", "-/-", "deficient"]),
        "lung_or_pneumonia": any(x in blob for x in ["lung", "pneumonia", "pulmonary", "alveolar"]),
        "pubmed_ids": " | ".join(map(str, pmids)), "relations": " | ".join(map(str, relations)),
        "soft_path": parsed["path"], "soft_sha256": parsed["sha256"], "runinfo_path": str(run_path),
    }
    sample_rows = []
    for x in samples:
        def j(k: str) -> str:
            v = x.get(k, "")
            return " | ".join(map(str, v)) if isinstance(v, list) else str(v)
        sample_rows.append({"series": spec["accession"], "sample": x.get("accession", ""), "title": j("title"),
                            "source": j("source_name_ch1"), "organism": j("organism_ch1"),
                            "characteristics": j("characteristics_ch1"), "treatment": j("treatment_protocol_ch1"),
                            "relation": j("relation")})
    return row, pd.DataFrame(sample_rows)


def ena_holdout(spec: dict[str, Any]) -> tuple[dict[str, Any], pd.DataFrame]:
    acc = spec["accession"]
    fields = "run_accession,experiment_accession,sample_accession,secondary_sample_accession,scientific_name,study_title,experiment_title,sample_title,sample_description,library_strategy,library_source,library_selection,library_layout,instrument_model,fastq_bytes,fastq_ftp"
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields={fields}&format=tsv"
    path = RAW / f"{acc}_ENA_run_report.tsv"
    df = pd.read_csv(io.BytesIO(get(url, path)), sep="\t")
    size = float(pd.to_numeric(df.get("fastq_bytes", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) / 1e6
    title = " | ".join(sorted(set(df.get("study_title", pd.Series(dtype=str)).dropna().astype(str))))
    return {"accession": acc, "role": spec["role"], "model_planned": spec["model"], "title": title,
            "sample_count": int(df.sample_accession.nunique()), "run_count": len(df), "estimated_raw_mb": round(size, 2),
            "organisms": " | ".join(sorted(set(df.scientific_name.dropna().astype(str)))), "path": str(path), "sha256": sha(path)}, df


def taxids(x: str) -> set[int]:
    return {int(v) for v in re.findall(r"taxid:(\d+)", x)}


def uniprot(x: str) -> str:
    for part in x.split("|"):
        if part.startswith("uniprotkb:"):
            return part.split(":", 1)[1].split("(", 1)[0]
    return x.split("|", 1)[0]


def audit_intact(path: Path) -> tuple[dict[str, Any], pd.DataFrame]:
    targets = set(CFG["all_target_taxids"])
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\r\n").split("\t")
            if len(c) < 15:
                continue
            ta, tb = taxids(c[9]), taxids(c[10])
            at, bt = bool(ta & targets), bool(tb & targets)
            ah, bh = 9606 in ta, 9606 in tb
            if not ((at and bh) or (bt and ah)):
                continue
            typ, method, source = c[11], c[6], c[12]
            negative = len(c) > 35 and c[35].lower() == "true"
            physical = (any(k in typ.lower() for k in ["direct interaction", "physical association", "association"])
                        and "genetic" not in typ.lower() and not negative)
            rows.append({"bacterial_accession": uniprot(c[0] if at else c[1]),
                         "human_accession": uniprot(c[1] if at else c[0]), "method": method, "type": typ,
                         "source": source, "pmid": "|".join(re.findall(r"pubmed:(\d+)", c[8])),
                         "negative": negative, "physical_or_direct": physical})
    df = pd.DataFrame(rows)
    pos = df[~df.negative].copy() if not df.empty else df
    phy = pos[pos.physical_or_direct] if not pos.empty else pos
    counts = pos.source.value_counts() if not pos.empty else pd.Series(dtype=int)
    pids = {p for x in pos.get("pmid", []) for p in str(x).split("|") if p}
    summary = {"records": len(pos), "unique_pairs": len(pos.drop_duplicates(["bacterial_accession", "human_accession"])) if not pos.empty else 0,
               "unique_bacterial_proteins": int(pos.bacterial_accession.nunique()) if not pos.empty else 0,
               "unique_human_targets": int(pos.human_accession.nunique()) if not pos.empty else 0,
               "unique_physical_human_targets": int(phy.human_accession.nunique()) if not phy.empty else 0,
               "unique_pmids": len(pids), "source_count": int(pos.source.nunique()) if not pos.empty else 0,
               "largest_source": str(counts.index[0]) if len(counts) else "",
               "largest_source_fraction": float(counts.iloc[0] / counts.sum()) if len(counts) else 1.0,
               "mitab_path": str(path), "mitab_sha256": sha(path)}
    return summary, df


def verify_pubmed(pmid: str) -> dict[str, Any]:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    item = json.loads(get(url).decode())["result"].get(pmid, {})
    return {"pmid": pmid, "title": item.get("title", ""), "pubdate": item.get("pubdate", ""),
            "source": item.get("source", ""), "verified": bool(item.get("title"))}


def novelty_search() -> pd.DataFrame:
    queries = ['"Legionella pneumophila" AND host AND interactome AND transcriptome',
               '"Legionella pneumophila" AND CRISPR AND transcriptome AND host',
               '"Legionella pneumophila" AND effector AND "ribosome profiling"',
               '"Legionella pneumophila" AND "multi-omics" AND host']
    rows = []
    for q in queries:
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + requests.utils.quote(q) + "&format=json&pageSize=25"
        for x in json.loads(get(url).decode()).get("resultList", {}).get("result", []):
            rows.append({"query": q, "title": x.get("title", ""), "authors": x.get("authorString", ""),
                         "journal": x.get("journalTitle", ""), "year": x.get("pubYear", ""),
                         "pmid": x.get("pmid", ""), "pmcid": x.get("pmcid", ""), "doi": x.get("doi", "")})
    return pd.DataFrame(rows).drop_duplicates(["title", "doi", "pmid"]) if rows else pd.DataFrame()


def independence(data: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    accs = data.accession.tolist()
    sets = {a: set(samples.loc[samples.series == a, "sample"].astype(str)) for a in accs if a != "PRJEB108343"}
    pids = {r.accession: {x for x in str(r.pubmed_ids).split(" | ") if x and x != "nan"} for _, r in data.iterrows()}
    rows = []
    for i, a in enumerate(accs):
        for b in accs[i + 1:]:
            so, po = sets.get(a, set()) & sets.get(b, set()), pids.get(a, set()) & pids.get(b, set())
            rows.append({"dataset_a": a, "dataset_b": b, "independent_preliminary": not so and not po,
                         "shared_samples": len(so), "shared_pmids": "|".join(sorted(po))})
    return pd.DataFrame(rows)


def main() -> int:
    for p in [RAW, OUT, CP]:
        p.mkdir(parents=True, exist_ok=True)
    failures, drows, sframes = [], [], []
    for spec in CFG["stage0"]["geo_datasets"]:
        try:
            row, smp = summarize_geo(spec); drows.append(row); sframes.append(smp)
        except Exception as e:  # noqa: BLE001
            failures.append({"resource": spec["accession"], "error": f"{type(e).__name__}: {e}"})
    ddf = pd.DataFrame(drows); sdf = pd.concat(sframes, ignore_index=True) if sframes else pd.DataFrame()
    hold, hdf = {}, pd.DataFrame()
    try:
        hold, hdf = ena_holdout(CFG["stage0"]["ena_holdout"])
        ddf = pd.concat([ddf, pd.DataFrame([{"accession": hold["accession"], "role": hold["role"], "model_planned": hold["model_planned"],
            "contrast_family": "multi_effector_external_holdout", "heldout": True, "title": hold["title"], "status": "Public ENA",
            "sample_count": hold["sample_count"], "sra_runs": hold["run_count"], "estimated_raw_mb": hold["estimated_raw_mb"],
            "organisms": hold["organisms"], "sources": "THP-1 macrophages", "human_model": True, "mouse_model": False,
            "mutant_or_perturbation": True, "lung_or_pneumonia": False, "pubmed_ids": "", "relations": "ENA BioProject"}])], ignore_index=True)
    except Exception as e:  # noqa: BLE001
        failures.append({"resource": "PRJEB108343", "error": f"{type(e).__name__}: {e}"})
    isum, idf = {}, pd.DataFrame()
    try:
        isum, idf = audit_intact(Path(CFG["stage0"]["cached_intact_mitab"]))
    except Exception as e:  # noqa: BLE001
        failures.append({"resource": "IntAct", "error": f"{type(e).__name__}: {e}"})
    crispr = {}
    try:
        crispr = verify_pubmed(CFG["stage0"]["crispr"]["pmid"])
    except Exception as e:  # noqa: BLE001
        failures.append({"resource": "CRISPR PMID", "error": f"{type(e).__name__}: {e}"})
    try:
        nov = novelty_search()
    except Exception as e:  # noqa: BLE001
        nov = pd.DataFrame(); failures.append({"resource": "Europe PMC", "error": f"{type(e).__name__}: {e}"})

    ddf.to_csv(OUT / "dataset_metadata_audit.tsv", sep="\t", index=False)
    sdf.to_csv(OUT / "geo_sample_manifest.tsv", sep="\t", index=False)
    hdf.to_csv(OUT / "holdout_run_manifest.tsv", sep="\t", index=False)
    idf.to_csv(OUT / "intact_legionella_human_records.tsv", sep="\t", index=False)
    nov.to_csv(OUT / "novelty_search_screen.tsv", sep="\t", index=False)
    pd.DataFrame([crispr]).to_csv(OUT / "crispr_source_verification.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(OUT / "stage0_technical_failures.tsv", sep="\t", index=False)
    independence(ddf, sdf).to_csv(OUT / "dataset_independence_matrix.tsv", sep="\t", index=False)

    g = CFG["gates"]; disc = ddf[~ddf.heldout.astype(bool)]
    n, nh = len(disc), int(disc.human_model.astype(bool).sum())
    nl, nm = int(disc.lung_or_pneumonia.astype(bool).sum()), int(disc.mutant_or_perturbation.astype(bool).sum())
    g0 = n >= g["minimum_independent_datasets"] and nh >= g["minimum_human_datasets"] and nl >= g["minimum_in_vivo_lung_datasets"] and nm >= g["minimum_mutant_datasets"] and bool(hold)
    pt, pm, dom = int(isum.get("unique_physical_human_targets", 0)), int(isum.get("unique_pmids", 0)), float(isum.get("largest_source_fraction", 1))
    g1 = pt >= g["minimum_physical_host_targets"] and pm >= g["minimum_interaction_pmids"] and dom <= g["maximum_single_source_fraction"]
    g2 = bool(crispr.get("verified")); decision = "PASS" if g0 and g1 and g2 else ("REVISE" if g0 or g1 or g2 else "NO_GO")
    gate = {"generated_at": now(), "overall_decision": decision, "gates": {
        "gate0_data": {"status": "PASS" if g0 else "FAIL", "reason": f"datasets={n}; human={nh}; lung={nl}; mutant={nm}; holdout={bool(hold)}"},
        "gate1_interactions": {"status": "PASS" if g1 else "FAIL", "reason": f"physical_targets={pt}; pmids={pm}; source_dominance={dom:.3f}"},
        "gate2_crispr": {"status": "PASS" if g2 else "FAIL", "reason": f"publication_verified={g2}"}},
        "technical_failures": failures,
        "authorization": "Stage 1-2 authorized; holdout locked." if decision == "PASS" else "Downstream analysis not authorized until repair/downgrade."}
    (OUT / "stage0_gate_decision.json").write_text(json.dumps(gate, indent=2, ensure_ascii=False), encoding="utf-8")
    report = ["# Stage 0 Feasibility Audit", "", f"Generated: {now()}", f"Decision: **{decision}**", "", "## Dataset census", "", ddf.to_markdown(index=False), "",
              "## Interaction foundation", ""] + [f"- {k}: {v}" for k, v in isum.items()] + ["", "## CRISPR evidence", "", f"- PMID: {crispr.get('pmid')}", f"- Verified: {crispr.get('verified')}", f"- Title: {crispr.get('title')}", "", "## Gate decisions", ""]
    report += [f"- **{k}: {v['status']}** — {v['reason']}" for k, v in gate["gates"].items()]
    report += ["", "## Authorization", "", gate["authorization"], "", "No differential expression, candidate ranking, or holdout outcome analysis was performed in Stage 0."]
    (OUT / "STAGE0_FEASIBILITY_AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    (CP / ("stage0.done" if decision == "PASS" else "stage0.failed")).write_text(json.dumps(gate, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(gate, indent=2, ensure_ascii=False))
    return 0 if decision == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
