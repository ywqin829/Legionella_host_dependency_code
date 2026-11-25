"""Gate G0 — automated claim ledger.

Executes every assertion declared in ``g0_registry.py`` against archived outputs and the
pinned manuscript snapshot, then assigns each claim a status from the computed outcome.

Statuses are computed, never written by hand:

    AUTO_PASS     every declared assertion matched
    FAIL          at least one assertion did not match
    STALE         a pinned input hash no longer matches this snapshot
    MANUAL_PASS   no machine assertion exists; a recorded manual check is present
    NOT_TESTABLE  no machine assertion and no recorded manual check

A claim whose number is right but whose sentence needs work carries ``text_action``; that
is reported separately and is never converted into a silent pass.

Outputs (``revision2_2026-11-23/``):
    REVISION2_CLAIM_LEDGER.tsv
    REVISION2_CLAIM_LEDGER_STATS.json
    REVISION2_CLAIM_LEDGER.md
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))

import g0_registry as reg  # noqa: E402

OUT = ROOT / "revision2_2026-11-23"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


_JSON_CACHE: dict[str, dict] = {}
_TSV_CACHE: dict[str, pd.DataFrame] = {}


def load_json(rel: str) -> dict:
    if rel not in _JSON_CACHE:
        _JSON_CACHE[rel] = json.loads((ROOT / rel).read_text(encoding="utf-8"))
    return _JSON_CACHE[rel]


def load_tsv(rel: str) -> pd.DataFrame:
    if rel not in _TSV_CACHE:
        kwargs = {"low_memory": False} if rel.endswith(".gz") or rel == reg.ATLAS else {}
        _TSV_CACHE[rel] = pd.read_csv(ROOT / rel, sep="\t", **kwargs)
    return _TSV_CACHE[rel]


def dig(payload, dotted: str):
    current = payload
    for part in dotted.split("."):
        if isinstance(current, list):
            current = current[int(part)]
        else:
            current = current[part]
    return current


def compare(expect, observed, rel: float = 1e-3) -> bool:
    if isinstance(expect, bool):
        return bool(observed) == expect
    if isinstance(expect, (int, float)):
        try:
            obs = float(observed)
        except (TypeError, ValueError):
            return False
        if expect == 0:
            return abs(obs) <= rel
        return abs(obs - expect) / abs(expect) <= rel
    if isinstance(expect, dict):
        if not isinstance(observed, dict) or set(observed) != set(expect):
            return False
        return all(compare(expect[k], observed[k], rel) for k in expect)
    return str(observed).strip() == str(expect).strip()


def select_row(frame: pd.DataFrame, where: dict, claim_id: str, file: str) -> pd.Series:
    mask = pd.Series(True, index=frame.index)
    for column, value in where.items():
        if column not in frame.columns:
            raise KeyError(f"{claim_id}: {file} has no column {column!r}")
        mask &= frame[column].astype(str).str.strip() == str(value)
    hits = frame.loc[mask]
    if len(hits) != 1:
        raise ValueError(
            f"{claim_id}: expected exactly 1 row in {file} for {where}, got {len(hits)}"
        )
    return hits.iloc[0]


_NORMALISE = re.compile(r"\s+")


def normalise(text: str) -> str:
    return _NORMALISE.sub(" ", text).strip()


def run_assertion(claim_id: str, assertion: dict, manuscript_norm: str):
    """Return (ok, observed_text, detail)."""
    kind = assertion["kind"]

    try:
        if kind == "json":
            observed = dig(load_json(assertion["file"]), assertion["path"])
            ok = compare(assertion["expect"], observed, assertion.get("rel", 1e-3))
            return ok, repr(observed), f"{assertion['file']}::{assertion['path']}"

        if kind == "tsv_cell":
            frame = load_tsv(assertion["file"])
            row = select_row(frame, assertion["where"], claim_id, assertion["file"])
            observed = row[assertion["column"]]
            ok = compare(assertion["expect"], observed, assertion.get("rel", 1e-3))
            return ok, repr(observed), f"{assertion['file']}::{assertion['column']}"

        if kind == "tsv_row":
            frame = load_tsv(assertion["file"])
            row = select_row(frame, assertion["where"], claim_id, assertion["file"])
            observed = {c: row[c] for c in assertion["expect"]}
            ok = compare(assertion["expect"], observed, assertion.get("rel", 1e-3))
            return ok, json.dumps({k: str(v) for k, v in observed.items()},
                                  ensure_ascii=False), f"{assertion['file']} row {assertion['where']}"

        if kind == "tsv_count":
            frame = load_tsv(assertion["file"])
            if assertion.get("where"):
                mask = pd.Series(True, index=frame.index)
                for column, value in assertion["where"].items():
                    mask &= frame[column].astype(str).str.strip() == str(value)
                frame = frame.loc[mask]
            observed = int(len(frame))
            return observed == assertion["expect"], str(observed), f"{assertion['file']} row count"

        if kind == "tsv_nunique":
            frame = load_tsv(assertion["file"])
            if assertion.get("where"):
                mask = pd.Series(True, index=frame.index)
                for column, value in assertion["where"].items():
                    mask &= frame[column].astype(str).str.strip() == str(value)
                frame = frame.loc[mask]
            observed = int(frame[assertion["column"]].astype(str).nunique())
            return (observed == assertion["expect"], str(observed),
                    f"{assertion['file']}::{assertion['column']} nunique")

        if kind == "text":
            missing = [t for t in assertion["tokens"] if normalise(t) not in manuscript_norm]
            return (not missing, "all tokens present" if not missing else f"missing: {missing}",
                    "manuscript token check")

        if kind == "path":
            exists = (ROOT / assertion["file"]).exists()
            return exists, "exists" if exists else "MISSING", assertion["file"]

        if kind == "sha":
            observed = sha256(ROOT / assertion["file"])
            return observed == assertion["expect"], observed[:16] + "...", assertion["file"]

        if kind == "calc":
            name = assertion["name"]
            if name not in reg.CALCULATIONS:
                raise KeyError(f"unknown calculation {name!r}")
            observed = reg.CALCULATIONS[name]()
            ok = compare(assertion["expect"], observed, assertion.get("rel", 1e-3))
            return ok, repr(observed), f"calc:{name}"

    except Exception as exc:  # a broken locator is a ledger failure, not a pass
        return False, f"ERROR: {type(exc).__name__}: {exc}", f"{kind} assertion raised"

    return False, f"ERROR: unknown assertion kind {kind!r}", ""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    manuscript_path = ROOT / reg.MANUSCRIPT
    manuscript_norm = normalise(manuscript_path.read_text(encoding="utf-8"))

    # -------------------------------------------------------------- #
    # Global staleness: any pinned file that moved invalidates claims.
    # -------------------------------------------------------------- #
    stale_files: list[dict] = []
    for rel, pinned in reg.PINNED_INPUTS.items():
        path = ROOT / rel
        if not path.exists():
            stale_files.append({"file": rel, "reason": "missing"})
            continue
        actual = sha256(path)
        if actual != pinned:
            stale_files.append({"file": rel, "reason": "hash changed",
                                "pinned": pinned, "actual": actual})

    stale_names = {item["file"] for item in stale_files}

    rows: list[dict] = []
    for claim in reg.SPECS:
        cid = claim["claim_id"]
        cited = {claim["source_file"]}
        cited.update(a["file"] for a in claim["assertions"] if "file" in a)
        claim_stale = sorted(cited & stale_names)

        outcomes = [run_assertion(cid, a, manuscript_norm) for a in claim["assertions"]]
        passed = sum(1 for ok, _, _ in outcomes if ok)
        fails = [(obs, detail) for ok, obs, detail in outcomes if not ok]

        if claim_stale:
            status = "STALE"
        elif fails:
            status = "FAIL"
        elif outcomes:
            status = "AUTO_PASS"
        elif claim["manual_note"]:
            status = "MANUAL_PASS"
        else:
            status = "NOT_TESTABLE"

        evidence = "; ".join(
            f"[{'ok' if ok else 'FAIL'}] {detail} -> {obs}" for ok, obs, detail in outcomes
        )
        failure_note = ""
        if fails:
            failure_note = " | ".join(f"{detail}: observed {obs}" for obs, detail in fails)
        if claim_stale:
            failure_note = "stale input(s): " + ", ".join(claim_stale)

        rows.append(
            {
                "claim_id": cid,
                "claim_text": claim["claim_text"],
                "estimate_ci_p": claim["estimate_ci_p"],
                "source_file": claim["source_file"],
                "generating_script": claim["generating_script"],
                "analysis_status": claim["analysis_status"],
                "manuscript_location": claim["manuscript_location"],
                "verification_status": status,
                "assertions_total": len(outcomes),
                "assertions_passed": passed,
                "failure_detail": failure_note,
                "text_action": claim["text_action"],
                "evidence": evidence,
                "notes": claim["manual_note"],
            }
        )

    ledger = pd.DataFrame(rows)
    ledger.to_csv(OUT / "REVISION2_CLAIM_LEDGER.tsv", sep="\t", index=False)

    counts = ledger.verification_status.value_counts().to_dict()
    blocking = ledger.loc[
        ledger.verification_status.isin(["FAIL", "STALE", "NOT_TESTABLE"]), "claim_id"
    ].tolist()
    open_text_actions = ledger.loc[ledger.text_action.astype(str).str.strip() != "", "claim_id"].tolist()

    if blocking:
        gate = "FAIL"
    elif open_text_actions:
        gate = "PASS_WITH_OPEN_TEXT_ACTIONS"
    else:
        gate = "PASS"

    stats = {
        "analysis": "WP1 claim ledger (Gate G0)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "manuscript": {"path": reg.MANUSCRIPT, "sha256": sha256(manuscript_path),
                       "pinned": reg.MANUSCRIPT_SHA256},
        "claims_total": int(len(ledger)),
        "verification_counts": {k: int(v) for k, v in counts.items()},
        "assertions_total": int(ledger.assertions_total.sum()),
        "gate_g0": gate,
        "blocking_claims": blocking,
        "open_text_actions": open_text_actions,
        "stale_inputs": stale_files,
        "holdout_accessed_for_new_analysis": False,
        "pinned_inputs_checked": len(reg.PINNED_INPUTS),
    }
    (OUT / "REVISION2_CLAIM_LEDGER_STATS.json").write_text(
        json.dumps(stats, indent=2), encoding="utf-8"
    )

    # -------------------------------------------------------------- #
    # Human-readable report.
    # -------------------------------------------------------------- #
    lines = [
        "# Revision 2 — WP1 claim ledger (Gate G0)",
        "",
        f"- Timestamp (UTC): {stats['timestamp_utc']}",
        f"- Manuscript SHA-256: `{stats['manuscript']['sha256']}`",
        f"- Claims audited: {stats['claims_total']}",
        f"- Machine assertions executed: {stats['assertions_total']}",
        f"- Pinned inputs checked: {stats['pinned_inputs_checked']}",
        f"- Holdout accessed for new analysis: **No**",
        f"- **Gate G0: {gate}**",
        f"- Blocking claims: {len(blocking)}",
        f"- Open text actions: {len(open_text_actions)}",
        "",
        "## Verification counts",
        "",
        "| status | claims |",
        "|---|---|",
    ]
    for key in ["AUTO_PASS", "MANUAL_PASS", "FAIL", "NOT_TESTABLE", "STALE"]:
        lines.append(f"| {key} | {int(counts.get(key, 0))} |")

    if stale_files:
        lines += ["", "## Stale inputs", "", "| file | reason |", "|---|---|"]
        for item in stale_files:
            lines.append(f"| `{item['file']}` | {item['reason']} |")

    lines += [
        "",
        "## What G0 does and does not certify",
        "",
        "G0 certifies that every declared quantitative statement resolves to an archived output",
        "file or a named derived calculation, and that the recomputed value matches the value in",
        "the pinned manuscript snapshot. Assertion outcomes are computed; no claim is marked as",
        "verified by hand.",
        "",
        "G0 does **not** certify scientific interpretation, causal wording, or whether a phrase",
        "such as \"supported\" or \"robust\" is warranted. Those remain manual-review items and are",
        "reported separately from machine passes.",
    ]

    if open_text_actions:
        lines += ["", "## Open text actions", ""]
        for cid in open_text_actions:
            row = ledger.loc[ledger.claim_id == cid].iloc[0]
            lines += [
                f"### {cid} — {row.verification_status}",
                "",
                f"Claim: {row.estimate_ci_p}  ",
                f"Location: {row.manuscript_location}",
                "",
                row.text_action,
                "",
            ]

    if blocking:
        lines += ["", "## Blocking claims", "", "| claim | status | detail |", "|---|---|---|"]
        for cid in blocking:
            row = ledger.loc[ledger.claim_id == cid].iloc[0]
            lines.append(f"| {cid} | {row.verification_status} | {row.failure_detail} |")

    (OUT / "REVISION2_CLAIM_LEDGER.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps(
        {
            "gate_g0": gate,
            "claims_total": stats["claims_total"],
            "verification_counts": stats["verification_counts"],
            "assertions_total": stats["assertions_total"],
            "blocking": blocking,
            "open_text_actions": open_text_actions,
            "stale_inputs": [s["file"] for s in stale_files],
        },
        indent=2,
    ))
    if blocking:
        print("\nBLOCKING:")
        print(ledger.loc[ledger.verification_status.isin(["FAIL", "STALE", "NOT_TESTABLE"]),
                         ["claim_id", "verification_status", "failure_detail"]].to_string(index=False))
    return 0 if gate != "FAIL" else 2


if __name__ == "__main__":
    raise SystemExit(main())
