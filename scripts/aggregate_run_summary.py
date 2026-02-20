#!/usr/bin/env python3
"""Build standardized per-run summary JSON from case summaries."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_case_summary(path: Path) -> bool:
    if path.name != "summary.json":
        return False
    if path.parent.name == "pod":
        return False
    return True


def _extract_case_record(path: Path) -> dict[str, Any]:
    payload = _load_json(path)
    meta = payload.get("metadata", {})
    if not isinstance(meta, dict):
        meta = {}
    case = str(payload.get("case", path.parent.name))
    rec = {
        "case": case,
        "case_dir": str(path.parent),
        "summary_path": str(path),
        "type": payload.get("type"),
        "metadata": {
            "git_hash": meta.get("git_hash"),
            "scenario_id": meta.get("scenario_id"),
            "seed": meta.get("seed"),
            "toggles": meta.get("toggles") if isinstance(meta.get("toggles"), dict) else {},
        },
        "metrics": {
            "ut_mean_err": payload.get("ut_mean_err"),
            "stm_mean_err": payload.get("stm_mean_err"),
            "ut_cov_rel_err": payload.get("ut_cov_rel_err"),
            "stm_cov_rel_err": payload.get("stm_cov_rel_err"),
        },
        "has_pod_summary": bool((path.parent / "pod" / "summary.json").exists()),
        "has_payload_impact": isinstance(payload.get("payload_impact"), dict),
        "has_model_discrepancy": isinstance(payload.get("model_discrepancy"), dict),
    }
    return rec


def build_standardized_run_summary(run_dir: Path) -> dict[str, Any]:
    summaries = [p for p in sorted(run_dir.rglob("summary.json")) if _is_case_summary(p)]
    cases = []
    missing_meta = []
    for path in summaries:
        rec = _extract_case_record(path)
        cases.append(rec)
        meta = rec["metadata"]
        missing = [
            k
            for k in ("git_hash", "scenario_id", "seed", "toggles")
            if meta.get(k) is None
            or (k == "toggles" and not isinstance(meta.get("toggles"), dict))
        ]
        if missing:
            missing_meta.append({"case": rec["case"], "missing": missing})

    out = {
        "schema_version": 1,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "case_count": int(len(cases)),
        "cases": cases,
        "metadata_validation": {
            "missing_metadata_count": int(len(missing_meta)),
            "missing_metadata_cases": missing_meta,
            "pass": bool(len(missing_meta) == 0),
        },
    }
    out["pass"] = bool(out["metadata_validation"]["pass"])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate standardized per-run summary JSON.")
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    summary = build_standardized_run_summary(run_dir)
    out_path = Path(args.out) if args.out else (run_dir / "run_summary.json")
    out_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[run-summary] wrote {out_path}")
    raise SystemExit(0 if bool(summary.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
