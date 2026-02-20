#!/usr/bin/env python3
"""Compare runs with model discrepancy OFF vs ON and report deltas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _as_float(value: Any) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    return out


def _metric_payload_p95(summary: dict[str, Any]) -> float | None:
    payload = summary.get("payload_impact")
    if isinstance(payload, dict):
        return _as_float(payload.get("final_p95_geolocation_error_m"))
    return None


def _collect_case_metrics(run_dir: Path) -> dict[str, dict[str, float]]:
    out = {}
    for path in sorted(run_dir.rglob("summary.json")):
        if path.parent.name == "pod":
            continue
        try:
            s = _load_json(path)
        except Exception:
            continue
        case = str(s.get("case", path.parent.name))
        vals = {}
        for key in ("ut_mean_err", "stm_mean_err", "ut_cov_rel_err", "stm_cov_rel_err"):
            v = _as_float(s.get(key))
            if v is not None:
                vals[key] = float(v)
        p95 = _metric_payload_p95(s)
        if p95 is not None:
            vals["payload_p95_geolocation_error_m"] = float(p95)
        if vals:
            out[case] = vals
    return out


def run_study(run_off: Path, run_on: Path, cases: list[str] | None = None) -> dict[str, Any]:
    off = _collect_case_metrics(run_off)
    on = _collect_case_metrics(run_on)
    keys = sorted(set(off.keys()) & set(on.keys()))
    if cases:
        allowed = {c.strip() for c in cases if c.strip()}
        keys = [k for k in keys if k in allowed]
    rows = []
    for case in keys:
        m_off = off.get(case, {})
        m_on = on.get(case, {})
        metrics = sorted(set(m_off.keys()) & set(m_on.keys()))
        delta = {}
        for k in metrics:
            delta[k] = float(m_on[k] - m_off[k])
        rows.append({"case": case, "off": m_off, "on": m_on, "delta": delta})
    return {
        "run_off": str(run_off),
        "run_on": str(run_on),
        "paired_case_count": int(len(rows)),
        "cases": rows,
        "pass": bool(len(rows) > 0),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Model discrepancy toggle study report.")
    parser.add_argument("--run_off", required=True, help="Run directory with model discrepancy disabled.")
    parser.add_argument("--run_on", required=True, help="Run directory with model discrepancy enabled.")
    parser.add_argument("--cases", default=None, help="Optional comma-separated case names.")
    parser.add_argument("--out", default="results/reports/model_discrepancy_toggle_study.json")
    args = parser.parse_args()

    run_off = Path(args.run_off)
    run_on = Path(args.run_on)
    if not run_off.is_dir():
        raise FileNotFoundError(run_off)
    if not run_on.is_dir():
        raise FileNotFoundError(run_on)
    cases = None
    if args.cases:
        cases = [x.strip() for x in str(args.cases).split(",") if x.strip()]
    report = run_study(run_off, run_on, cases=cases)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[model-discrepancy-study] wrote {out_path}")
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
