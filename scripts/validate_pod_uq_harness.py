#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _scenario_case_paths(scenario: dict, results_dir: Path) -> list[Path]:
    name = str(scenario.get("name", "")).strip()
    if not name:
        return []
    if str(scenario.get("type", "")).strip().lower() == "formation":
        offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
        return [results_dir / name / f"obj_{idx+1:02d}" for idx, _ in enumerate(offsets)]
    return [results_dir / name]


def _validate_arc_summary(arc: dict) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    for key in ("radial_rms_m", "coverage_sigma", "coverage_3sigma"):
        if key not in arc:
            reasons.append(f"missing_{key}")
    cov = arc.get("coverage_sigma", {})
    if not isinstance(cov, dict):
        reasons.append("invalid_coverage_sigma")
    else:
        for band in ("1sigma", "2sigma", "3sigma"):
            if band not in cov:
                reasons.append(f"missing_coverage_{band}")
    return (len(reasons) == 0), reasons


def _validate_case_pod_summary(path: Path) -> dict:
    out = {
        "case": str(path.parent),
        "summary_path": str(path),
        "pass": True,
        "reasons": [],
        "n_arcs": 0,
    }
    if not path.exists():
        out["pass"] = False
        out["reasons"].append("missing_pod_summary")
        return out
    payload = _load_json(path)
    arcs = payload.get("pod_summaries")
    if not isinstance(arcs, list) or len(arcs) == 0:
        out["pass"] = False
        out["reasons"].append("empty_pod_summaries")
        return out
    out["n_arcs"] = int(len(arcs))
    for idx, arc in enumerate(arcs):
        if not isinstance(arc, dict):
            out["pass"] = False
            out["reasons"].append(f"arc_{idx}_invalid_type")
            continue
        status = str(arc.get("status", "ok"))
        if status.startswith("skipped"):
            continue
        ok, reasons = _validate_arc_summary(arc)
        if not ok:
            out["pass"] = False
            out["reasons"].extend([f"arc_{idx}_{r}" for r in reasons])
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Validate POD UQ harness outputs.")
    p.add_argument("--config", required=True)
    p.add_argument("--results_dir", required=True)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    cfg = _load_json(Path(args.config))
    scenarios = cfg.get("scenarios", [])
    if not isinstance(scenarios, list):
        raise ValueError("config.scenarios must be a list")

    results_dir = Path(args.results_dir)
    checks = []
    required_cases = 0
    for scenario in scenarios:
        if not bool(scenario.get("pod_uq", False)):
            continue
        for case_dir in _scenario_case_paths(scenario, results_dir):
            required_cases += 1
            checks.append(_validate_case_pod_summary(case_dir / "pod" / "summary.json"))

    overall_pass = all(bool(c.get("pass", False)) for c in checks) if checks else True
    out = {
        "overall_pass": bool(overall_pass),
        "required_cases": int(required_cases),
        "checked_cases": int(len(checks)),
        "checks": checks,
    }
    out_path = Path(args.out) if args.out else (results_dir / "pod_uq_harness_validation.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[pod-harness] wrote {out_path}")
    if required_cases == 0:
        print("[pod-harness] no pod_uq scenarios in config; treated as pass")
    if not overall_pass:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
