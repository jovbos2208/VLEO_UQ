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


def _check_stat(
    label: str,
    stat: dict,
    *,
    min_count: int,
    max_abs_norm_mean: float,
    min_norm_std: float,
    max_norm_std: float,
    max_norm_rms: float,
) -> tuple[bool, list[str], bool]:
    count = int(stat.get("count", 0))
    if count < min_count:
        return True, [], False
    reasons: list[str] = []
    mean_norm = float(stat.get("mean_norm", 0.0))
    std_norm = float(stat.get("std_norm", 0.0))
    rms_norm = float(stat.get("rms_norm", 0.0))
    if abs(mean_norm) > max_abs_norm_mean:
        reasons.append(f"{label}_mean_norm_exceeds")
    if std_norm < min_norm_std or std_norm > max_norm_std:
        reasons.append(f"{label}_std_norm_out_of_range")
    if rms_norm > max_norm_rms:
        reasons.append(f"{label}_rms_norm_exceeds")
    return (len(reasons) == 0), reasons, True


def _validate_case(
    pod_summary_path: Path,
    *,
    min_count: int,
    max_abs_norm_mean: float,
    min_norm_std: float,
    max_norm_std: float,
    max_norm_rms: float,
) -> dict:
    out = {
        "case": str(pod_summary_path.parent.parent),
        "summary_path": str(pod_summary_path),
        "pass": True,
        "reasons": [],
        "n_arcs": 0,
        "n_checked_stats": 0,
    }
    if not pod_summary_path.exists():
        out["pass"] = False
        out["reasons"].append("missing_pod_summary")
        return out
    payload = _load_json(pod_summary_path)
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
        residuals = arc.get("measurement_residuals")
        if not isinstance(residuals, dict) or len(residuals) == 0:
            out["pass"] = False
            out["reasons"].append(f"arc_{idx}_missing_measurement_residuals")
            continue
        checked_in_arc = 0
        for label, stat in residuals.items():
            if not isinstance(stat, dict):
                continue
            ok, reasons, checked = _check_stat(
                str(label),
                stat,
                min_count=min_count,
                max_abs_norm_mean=max_abs_norm_mean,
                min_norm_std=min_norm_std,
                max_norm_std=max_norm_std,
                max_norm_rms=max_norm_rms,
            )
            if checked:
                checked_in_arc += 1
                out["n_checked_stats"] = int(out["n_checked_stats"]) + 1
            if not ok:
                out["pass"] = False
                out["reasons"].extend([f"arc_{idx}_{r}" for r in reasons])
        if checked_in_arc == 0:
            out["pass"] = False
            out["reasons"].append(f"arc_{idx}_no_stats_above_min_count")
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Validate POD measurement realism residual gates.")
    p.add_argument("--config", required=True)
    p.add_argument("--results_dir", required=True)
    p.add_argument("--out", default=None)
    p.add_argument("--min_count", type=int, default=20)
    p.add_argument("--max_abs_norm_mean", type=float, default=0.35)
    p.add_argument("--min_norm_std", type=float, default=0.5)
    p.add_argument("--max_norm_std", type=float, default=1.8)
    p.add_argument("--max_norm_rms", type=float, default=2.0)
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
            checks.append(
                _validate_case(
                    case_dir / "pod" / "summary.json",
                    min_count=int(args.min_count),
                    max_abs_norm_mean=float(args.max_abs_norm_mean),
                    min_norm_std=float(args.min_norm_std),
                    max_norm_std=float(args.max_norm_std),
                    max_norm_rms=float(args.max_norm_rms),
                )
            )

    overall_pass = all(bool(c.get("pass", False)) for c in checks) if checks else True
    out = {
        "overall_pass": bool(overall_pass),
        "required_cases": int(required_cases),
        "checked_cases": int(len(checks)),
        "thresholds": {
            "min_count": int(args.min_count),
            "max_abs_norm_mean": float(args.max_abs_norm_mean),
            "min_norm_std": float(args.min_norm_std),
            "max_norm_std": float(args.max_norm_std),
            "max_norm_rms": float(args.max_norm_rms),
        },
        "checks": checks,
    }
    out_path = Path(args.out) if args.out else (results_dir / "measurement_realism_gate.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[vg5] wrote {out_path}")
    if required_cases == 0:
        print("[vg5] no pod_uq scenarios in config; treated as pass")
    if not overall_pass:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
