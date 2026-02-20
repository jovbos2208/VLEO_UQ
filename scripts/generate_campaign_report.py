#!/usr/bin/env python3
"""Generate campaign completion report from run outputs.

This script supports both:
- legacy phase-based runs (`phase_a`/`phase_b`/`phase_c`), and
- catalog-area runs (`att`/`mis`/`orb`/`aero`) produced by newer workflows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


LEGACY_PHASE_SCENARIOS = {
    "phase_a": [
        "att_a1_detumble_aero_only",
        "att_a2_nadir_pointing",
        "od_c1_gnss_continuous_bestcase",
        "od_c2_gnss_realistic_outages",
        "od_c3_gnss_plus_accel",
    ],
    "phase_b": [
        "att_a7_eclipse_thermal_transition",
        "od_c9_meas_error_stress_tests",
        "od_c8_station_weather_availability",
        "form_b6_robust_formation_storm",
    ],
    "phase_c": [
        "form_b2_alongtrack_2sat_dd_keeping",
        "form_b3_lvlh_box_constraint_3sat",
        "form_b5_close_approach_risk",
        "form_b7_cluster_sparse_ground_contact",
    ],
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _case_present(run_dir: Path, case_name: str) -> tuple[bool, dict[str, Any] | None]:
    summary_path = run_dir / case_name / "summary.json"
    if not summary_path.exists():
        return False, None
    try:
        payload = _load_json(summary_path)
    except Exception:
        return False, None
    return True, payload


def _collect_groups_from_configs(run_dir: Path) -> dict[str, list[str]]:
    cfg_dir = run_dir / "configs"
    if not cfg_dir.is_dir():
        return {}
    groups: dict[str, list[str]] = {}
    for cfg_path in sorted(cfg_dir.glob("*.json")):
        try:
            payload = _load_json(cfg_path)
        except Exception:
            continue
        scenarios = payload.get("scenarios", [])
        if not isinstance(scenarios, list):
            continue
        for s in scenarios:
            if not isinstance(s, dict):
                continue
            case_name = str(s.get("name", "")).strip()
            if not case_name:
                continue
            area = str(s.get("catalog_area", "")).strip().upper()
            if area:
                group = area.lower()
            else:
                sid = str(s.get("catalog_scenario_id", "")).strip().upper()
                prefix = sid.split("_", 1)[0].lower() if sid else "unknown"
                group = prefix
            groups.setdefault(group, [])
            if case_name not in groups[group]:
                groups[group].append(case_name)
    return groups


def _build_group_report(run_dir: Path, groups: dict[str, list[str]]) -> tuple[dict[str, Any], bool]:
    out: dict[str, Any] = {}
    all_ok = True
    for group, cases in groups.items():
        rows = []
        ok_count = 0
        for case in cases:
            present, payload = _case_present(run_dir, case)
            row = {
                "case": case,
                "present": bool(present),
                "has_pod": bool((run_dir / case / "pod" / "summary.json").exists()),
            }
            if isinstance(payload, dict):
                row["ut_mean_err"] = payload.get("ut_mean_err")
                row["stm_mean_err"] = payload.get("stm_mean_err")
                row["ut_cov_rel_err"] = payload.get("ut_cov_rel_err")
                row["stm_cov_rel_err"] = payload.get("stm_cov_rel_err")
            rows.append(row)
            if present:
                ok_count += 1
        group_pass = bool(ok_count == len(cases))
        all_ok = bool(all_ok and group_pass)
        out[group] = {
            "required_cases": int(len(cases)),
            "present_cases": int(ok_count),
            "pass": group_pass,
            "cases": rows,
        }
    return out, all_ok


def build_campaign_report(run_dir: Path) -> dict[str, Any]:
    groups = _collect_groups_from_configs(run_dir)
    mode = "config_groups"
    if not groups:
        mode = "legacy_phases"
        groups = LEGACY_PHASE_SCENARIOS
    groups_out, all_ok = _build_group_report(run_dir, groups)
    report: dict[str, Any] = {
        "run_dir": str(run_dir),
        "mode": mode,
        "groups": groups_out,
        "pass": bool(all_ok),
    }
    if mode == "legacy_phases":
        report["phases"] = groups_out
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate campaign report.")
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    report = build_campaign_report(run_dir)
    out_path = Path(args.out) if args.out else (run_dir / "campaign_report.json")
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[campaign-report] wrote {out_path}")
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
