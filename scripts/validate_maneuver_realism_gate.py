#!/usr/bin/env python3
"""VG-4 maneuver realism gate: static authority precheck + post-run realism checks."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _as_float(v: Any, default: float | None = None) -> float | None:
    try:
        x = float(v)
    except Exception:
        return default
    if not np.isfinite(x):
        return default
    return x


def _axis_index(axis: str) -> int:
    key = str(axis).strip().lower()
    return {"roll": 0, "pitch": 1, "yaw": 2}[key]


def _axis_eta_cmd(axis: str, wing_deg: float) -> tuple[float, float]:
    key = str(axis).strip().lower()
    d = np.deg2rad(float(wing_deg))
    if key == "roll":
        return float(d), float(-d)
    if key == "pitch":
        return float(d), float(d)
    if key == "yaw":
        return float(d), 0.0
    raise ValueError(f"unsupported maneuver axis '{axis}'")


def _inertia_from_scenario(s: dict) -> np.ndarray:
    vals = s.get("spacecraft_inertia_kgm2")
    if isinstance(vals, (list, tuple)):
        arr = np.asarray(vals, dtype=float).reshape(-1)
        if arr.size == 3 and np.isfinite(arr).all():
            return np.diag(arr)
        if arr.size == 9 and np.isfinite(arr).all():
            return arr.reshape(3, 3)
    return np.diag(np.array([0.15, 0.12, 0.20], dtype=float))


def _is_maneuver_attitude_case(s: dict) -> bool:
    if str(s.get("type", "")).lower() != "attitude":
        return False
    axis = str(s.get("maneuver_axis", "")).strip().lower()
    wing = _as_float(s.get("wing_constant_deg"))
    return axis in {"roll", "pitch", "yaw"} and wing is not None and wing != 0.0


def _authority_eval(
    *,
    max_abs_torque_axis_nm: float,
    inertia_axis_kgm2: float,
    full_turn_deg: float,
    duration_s: float,
    min_authority_ratio: float,
) -> dict:
    theta = np.deg2rad(max(0.0, float(full_turn_deg)))
    T = max(1.0, float(duration_s))
    alpha_req = 4.0 * theta / (T * T)
    torque_req = abs(float(inertia_axis_kgm2)) * alpha_req
    ratio = float(max_abs_torque_axis_nm / max(torque_req, 1e-15))
    passed = bool(ratio >= float(min_authority_ratio))
    return {
        "required_torque_axis_nm": float(torque_req),
        "max_abs_torque_axis_nm": float(max_abs_torque_axis_nm),
        "authority_ratio": ratio,
        "min_authority_ratio": float(min_authority_ratio),
        "pass": passed,
    }


def _evaluate_post_stats(
    *,
    stats: dict,
    max_saturation_fraction: float,
    max_rate_limit_ratio: float,
    min_reached_fraction_180: float,
    require_full_turn: bool,
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    sat = _as_float(stats.get("saturation_time_fraction"))
    if sat is not None and sat > max_saturation_fraction:
        reasons.append(
            f"saturation_time_fraction={sat:.3f} > max={max_saturation_fraction:.3f}"
        )

    cmd_rate = _as_float(stats.get("command_rate_max_abs_deg_s"))
    rate_max = _as_float(stats.get("eta_rate_max_deg_s"))
    if cmd_rate is not None and rate_max is not None and rate_max > 0.0:
        if cmd_rate > rate_max * max_rate_limit_ratio:
            reasons.append(
                f"command_rate_max_abs_deg_s={cmd_rate:.3f} exceeds "
                f"limit={rate_max:.3f} * {max_rate_limit_ratio:.3f}"
            )

    frac180 = _as_float(stats.get("reached_fraction_180"))
    if frac180 is not None and frac180 < min_reached_fraction_180:
        reasons.append(
            f"reached_fraction_180={frac180:.3f} < min={min_reached_fraction_180:.3f}"
        )

    if require_full_turn:
        reached = stats.get("full_turn_reached")
        if reached is False:
            reasons.append("full_turn_reached=false")

    return (len(reasons) == 0), reasons


def _run_pre_gate_for_scenario(s: dict, *, min_authority_ratio: float, map_points: int) -> dict:
    from vleo_uq import AeroAdapter, PropagatorConfig, default_geometry
    from scripts.run_case_studies import build_initial_state, scenario_initial_w0_rad_s

    axis = str(s["maneuver_axis"]).strip().lower()
    wing = float(s["wing_constant_deg"])
    eta_max = _as_float(s.get("eta_max_deg"), default=110.0)
    full_turn_deg = _as_float(s.get("maneuver_full_turn_deg"), default=360.0)
    duration_s = _as_float(s.get("duration_s"), default=3600.0)
    density_scale = _as_float(s.get("density_scale"), default=1.0)

    cfg = PropagatorConfig()
    x0 = build_initial_state(
        cfg.mu_earth_m3_s2,
        18,
        w0=scenario_initial_w0_rad_s(s),
        orbit=s.get("orbit"),
        attitude_reference=str(s.get("attitude_reference", "flow_x")),
        initial_attitude_euler_deg_zyx=s.get("initial_attitude_euler_deg_zyx"),
    )
    q_wxyz = np.asarray(x0[6:10], dtype=float)
    w_b = np.asarray(x0[10:13], dtype=float)
    v_i = np.asarray(x0[3:6], dtype=float)

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)

    n_points = int(max(5, map_points))
    cmds = np.linspace(-abs(float(eta_max)), abs(float(eta_max)), n_points)
    torque_axis = np.zeros(n_points, dtype=float)
    for i, cmd in enumerate(cmds):
        eta1, eta2 = _axis_eta_cmd(axis, float(cmd))
        _force, tau = aero.compute_ft(
            q_wxyz,
            w_b,
            v_i,
            np.zeros(3),
            float(1e-12 * float(density_scale)),
            1000.0,
            28.0 * 1.6605390689252e-27,
            float(eta1),
            float(eta2),
            1,
        )
        torque_axis[i] = float(np.asarray(tau, dtype=float)[_axis_index(axis)])

    max_abs_tau = float(np.max(np.abs(torque_axis)))
    inertia_axis = float(_inertia_from_scenario(s)[_axis_index(axis), _axis_index(axis)])
    eval_out = _authority_eval(
        max_abs_torque_axis_nm=max_abs_tau,
        inertia_axis_kgm2=inertia_axis,
        full_turn_deg=float(full_turn_deg),
        duration_s=float(duration_s),
        min_authority_ratio=float(min_authority_ratio),
    )
    has_pos = bool(np.any(torque_axis > 0.0))
    has_neg = bool(np.any(torque_axis < 0.0))
    return {
        "scenario": str(s.get("name", "unknown")),
        "axis": axis,
        "wing_constant_deg": float(wing),
        "eta_max_deg": float(eta_max),
        "map_cmd_deg": [float(x) for x in cmds],
        "map_torque_axis_nm": [float(x) for x in torque_axis],
        "inertia_axis_kgm2": inertia_axis,
        "has_positive_torque": has_pos,
        "has_negative_torque": has_neg,
        **eval_out,
    }


def _run_post_gate_for_scenario(
    s: dict,
    *,
    results_dir: Path,
    pre_result_by_scenario: dict[str, dict] | None,
    max_saturation_fraction: float,
    max_rate_limit_ratio: float,
    min_reached_fraction_180: float,
    require_full_turn: bool,
    require_pre_authority_pass: bool,
) -> dict:
    name = str(s.get("name", "unknown"))
    summary_path = results_dir / name / "summary.json"
    if not summary_path.exists():
        return {"scenario": name, "pass": False, "reasons": ["missing_summary_json"]}

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    stats = summary.get("maneuver_actuator_stats")
    if not isinstance(stats, dict):
        return {"scenario": name, "pass": False, "reasons": ["missing_maneuver_actuator_stats"]}

    ok, reasons = _evaluate_post_stats(
        stats=stats,
        max_saturation_fraction=max_saturation_fraction,
        max_rate_limit_ratio=max_rate_limit_ratio,
        min_reached_fraction_180=min_reached_fraction_180,
        require_full_turn=require_full_turn,
    )

    if require_pre_authority_pass:
        pre = (pre_result_by_scenario or {}).get(name)
        if not isinstance(pre, dict):
            ok = False
            reasons.append("missing_pre_authority_result")
        elif not bool(pre.get("pass", False)):
            ok = False
            reasons.append("pre_authority_gate_failed")

    return {"scenario": name, "pass": bool(ok), "reasons": reasons, "stats": stats}


def _filter_scenarios(config: dict, scenario_names: set[str] | None) -> list[dict]:
    scenarios = [s for s in config.get("scenarios", []) if isinstance(s, dict)]
    if scenario_names:
        scenarios = [s for s in scenarios if str(s.get("name", "")) in scenario_names]
    return [s for s in scenarios if _is_maneuver_attitude_case(s)]


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate maneuver realism gate (VG-4).")
    parser.add_argument("--config", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--stage", choices=("pre", "post", "both"), default="both")
    parser.add_argument("--scenarios", default=None, help="Optional comma-separated scenario names.")
    parser.add_argument("--map_points", type=int, default=25)
    parser.add_argument("--min_authority_ratio", type=float, default=0.0)
    parser.add_argument("--min_static_torque_nm", type=float, default=1e-12)
    parser.add_argument("--max_saturation_fraction", type=float, default=0.85)
    parser.add_argument("--max_rate_limit_ratio", type=float, default=1.05)
    parser.add_argument("--min_reached_fraction_180", type=float, default=0.5)
    parser.add_argument("--no_require_full_turn", action="store_true")
    parser.add_argument("--no_require_pre_authority_pass", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenario_names = None
    if args.scenarios:
        scenario_names = {x.strip() for x in args.scenarios.split(",") if x.strip()}
    scenarios = _filter_scenarios(cfg, scenario_names)

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    pre_report_path = results_dir / "maneuver_gate_pre.json"
    post_report_path = results_dir / "maneuver_gate_post.json"

    pre_rows: list[dict] = []
    if args.stage in {"pre", "both"}:
        for s in scenarios:
            row = _run_pre_gate_for_scenario(
                s,
                min_authority_ratio=float(args.min_authority_ratio),
                map_points=int(args.map_points),
            )
            row["pass"] = bool(
                row.get("pass", False)
                and bool(row.get("has_positive_torque", False))
                and bool(row.get("has_negative_torque", False))
                and float(row.get("max_abs_torque_axis_nm", 0.0)) >= float(args.min_static_torque_nm)
            )
            pre_rows.append(row)
        pre_payload = {
            "stage": "pre",
            "count": len(pre_rows),
            "all_pass": bool(all(bool(r.get("pass", False)) for r in pre_rows)),
            "results": pre_rows,
        }
        pre_report_path.write_text(json.dumps(pre_payload, indent=2), encoding="utf-8")
        print(f"[maneuver-gate] wrote {pre_report_path}")
    else:
        if pre_report_path.exists():
            pre_payload = json.loads(pre_report_path.read_text(encoding="utf-8"))
            pre_rows = list(pre_payload.get("results", []))

    pre_by_name = {str(r.get("scenario")): r for r in pre_rows if isinstance(r, dict)}

    post_rows: list[dict] = []
    if args.stage in {"post", "both"}:
        for s in scenarios:
            row = _run_post_gate_for_scenario(
                s,
                results_dir=results_dir,
                pre_result_by_scenario=pre_by_name,
                max_saturation_fraction=float(args.max_saturation_fraction),
                max_rate_limit_ratio=float(args.max_rate_limit_ratio),
                min_reached_fraction_180=float(args.min_reached_fraction_180),
                require_full_turn=not bool(args.no_require_full_turn),
                require_pre_authority_pass=not bool(args.no_require_pre_authority_pass),
            )
            post_rows.append(row)
        post_payload = {
            "stage": "post",
            "count": len(post_rows),
            "all_pass": bool(all(bool(r.get("pass", False)) for r in post_rows)),
            "results": post_rows,
        }
        post_report_path.write_text(json.dumps(post_payload, indent=2), encoding="utf-8")
        print(f"[maneuver-gate] wrote {post_report_path}")

    payload = {
        "stage": args.stage,
        "scenarios_checked": len(scenarios),
        "pre_report": str(pre_report_path),
        "post_report": str(post_report_path),
        "pre_all_pass": bool(all(bool(r.get("pass", False)) for r in pre_rows)) if pre_rows else True,
        "post_all_pass": bool(all(bool(r.get("pass", False)) for r in post_rows)) if post_rows else True,
    }
    payload["pass"] = bool(payload["pre_all_pass"] and payload["post_all_pass"])

    out_path = Path(args.out) if args.out else (results_dir / "maneuver_gate_summary.json")
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[maneuver-gate] wrote {out_path}")
    print(json.dumps(payload, indent=2))
    raise SystemExit(0 if payload["pass"] else 1)


if __name__ == "__main__":
    main()
