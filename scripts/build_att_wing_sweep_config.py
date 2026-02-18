#!/usr/bin/env python3
"""Build an ATT wing-angle sweep config from a base catalog config."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

RECOMMENDED_WING_ABS_MAX_DEG = 60.0


def parse_angles(text: str) -> list[float]:
    vals: list[float] = []
    for tok in str(text).split(","):
        s = tok.strip()
        if not s:
            continue
        vals.append(float(s))
    if not vals:
        raise ValueError("no valid wing angles provided")
    return vals


def enforce_recommended_range(
    angles_deg: list[float],
    recommended_abs_max_deg: float,
    allow_out_of_recommended_range: bool,
) -> list[float]:
    if allow_out_of_recommended_range:
        return list(angles_deg)
    bad = [a for a in angles_deg if abs(float(a)) > float(recommended_abs_max_deg)]
    if bad:
        bad_s = ",".join(f"{a:g}" for a in bad)
        raise ValueError(
            f"angles exceed recommended range |wing_deg| <= {recommended_abs_max_deg:g}: {bad_s}. "
            "Pass --allow_out_of_recommended_range to override."
        )
    return list(angles_deg)


def sanitize_angle_tag(angle_deg: float) -> str:
    s = f"{angle_deg:.3f}".rstrip("0").rstrip(".")
    return s.replace("-", "m").replace(".", "p")


def scenario_axis(s: dict) -> str:
    axis = str(s.get("maneuver_axis", "")).strip().lower()
    if axis in {"roll", "pitch", "yaw"}:
        return axis
    name = str(s.get("name", "")).lower()
    if "roll" in name:
        return "roll"
    if "pitch" in name:
        return "pitch"
    if "yaw" in name:
        return "yaw"
    raise ValueError(f"cannot infer maneuver axis for scenario '{name}'")


def build_metrics_for_axis(axis: str) -> list[str]:
    return [
        f"time_to_{axis}_90_deg",
        f"time_to_{axis}_180_deg",
        f"peak_{axis}_rate_deg_s",
    ]


def main() -> None:
    default_angles_deg = ",".join(str(a) for a in range(0, 61, 5))
    parser = argparse.ArgumentParser(description="Create A8/A10 wing-angle sweep scenarios.")
    parser.add_argument("--in_config", required=True, help="input JSON config from catalog_to_case_config.py")
    parser.add_argument("--out", required=True, help="output JSON config path")
    parser.add_argument("--angles_deg", default=default_angles_deg, help="comma-separated wing angles in degrees")
    parser.add_argument(
        "--recommended_abs_max_deg",
        type=float,
        default=RECOMMENDED_WING_ABS_MAX_DEG,
        help="recommended absolute wing-angle bound used by default policy",
    )
    parser.add_argument(
        "--allow_out_of_recommended_range",
        action="store_true",
        help="allow angles outside recommended_abs_max_deg",
    )
    parser.add_argument("--duration_s", type=float, default=120.0, help="scenario duration in seconds")
    parser.add_argument("--particles", type=int, default=1000, help="MC particles per scenario")
    parser.add_argument("--dt_s", type=float, default=1.0, help="integration step in seconds")
    parser.add_argument("--eta_max_deg", type=float, default=110.0, help="maximum wing deflection magnitude in degrees")
    parser.add_argument("--eta_rate_max_deg_s", type=float, default=5.0, help="maximum wing rate in deg/s")
    parser.add_argument("--eta_accel_max_deg_s2", type=float, default=1.0, help="maximum wing acceleration in deg/s^2")
    parser.add_argument(
        "--include_names",
        default="att_a8_roll_rotation_baseline,att_a10_yaw_rotation_baseline",
        help="comma-separated base scenario names to include from input config",
    )
    args = parser.parse_args()

    in_path = Path(args.in_config)
    out_path = Path(args.out)
    payload = json.loads(in_path.read_text(encoding="utf-8"))
    base_scenarios = payload.get("scenarios", [])
    if not isinstance(base_scenarios, list) or not base_scenarios:
        raise ValueError(f"no scenarios found in {in_path}")

    include = {s.strip().lower() for s in str(args.include_names).split(",") if s.strip()}
    selected = [s for s in base_scenarios if str(s.get("name", "")).lower() in include]
    if not selected:
        raise ValueError("no matching scenarios found for include_names")

    angles = enforce_recommended_range(
        parse_angles(args.angles_deg),
        recommended_abs_max_deg=float(args.recommended_abs_max_deg),
        allow_out_of_recommended_range=bool(args.allow_out_of_recommended_range),
    )
    out_scenarios: list[dict] = []
    for base in selected:
        axis = scenario_axis(base)
        metrics = build_metrics_for_axis(axis)
        for angle_deg in angles:
            s = copy.deepcopy(base)
            tag = sanitize_angle_tag(angle_deg)
            s["name"] = f"{str(base['name']).lower()}_wing_{tag}deg"
            s["duration_s"] = float(args.duration_s)
            s["particles"] = int(args.particles)
            s["dt_s"] = float(args.dt_s)
            s["wing_constant_deg"] = float(angle_deg)
            s["catalog_output_metrics"] = metrics
            # Keep runtime maneuver command semantics explicit.
            s["maneuver_axis"] = axis
            s["maneuver_half_turn_deg"] = 180.0
            s["maneuver_full_turn_deg"] = 360.0
            s["eta_max_deg"] = float(args.eta_max_deg)
            s["eta_rate_max_deg_s"] = float(args.eta_rate_max_deg_s)
            s["eta_accel_max_deg_s2"] = float(args.eta_accel_max_deg_s2)
            out_scenarios.append(s)

    out = {"scenarios": out_scenarios}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[att-sweep] wrote {len(out_scenarios)} scenarios -> {out_path}")


if __name__ == "__main__":
    main()
