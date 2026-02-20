#!/usr/bin/env python3
"""Convert scenario_catalog.md YAML blocks into run_case_studies JSON config."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def extract_yaml_blocks(text: str) -> list[str]:
    pattern = re.compile(r"```yaml\s*\n(.*?)```", re.DOTALL)
    return [m.group(1).strip() for m in pattern.finditer(text)]


def map_type(area: str) -> str:
    area = (area or "").strip().upper()
    if area == "ATT":
        return "attitude"
    if area in {"FORM", "MIS"}:
        return "formation"
    if area in {"OD", "ORB", "AERO"}:
        return "mission"
    raise ValueError(f"unsupported area '{area}'")


def infer_type_from_controls(area: str, control_toggles: set[str]) -> str:
    keys = {str(k).strip().upper() for k in control_toggles if str(k).strip()}
    if any(k.startswith("CTL_DD_") for k in keys):
        return "formation"
    if any(k.startswith("CTL_ATT_") for k in keys):
        return "attitude"
    return map_type(area)


def extract_toggles(section: dict | None, *keys: str) -> set[str]:
    if not isinstance(section, dict):
        return set()
    out: set[str] = set()
    for key in keys:
        raw = section.get(key)
        if raw is None:
            continue
        if isinstance(raw, str):
            s = raw.strip()
            if s:
                out.add(s)
            continue
        if isinstance(raw, (list, tuple, set)):
            for item in raw:
                s = str(item).strip()
                if s:
                    out.add(s)
    return out


def choose_primary_control_toggle(control_toggles: set[str]) -> str:
    if not control_toggles:
        return ""
    ordered = [
        "CTL_DD_REPHASE",
        "CTL_DD_ALONGTRACK",
        "CTL_ATT_DRAG_MIN",
        "CTL_ATT_AERO_POINT",
        "CTL_ATT_AERO_RATE",
        "CTL_ORB_THRUST_SMA",
        "CTL_ORB_DRAG_COMP",
        "CTL_NONE",
    ]
    keys = {str(k).strip() for k in control_toggles if str(k).strip()}
    for k in ordered:
        if k in keys:
            return k
    return sorted(keys)[0]


def choose_dt_seconds(block: dict, scenario_type: str, dt_floor: float) -> float:
    vals = []
    for key in ("dt_control_s", "dt_meas_s"):
        val = block.get(key)
        if val is None:
            continue
        try:
            fv = float(val)
        except Exception:
            continue
        if fv > 0.0:
            vals.append(fv)
    if vals:
        dt = min(vals)
    else:
        dt = 5.0 if scenario_type == "attitude" else 60.0
    return max(float(dt_floor), float(dt))


def map_env_to_overrides(env_toggles: set[str]) -> tuple[dict, bool, float]:
    overrides = {}
    augment = False
    density_scale = 1.0
    if "ENV_Q2_SCALE_UP" in env_toggles:
        density_scale = 1.5
    if "ENV_Q3_SCALE_DOWN" in env_toggles:
        density_scale = 0.7
    if "ENV_S1_OU_LOGRHO" in env_toggles:
        overrides["rho_fast_tau_s"] = 900.0
        overrides["rho_fast_sigma"] = 0.15
        augment = True
    if "ENV_W2_WIND_OU_BIAS" in env_toggles:
        overrides["wind_tau_s"] = 900.0
        overrides["wind_sigma"] = 10.0
        augment = True
    if "ENV_S2_STORM_PULSE" in env_toggles:
        overrides["rho_bias_tau_s"] = 7200.0
        overrides["rho_bias_sigma"] = 0.08
        augment = True
    return overrides, augment, density_scale


def default_formation_offsets(scenario_id: str) -> list[list[float]]:
    sid = scenario_id.upper()
    if "SINGLE" in sid:
        return [[0.0, 0.0, 0.0]]
    if "2SAT" in sid:
        return [[0.0, 0.0, 0.0], [120.0, 0.0, 0.0]]
    if "CLUSTER" in sid:
        return [[0.0, 0.0, 0.0], [120.0, 0.0, 0.0], [240.0, 0.0, 0.0], [360.0, 0.0, 0.0]]
    return [[0.0, 0.0, 0.0], [120.0, 0.0, 0.0], [240.0, 0.0, 0.0]]


def maybe_rephase_event(control_toggles: set[str], duration_s: float) -> list[dict]:
    if "CTL_DD_REPHASE" not in {str(k).strip() for k in control_toggles}:
        return []
    duration = float(duration_s)
    if duration <= 0.0:
        return []
    return [
        {
            "t_s": 0.5 * duration,
            "kind": "delta_v",
            "params": {"dv_eci_m_s": [0.0, 0.02, 0.0], "sigma_dv_m_s": 0.002},
        }
    ]


def convert_block(
    block: dict,
    *,
    default_particles: int,
    duration_scale: float,
    start_utc: str | None,
    use_env_sources: bool,
    dt_floor: float,
) -> dict:
    scenario_id = str(block["scenario_id"]).strip()
    area = str(block.get("area", "")).strip()
    control = block.get("control", {}) or {}
    control_toggles = extract_toggles(control, "toggles", "control_toggle")
    scenario_type = infer_type_from_controls(area, control_toggles)
    duration_s = float(block.get("duration_s", 3600.0)) * float(duration_scale)
    dt_s = choose_dt_seconds(block, scenario_type, dt_floor)

    uq = block.get("uq", {}) or {}
    method = str(uq.get("propagator", uq.get("method", "PROP_MC"))).upper()
    particles = int(uq.get("n_mc", default_particles))
    if particles < 2:
        particles = max(2, default_particles)

    env = block.get("environment", {}) or {}
    env_toggles = extract_toggles(env, "toggles", "env_toggle")
    overrides, augment_process_noise, density_scale = map_env_to_overrides(env_toggles)
    aero = block.get("aero", {}) or {}
    gsi_toggles = extract_toggles(aero, "toggles", "gsi_toggle")

    sensors = block.get("sensors", {}) or {}
    sensor_toggles = extract_toggles(sensors, "toggles", "sensor_toggle")
    ground = block.get("ground", {}) or {}
    ground_toggles = extract_toggles(ground, "toggles", "ground_toggle")
    area_upper = area.strip().upper()
    has_gnss = ("SENS_GNSS_RAW_DUAL" in sensor_toggles) or ("SENS_GNSS_RAW_SINGLE" in sensor_toggles)
    has_slr = "SENS_SLR" in sensor_toggles
    pod_uq = bool(
        (area_upper == "OD")
        or has_gnss
        or has_slr
    )
    pod_skip_slr = not has_slr
    pod_disable_gnss = has_slr and (not has_gnss)
    estimation = block.get("estimation", {}) or {}
    est_toggles = extract_toggles(estimation, "toggles", "est_toggle")
    pod_estimator = "enkf" if "EST_ENKF" in est_toggles else "batch"
    control_toggle = choose_primary_control_toggle(control_toggles)

    out = {
        "name": scenario_id.lower(),
        "catalog_scenario_id": scenario_id,
        "type": scenario_type,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "particles": particles,
        "augment_process_noise": bool(augment_process_noise),
        "propagator_overrides": overrides,
        "density_scale": float(density_scale),
        "pod_uq": bool(pod_uq),
        "pod_skip_slr": bool(pod_skip_slr),
        "pod_disable_gnss": bool(pod_disable_gnss),
        "pod_estimator": pod_estimator,
        "use_env_sources": bool(use_env_sources),
        "catalog_area": area,
        "catalog_purpose": block.get("purpose"),
        "catalog_uq_method": method,
        "catalog_env_toggles": sorted(str(t) for t in env_toggles),
        "catalog_gsi_toggles": sorted(str(t) for t in gsi_toggles),
        "catalog_sensor_toggles": sorted(str(t) for t in sensor_toggles),
        "catalog_ground_toggles": sorted(str(t) for t in ground_toggles),
        "catalog_est_toggles": sorted(str(t) for t in est_toggles),
        "catalog_control_toggles": [control_toggle] if control_toggle else [],
    }

    # Optional extended stochastic channels from catalog env toggles.
    if ("ENV_C3_COMPOSITION_OU" in env_toggles) or ("ENV_COMP_O_N2" in env_toggles):
        out["composition_discrepancy_on"] = True
        out["composition_sigma_rel"] = 0.15
        out["composition_tau_s"] = 1800.0
        out["composition_df"] = 4.0
        out["composition_clip_rel"] = 0.5

    if "ENV_S3_STORM_JUMP" in env_toggles:
        out["storm_jump_on"] = True
        out["storm_jump_rate_per_day"] = 2.0
        out["storm_jump_duration_s"] = 10800.0
        out["storm_jump_sigma_rel"] = 0.35
        out["storm_jump_mean_rel"] = 0.15
        out["storm_jump_df"] = 4.0
        out["storm_jump_clip_rel"] = 0.9

    weather_on = ("GRD_WEATHER_ON" in ground_toggles) and ("GRD_WEATHER_OFF" not in ground_toggles)
    if weather_on:
        out["pod_slr_weather_on"] = True
        out["pod_slr_weather_clear_prob"] = 0.65 if ("GRD_REG_EU" in ground_toggles) else 0.8
        out["pod_slr_weather_p_stay_clear"] = 0.97
        out["pod_slr_weather_p_stay_blocked"] = 0.92
    elif "GRD_WEATHER_OFF" in ground_toggles:
        out["pod_slr_weather_on"] = False

    if scenario_id.upper().startswith("OD_C9_"):
        out["pod_gnss_cycle_slip_gap_s"] = 20.0
        out["pod_gnss_cycle_slip_prob_per_min"] = 0.02
        out["pod_gnss_code_outlier_prob"] = 0.02
        out["pod_gnss_code_outlier_sigma_scale"] = 50.0
        out["pod_gnss_carrier_outlier_prob"] = 0.01
        out["pod_gnss_carrier_outlier_sigma_scale"] = 40.0

    outputs = block.get("outputs")
    if isinstance(outputs, dict):
        metrics = outputs.get("metrics")
        if isinstance(metrics, (list, tuple)):
            out["catalog_output_metrics"] = [str(m) for m in metrics]

    maneuver = block.get("maneuver")
    if isinstance(maneuver, dict):
        axis = str(maneuver.get("axis", "")).strip().lower()
        if axis in {"roll", "pitch", "yaw"}:
            out["maneuver_axis"] = axis
        wing_initial = maneuver.get("wing_initial_deg")
        if wing_initial is not None:
            try:
                wing_val = float(wing_initial)
                if wing_val >= 0.0:
                    out["wing_constant_deg"] = wing_val
            except Exception:
                pass
        for src_key, dst_key in (
            ("half_turn_deg", "maneuver_half_turn_deg"),
            ("full_turn_deg", "maneuver_full_turn_deg"),
        ):
            v = maneuver.get(src_key)
            if v is not None:
                try:
                    out[dst_key] = float(v)
                except Exception:
                    pass
        att0 = maneuver.get("initial_attitude_euler_deg_zyx")
        if isinstance(att0, (list, tuple)) and len(att0) == 3:
            vals = []
            for a in att0:
                try:
                    vals.append(float(a))
                except Exception:
                    vals = []
                    break
            if len(vals) == 3:
                out["initial_attitude_euler_deg_zyx"] = vals

    orbit = block.get("orbit")
    if isinstance(orbit, dict):
        o = dict(orbit)
        # Keep only fields consumed by runtime initialization.
        keep = {
            "type",
            "a_km",
            "alt_km",
            "e",
            "i_deg",
            "raan_deg",
            "argp_deg",
            "M_deg",
            "nu_deg",
            "r_eci_m",
            "v_eci_m_s",
        }
        out["orbit"] = {k: o[k] for k in o.keys() if k in keep}

    spacecraft = block.get("spacecraft")
    if isinstance(spacecraft, dict):
        mass_kg = spacecraft.get("mass_kg")
        if mass_kg is not None:
            try:
                out["spacecraft_mass_kg"] = float(mass_kg)
            except Exception:
                pass
        inertia = spacecraft.get("inertia_kgm2")
        if isinstance(inertia, (list, tuple)):
            vals = []
            for v in inertia:
                try:
                    vals.append(float(v))
                except Exception:
                    vals = []
                    break
            if len(vals) in (3, 9):
                out["spacecraft_inertia_kgm2"] = vals

    if start_utc is not None:
        out["start_utc"] = start_utc

    if scenario_type == "formation":
        out["formation_offsets_m"] = default_formation_offsets(scenario_id)

    if scenario_type == "attitude":
        toggle = control_toggle
        if toggle == "CTL_ATT_AERO_RATE":
            out["wing_amp_deg"] = 15.0
            out["wing_period_s"] = 90.0
        else:
            out["wing_amp_deg"] = 10.0
            out["wing_period_s"] = 120.0
        out["eta_max_deg"] = 110.0
        out["eta_rate_max_deg_s"] = 5.0
        out["eta_accel_max_deg_s2"] = 1.0
        # Detumble case should start with a clear high-rate initial condition.
        if scenario_id.upper() in {"ATT_A1_DETUMBLE_AERO_ONLY", "ATT_01_DETUMBLE_AERO_ONLY"}:
            out["initial_w_BI_B_deg_s"] = [8.0, -6.0, 10.0]
        # Baseline axis-rotation maneuvers should start from zero body rates by default.
        if "maneuver_axis" in out and ("initial_w_BI_B_deg_s" not in out) and ("initial_w_BI_B_rad_s" not in out):
            out["initial_w_BI_B_deg_s"] = [0.0, 0.0, 0.0]

    if "ut_alpha" in uq:
        out["ut_alpha"] = float(uq["ut_alpha"])
    if "ut_beta" in uq:
        out["ut_beta"] = float(uq["ut_beta"])
    if "ut_kappa" in uq:
        out["ut_kappa"] = float(uq["ut_kappa"])

    events = maybe_rephase_event(control_toggles, duration_s)
    if events:
        out["events"] = events

    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build run_case_studies config from scenario catalog markdown.")
    parser.add_argument("--catalog", default="scenario_catalog.md")
    parser.add_argument("--out", default="configs/case_studies_catalog.json")
    parser.add_argument("--ids", default=None, help="Comma-separated scenario IDs to include.")
    parser.add_argument("--default_particles", type=int, default=256)
    parser.add_argument("--duration_scale", type=float, default=1.0)
    parser.add_argument("--start_utc", default=None, help="Optional UTC start passed to all scenarios.")
    parser.add_argument("--no_env_sources", action="store_true")
    parser.add_argument("--dt_floor", type=float, default=1.0)
    args = parser.parse_args()

    try:
        import yaml  # type: ignore
    except Exception as exc:
        raise SystemExit("PyYAML is required for this converter. Install with `pip install pyyaml`.") from exc

    include = None
    if args.ids:
        include = {x.strip() for x in args.ids.split(",") if x.strip()}

    text = Path(args.catalog).read_text(encoding="utf-8")
    blocks = extract_yaml_blocks(text)
    scenarios = []
    for raw in blocks:
        obj = yaml.safe_load(raw)
        if not isinstance(obj, dict):
            continue
        sid = obj.get("scenario_id")
        if sid is None or str(sid).strip() in {"<string>", ""}:
            continue
        if include is not None and str(sid).strip() not in include:
            continue
        scenarios.append(
            convert_block(
                obj,
                default_particles=args.default_particles,
                duration_scale=args.duration_scale,
                start_utc=args.start_utc,
                use_env_sources=not args.no_env_sources,
                dt_floor=args.dt_floor,
            )
        )

    out = {"scenarios": scenarios}
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[catalog] wrote {len(scenarios)} scenarios -> {out_path}")


if __name__ == "__main__":
    main()
