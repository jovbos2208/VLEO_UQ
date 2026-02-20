from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path

import numpy as np

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_case_studies import (
    build_env,
    build_env_from_sources,
    build_initial_state,
    clone_propagator_config,
    scenario_initial_w0_rad_s,
    wing_profile_for_scenario,
)
from scripts.model_discrepancy import apply_density_model_discrepancy
from scripts.uq_parameter_channels import apply_uq_parameter_channels


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 1: run DET+STM and store outputs.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--scenario", default=None, help="Optional single scenario name to run.")
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    parser.add_argument("--disable_stm", action="store_true", help="Skip STM propagation and write NaN placeholders.")
    args = parser.parse_args()

    from vleo_uq import (
        AeroAdapter,
        DeterministicPropagator,
        EnvInputs,
        PropagatorConfig,
        StmPropagator,
        VehicleParams,
        default_geometry,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenarios = raw.get("scenarios", [])
    if not scenarios:
        raise ValueError("config must include scenarios")
    if args.scenario:
        scenarios = [s for s in scenarios if s.get("name") == args.scenario]
        if not scenarios:
            raise ValueError(f"scenario '{args.scenario}' not found in config")

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)

    def make_config(seed: int, freeze_attitude: bool) -> PropagatorConfig:
        cfg = PropagatorConfig()
        cfg.rtol = 1e-6
        cfg.atol = 1e-6
        cfg.max_step_s = 120.0
        cfg.rng_seed = seed
        cfg.freeze_attitude = freeze_attitude
        return cfg

    config_mission = make_config(42, True)
    config_att = make_config(42, False)

    def apply_vehicle_from_scenario(vehicle, scenario: dict) -> None:
        mass_kg = scenario.get("spacecraft_mass_kg")
        if mass_kg is not None:
            try:
                mass = float(mass_kg)
                if np.isfinite(mass) and mass > 0.0:
                    vehicle.mass_kg = mass
            except Exception:
                pass
        inertia = scenario.get("spacecraft_inertia_kgm2")
        if isinstance(inertia, (list, tuple)):
            try:
                arr = np.asarray(inertia, dtype=float)
            except Exception:
                arr = np.array([], dtype=float)
            if arr.size == 3:
                vehicle.inertia_B = np.diag(arr)
            elif arr.size == 9:
                vehicle.inertia_B = arr.reshape(3, 3)

    def make_props(scenario: dict, freeze_attitude: bool):
        cfg_base = config_mission if freeze_attitude else config_att
        cfg = clone_propagator_config(
            PropagatorConfig=PropagatorConfig,
            base_config=cfg_base,
            overrides=scenario.get("propagator_overrides", {}),
            freeze_attitude=freeze_attitude,
        )
        cfg.rng_seed = int(os.environ.get("VLEO_SEED", "42"))
        vehicle = VehicleParams()
        apply_vehicle_from_scenario(vehicle, scenario)
        prop_det = DeterministicPropagator(aero, vehicle, cfg)
        prop_stm = StmPropagator(aero, vehicle, cfg)
        return prop_det, prop_stm

    disable_stm = bool(args.disable_stm) or os.environ.get("VLEO_DISABLE_STM", "0").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }
    if disable_stm:
        print("[phase1] STM disabled: writing NaN stm_mean/stm_cov placeholders")

    start_utc = None
    if args.start_utc:
        start_utc = dt.datetime.fromisoformat(args.start_utc)

    def build_env_case(
        t_grid,
        x0,
        prop_det,
        eta1_rad,
        eta2_rad,
        scenario: dict,
        scenario_seed: int,
        density_scale: float = 1.0,
    ):
        if not args.use_env_sources:
            env = build_env(
                t_grid,
                EnvInputs,
                density=1e-12,
                temperature_K=1000.0,
                particle_mass_kg=28.0 * 1.6605390689252e-27,
                eta1_rad=eta1_rad,
                eta2_rad=eta2_rad,
            )
            if density_scale != 1.0:
                for e in env:
                    e.density *= float(density_scale)
            apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(scenario_seed))
            return env
        if start_utc is None:
            raise ValueError("start_utc required when use_env_sources is set")
        strict = os.environ.get("VLEO_ENV_SOURCE_STRICT", "0").strip().lower() in {"1", "true", "yes", "on"}
        try:
            series = build_env_from_sources(
                t_grid_s=t_grid,
                t0_utc=start_utc,
                x0=x0,
                prop_det=prop_det,
                EnvInputs=EnvInputs,
                eta1_rad=eta1_rad,
                eta2_rad=eta2_rad,
                space_weather_source="omni",
                omni_path=args.omni_path,
                hwm14_lib=args.hwm14_lib,
                hwm14_data=args.hwm14_data,
                interpolation=args.env_interpolation,
            )
            env = series.env_inputs
            if density_scale != 1.0:
                for e in env:
                    e.density *= float(density_scale)
            apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(scenario_seed))
            return env
        except Exception as exc:
            if strict:
                raise
            print(
                f"[phase1] warning: env source build failed ({type(exc).__name__}: {exc}); "
                "falling back to constant environment.",
                flush=True,
            )
            env = build_env(
                t_grid,
                EnvInputs,
                density=1e-12,
                temperature_K=1000.0,
                particle_mass_kg=28.0 * 1.6605390689252e-27,
                eta1_rad=eta1_rad,
                eta2_rad=eta2_rad,
            )
            if density_scale != 1.0:
                for e in env:
                    e.density *= float(density_scale)
            apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(scenario_seed))
            return env

    def require_finite(name: str, label: str, arr: np.ndarray, dt_s: float) -> None:
        if np.isfinite(arr).all():
            return
        bad = int(np.size(arr) - np.count_nonzero(np.isfinite(arr)))
        if arr.ndim >= 2:
            flat = arr.reshape(arr.shape[0], -1)
            bad_rows = int(np.sum(~np.isfinite(flat).all(axis=1)))
            raise RuntimeError(
                f"[phase1] {name}: {label} contains non-finite values "
                f"(bad_rows={bad_rows}/{arr.shape[0]}, bad_values={bad}). "
                f"Reduce scenario dt_s (current {dt_s:g}s; ATT cases typically need <=5s)."
            )
        raise RuntimeError(
            f"[phase1] {name}: {label} contains non-finite values "
            f"(bad_values={bad}). Reduce scenario dt_s (current {dt_s:g}s)."
        )

    for scenario in scenarios:
        scenario_seed = int(os.environ.get("VLEO_SEED", "42"))
        scenario, _ = apply_uq_parameter_channels(scenario, seed=scenario_seed)
        name = scenario["name"]
        typ = scenario["type"]
        duration_s = float(scenario.get("duration_s", 3600.0))
        dt_s = float(scenario.get("dt_s", 60.0))
        t_grid = np.arange(0.0, duration_s + dt_s, dt_s)

        if typ == "mission":
            prop_det, prop_stm = make_props(scenario, True)
            x0 = build_initial_state(3.986004418e14, prop_det.state_size, orbit=scenario.get("orbit"))
            env = build_env_case(
                t_grid,
                x0,
                prop_det,
                None,
                None,
                scenario,
                scenario_seed,
                density_scale=float(scenario.get("density_scale", 1.0)),
            )
            det_states = prop_det.propagate(x0, t_grid, env)
            require_finite(name, "det_states", det_states, dt_s)
            P0 = np.diag([10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [0.0] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3)
            if disable_stm:
                n_state = det_states.shape[1]
                stm_mean = np.full_like(det_states, np.nan)
                stm_cov = np.full((det_states.shape[0], n_state, n_state), np.nan, dtype=float)
            else:
                stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)
                require_finite(name, "stm_mean", stm_mean, dt_s)
                require_finite(name, "stm_cov", stm_cov, dt_s)
            case_dir = outdir / name
            case_dir.mkdir(parents=True, exist_ok=True)
            np.savez(
                case_dir / "det_stm.npz",
                t_grid=t_grid,
                det_states=det_states,
                stm_mean=stm_mean,
                stm_cov=stm_cov,
            )
            print(f"[phase1] saved {case_dir}/det_stm.npz")
        elif typ == "formation":
            prop_det, prop_stm = make_props(scenario, True)
            offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
            for idx, offset in enumerate(offsets):
                obj_name = f"{name}/obj_{idx+1:02d}"
                x0 = build_initial_state(
                    3.986004418e14,
                    prop_det.state_size,
                    r_offset=np.array(offset, dtype=float),
                    orbit=scenario.get("orbit"),
                )
                env = build_env_case(
                    t_grid,
                    x0,
                    prop_det,
                    None,
                    None,
                    scenario,
                    scenario_seed,
                    density_scale=float(scenario.get("density_scale", 1.0)),
                )
                det_states = prop_det.propagate(x0, t_grid, env)
                require_finite(obj_name, "det_states", det_states, dt_s)
                P0 = np.diag([10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [0.0] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3)
                if disable_stm:
                    n_state = det_states.shape[1]
                    stm_mean = np.full_like(det_states, np.nan)
                    stm_cov = np.full((det_states.shape[0], n_state, n_state), np.nan, dtype=float)
                else:
                    stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)
                    require_finite(obj_name, "stm_mean", stm_mean, dt_s)
                    require_finite(obj_name, "stm_cov", stm_cov, dt_s)
                case_dir = outdir / obj_name
                case_dir.mkdir(parents=True, exist_ok=True)
                np.savez(
                    case_dir / "det_stm.npz",
                    t_grid=t_grid,
                    det_states=det_states,
                    stm_mean=stm_mean,
                    stm_cov=stm_cov,
                )
                print(f"[phase1] saved {case_dir}/det_stm.npz")
        elif typ == "attitude":
            prop_det, prop_stm = make_props(scenario, False)
            eta1_rad, eta2_rad = wing_profile_for_scenario(
                t_grid,
                scenario,
                default_amp_deg=float(scenario.get("wing_amp_deg", 10.0)),
                default_period_s=float(scenario.get("wing_period_s", 120.0)),
            )
            x0 = build_initial_state(
                3.986004418e14,
                prop_det.state_size,
                orbit=scenario.get("orbit"),
                w0=scenario_initial_w0_rad_s(scenario),
            )
            env = build_env_case(
                t_grid,
                x0,
                prop_det,
                eta1_rad,
                eta2_rad,
                scenario,
                scenario_seed,
                density_scale=float(scenario.get("density_scale", 1.0)),
            )
            det_states = prop_det.propagate(x0, t_grid, env)
            require_finite(name, "det_states", det_states, dt_s)
            P0 = np.diag([10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [1e-4] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3)
            if disable_stm:
                n_state = det_states.shape[1]
                stm_mean = np.full_like(det_states, np.nan)
                stm_cov = np.full((det_states.shape[0], n_state, n_state), np.nan, dtype=float)
            else:
                stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)
                require_finite(name, "stm_mean", stm_mean, dt_s)
                require_finite(name, "stm_cov", stm_cov, dt_s)
            case_dir = outdir / name
            case_dir.mkdir(parents=True, exist_ok=True)
            np.savez(
                case_dir / "det_stm.npz",
                t_grid=t_grid,
                det_states=det_states,
                stm_mean=stm_mean,
                stm_cov=stm_cov,
            )
            print(f"[phase1] saved {case_dir}/det_stm.npz")
        else:
            raise ValueError(f"unknown scenario type: {typ}")


if __name__ == "__main__":
    main()
