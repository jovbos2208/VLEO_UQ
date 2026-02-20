from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnsemblePropagatorMC,
    EnvInputs,
    PropagatorConfig,
    VehicleParams,
    default_geometry,
)

from scripts.run_case_studies import (
    build_env,
    build_env_from_sources,
    build_initial_state,
    expand_covariance,
    scenario_initial_w0_rad_s,
    wing_profile_for_scenario,
)
from scripts.model_discrepancy import apply_density_model_discrepancy
from scripts.uq_parameter_channels import apply_propagator_overrides, apply_uq_parameter_channels


def _make_config(seed: int, freeze_attitude: bool) -> PropagatorConfig:
    cfg = PropagatorConfig()
    cfg.rtol = float(os.environ.get("VLEO_RTOL", "1e-6"))
    cfg.atol = float(os.environ.get("VLEO_ATOL", "1e-6"))
    cfg.max_step_s = float(os.environ.get("VLEO_MAX_STEP_S", "120"))
    cfg.rng_seed = seed
    cfg.freeze_attitude = freeze_attitude
    return cfg


def _load_case_config() -> dict:
    cfg_path = Path(os.environ.get("VLEO_CASE_CONFIG", "configs/case_studies_1000.json"))
    return json.loads(cfg_path.read_text(encoding="utf-8"))


def _select_scenario(name: str) -> dict:
    raw = _load_case_config()
    for scenario in raw.get("scenarios", []):
        if scenario.get("name") == name:
            return scenario
    raise ValueError(f"scenario '{name}' not found in config")


def _apply_vehicle_from_scenario(vehicle, scenario: dict) -> None:
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


def _make_props(seed: int, freeze_attitude: bool, scenario: dict):
    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()
    _apply_vehicle_from_scenario(vehicle, scenario)
    cfg = _make_config(seed, freeze_attitude)
    apply_propagator_overrides(cfg, scenario.get("propagator_overrides"))
    prop_det = DeterministicPropagator(aero, vehicle, cfg)
    prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
    return prop_det, prop_mc, cfg


def _p0_mission() -> np.ndarray:
    return np.diag(
        [10.0 ** 2] * 3
        + [0.01 ** 2] * 3
        + [1e-6] * 4
        + [0.0] * 3
        + [0.05 ** 2] * 2
        + [1.0 ** 2] * 3
    )


def _p0_attitude() -> np.ndarray:
    return np.diag(
        [10.0 ** 2] * 3
        + [0.01 ** 2] * 3
        + [1e-6] * 4
        + [1e-4] * 3
        + [0.05 ** 2] * 2
        + [1.0 ** 2] * 3
    )


def _p0_attitude_for_scenario(scenario: dict) -> np.ndarray:
    p0 = _p0_attitude()
    axis = str(scenario.get("maneuver_axis", "")).strip().lower()
    if axis in {"roll", "pitch", "yaw"}:
        sigma_rate_deg_s = float(scenario.get("initial_rate_sigma_deg_s", 0.5))
        sigma_rate_deg_s = max(0.05, abs(sigma_rate_deg_s))
        sigma_rate_rad_s = np.deg2rad(sigma_rate_deg_s)
        p0[10:13, 10:13] = np.eye(3) * (sigma_rate_rad_s ** 2)
    return p0


def _build_env_for_case(
    t_grid: np.ndarray,
    x0: np.ndarray,
    prop_det,
    eta1_rad: np.ndarray | None,
    eta2_rad: np.ndarray | None,
    scenario: dict,
    seed: int,
    density_scale: float = 1.0,
) -> list:
    use_env_sources = os.environ.get("VLEO_USE_ENV_SOURCES", "0") == "1"
    if not use_env_sources:
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
        apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(seed))
        return env

    start_utc = os.environ.get("VLEO_START_UTC")
    if not start_utc:
        raise RuntimeError("VLEO_START_UTC must be set when VLEO_USE_ENV_SOURCES=1")
    omni_path = os.environ.get("VLEO_OMNI_PATH", "data/space_weather/omni/omni2_all_years.dat")
    hwm14_lib = os.environ.get("VLEO_HWM14_LIB", "data/hwm14/libhwm14.so")
    hwm14_data = os.environ.get("VLEO_HWM14_DATA", "data/hwm14")
    interpolation = os.environ.get("VLEO_ENV_INTERP", "nearest")

    series = build_env_from_sources(
        t_grid_s=t_grid,
        t0_utc=dt.datetime.fromisoformat(start_utc),
        x0=x0,
        prop_det=prop_det,
        EnvInputs=EnvInputs,
        eta1_rad=eta1_rad,
        eta2_rad=eta2_rad,
        space_weather_source="omni",
        omni_path=omni_path,
        hwm14_lib=hwm14_lib,
        hwm14_data=hwm14_data,
        interpolation=interpolation,
    )
    env = series.env_inputs
    if density_scale != 1.0:
        for e in env:
            e.density *= float(density_scale)
    apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(seed))
    return env


def mission_ut_scenario() -> dict:
    name = os.environ.get("VLEO_CASE_NAME", "mission")
    seed = int(os.environ.get("VLEO_SEED", "42"))
    scenario = _select_scenario(name)
    scenario, _ = apply_uq_parameter_channels(scenario, seed=seed)
    prop_det, prop_mc, cfg = _make_props(seed, freeze_attitude=True, scenario=scenario)
    duration_s = float(scenario.get("duration_s", 3600.0))
    dt_s = float(scenario.get("dt_s", 60.0))
    t_grid = np.arange(0.0, duration_s + dt_s, dt_s)
    x0 = build_initial_state(3.986004418e14, prop_det.state_size, orbit=scenario.get("orbit"))
    env = _build_env_for_case(
        t_grid,
        x0,
        prop_det,
        None,
        None,
        scenario,
        seed,
        density_scale=float(scenario.get("density_scale", 1.0)),
    )
    P0 = expand_covariance(_p0_mission(), prop_det.state_size)
    return {
        "propagator": prop_mc,
        "t_grid": t_grid,
        "env": env,
        "x0": x0,
        "P0": P0,
        "freeze_attitude": True,
    }


def formation_ut_scenario() -> dict:
    name = os.environ.get("VLEO_CASE_NAME", "formation")
    idx = int(os.environ.get("VLEO_FORMATION_INDEX", "0"))
    scenario = _select_scenario(name)
    offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
    if idx < 0 or idx >= len(offsets):
        raise ValueError("VLEO_FORMATION_INDEX out of range")
    seed = int(os.environ.get("VLEO_SEED", "42")) + idx
    scenario, _ = apply_uq_parameter_channels(scenario, seed=seed)
    prop_det, prop_mc, cfg = _make_props(seed, freeze_attitude=True, scenario=scenario)
    duration_s = float(scenario.get("duration_s", 3600.0))
    dt_s = float(scenario.get("dt_s", 60.0))
    t_grid = np.arange(0.0, duration_s + dt_s, dt_s)
    x0 = build_initial_state(
        3.986004418e14,
        prop_det.state_size,
        r_offset=np.array(offsets[idx]),
        orbit=scenario.get("orbit"),
    )
    env = _build_env_for_case(
        t_grid,
        x0,
        prop_det,
        None,
        None,
        scenario,
        seed,
        density_scale=float(scenario.get("density_scale", 1.0)),
    )
    P0 = expand_covariance(_p0_mission(), prop_det.state_size)
    return {
        "propagator": prop_mc,
        "t_grid": t_grid,
        "env": env,
        "x0": x0,
        "P0": P0,
        "freeze_attitude": True,
    }


def attitude_ut_scenario() -> dict:
    name = os.environ.get("VLEO_CASE_NAME", "attitude")
    seed = int(os.environ.get("VLEO_SEED", "42"))
    scenario = _select_scenario(name)
    scenario, _ = apply_uq_parameter_channels(scenario, seed=seed)
    prop_det, prop_mc, cfg = _make_props(seed, freeze_attitude=False, scenario=scenario)
    duration_s = float(scenario.get("duration_s", 600.0))
    dt_s = float(scenario.get("dt_s", 5.0))
    t_grid = np.arange(0.0, duration_s + dt_s, dt_s)
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
    env = _build_env_for_case(
        t_grid,
        x0,
        prop_det,
        eta1_rad,
        eta2_rad,
        scenario,
        seed,
        density_scale=float(scenario.get("density_scale", 1.0)),
    )
    P0 = expand_covariance(_p0_attitude_for_scenario(scenario), prop_det.state_size)
    return {
        "propagator": prop_mc,
        "t_grid": t_grid,
        "env": env,
        "x0": x0,
        "P0": P0,
        "freeze_attitude": False,
    }
