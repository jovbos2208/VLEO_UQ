from __future__ import annotations

import copy
import hashlib
from typing import Any

import numpy as np

_SLOW_PARAM_DEFAULTS = {
    "rho_fast_tau_s": 900.0,
    "rho_fast_sigma": 0.15,
    "rho_bias_tau_s": 7200.0,
    "rho_bias_sigma": 0.08,
    "wind_tau_s": 900.0,
    "wind_sigma": 10.0,
    "j2_earth": 1.08262668e-3,
    "j3_earth": -2.53215306e-6,
    "j4_earth": -1.61098761e-6,
    "sun_ephemeris_scale": 1.0,
    "moon_ephemeris_scale": 1.0,
    "srp_cr": 1.2,
    "srp_area_m2": 0.2,
    "albedo_ir_cr": 1.0,
    "albedo_ir_area_m2": 0.2,
    "albedo_pressure_n_m2": 1.2e-6,
    "earth_ir_pressure_n_m2": 1.0e-6,
    "tide_loading_scale": 1.0,
    "magnetic_field_scale": 1.0,
    "residual_dipole_scale": 1.0,
}

_SLOW_PARAM_SIGMA_KEYS = {
    "rho_fast_tau_s": "uq_rho_fast_tau_rel_sigma",
    "rho_fast_sigma": "uq_rho_fast_sigma_rel_sigma",
    "rho_bias_tau_s": "uq_rho_bias_tau_rel_sigma",
    "rho_bias_sigma": "uq_rho_bias_sigma_rel_sigma",
    "wind_tau_s": "uq_wind_tau_rel_sigma",
    "wind_sigma": "uq_wind_sigma_rel_sigma",
    "j2_earth": "uq_j2_rel_sigma",
    "j3_earth": "uq_j3_rel_sigma",
    "j4_earth": "uq_j4_rel_sigma",
    "sun_ephemeris_scale": "uq_sun_ephemeris_scale_rel_sigma",
    "moon_ephemeris_scale": "uq_moon_ephemeris_scale_rel_sigma",
    "srp_cr": "uq_srp_cr_rel_sigma",
    "srp_area_m2": "uq_srp_area_rel_sigma",
    "albedo_ir_cr": "uq_albedo_ir_cr_rel_sigma",
    "albedo_ir_area_m2": "uq_albedo_ir_area_rel_sigma",
    "albedo_pressure_n_m2": "uq_albedo_pressure_rel_sigma",
    "earth_ir_pressure_n_m2": "uq_earth_ir_pressure_rel_sigma",
    "tide_loading_scale": "uq_tide_loading_scale_rel_sigma",
    "magnetic_field_scale": "uq_magnetic_field_scale_rel_sigma",
    "residual_dipole_scale": "uq_residual_dipole_scale_rel_sigma",
}

_SIGNED_OVERRIDE_KEYS = {
    "j3_earth",
    "j4_earth",
}


def _as_finite_float(v: Any) -> float | None:
    try:
        out = float(v)
    except Exception:
        return None
    if not np.isfinite(out):
        return None
    return out


def _stable_mix_from_name(name: str) -> int:
    digest = hashlib.sha256(name.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def _draw_lognormal_scale(rng: np.random.Generator, rel_sigma: float) -> float:
    rel_sigma = float(abs(rel_sigma))
    if rel_sigma <= 0.0:
        return 1.0
    return float(np.exp(rng.normal(0.0, rel_sigma)))


def apply_uq_parameter_channels(
    scenario: dict,
    *,
    seed: int,
) -> tuple[dict, dict]:
    if not isinstance(scenario, dict):
        return {}, {}
    out = copy.deepcopy(scenario)
    name = str(out.get("name", out.get("catalog_scenario_id", "scenario")))
    seed_offset = int(_as_finite_float(out.get("uq_parameter_seed_offset")) or 0)
    mix = _stable_mix_from_name(name)
    rng_seed = int((np.uint64(seed) ^ np.uint64(mix) ^ np.uint64(seed_offset)).item())
    rng = np.random.default_rng(rng_seed)

    draw: dict[str, Any] = {"seed": int(rng_seed)}

    mass_base = _as_finite_float(out.get("spacecraft_mass_kg"))
    if mass_base is not None and mass_base > 0.0:
        mass_sigma_abs = _as_finite_float(out.get("uq_mass_sigma_kg")) or 0.0
        mass_rel_sigma = _as_finite_float(out.get("uq_mass_rel_sigma")) or 0.0
        fuel_sigma_abs = _as_finite_float(out.get("uq_fuel_gauging_sigma_kg")) or 0.0
        fuel_rel_sigma = _as_finite_float(out.get("uq_fuel_gauging_rel_sigma")) or 0.0
        mass_sigma = (
            abs(float(mass_sigma_abs))
            + abs(float(mass_rel_sigma)) * mass_base
            + abs(float(fuel_sigma_abs))
            + abs(float(fuel_rel_sigma)) * mass_base
        )
        if mass_sigma > 0.0:
            mass_draw = max(1e-6, float(rng.normal(mass_base, mass_sigma)))
            out["spacecraft_mass_kg"] = mass_draw
            draw["spacecraft_mass_kg"] = {
                "base": float(mass_base),
                "value": float(mass_draw),
                "sigma_abs": float(mass_sigma_abs),
                "sigma_rel": float(mass_rel_sigma),
                "fuel_sigma_abs": float(fuel_sigma_abs),
                "fuel_sigma_rel": float(fuel_rel_sigma),
            }

    inertia_rel_sigma = abs(float(_as_finite_float(out.get("uq_inertia_rel_sigma")) or 0.0))
    inertia = out.get("spacecraft_inertia_kgm2")
    if inertia_rel_sigma > 0.0 and isinstance(inertia, (list, tuple)):
        arr = np.asarray(inertia, dtype=float).reshape(-1)
        if arr.size in {3, 9} and np.isfinite(arr).all():
            scale = _draw_lognormal_scale(rng, inertia_rel_sigma)
            out["spacecraft_inertia_kgm2"] = (arr * scale).tolist()
            draw["spacecraft_inertia_scale"] = {
                "base_scale": 1.0,
                "scale": float(scale),
                "sigma_rel": float(inertia_rel_sigma),
            }

    density_base = _as_finite_float(out.get("density_scale"))
    if density_base is None or density_base <= 0.0:
        density_base = 1.0
    density_rel_sigma = abs(float(_as_finite_float(out.get("uq_density_scale_rel_sigma")) or 0.0))
    if density_rel_sigma > 0.0:
        density_scale = _draw_lognormal_scale(rng, density_rel_sigma)
        density_draw = max(1e-8, float(density_base * density_scale))
        out["density_scale"] = density_draw
        draw["density_scale"] = {
            "base": float(density_base),
            "value": float(density_draw),
            "sigma_rel": float(density_rel_sigma),
        }

    overrides = dict(out.get("propagator_overrides", {}) or {})
    for key, sigma_key in _SLOW_PARAM_SIGMA_KEYS.items():
        rel_sigma = abs(float(_as_finite_float(out.get(sigma_key)) or 0.0))
        if rel_sigma <= 0.0:
            continue
        base = _as_finite_float(overrides.get(key))
        if base is None:
            base = _SLOW_PARAM_DEFAULTS[key]
        scale = _draw_lognormal_scale(rng, rel_sigma)
        val = float(base * scale)
        if key not in _SIGNED_OVERRIDE_KEYS:
            val = max(1e-12, val)
        overrides[key] = val
        draw[key] = {
            "base": float(base),
            "value": float(val),
            "sigma_rel": float(rel_sigma),
        }
    if overrides:
        out["propagator_overrides"] = overrides

    if len(draw) == 1:
        return out, {}
    out["uq_parameter_draw"] = draw
    return out, draw


def apply_propagator_overrides(cfg, overrides: dict | None) -> None:
    if not isinstance(overrides, dict):
        return
    for key, value in overrides.items():
        if not hasattr(cfg, key):
            continue
        if isinstance(value, (list, tuple, np.ndarray)):
            try:
                setattr(cfg, key, np.asarray(value, dtype=float))
                continue
            except Exception:
                continue
        v = _as_finite_float(value)
        if v is None:
            continue
        setattr(cfg, key, float(v))
