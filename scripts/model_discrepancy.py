from __future__ import annotations

import hashlib
import os
from typing import Any, Mapping

import numpy as np


def _as_finite_float(value: Any, default: float) -> float:
    try:
        out = float(value)
    except Exception:
        return float(default)
    if not np.isfinite(out):
        return float(default)
    return float(out)


def _as_nonnegative_float(value: Any, default: float) -> float:
    return max(0.0, _as_finite_float(value, default))


def _as_int(value: Any, default: int) -> int:
    try:
        return int(value)
    except Exception:
        return int(default)


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(int(value) != 0)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def _stable_u64(text: str) -> int:
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="little", signed=False)


def _scenario_name(scenario: Mapping[str, Any] | None) -> str:
    if not isinstance(scenario, Mapping):
        return "scenario"
    return str(scenario.get("name", scenario.get("catalog_scenario_id", "scenario")))


def _scenario_get(scenario: Mapping[str, Any] | None, key: str, default: Any) -> Any:
    if not isinstance(scenario, Mapping):
        return default
    return scenario.get(key, default)


def _env_bool(env_key: str, fallback: bool) -> bool:
    raw = os.environ.get(env_key)
    if raw is None:
        return bool(fallback)
    return _as_bool(raw, fallback)


def _env_float(env_key: str, fallback: float) -> float:
    raw = os.environ.get(env_key)
    if raw is None:
        return float(fallback)
    return _as_finite_float(raw, fallback)


def _env_int(env_key: str, fallback: int) -> int:
    raw = os.environ.get(env_key)
    if raw is None:
        return int(fallback)
    return _as_int(raw, fallback)


def _student_t_unit_variance(rng: np.random.Generator, dof: float) -> float:
    dof = float(dof)
    if dof > 2.0:
        # Var(t_nu) = nu/(nu-2); apply scale so output has unit variance.
        scale = np.sqrt((dof - 2.0) / dof)
        return float(rng.standard_t(dof) * scale)
    return float(rng.normal())


def _bounded_log_from_rel(rel: float, clip_rel: float) -> float:
    rel = float(rel)
    clip_rel = max(0.0, float(clip_rel))
    if clip_rel > 0.0:
        rel = float(np.clip(rel, -clip_rel, clip_rel))
    rel = max(rel, -0.999999)
    return float(np.log1p(rel))


def density_discrepancy_scale_series(
    t_grid_s: np.ndarray,
    *,
    rng: np.random.Generator,
    sigma_rel: float,
    tau_s: float,
    dof: float,
    clip_rel: float,
) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t_grid_s, dtype=float).reshape(-1)
    n = int(t.size)
    if n == 0:
        return np.zeros(0, dtype=float), np.zeros(0, dtype=float)
    if n == 1:
        return np.ones(1, dtype=float), np.zeros(1, dtype=float)

    sigma_rel = max(0.0, float(sigma_rel))
    tau_s = max(0.0, float(tau_s))
    dof = float(dof)
    clip_rel = max(0.0, float(clip_rel))

    sigma_log = float(np.log1p(sigma_rel))
    x = np.zeros(n, dtype=float)

    for i in range(1, n):
        dt = max(0.0, float(t[i] - t[i - 1]))
        if tau_s > 0.0:
            phi = float(np.exp(-dt / tau_s))
            step_std = sigma_log * float(np.sqrt(max(0.0, 1.0 - phi * phi)))
            x[i] = phi * x[i - 1] + step_std * _student_t_unit_variance(rng, dof)
        else:
            step_std = sigma_log * float(np.sqrt(dt))
            x[i] = x[i - 1] + step_std * _student_t_unit_variance(rng, dof)

    scale = np.exp(x)
    if clip_rel > 0.0:
        lo = max(1e-6, 1.0 - clip_rel)
        hi = 1.0 + clip_rel
        scale = np.clip(scale, lo, hi)
    return scale, x


def storm_jump_log_series(
    t_grid_s: np.ndarray,
    *,
    rng: np.random.Generator,
    enabled: bool,
    rate_per_day: float,
    duration_s: float,
    jump_sigma_rel: float,
    jump_mean_rel: float,
    dof: float,
    clip_rel: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    t = np.asarray(t_grid_s, dtype=float).reshape(-1)
    n = int(t.size)
    out = np.zeros(n, dtype=float)
    meta = {
        "enabled": bool(enabled),
        "rate_per_day": float(max(0.0, float(rate_per_day))),
        "duration_s": float(max(0.0, float(duration_s))),
        "jump_sigma_rel": float(max(0.0, float(jump_sigma_rel))),
        "jump_mean_rel": float(float(jump_mean_rel)),
        "df": float(dof),
        "clip_rel": float(max(0.0, float(clip_rel))),
        "event_count": 0,
    }
    if (not enabled) or n <= 1 or meta["rate_per_day"] <= 0.0 or meta["duration_s"] <= 0.0:
        return out, meta

    jump_std_log = float(np.log1p(meta["jump_sigma_rel"]))
    jump_mean_log = _bounded_log_from_rel(meta["jump_mean_rel"], clip_rel=meta["clip_rel"])
    active_until = -np.inf
    active_jump_log = 0.0
    event_count = 0

    for i in range(1, n):
        ti = float(t[i])
        dt = max(0.0, float(t[i] - t[i - 1]))
        if ti >= active_until:
            active_jump_log = 0.0
            p_start = 1.0 - float(np.exp(-meta["rate_per_day"] * dt / 86400.0))
            if rng.random() < p_start:
                sample = jump_mean_log + jump_std_log * _student_t_unit_variance(rng, dof)
                if meta["clip_rel"] > 0.0:
                    lo = _bounded_log_from_rel(-meta["clip_rel"], clip_rel=0.0)
                    hi = _bounded_log_from_rel(+meta["clip_rel"], clip_rel=0.0)
                    sample = float(np.clip(sample, lo, hi))
                active_jump_log = float(sample)
                active_until = ti + meta["duration_s"]
                event_count += 1
        out[i] = active_jump_log

    meta["event_count"] = int(event_count)
    if n > 0:
        scale = np.exp(out)
        meta["scale_mean"] = float(np.mean(scale))
        meta["scale_std"] = float(np.std(scale))
        meta["scale_min"] = float(np.min(scale))
        meta["scale_max"] = float(np.max(scale))
    return out, meta


def _sample_switch_active_mask(
    t_grid_s: np.ndarray,
    *,
    rng: np.random.Generator,
    start_rate_per_hour: float,
    mean_duration_s: float,
) -> np.ndarray:
    t = np.asarray(t_grid_s, dtype=float).reshape(-1)
    n = int(t.size)
    active = np.zeros(n, dtype=bool)
    if n == 0 or start_rate_per_hour <= 0.0 or mean_duration_s <= 0.0:
        return active
    active_until = -np.inf
    for i in range(n):
        ti = float(t[i])
        if ti < active_until:
            active[i] = True
            continue
        if i == 0:
            dt = max(1.0, float(np.median(np.diff(t))) if n > 1 else 1.0)
        else:
            dt = max(0.0, float(t[i] - t[i - 1]))
        p_start = 1.0 - float(np.exp(-start_rate_per_hour * dt / 3600.0))
        if rng.random() < p_start:
            dur = float(rng.exponential(mean_duration_s))
            active_until = ti + max(1.0, dur)
            active[i] = True
    return active


def apply_density_model_discrepancy(
    env: list,
    t_grid_s: np.ndarray,
    *,
    scenario: Mapping[str, Any] | None,
    seed: int,
) -> dict[str, Any]:
    t = np.asarray(t_grid_s, dtype=float).reshape(-1)
    if len(env) != int(t.size):
        raise ValueError("env length must match t_grid size for model discrepancy")

    enabled_default = _as_bool(_scenario_get(scenario, "model_discrepancy_on", False), False)
    enabled = _env_bool("VLEO_MODEL_DISCREPANCY_ON", enabled_default)

    sigma_rel_default = _as_nonnegative_float(_scenario_get(scenario, "model_discrepancy_sigma_rel", 0.0), 0.0)
    sigma_rel = _as_nonnegative_float(_env_float("VLEO_MODEL_DISCREPANCY_SIGMA_REL", sigma_rel_default), 0.0)

    tau_s_default = _as_nonnegative_float(_scenario_get(scenario, "model_discrepancy_tau_s", 1800.0), 1800.0)
    tau_s = _as_nonnegative_float(_env_float("VLEO_MODEL_DISCREPANCY_TAU_S", tau_s_default), 1800.0)

    dof_default = _as_finite_float(_scenario_get(scenario, "model_discrepancy_df", 4.0), 4.0)
    dof = _as_finite_float(_env_float("VLEO_MODEL_DISCREPANCY_DF", dof_default), 4.0)

    clip_rel_default = _as_nonnegative_float(_scenario_get(scenario, "model_discrepancy_clip_rel", 0.8), 0.8)
    clip_rel = _as_nonnegative_float(_env_float("VLEO_MODEL_DISCREPANCY_CLIP_REL", clip_rel_default), 0.8)

    seed_offset_default = _as_int(_scenario_get(scenario, "model_discrepancy_seed_offset", 0), 0)
    seed_offset = _env_int("VLEO_MODEL_DISCREPANCY_SEED_OFFSET", seed_offset_default)

    comp_on_default = _as_bool(_scenario_get(scenario, "composition_discrepancy_on", False), False)
    comp_on = _env_bool("VLEO_COMPOSITION_DISCREPANCY_ON", comp_on_default)
    comp_sigma_default = _as_nonnegative_float(_scenario_get(scenario, "composition_sigma_rel", 0.0), 0.0)
    comp_sigma_rel = _as_nonnegative_float(
        _env_float("VLEO_COMPOSITION_DISCREPANCY_SIGMA_REL", comp_sigma_default), 0.0
    )
    comp_tau_default = _as_nonnegative_float(_scenario_get(scenario, "composition_tau_s", 1800.0), 1800.0)
    comp_tau_s = _as_nonnegative_float(_env_float("VLEO_COMPOSITION_DISCREPANCY_TAU_S", comp_tau_default), 1800.0)
    comp_df_default = _as_finite_float(_scenario_get(scenario, "composition_df", 4.0), 4.0)
    comp_df = _as_finite_float(_env_float("VLEO_COMPOSITION_DISCREPANCY_DF", comp_df_default), 4.0)
    comp_clip_default = _as_nonnegative_float(_scenario_get(scenario, "composition_clip_rel", 0.5), 0.5)
    comp_clip_rel = _as_nonnegative_float(
        _env_float("VLEO_COMPOSITION_DISCREPANCY_CLIP_REL", comp_clip_default), 0.5
    )
    comp_seed_offset_default = _as_int(_scenario_get(scenario, "composition_seed_offset", 101), 101)
    comp_seed_offset = _env_int("VLEO_COMPOSITION_DISCREPANCY_SEED_OFFSET", comp_seed_offset_default)

    storm_on_default = _as_bool(_scenario_get(scenario, "storm_jump_on", False), False)
    storm_on = _env_bool("VLEO_STORM_JUMP_ON", storm_on_default)
    storm_rate_default = _as_nonnegative_float(_scenario_get(scenario, "storm_jump_rate_per_day", 0.0), 0.0)
    storm_rate_per_day = _as_nonnegative_float(
        _env_float("VLEO_STORM_JUMP_RATE_PER_DAY", storm_rate_default), 0.0
    )
    storm_duration_default = _as_nonnegative_float(_scenario_get(scenario, "storm_jump_duration_s", 10800.0), 10800.0)
    storm_duration_s = _as_nonnegative_float(
        _env_float("VLEO_STORM_JUMP_DURATION_S", storm_duration_default), 10800.0
    )
    storm_sigma_default = _as_nonnegative_float(_scenario_get(scenario, "storm_jump_sigma_rel", 0.35), 0.35)
    storm_sigma_rel = _as_nonnegative_float(
        _env_float("VLEO_STORM_JUMP_SIGMA_REL", storm_sigma_default), 0.35
    )
    storm_mean_default = _as_finite_float(_scenario_get(scenario, "storm_jump_mean_rel", 0.0), 0.0)
    storm_mean_rel = _as_finite_float(_env_float("VLEO_STORM_JUMP_MEAN_REL", storm_mean_default), 0.0)
    storm_df_default = _as_finite_float(_scenario_get(scenario, "storm_jump_df", 4.0), 4.0)
    storm_df = _as_finite_float(_env_float("VLEO_STORM_JUMP_DF", storm_df_default), 4.0)
    storm_clip_default = _as_nonnegative_float(_scenario_get(scenario, "storm_jump_clip_rel", 0.9), 0.9)
    storm_clip_rel = _as_nonnegative_float(_env_float("VLEO_STORM_JUMP_CLIP_REL", storm_clip_default), 0.9)
    storm_seed_offset_default = _as_int(_scenario_get(scenario, "storm_jump_seed_offset", 202), 202)
    storm_seed_offset = _env_int("VLEO_STORM_JUMP_SEED_OFFSET", storm_seed_offset_default)

    ao_on_default = _as_bool(_scenario_get(scenario, "ao_erosion_on", False), False)
    ao_on = _env_bool("VLEO_AO_EROSION_ON", ao_on_default)
    ao_sigma_default = _as_nonnegative_float(_scenario_get(scenario, "ao_erosion_sigma_rel", 0.0), 0.0)
    ao_sigma_rel = _as_nonnegative_float(_env_float("VLEO_AO_EROSION_SIGMA_REL", ao_sigma_default), 0.0)
    ao_tau_default = _as_nonnegative_float(_scenario_get(scenario, "ao_erosion_tau_s", 7200.0), 7200.0)
    ao_tau_s = _as_nonnegative_float(_env_float("VLEO_AO_EROSION_TAU_S", ao_tau_default), 7200.0)
    ao_df_default = _as_finite_float(_scenario_get(scenario, "ao_erosion_df", 4.0), 4.0)
    ao_df = _as_finite_float(_env_float("VLEO_AO_EROSION_DF", ao_df_default), 4.0)
    ao_clip_default = _as_nonnegative_float(_scenario_get(scenario, "ao_erosion_clip_rel", 0.5), 0.5)
    ao_clip_rel = _as_nonnegative_float(_env_float("VLEO_AO_EROSION_CLIP_REL", ao_clip_default), 0.5)
    ao_rate_default = _as_finite_float(_scenario_get(scenario, "ao_erosion_rate_rel_per_day", 0.0), 0.0)
    ao_rate_rel_per_day = _as_finite_float(
        _env_float("VLEO_AO_EROSION_RATE_REL_PER_DAY", ao_rate_default),
        0.0,
    )
    ao_seed_offset_default = _as_int(_scenario_get(scenario, "ao_erosion_seed_offset", 303), 303)
    ao_seed_offset = _env_int("VLEO_AO_EROSION_SEED_OFFSET", ao_seed_offset_default)

    plasma_on_default = _as_bool(_scenario_get(scenario, "plasma_charging_on", False), False)
    plasma_on = _env_bool("VLEO_PLASMA_CHARGING_ON", plasma_on_default)
    plasma_sigma_default = _as_nonnegative_float(_scenario_get(scenario, "plasma_charging_sigma_rel", 0.0), 0.0)
    plasma_sigma_rel = _as_nonnegative_float(
        _env_float("VLEO_PLASMA_CHARGING_SIGMA_REL", plasma_sigma_default),
        0.0,
    )
    plasma_tau_default = _as_nonnegative_float(_scenario_get(scenario, "plasma_charging_tau_s", 1200.0), 1200.0)
    plasma_tau_s = _as_nonnegative_float(_env_float("VLEO_PLASMA_CHARGING_TAU_S", plasma_tau_default), 1200.0)
    plasma_df_default = _as_finite_float(_scenario_get(scenario, "plasma_charging_df", 4.0), 4.0)
    plasma_df = _as_finite_float(_env_float("VLEO_PLASMA_CHARGING_DF", plasma_df_default), 4.0)
    plasma_clip_default = _as_nonnegative_float(_scenario_get(scenario, "plasma_charging_clip_rel", 0.8), 0.8)
    plasma_clip_rel = _as_nonnegative_float(
        _env_float("VLEO_PLASMA_CHARGING_CLIP_REL", plasma_clip_default),
        0.8,
    )
    plasma_seed_offset_default = _as_int(_scenario_get(scenario, "plasma_charging_seed_offset", 404), 404)
    plasma_seed_offset = _env_int("VLEO_PLASMA_CHARGING_SEED_OFFSET", plasma_seed_offset_default)

    structural_on_default = _as_bool(_scenario_get(scenario, "structural_flex_on", False), False)
    structural_on = _env_bool("VLEO_STRUCTURAL_FLEX_ON", structural_on_default)
    structural_sigma_default = _as_nonnegative_float(_scenario_get(scenario, "structural_flex_sigma_rel", 0.0), 0.0)
    structural_sigma_rel = _as_nonnegative_float(
        _env_float("VLEO_STRUCTURAL_FLEX_SIGMA_REL", structural_sigma_default),
        0.0,
    )
    structural_tau_default = _as_nonnegative_float(_scenario_get(scenario, "structural_flex_tau_s", 1800.0), 1800.0)
    structural_tau_s = _as_nonnegative_float(
        _env_float("VLEO_STRUCTURAL_FLEX_TAU_S", structural_tau_default),
        1800.0,
    )
    structural_df_default = _as_finite_float(_scenario_get(scenario, "structural_flex_df", 4.0), 4.0)
    structural_df = _as_finite_float(_env_float("VLEO_STRUCTURAL_FLEX_DF", structural_df_default), 4.0)
    structural_clip_default = _as_nonnegative_float(_scenario_get(scenario, "structural_flex_clip_rel", 0.3), 0.3)
    structural_clip_rel = _as_nonnegative_float(
        _env_float("VLEO_STRUCTURAL_FLEX_CLIP_REL", structural_clip_default),
        0.3,
    )
    structural_seed_offset_default = _as_int(_scenario_get(scenario, "structural_flex_seed_offset", 505), 505)
    structural_seed_offset = _env_int("VLEO_STRUCTURAL_FLEX_SEED_OFFSET", structural_seed_offset_default)

    rarefied_on_default = _as_bool(_scenario_get(scenario, "rarefied_regime_switch_on", False), False)
    rarefied_on = _env_bool("VLEO_RAREFIED_REGIME_SWITCH_ON", rarefied_on_default)
    rarefied_rate_default = _as_nonnegative_float(
        _scenario_get(scenario, "rarefied_switch_rate_per_hour", 0.0),
        0.0,
    )
    rarefied_rate_per_hour = _as_nonnegative_float(
        _env_float("VLEO_RAREFIED_SWITCH_RATE_PER_HOUR", rarefied_rate_default),
        0.0,
    )
    rarefied_duration_default = _as_nonnegative_float(
        _scenario_get(scenario, "rarefied_switch_mean_duration_s", 0.0),
        0.0,
    )
    rarefied_mean_duration_s = _as_nonnegative_float(
        _env_float("VLEO_RAREFIED_SWITCH_MEAN_DURATION_S", rarefied_duration_default),
        0.0,
    )
    rarefied_alt_method_default = _as_int(
        _scenario_get(scenario, "rarefied_switch_alt_temperature_ratio_method", 2),
        2,
    )
    rarefied_alt_method = _env_int(
        "VLEO_RAREFIED_SWITCH_ALT_TEMPERATURE_RATIO_METHOD",
        rarefied_alt_method_default,
    )
    rarefied_seed_offset_default = _as_int(_scenario_get(scenario, "rarefied_switch_seed_offset", 606), 606)
    rarefied_seed_offset = _env_int("VLEO_RAREFIED_SWITCH_SEED_OFFSET", rarefied_seed_offset_default)

    meta = {
        "enabled": bool(enabled),
        "sigma_rel": float(sigma_rel),
        "tau_s": float(tau_s),
        "df": float(dof),
        "clip_rel": float(clip_rel),
        "seed_offset": int(seed_offset),
        "applied_density": False,
        "applied_composition": False,
        "applied_storm_jump": False,
        "composition": {
            "enabled": bool(comp_on),
            "sigma_rel": float(comp_sigma_rel),
            "tau_s": float(comp_tau_s),
            "df": float(comp_df),
            "clip_rel": float(comp_clip_rel),
            "seed_offset": int(comp_seed_offset),
        },
        "storm_jump": {
            "enabled": bool(storm_on),
            "rate_per_day": float(storm_rate_per_day),
            "duration_s": float(storm_duration_s),
            "jump_sigma_rel": float(storm_sigma_rel),
            "jump_mean_rel": float(storm_mean_rel),
            "df": float(storm_df),
            "clip_rel": float(storm_clip_rel),
            "seed_offset": int(storm_seed_offset),
        },
        "applied_ao_erosion": False,
        "ao_erosion": {
            "enabled": bool(ao_on),
            "sigma_rel": float(ao_sigma_rel),
            "tau_s": float(ao_tau_s),
            "df": float(ao_df),
            "clip_rel": float(ao_clip_rel),
            "rate_rel_per_day": float(ao_rate_rel_per_day),
            "seed_offset": int(ao_seed_offset),
        },
        "applied_plasma_charging": False,
        "plasma_charging": {
            "enabled": bool(plasma_on),
            "sigma_rel": float(plasma_sigma_rel),
            "tau_s": float(plasma_tau_s),
            "df": float(plasma_df),
            "clip_rel": float(plasma_clip_rel),
            "seed_offset": int(plasma_seed_offset),
        },
        "applied_structural_flex": False,
        "structural_flex": {
            "enabled": bool(structural_on),
            "sigma_rel": float(structural_sigma_rel),
            "tau_s": float(structural_tau_s),
            "df": float(structural_df),
            "clip_rel": float(structural_clip_rel),
            "seed_offset": int(structural_seed_offset),
        },
        "applied_rarefied_regime_switch": False,
        "rarefied_regime_switch": {
            "enabled": bool(rarefied_on),
            "rate_per_hour": float(rarefied_rate_per_hour),
            "mean_duration_s": float(rarefied_mean_duration_s),
            "alt_temperature_ratio_method": int(rarefied_alt_method),
            "seed_offset": int(rarefied_seed_offset),
        },
    }
    name = _scenario_name(scenario)
    base = np.uint64(seed) ^ np.uint64(_stable_u64(name))
    density_scale = np.ones_like(t, dtype=float)
    density_log = np.zeros_like(t, dtype=float)

    if enabled and sigma_rel > 0.0:
        seed_mix = base ^ np.uint64(seed_offset) ^ np.uint64(0x9E3779B97F4A7C15)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        scale, log_scale = density_discrepancy_scale_series(
            t,
            rng=rng,
            sigma_rel=sigma_rel,
            tau_s=tau_s,
            dof=dof,
            clip_rel=clip_rel,
        )
        density_scale *= scale
        density_log += log_scale
        meta["applied_density"] = True
        meta["seed"] = int(rng_seed)
        meta["scale_mean"] = float(np.mean(scale))
        meta["scale_std"] = float(np.std(scale))
        meta["scale_min"] = float(np.min(scale))
        meta["scale_max"] = float(np.max(scale))
        meta["log_scale_std"] = float(np.std(log_scale))

    storm_log = np.zeros_like(t, dtype=float)
    if storm_on and storm_rate_per_day > 0.0 and storm_duration_s > 0.0:
        seed_mix = base ^ np.uint64(storm_seed_offset) ^ np.uint64(0xD1B54A32D192ED03)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        storm_log, storm_meta = storm_jump_log_series(
            t,
            rng=rng,
            enabled=True,
            rate_per_day=storm_rate_per_day,
            duration_s=storm_duration_s,
            jump_sigma_rel=storm_sigma_rel,
            jump_mean_rel=storm_mean_rel,
            dof=storm_df,
            clip_rel=storm_clip_rel,
        )
        storm_scale = np.exp(storm_log)
        density_scale *= storm_scale
        density_log += storm_log
        meta["applied_storm_jump"] = bool(storm_meta.get("event_count", 0) > 0)
        storm_meta["seed"] = int(rng_seed)
        meta["storm_jump"] = storm_meta

    if meta["applied_density"] or meta["applied_storm_jump"]:
        for i, entry in enumerate(env):
            entry.density = float(entry.density) * float(density_scale[i])
        meta["density_total_scale_mean"] = float(np.mean(density_scale))
        meta["density_total_scale_std"] = float(np.std(density_scale))
        meta["density_total_scale_min"] = float(np.min(density_scale))
        meta["density_total_scale_max"] = float(np.max(density_scale))
        meta["density_total_log_std"] = float(np.std(density_log))

    if comp_on and comp_sigma_rel > 0.0:
        seed_mix = base ^ np.uint64(comp_seed_offset) ^ np.uint64(0x94D049BB133111EB)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        comp_scale, comp_log = density_discrepancy_scale_series(
            t,
            rng=rng,
            sigma_rel=comp_sigma_rel,
            tau_s=comp_tau_s,
            dof=comp_df,
            clip_rel=comp_clip_rel,
        )
        missing_attr = 0
        for i, entry in enumerate(env):
            if hasattr(entry, "particles_mass_kg"):
                entry.particles_mass_kg = float(entry.particles_mass_kg) * float(comp_scale[i])
            else:
                missing_attr += 1
        meta["applied_composition"] = True
        meta["composition"].update(
            {
                "seed": int(rng_seed),
                "scale_mean": float(np.mean(comp_scale)),
                "scale_std": float(np.std(comp_scale)),
                "scale_min": float(np.min(comp_scale)),
                "scale_max": float(np.max(comp_scale)),
                "log_scale_std": float(np.std(comp_log)),
                "missing_particles_mass_attr": int(missing_attr),
            }
        )

    if ao_on and (ao_sigma_rel > 0.0 or abs(ao_rate_rel_per_day) > 0.0):
        seed_mix = base ^ np.uint64(ao_seed_offset) ^ np.uint64(0xA24BAED4963EE407)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        ao_scale, ao_log = density_discrepancy_scale_series(
            t,
            rng=rng,
            sigma_rel=ao_sigma_rel,
            tau_s=ao_tau_s,
            dof=ao_df,
            clip_rel=ao_clip_rel,
        )
        if abs(ao_rate_rel_per_day) > 0.0:
            drift = 1.0 + float(ao_rate_rel_per_day) * (t / 86400.0)
            drift = np.maximum(drift, 1e-6)
            if ao_clip_rel > 0.0:
                drift = np.clip(drift, max(1e-6, 1.0 - ao_clip_rel), 1.0 + ao_clip_rel)
            ao_scale = ao_scale * drift
        for i, entry in enumerate(env):
            scale_i = float(max(1e-6, ao_scale[i]))
            if hasattr(entry, "eta1_rad"):
                entry.eta1_rad = float(np.clip(float(entry.eta1_rad) * scale_i, 0.0, 0.5 * np.pi))
            if hasattr(entry, "eta2_rad"):
                entry.eta2_rad = float(np.clip(float(entry.eta2_rad) * scale_i, 0.0, 0.5 * np.pi))
        meta["applied_ao_erosion"] = True
        meta["ao_erosion"].update(
            {
                "seed": int(rng_seed),
                "scale_mean": float(np.mean(ao_scale)),
                "scale_std": float(np.std(ao_scale)),
                "scale_min": float(np.min(ao_scale)),
                "scale_max": float(np.max(ao_scale)),
                "log_scale_std": float(np.std(ao_log)),
            }
        )

    if plasma_on and plasma_sigma_rel > 0.0:
        seed_mix = base ^ np.uint64(plasma_seed_offset) ^ np.uint64(0x9E6C63D0676A9A99)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        plasma_scale, plasma_log = density_discrepancy_scale_series(
            t,
            rng=rng,
            sigma_rel=plasma_sigma_rel,
            tau_s=plasma_tau_s,
            dof=plasma_df,
            clip_rel=plasma_clip_rel,
        )
        for i, entry in enumerate(env):
            if hasattr(entry, "magnetic_field_I_T"):
                b = np.asarray(entry.magnetic_field_I_T, dtype=float).reshape(3)
                entry.magnetic_field_I_T = b * float(plasma_scale[i])
        meta["applied_plasma_charging"] = True
        meta["plasma_charging"].update(
            {
                "seed": int(rng_seed),
                "scale_mean": float(np.mean(plasma_scale)),
                "scale_std": float(np.std(plasma_scale)),
                "scale_min": float(np.min(plasma_scale)),
                "scale_max": float(np.max(plasma_scale)),
                "log_scale_std": float(np.std(plasma_log)),
            }
        )

    if structural_on and structural_sigma_rel > 0.0:
        seed_mix = base ^ np.uint64(structural_seed_offset) ^ np.uint64(0xC2B2AE3D27D4EB4F)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        structural_scale, structural_log = density_discrepancy_scale_series(
            t,
            rng=rng,
            sigma_rel=structural_sigma_rel,
            tau_s=structural_tau_s,
            dof=structural_df,
            clip_rel=structural_clip_rel,
        )
        for i, entry in enumerate(env):
            scale_i = float(max(1e-6, structural_scale[i]))
            if hasattr(entry, "srp_scale"):
                entry.srp_scale = float(max(0.0, float(entry.srp_scale) * scale_i))
            if hasattr(entry, "eta1_rad"):
                entry.eta1_rad = float(np.clip(float(entry.eta1_rad) * np.sqrt(scale_i), 0.0, 0.5 * np.pi))
            if hasattr(entry, "eta2_rad"):
                entry.eta2_rad = float(np.clip(float(entry.eta2_rad) * np.sqrt(scale_i), 0.0, 0.5 * np.pi))
        meta["applied_structural_flex"] = True
        meta["structural_flex"].update(
            {
                "seed": int(rng_seed),
                "scale_mean": float(np.mean(structural_scale)),
                "scale_std": float(np.std(structural_scale)),
                "scale_min": float(np.min(structural_scale)),
                "scale_max": float(np.max(structural_scale)),
                "log_scale_std": float(np.std(structural_log)),
            }
        )

    if rarefied_on and rarefied_rate_per_hour > 0.0 and rarefied_mean_duration_s > 0.0:
        seed_mix = base ^ np.uint64(rarefied_seed_offset) ^ np.uint64(0xC6BC279692B5CC83)
        rng_seed = int(seed_mix.item())
        rng = np.random.default_rng(rng_seed)
        mask = _sample_switch_active_mask(
            t,
            rng=rng,
            start_rate_per_hour=float(rarefied_rate_per_hour),
            mean_duration_s=float(rarefied_mean_duration_s),
        )
        switched = 0
        for i, entry in enumerate(env):
            if i < mask.size and bool(mask[i]) and hasattr(entry, "temperature_ratio_method"):
                entry.temperature_ratio_method = int(rarefied_alt_method)
                switched += 1
        meta["applied_rarefied_regime_switch"] = bool(switched > 0)
        meta["rarefied_regime_switch"].update(
            {
                "seed": int(rng_seed),
                "switched_count": int(switched),
                "switched_fraction": float(switched / max(len(env), 1)),
            }
        )

    meta["name"] = name
    meta["applied"] = bool(
        meta["applied_density"]
        or meta["applied_storm_jump"]
        or meta["applied_composition"]
        or meta["applied_ao_erosion"]
        or meta["applied_plasma_charging"]
        or meta["applied_structural_flex"]
        or meta["applied_rarefied_regime_switch"]
    )
    return meta
