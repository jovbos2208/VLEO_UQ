"""Discrete propagation event utilities."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class PropagationEvent:
    t_s: float
    kind: str
    params: dict


def sample_event_realizations(
    events: Sequence[PropagationEvent],
    rng: np.random.Generator,
) -> list[PropagationEvent]:
    """Sample uncertain event timing/parameters.

    Supported stochastic keys:
    - ``sigma_time_s`` for Gaussian event-time jitter.
    - ``sigma_scale`` for multiplicative Gaussian scaling of ``dv_eci_m_s`` or
      ``accel_eci_m_s2``.
    - ``sigma_duration_s`` for Gaussian duration jitter of ``finite_burn`` events.
    """
    sampled: list[PropagationEvent] = []
    for ev in events:
        t_s = float(ev.t_s)
        params = dict(ev.params)
        sigma_time = float(params.pop("sigma_time_s", 0.0) or 0.0)
        if sigma_time > 0.0:
            t_s = t_s + sigma_time * rng.normal()

        if ev.kind == "delta_v" and "dv_eci_m_s" in params:
            sigma_scale = float(params.pop("sigma_scale", 0.0) or 0.0)
            dv = np.array(params["dv_eci_m_s"], dtype=float).reshape(3)
            if sigma_scale > 0.0:
                dv = dv * (1.0 + sigma_scale * rng.normal())
            params["dv_eci_m_s"] = dv
        elif ev.kind == "finite_burn":
            sigma_scale = float(params.pop("sigma_scale", 0.0) or 0.0)
            sigma_duration = float(params.pop("sigma_duration_s", 0.0) or 0.0)
            if "accel_eci_m_s2" in params:
                accel = np.array(params["accel_eci_m_s2"], dtype=float).reshape(3)
                if sigma_scale > 0.0:
                    accel = accel * (1.0 + sigma_scale * rng.normal())
                params["accel_eci_m_s2"] = accel
            elif "dv_eci_m_s" in params:
                dv = np.array(params["dv_eci_m_s"], dtype=float).reshape(3)
                if sigma_scale > 0.0:
                    dv = dv * (1.0 + sigma_scale * rng.normal())
                params["dv_eci_m_s"] = dv
            duration = float(params.get("duration_s", 0.0) or 0.0)
            if sigma_duration > 0.0:
                duration = max(0.0, duration + sigma_duration * rng.normal())
            params["duration_s"] = duration

        sampled.append(replace(ev, t_s=t_s, params=params))
    return sampled


def _finite_burn_pulses(
    t_grid: np.ndarray,
    event: PropagationEvent,
) -> list[tuple[int, PropagationEvent]]:
    params = dict(event.params)
    duration_s = float(params.get("duration_s", 0.0) or 0.0)
    if duration_s <= 0.0:
        return []
    start_t = float(event.t_s)
    end_t = start_t + duration_s
    if end_t <= float(t_grid[0]) or start_t >= float(t_grid[-1]):
        return []

    if "accel_eci_m_s2" in params:
        accel = np.array(params["accel_eci_m_s2"], dtype=float).reshape(3)
    elif "dv_eci_m_s" in params:
        dv = np.array(params["dv_eci_m_s"], dtype=float).reshape(3)
        accel = dv / duration_s
    else:
        raise ValueError("finite_burn requires accel_eci_m_s2 or dv_eci_m_s")

    pulses: list[tuple[int, PropagationEvent]] = []
    for i in range(t_grid.size - 1):
        t0 = float(t_grid[i])
        t1 = float(t_grid[i + 1])
        dt_eff = max(0.0, min(t1, end_t) - max(t0, start_t))
        if dt_eff <= 0.0:
            continue
        dv_seg = accel * dt_eff
        pulses.append(
            (
                i + 1,
                PropagationEvent(
                    t_s=t1,
                    kind="delta_v",
                    params={"dv_eci_m_s": dv_seg},
                ),
            )
        )
    return pulses


def event_indices(
    t_grid: np.ndarray,
    events: Sequence[PropagationEvent],
    *,
    time_tol_s: float = 1e-9,
) -> dict[int, list[PropagationEvent]]:
    if time_tol_s < 0.0:
        raise ValueError("time_tol_s must be non-negative")
    t_grid = np.asarray(t_grid, dtype=float)
    out: dict[int, list[PropagationEvent]] = {}
    for ev in events:
        if ev.kind == "finite_burn":
            for idx, pulse in _finite_burn_pulses(t_grid, ev):
                out.setdefault(int(idx), []).append(pulse)
            continue
        idx = int(np.argmin(np.abs(t_grid - ev.t_s)))
        if abs(float(t_grid[idx]) - float(ev.t_s)) > time_tol_s:
            raise ValueError(
                f"event at t={ev.t_s} does not align with t_grid within tolerance {time_tol_s}"
            )
        out.setdefault(idx, []).append(ev)
    return out


def apply_event_to_state(x: np.ndarray, event: PropagationEvent) -> np.ndarray:
    x_new = np.array(x, dtype=float, copy=True)
    if event.kind == "delta_v":
        dv = np.array(event.params.get("dv_eci_m_s"), dtype=float).reshape(3)
        x_new[3:6] += dv
        return x_new
    if event.kind == "finite_burn":
        params = dict(event.params)
        duration_s = float(params.get("duration_s", 0.0) or 0.0)
        if duration_s <= 0.0:
            return x_new
        if "dv_eci_m_s" in params:
            dv = np.array(params.get("dv_eci_m_s"), dtype=float).reshape(3)
        elif "accel_eci_m_s2" in params:
            accel = np.array(params.get("accel_eci_m_s2"), dtype=float).reshape(3)
            dv = accel * duration_s
        else:
            raise ValueError("finite_burn requires dv_eci_m_s or accel_eci_m_s2")
        x_new[3:6] += dv
        return x_new
    if event.kind == "attitude_reset":
        q = np.array(event.params.get("q_wxyz"), dtype=float).reshape(4)
        qn = np.linalg.norm(q)
        if qn == 0.0:
            raise ValueError("attitude_reset requires non-zero q_wxyz")
        x_new[6:10] = q / qn
        if "w_BI_B" in event.params:
            x_new[10:13] = np.array(event.params["w_BI_B"], dtype=float).reshape(3)
        return x_new
    raise ValueError(f"unsupported event kind '{event.kind}'")


def propagate_with_events(
    prop_det,
    x0: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    events: Iterable[PropagationEvent],
    *,
    time_tol_s: float = 1e-9,
) -> np.ndarray:
    """Propagate deterministically with state jumps at event epochs."""
    t_grid = np.asarray(t_grid, dtype=float)
    if t_grid.ndim != 1 or t_grid.size < 2:
        raise ValueError("t_grid must be 1D with at least 2 entries")
    if len(env) != t_grid.size:
        raise ValueError("env length must match t_grid")

    ev_map = event_indices(t_grid, list(events), time_tol_s=time_tol_s)
    x = np.array(x0, dtype=float, copy=True)
    out = np.zeros((t_grid.size, x.size), dtype=float)
    out[0] = x

    for i in range(t_grid.size - 1):
        seg = prop_det.propagate(x, t_grid[i : i + 2], env[i : i + 2])
        x = np.array(seg[-1], dtype=float, copy=True)
        for ev in ev_map.get(i + 1, []):
            x = apply_event_to_state(x, ev)
        out[i + 1] = x
    return out
