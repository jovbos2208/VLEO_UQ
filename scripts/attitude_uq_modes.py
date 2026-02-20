from __future__ import annotations

from typing import Any

import numpy as np


def _normalize_quat(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    n = np.linalg.norm(q, axis=-1, keepdims=True)
    n = np.where(n > 0.0, n, 1.0)
    return q / n


def _quat_angle_error_deg(q_ref: np.ndarray, q_cmp: np.ndarray) -> np.ndarray:
    qa = _normalize_quat(q_ref)
    qb = _normalize_quat(q_cmp)
    dots = np.clip(np.abs(np.sum(qa * qb, axis=-1)), -1.0, 1.0)
    return np.rad2deg(2.0 * np.arccos(dots))


def _expand3(value: Any, default: float) -> np.ndarray:
    if value is None:
        return np.full(3, float(default), dtype=float)
    if np.isscalar(value):
        return np.full(3, float(value), dtype=float)
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size == 3:
        return arr
    return np.full(3, float(default), dtype=float)


def _simulate_attitude_fast_ou(
    t_grid: np.ndarray,
    *,
    tau_s: np.ndarray,
    sigma_rad: np.ndarray,
    seed: int,
) -> np.ndarray:
    t = np.asarray(t_grid, dtype=float).reshape(-1)
    n = t.size
    if n == 0:
        return np.zeros((0, 3), dtype=float)
    rng = np.random.default_rng(int(seed))
    x = np.zeros((n, 3), dtype=float)
    for i in range(1, n):
        dt = max(0.0, float(t[i] - t[i - 1]))
        for k in range(3):
            tau = max(0.0, float(tau_s[k]))
            sig = max(0.0, float(sigma_rad[k]))
            if tau > 0.0:
                phi = np.exp(-dt / tau)
                var = max(0.0, 1.0 - phi * phi)
                x[i, k] = phi * x[i - 1, k] + sig * np.sqrt(var) * rng.normal()
            else:
                x[i, k] = x[i - 1, k] + sig * np.sqrt(dt) * rng.normal()
    return x


def fast_attitude_mode_metrics(
    *,
    t_grid: np.ndarray,
    seed: int,
    tau_s: Any = 200.0,
    sigma_rad: Any = 1e-4,
) -> dict[str, Any]:
    t = np.asarray(t_grid, dtype=float).reshape(-1)
    if t.size < 2:
        return {"available": False, "reason": "short_time_grid"}
    tau = _expand3(tau_s, 200.0)
    sig = _expand3(sigma_rad, 1e-4)
    err = _simulate_attitude_fast_ou(t, tau_s=tau, sigma_rad=sig, seed=int(seed))
    norm = np.linalg.norm(err, axis=1)
    return {
        "available": True,
        "mode": "fast_ou_error_state",
        "seed": int(seed),
        "tau_s": [float(x) for x in tau],
        "sigma_rad": [float(x) for x in sig],
        "rms_error_deg": float(np.rad2deg(np.sqrt(np.mean(norm * norm)))),
        "p95_error_deg": float(np.rad2deg(np.percentile(norm, 95.0))),
        "final_error_deg": float(np.rad2deg(norm[-1])),
    }


def full_attitude_mode_metrics(
    *,
    det_states: np.ndarray,
    mc_states: np.ndarray,
) -> dict[str, Any]:
    det = np.asarray(det_states, dtype=float)
    mc = np.asarray(mc_states, dtype=float)
    if det.ndim != 2 or mc.ndim != 3:
        return {"available": False, "reason": "invalid_shapes"}
    if det.shape[0] != mc.shape[0]:
        return {"available": False, "reason": "time_mismatch"}
    if det.shape[1] < 10 or mc.shape[2] < 10:
        return {"available": False, "reason": "state_size_lt_10"}
    q_det = np.asarray(det[:, 6:10], dtype=float)
    q_mc = np.asarray(mc[:, :, 6:10], dtype=float)
    n_t, n_p, _ = q_mc.shape
    err = np.zeros((n_t, n_p), dtype=float)
    for i in range(n_t):
        err[i] = _quat_angle_error_deg(np.repeat(q_det[i : i + 1], n_p, axis=0), q_mc[i])
    final = err[-1]
    return {
        "available": True,
        "mode": "full_fidelity_ensemble",
        "particles": int(n_p),
        "rms_error_deg": float(np.sqrt(np.mean(err * err))),
        "p95_error_deg": float(np.percentile(err, 95.0)),
        "final_mean_error_deg": float(np.mean(final)),
        "final_p95_error_deg": float(np.percentile(final, 95.0)),
        "final_max_error_deg": float(np.max(final)),
    }


def summarize_attitude_uq_modes(
    *,
    t_grid: np.ndarray,
    det_states: np.ndarray,
    mc_states: np.ndarray,
    seed: int,
    fast_tau_s: Any = 200.0,
    fast_sigma_rad: Any = 1e-4,
) -> dict[str, Any]:
    fast = fast_attitude_mode_metrics(
        t_grid=t_grid,
        seed=int(seed),
        tau_s=fast_tau_s,
        sigma_rad=fast_sigma_rad,
    )
    full = full_attitude_mode_metrics(det_states=det_states, mc_states=mc_states)
    return {
        "available": bool(fast.get("available", False) and full.get("available", False)),
        "fast_mode": fast,
        "full_mode": full,
    }
