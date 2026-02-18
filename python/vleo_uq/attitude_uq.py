"""Attitude UQ: fast error-state model and payload impact mapping utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np


@dataclass
class AttitudeErrorSeries:
    t_grid: np.ndarray
    errors_rad: np.ndarray


def _expand_vec3(value: Sequence[float] | float) -> np.ndarray:
    if np.isscalar(value):
        return np.full(3, float(value))
    arr = np.array(value, dtype=float)
    if arr.shape != (3,):
        raise ValueError("value must be scalar or length-3")
    return arr


def _skew(v: np.ndarray) -> np.ndarray:
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ]
    )


def rotvec_to_matrix(rotvec: np.ndarray) -> np.ndarray:
    theta = np.linalg.norm(rotvec)
    if theta < 1e-12:
        return np.eye(3) + _skew(rotvec)
    k = rotvec / theta
    K = _skew(k)
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def quat_wxyz_to_dcm(q_wxyz: np.ndarray) -> np.ndarray:
    q = np.array(q_wxyz, dtype=float)
    if q.shape != (4,):
        raise ValueError("q_wxyz must be length-4")
    norm = np.linalg.norm(q)
    if norm == 0.0:
        raise ValueError("quaternion must be non-zero")
    q = q / norm
    w, x, y, z = q
    return np.array(
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
            [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
            [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)],
        ]
    )


def dcm_to_quat_wxyz(R: np.ndarray) -> np.ndarray:
    R = np.array(R, dtype=float)
    if R.shape != (3, 3):
        raise ValueError("R must be 3x3")
    trace = np.trace(R)
    if trace > 0.0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (R[2, 1] - R[1, 2]) * s
        y = (R[0, 2] - R[2, 0]) * s
        z = (R[1, 0] - R[0, 1]) * s
    else:
        if R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
            w = (R[2, 1] - R[1, 2]) / s
            x = 0.25 * s
            y = (R[0, 1] + R[1, 0]) / s
            z = (R[0, 2] + R[2, 0]) / s
        elif R[1, 1] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
            w = (R[0, 2] - R[2, 0]) / s
            x = (R[0, 1] + R[1, 0]) / s
            y = 0.25 * s
            z = (R[1, 2] + R[2, 1]) / s
        else:
            s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
            w = (R[1, 0] - R[0, 1]) / s
            x = (R[0, 2] + R[2, 0]) / s
            y = (R[1, 2] + R[2, 1]) / s
            z = 0.25 * s
    q = np.array([w, x, y, z], dtype=float)
    return q / np.linalg.norm(q)


def simulate_attitude_ou(
    t_grid: np.ndarray,
    tau_s: Sequence[float] | float,
    sigma_rad: Sequence[float] | float,
    rng: Optional[np.random.Generator] = None,
    x0: Optional[np.ndarray] = None,
) -> AttitudeErrorSeries:
    if rng is None:
        rng = np.random.default_rng(0)
    tau = _expand_vec3(tau_s)
    sigma = _expand_vec3(sigma_rad)
    if np.any(tau < 0.0):
        raise ValueError("tau_s must be non-negative")
    if np.any(sigma < 0.0):
        raise ValueError("sigma_rad must be non-negative")
    t_grid = np.array(t_grid, dtype=float)
    if t_grid.ndim != 1 or t_grid.size < 2:
        raise ValueError("t_grid must be 1D with at least 2 entries")
    if np.any(np.diff(t_grid) <= 0.0):
        raise ValueError("t_grid must be strictly increasing")

    errors = np.zeros((t_grid.size, 3))
    if x0 is not None:
        x0 = np.array(x0, dtype=float)
        if x0.shape != (3,):
            raise ValueError("x0 must be length-3")
        errors[0] = x0

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]
        for k in range(3):
            if tau[k] > 0.0:
                phi = np.exp(-dt / tau[k])
                var = 1.0 - phi * phi
                errors[i, k] = errors[i - 1, k] * phi
                if sigma[k] > 0.0 and var > 0.0:
                    errors[i, k] += sigma[k] * np.sqrt(var) * rng.normal()
            else:
                errors[i, k] = errors[i - 1, k]
                if sigma[k] > 0.0:
                    errors[i, k] += sigma[k] * np.sqrt(dt) * rng.normal()
    return AttitudeErrorSeries(t_grid=t_grid, errors_rad=errors)


def add_white_jitter(
    errors_rad: np.ndarray,
    sigma_rad: Sequence[float] | float,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng(0)
    sigma = _expand_vec3(sigma_rad)
    jitter = rng.normal(0.0, 1.0, size=errors_rad.shape) * sigma.reshape(1, 3)
    return errors_rad + jitter


def boresight_eci_from_states(
    states: np.ndarray,
    boresight_body: Optional[np.ndarray] = None,
) -> np.ndarray:
    if boresight_body is None:
        boresight_body = np.array([0.0, 0.0, 1.0])
    boresight_body = np.array(boresight_body, dtype=float)
    if boresight_body.shape != (3,):
        raise ValueError("boresight_body must be length-3")
    out = np.zeros((states.shape[0], 3))
    for i in range(states.shape[0]):
        q = states[i, 6:10]
        R_BI = quat_wxyz_to_dcm(q)
        out[i] = R_BI.T @ boresight_body
    return out


def boresight_eci_with_errors(
    states: np.ndarray,
    error_rad: np.ndarray,
    boresight_body: Optional[np.ndarray] = None,
) -> np.ndarray:
    if states.shape[0] != error_rad.shape[0]:
        raise ValueError("states and error_rad must have same length")
    if boresight_body is None:
        boresight_body = np.array([0.0, 0.0, 1.0])
    boresight_body = np.array(boresight_body, dtype=float)
    out = np.zeros((states.shape[0], 3))
    for i in range(states.shape[0]):
        q = states[i, 6:10]
        R_BI = quat_wxyz_to_dcm(q)
        R_err = rotvec_to_matrix(error_rad[i])
        b_body = R_err @ boresight_body
        out[i] = R_BI.T @ b_body
    return out


def ground_intersect_spherical(
    r_sc_eci: np.ndarray,
    los_eci: np.ndarray,
    earth_radius_m: float = 6378137.0,
) -> np.ndarray:
    r = np.array(r_sc_eci, dtype=float)
    l = np.array(los_eci, dtype=float)
    l_norm = np.linalg.norm(l)
    if l_norm == 0.0:
        raise ValueError("los_eci must be non-zero")
    l = l / l_norm
    b = np.dot(r, l)
    c = np.dot(r, r) - earth_radius_m * earth_radius_m
    disc = b * b - c
    if disc < 0.0:
        return np.full(3, np.nan)
    s = -b - np.sqrt(disc)
    if s <= 0.0:
        s = -b + np.sqrt(disc)
    if s <= 0.0:
        return np.full(3, np.nan)
    return r + s * l


def ground_points_from_los(
    r_sc_eci: np.ndarray,
    los_eci: np.ndarray,
    earth_radius_m: float = 6378137.0,
) -> np.ndarray:
    points = np.zeros_like(r_sc_eci)
    for i in range(r_sc_eci.shape[0]):
        points[i] = ground_intersect_spherical(r_sc_eci[i], los_eci[i], earth_radius_m)
    return points


def geolocation_error_m(
    ground_nominal: np.ndarray,
    ground_actual: np.ndarray,
    earth_radius_m: float = 6378137.0,
) -> np.ndarray:
    dots = np.sum(ground_nominal * ground_actual, axis=1)
    cos_angle = dots / (earth_radius_m * earth_radius_m)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    angle = np.arccos(cos_angle)
    return earth_radius_m * angle


def range_bias_m(
    r_sc_eci: np.ndarray,
    los_eci: np.ndarray,
    earth_radius_m: float = 6378137.0,
) -> np.ndarray:
    out = np.zeros(r_sc_eci.shape[0])
    for i in range(r_sc_eci.shape[0]):
        r = r_sc_eci[i]
        los = los_eci[i]
        ground = ground_intersect_spherical(r, los, earth_radius_m)
        if np.any(np.isnan(ground)):
            out[i] = np.nan
            continue
        range_actual = np.linalg.norm(ground - r)
        range_nominal = np.linalg.norm(r) - earth_radius_m
        out[i] = range_actual - range_nominal
    return out


def smear_over_exposure_m(
    t_grid: np.ndarray,
    ground_points: np.ndarray,
    exposure_s: float,
) -> np.ndarray:
    t_grid = np.array(t_grid, dtype=float)
    if exposure_s <= 0.0:
        raise ValueError("exposure_s must be positive")
    out = np.zeros(t_grid.size)
    for i, t in enumerate(t_grid):
        start = t - 0.5 * exposure_s
        end = t + 0.5 * exposure_s
        mask = (t_grid >= start) & (t_grid <= end)
        pts = ground_points[mask]
        if pts.shape[0] < 2:
            out[i] = 0.0
            continue
        center = np.mean(pts, axis=0)
        d = np.linalg.norm(pts - center, axis=1)
        out[i] = np.max(d)
    return out
