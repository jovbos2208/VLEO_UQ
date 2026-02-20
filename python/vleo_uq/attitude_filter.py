"""Attitude dynamics, sensor simulation, and a simple MEKF."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from .attitude_uq import (
    dcm_to_quat_wxyz,
    quat_wxyz_to_dcm,
    rotvec_to_matrix,
)
from ._vleo_uq import (
    AeroAdapter,
    EnsemblePropagatorMC,
    EnvInputs,
    SigmaPointPropagatorUT,
    StmPropagator,
    VehicleParams,
)


@dataclass
class AttitudeTruth:
    t_grid: np.ndarray
    q_wxyz: np.ndarray
    w_rad_s: np.ndarray
    gyro_bias_rad_s: np.ndarray


@dataclass
class AttitudeFilterResult:
    q_wxyz: np.ndarray
    gyro_bias_rad_s: np.ndarray
    cov: np.ndarray


@dataclass
class FullStateFilterResult:
    states: np.ndarray
    cov: np.ndarray


def _skew(v: np.ndarray) -> np.ndarray:
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ]
    )


def quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ],
        dtype=float,
    )


def quat_conjugate(q: np.ndarray) -> np.ndarray:
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=float)


def quat_normalize(q: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(q)
    if norm == 0.0:
        raise ValueError("quaternion must be non-zero")
    return q / norm


def small_angle_to_quat(delta: np.ndarray) -> np.ndarray:
    angle = np.linalg.norm(delta)
    if angle < 1e-12:
        return np.array([1.0, 0.5 * delta[0], 0.5 * delta[1], 0.5 * delta[2]])
    axis = delta / angle
    half = 0.5 * angle
    return np.array(
        [np.cos(half), *(np.sin(half) * axis)],
        dtype=float,
    )


def quat_error_small(q_meas: np.ndarray, q_pred: np.ndarray) -> np.ndarray:
    q_err = quat_multiply(q_meas, quat_conjugate(q_pred))
    if q_err[0] < 0.0:
        q_err = -q_err
    return 2.0 * q_err[1:4]


def integrate_quat(q: np.ndarray, w: np.ndarray, dt: float) -> np.ndarray:
    angle = np.linalg.norm(w) * dt
    if angle < 1e-12:
        dq = np.array([1.0, 0.5 * w[0] * dt, 0.5 * w[1] * dt, 0.5 * w[2] * dt])
    else:
        axis = w / np.linalg.norm(w)
        half = 0.5 * angle
        dq = np.array([np.cos(half), *(np.sin(half) * axis)], dtype=float)
    return quat_normalize(quat_multiply(dq, q))


def quat_derivative_BI(q: np.ndarray, w: np.ndarray) -> np.ndarray:
    w_quat = np.array([0.0, w[0], w[1], w[2]])
    qdot = quat_multiply(q, w_quat)
    return -0.5 * qdot


def simulate_attitude_truth(
    t_grid: np.ndarray,
    q0_wxyz: np.ndarray,
    w0_rad_s: np.ndarray,
    inertia_kg_m2: np.ndarray,
    torque_fn: Optional[Callable[[float, np.ndarray, np.ndarray], np.ndarray]] = None,
    bias_tau_s: float = 0.0,
    bias_sigma_rad_s: Sequence[float] | float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> AttitudeTruth:
    if rng is None:
        rng = np.random.default_rng(0)
    t_grid = np.array(t_grid, dtype=float)
    if np.any(np.diff(t_grid) <= 0.0):
        raise ValueError("t_grid must be strictly increasing")
    inertia = np.array(inertia_kg_m2, dtype=float)
    if inertia.shape != (3, 3):
        raise ValueError("inertia_kg_m2 must be 3x3")
    inv_inertia = np.linalg.inv(inertia)

    q = np.zeros((t_grid.size, 4))
    w = np.zeros((t_grid.size, 3))
    b = np.zeros((t_grid.size, 3))
    q[0] = quat_normalize(np.array(q0_wxyz, dtype=float))
    w[0] = np.array(w0_rad_s, dtype=float)
    b[0] = 0.0

    sigma = np.array(bias_sigma_rad_s if np.isscalar(bias_sigma_rad_s) else bias_sigma_rad_s, dtype=float)
    if sigma.shape == ():
        sigma = np.full(3, float(sigma))
    if sigma.shape != (3,):
        raise ValueError("bias_sigma_rad_s must be scalar or length-3")

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]
        torque = np.zeros(3)
        if torque_fn is not None:
            torque = np.array(torque_fn(t_grid[i - 1], q[i - 1], w[i - 1]), dtype=float)
        wdot = inv_inertia @ (torque - np.cross(w[i - 1], inertia @ w[i - 1]))
        w_mid = w[i - 1] + 0.5 * dt * wdot
        q[i] = integrate_quat(q[i - 1], w_mid, dt)
        w[i] = w[i - 1] + dt * wdot

        if bias_tau_s > 0.0:
            phi = np.exp(-dt / bias_tau_s)
            var = 1.0 - phi * phi
            b[i] = b[i - 1] * phi
            if np.any(sigma > 0.0) and var > 0.0:
                b[i] += sigma * np.sqrt(var) * rng.normal(size=3)
        else:
            b[i] = b[i - 1]
            if np.any(sigma > 0.0):
                b[i] += sigma * np.sqrt(dt) * rng.normal(size=3)

    return AttitudeTruth(t_grid=t_grid, q_wxyz=q, w_rad_s=w, gyro_bias_rad_s=b)


def _aero_torque(
    aero: AeroAdapter,
    q_wxyz: np.ndarray,
    w_rad_s: np.ndarray,
    v_I: np.ndarray,
    env: EnvInputs,
) -> np.ndarray:
    force, torque = aero.compute_ft(
        q_wxyz,
        w_rad_s,
        v_I,
        env.wind_I,
        env.density,
        env.temperature_K,
        env.particles_mass_kg,
        env.eta1_rad,
        env.eta2_rad,
        env.temperature_ratio_method,
    )
    return torque


def propagate_attitude_aero_step(
    q_wxyz: np.ndarray,
    w_rad_s: np.ndarray,
    v_I: np.ndarray,
    env: EnvInputs,
    vehicle: VehicleParams,
    aero: AeroAdapter,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    inertia = np.array(vehicle.inertia_B, dtype=float)
    inv_inertia = np.linalg.inv(inertia)

    def f(q, w):
        torque = _aero_torque(aero, q, w, v_I, env)
        wdot = inv_inertia @ (torque - np.cross(w, inertia @ w))
        qdot = quat_derivative_BI(q, w)
        return qdot, wdot

    q0 = quat_normalize(np.array(q_wxyz, dtype=float))
    w0 = np.array(w_rad_s, dtype=float)

    k1_q, k1_w = f(q0, w0)
    q1 = quat_normalize(q0 + 0.5 * dt * k1_q)
    w1 = w0 + 0.5 * dt * k1_w

    k2_q, k2_w = f(q1, w1)
    q2 = quat_normalize(q0 + 0.5 * dt * k2_q)
    w2 = w0 + 0.5 * dt * k2_w

    k3_q, k3_w = f(q2, w2)
    q3 = quat_normalize(q0 + dt * k3_q)
    w3 = w0 + dt * k3_w

    k4_q, k4_w = f(q3, w3)

    q_next = q0 + (dt / 6.0) * (k1_q + 2.0 * k2_q + 2.0 * k3_q + k4_q)
    w_next = w0 + (dt / 6.0) * (k1_w + 2.0 * k2_w + 2.0 * k3_w + k4_w)
    return quat_normalize(q_next), w_next


def simulate_attitude_truth_aero(
    t_grid: np.ndarray,
    q0_wxyz: np.ndarray,
    w0_rad_s: np.ndarray,
    v_I: np.ndarray,
    env: Sequence[EnvInputs],
    vehicle: VehicleParams,
    aero: AeroAdapter,
    bias_tau_s: float = 0.0,
    bias_sigma_rad_s: Sequence[float] | float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> AttitudeTruth:
    if rng is None:
        rng = np.random.default_rng(0)
    t_grid = np.array(t_grid, dtype=float)
    if np.any(np.diff(t_grid) <= 0.0):
        raise ValueError("t_grid must be strictly increasing")
    if v_I.shape[0] != t_grid.size:
        raise ValueError("v_I must align with t_grid")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")

    q = np.zeros((t_grid.size, 4))
    w = np.zeros((t_grid.size, 3))
    b = np.zeros((t_grid.size, 3))
    q[0] = quat_normalize(np.array(q0_wxyz, dtype=float))
    w[0] = np.array(w0_rad_s, dtype=float)
    b[0] = 0.0

    sigma = np.array(bias_sigma_rad_s if np.isscalar(bias_sigma_rad_s) else bias_sigma_rad_s, dtype=float)
    if sigma.shape == ():
        sigma = np.full(3, float(sigma))
    if sigma.shape != (3,):
        raise ValueError("bias_sigma_rad_s must be scalar or length-3")

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]
        q[i], w[i] = propagate_attitude_aero_step(
            q[i - 1],
            w[i - 1],
            v_I[i - 1],
            env[i - 1],
            vehicle,
            aero,
            dt,
        )

        if bias_tau_s > 0.0:
            phi = np.exp(-dt / bias_tau_s)
            var = 1.0 - phi * phi
            b[i] = b[i - 1] * phi
            if np.any(sigma > 0.0) and var > 0.0:
                b[i] += sigma * np.sqrt(var) * rng.normal(size=3)
        else:
            b[i] = b[i - 1]
            if np.any(sigma > 0.0):
                b[i] += sigma * np.sqrt(dt) * rng.normal(size=3)

    return AttitudeTruth(t_grid=t_grid, q_wxyz=q, w_rad_s=w, gyro_bias_rad_s=b)


def simulate_gyro_measurements(
    truth: AttitudeTruth,
    gyro_sigma_rad_s: Sequence[float] | float,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng(0)
    sigma = np.array(gyro_sigma_rad_s if np.isscalar(gyro_sigma_rad_s) else gyro_sigma_rad_s, dtype=float)
    if sigma.shape == ():
        sigma = np.full(3, float(sigma))
    if sigma.shape != (3,):
        raise ValueError("gyro_sigma_rad_s must be scalar or length-3")
    noise = rng.normal(0.0, 1.0, size=truth.w_rad_s.shape) * sigma.reshape(1, 3)
    return truth.w_rad_s + truth.gyro_bias_rad_s + noise


def _measurement_indices(t_grid: np.ndarray, cadence_s: Optional[float]) -> np.ndarray:
    if cadence_s is None:
        return np.arange(t_grid.size, dtype=int)
    indices = [0]
    next_t = t_grid[0] + cadence_s
    for i in range(1, t_grid.size):
        if t_grid[i] >= next_t - 1e-9:
            indices.append(i)
            next_t += cadence_s
    return np.array(indices, dtype=int)


def simulate_star_tracker_measurements(
    truth: AttitudeTruth,
    sigma_rad: float = 1e-4,
    cadence_s: Optional[float] = 1.0,
    rng: Optional[np.random.Generator] = None,
) -> tuple[np.ndarray, np.ndarray]:
    if rng is None:
        rng = np.random.default_rng(0)
    t_grid = truth.t_grid
    indices = _measurement_indices(t_grid, cadence_s)
    meas = np.zeros((indices.size, 4))
    for i, idx in enumerate(indices):
        delta = rng.normal(0.0, sigma_rad, size=3)
        dq = small_angle_to_quat(delta)
        meas[i] = quat_normalize(quat_multiply(dq, truth.q_wxyz[idx]))
    return indices, meas


def _predict_magnetometer_body(
    q_wxyz: np.ndarray,
    env: EnvInputs,
) -> np.ndarray:
    R_BI = quat_wxyz_to_dcm(quat_normalize(np.array(q_wxyz, dtype=float)))
    b_i = np.array(env.magnetic_field_I_T, dtype=float).reshape(3)
    return R_BI @ b_i


def simulate_magnetometer_measurements(
    truth: AttitudeTruth,
    env: Sequence[EnvInputs],
    sigma_T: Sequence[float] | float = 1e-7,
    bias_T: Sequence[float] | float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng(0)
    if len(env) != truth.t_grid.size:
        raise ValueError("env must align with truth timeline")
    sigma = np.array(sigma_T if np.isscalar(sigma_T) else sigma_T, dtype=float)
    if sigma.shape == ():
        sigma = np.full(3, float(sigma))
    if sigma.shape != (3,):
        raise ValueError("sigma_T must be scalar or length-3")
    bias = np.array(bias_T if np.isscalar(bias_T) else bias_T, dtype=float)
    if bias.shape == ():
        bias = np.full(3, float(bias))
    if bias.shape != (3,):
        raise ValueError("bias_T must be scalar or length-3")
    out = np.zeros((truth.t_grid.size, 3), dtype=float)
    for i in range(truth.t_grid.size):
        b_body = _predict_magnetometer_body(truth.q_wxyz[i], env[i])
        out[i] = b_body + bias + sigma * rng.normal(size=3)
    return out


def run_mekf(
    t_grid: np.ndarray,
    gyro_meas: np.ndarray,
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    q0_wxyz: np.ndarray,
    bias0: np.ndarray,
    P0: np.ndarray,
    gyro_noise_std: float = 1e-4,
    bias_rw_std: float = 1e-6,
    meas_noise_std: float = 1e-4,
) -> AttitudeFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    q = np.zeros((t_grid.size, 4))
    b = np.zeros((t_grid.size, 3))
    P = np.zeros((t_grid.size, 6, 6))

    q[0] = quat_normalize(np.array(q0_wxyz, dtype=float))
    b[0] = np.array(bias0, dtype=float)
    P[0] = np.array(P0, dtype=float)

    meas_map = {int(idx): i for i, idx in enumerate(star_indices)}
    I = np.eye(6)

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]
        w_meas = gyro_meas[i - 1] - b[i - 1]
        q[i] = integrate_quat(q[i - 1], w_meas, dt)
        b[i] = b[i - 1]

        F = np.eye(6)
        F[0:3, 0:3] -= _skew(w_meas) * dt
        F[0:3, 3:6] = -np.eye(3) * dt
        Q = np.zeros((6, 6))
        Q[0:3, 0:3] = np.eye(3) * (gyro_noise_std ** 2) * dt
        Q[3:6, 3:6] = np.eye(3) * (bias_rw_std ** 2) * dt
        P[i] = F @ P[i - 1] @ F.T + Q

        if i in meas_map:
            idx = meas_map[i]
            z = star_meas[idx]
            y = quat_error_small(z, q[i])
            H = np.zeros((3, 6))
            H[:, 0:3] = np.eye(3)
            R = np.eye(3) * (meas_noise_std ** 2)
            S = H @ P[i] @ H.T + R
            K = P[i] @ H.T @ np.linalg.inv(S)
            dx = K @ y
            dq = small_angle_to_quat(dx[0:3])
            q[i] = quat_normalize(quat_multiply(dq, q[i]))
            b[i] += dx[3:6]
            P[i] = (I - K @ H) @ P[i] @ (I - K @ H).T + K @ R @ K.T

    return AttitudeFilterResult(q_wxyz=q, gyro_bias_rad_s=b, cov=P)


def run_mekf_aero(
    t_grid: np.ndarray,
    v_I: np.ndarray,
    env: Sequence[EnvInputs],
    aero: AeroAdapter,
    vehicle: VehicleParams,
    gyro_meas: np.ndarray,
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    q0_wxyz: np.ndarray,
    w0_rad_s: np.ndarray,
    bias0: np.ndarray,
    P0: np.ndarray,
    gyro_noise_std: float = 1e-4,
    bias_rw_std: float = 1e-6,
    meas_noise_std: float = 1e-4,
) -> AttitudeFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    if v_I.shape != (t_grid.size, 3):
        raise ValueError("v_I must be (Nt, 3)")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")

    q = np.zeros((t_grid.size, 4))
    w = np.zeros((t_grid.size, 3))
    b = np.zeros((t_grid.size, 3))
    P = np.zeros((t_grid.size, 9, 9))

    q[0] = quat_normalize(np.array(q0_wxyz, dtype=float))
    w[0] = np.array(w0_rad_s, dtype=float)
    b[0] = np.array(bias0, dtype=float)
    P[0] = np.array(P0, dtype=float)

    meas_map = {int(idx): i for i, idx in enumerate(star_indices)}
    I = np.eye(9)

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]

        q_pred, w_pred = propagate_attitude_aero_step(
            q[i - 1],
            w[i - 1],
            v_I[i - 1],
            env[i - 1],
            vehicle,
            aero,
            dt,
        )
        b_pred = b[i - 1]

        F = np.zeros((9, 9))
        eps = 1e-6
        for j in range(9):
            delta = np.zeros(9)
            delta[j] = eps
            dq = small_angle_to_quat(delta[0:3])
            q_pert = quat_normalize(quat_multiply(dq, q[i - 1]))
            w_pert = w[i - 1] + delta[3:6]
            b_pert = b[i - 1] + delta[6:9]

            q_p, w_p = propagate_attitude_aero_step(
                q_pert,
                w_pert,
                v_I[i - 1],
                env[i - 1],
                vehicle,
                aero,
                dt,
            )
            b_p = b_pert

            dtheta = quat_error_small(q_p, q_pred)
            dw = w_p - w_pred
            db = b_p - b_pred
            F[:, j] = np.hstack((dtheta, dw, db)) / eps

        Q = np.zeros((9, 9))
        Q[3:6, 3:6] = np.eye(3) * (gyro_noise_std ** 2) * dt
        Q[6:9, 6:9] = np.eye(3) * (bias_rw_std ** 2) * dt

        P_pred = F @ P[i - 1] @ F.T + Q
        q[i] = q_pred
        w[i] = w_pred
        b[i] = b_pred
        P[i] = P_pred

        if i in meas_map:
            idx = meas_map[i]
            z = star_meas[idx]
            y = quat_error_small(z, q[i])
            H = np.zeros((3, 9))
            H[:, 0:3] = np.eye(3)
            R = np.eye(3) * (meas_noise_std ** 2)
            S = H @ P[i] @ H.T + R
            K = P[i] @ H.T @ np.linalg.inv(S)
            dx = K @ y
            dq = small_angle_to_quat(dx[0:3])
            q[i] = quat_normalize(quat_multiply(dq, q[i]))
            w[i] += dx[3:6]
            b[i] += dx[6:9]
            P[i] = (I - K @ H) @ P[i] @ (I - K @ H).T + K @ R @ K.T

    return AttitudeFilterResult(q_wxyz=q, gyro_bias_rad_s=b, cov=P)


def run_ukf_aero(
    t_grid: np.ndarray,
    v_I: np.ndarray,
    env: Sequence[EnvInputs],
    aero: AeroAdapter,
    vehicle: VehicleParams,
    gyro_meas: np.ndarray,
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    q0_wxyz: np.ndarray,
    w0_rad_s: np.ndarray,
    bias0: np.ndarray,
    P0: np.ndarray,
    alpha: float = 1e-3,
    beta: float = 2.0,
    kappa: float = 0.0,
    gyro_noise_std: float = 1e-4,
    bias_rw_std: float = 1e-6,
    meas_noise_std: float = 1e-4,
) -> AttitudeFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    if v_I.shape != (t_grid.size, 3):
        raise ValueError("v_I must be (Nt, 3)")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")

    q = np.zeros((t_grid.size, 4))
    w = np.zeros((t_grid.size, 3))
    b = np.zeros((t_grid.size, 3))
    P = np.zeros((t_grid.size, 9, 9))

    q[0] = quat_normalize(np.array(q0_wxyz, dtype=float))
    w[0] = np.array(w0_rad_s, dtype=float)
    b[0] = np.array(bias0, dtype=float)
    P[0] = np.array(P0, dtype=float)

    meas_map = {int(idx): i for i, idx in enumerate(star_indices)}

    n = 9
    lambda_ = alpha * alpha * (n + kappa) - n
    scale = n + lambda_
    w0m = lambda_ / scale
    w0c = w0m + (1.0 - alpha * alpha + beta)
    wi = 1.0 / (2.0 * scale)

    for i in range(1, t_grid.size):
        dt = t_grid[i] - t_grid[i - 1]

        P_sym = 0.5 * (P[i - 1] + P[i - 1].T)
        L = np.linalg.cholesky(scale * (P_sym + 1e-12 * np.eye(n)))

        sigma = np.zeros((2 * n + 1, n))
        for k in range(n):
            sigma[1 + k] = L[:, k]
            sigma[1 + n + k] = -L[:, k]

        q_sigma = []
        w_sigma = []
        b_sigma = []
        for s in range(2 * n + 1):
            dtheta = sigma[s, 0:3]
            dw = sigma[s, 3:6]
            db = sigma[s, 6:9]
            dq = small_angle_to_quat(dtheta)
            q_s = quat_normalize(quat_multiply(dq, q[i - 1]))
            w_s = w[i - 1] + dw
            b_s = b[i - 1] + db

            q_p, w_p = propagate_attitude_aero_step(
                q_s,
                w_s,
                v_I[i - 1],
                env[i - 1],
                vehicle,
                aero,
                dt,
            )
            q_sigma.append(q_p)
            w_sigma.append(w_p)
            b_sigma.append(b_s)

        q_mean = q_sigma[0] * w0m
        w_mean = w0m * w_sigma[0]
        b_mean = w0m * b_sigma[0]
        for s in range(1, 2 * n + 1):
            w_mean += wi * w_sigma[s]
            b_mean += wi * b_sigma[s]
            q_mean += wi * q_sigma[s]
        q_mean = quat_normalize(q_mean)

        P_pred = np.zeros((n, n))
        for s in range(2 * n + 1):
            weight = w0c if s == 0 else wi
            dtheta = quat_error_small(q_sigma[s], q_mean)
            dw = w_sigma[s] - w_mean
            db = b_sigma[s] - b_mean
            dx = np.hstack((dtheta, dw, db))
            P_pred += weight * np.outer(dx, dx)

        Q = np.zeros((n, n))
        Q[3:6, 3:6] = np.eye(3) * (gyro_noise_std ** 2) * dt
        Q[6:9, 6:9] = np.eye(3) * (bias_rw_std ** 2) * dt
        P_pred += Q

        q[i] = q_mean
        w[i] = w_mean
        b[i] = b_mean
        P[i] = P_pred

        if i in meas_map:
            idx = meas_map[i]
            z_meas = star_meas[idx]
            z_mean = quat_error_small(z_meas, q_mean)

            Z = np.zeros((2 * n + 1, 3))
            for s in range(2 * n + 1):
                Z[s] = quat_error_small(q_sigma[s], q_mean)

            z_hat = w0m * Z[0]
            for s in range(1, 2 * n + 1):
                z_hat += wi * Z[s]

            S = np.zeros((3, 3))
            Pxz = np.zeros((n, 3))
            for s in range(2 * n + 1):
                weight = w0c if s == 0 else wi
                dz = Z[s] - z_hat
                dtheta = quat_error_small(q_sigma[s], q_mean)
                dw = w_sigma[s] - w_mean
                db = b_sigma[s] - b_mean
                dx = np.hstack((dtheta, dw, db))
                S += weight * np.outer(dz, dz)
                Pxz += weight * np.outer(dx, dz)
            S += np.eye(3) * (meas_noise_std ** 2)

            K = Pxz @ np.linalg.inv(S)
            dx = K @ (z_mean - z_hat)
            dq = small_angle_to_quat(dx[0:3])
            q[i] = quat_normalize(quat_multiply(dq, q[i]))
            w[i] += dx[3:6]
            b[i] += dx[6:9]
            P[i] = P_pred - K @ S @ K.T

    return AttitudeFilterResult(q_wxyz=q, gyro_bias_rad_s=b, cov=P)


def _kalman_update(
    x: np.ndarray,
    P: np.ndarray,
    z: np.ndarray,
    H: np.ndarray,
    R: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    S = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(S)
    x_new = x + K @ (z - H @ x)
    I = np.eye(P.shape[0])
    P_new = (I - K @ H) @ P @ (I - K @ H).T + K @ R @ K.T
    return x_new, P_new


def _sigma_points(x: np.ndarray, P: np.ndarray, alpha: float, beta: float, kappa: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = x.size
    lambda_ = alpha * alpha * (n + kappa) - n
    scale = n + lambda_
    P_sym = 0.5 * (P + P.T)
    jitter = 1e-12 * np.eye(n)
    L = np.linalg.cholesky(scale * (P_sym + jitter))

    sigma = np.zeros((2 * n + 1, n))
    sigma[0] = x
    for i in range(n):
        sigma[1 + i] = x + L[:, i]
        sigma[1 + n + i] = x - L[:, i]

    w_m = np.full(2 * n + 1, 1.0 / (2.0 * scale))
    w_c = np.full(2 * n + 1, 1.0 / (2.0 * scale))
    w_m[0] = lambda_ / scale
    w_c[0] = w_m[0] + (1.0 - alpha * alpha + beta)
    return sigma, w_m, w_c


def _ukf_measurement_update(
    x: np.ndarray,
    P: np.ndarray,
    z: np.ndarray,
    h_fn: Callable[[np.ndarray], np.ndarray],
    R: np.ndarray,
    alpha: float,
    beta: float,
    kappa: float,
) -> tuple[np.ndarray, np.ndarray]:
    sigma, w_m, w_c = _sigma_points(x, P, alpha, beta, kappa)
    y0 = h_fn(x)
    y_sigma = np.zeros((sigma.shape[0], y0.size))
    for i in range(sigma.shape[0]):
        y_sigma[i] = h_fn(sigma[i])
        if y_sigma[i].size == 4 and np.dot(y_sigma[i], y0) < 0.0:
            y_sigma[i] = -y_sigma[i]

    y_mean = np.sum(y_sigma * w_m[:, None], axis=0)
    if y_mean.size == 4:
        y_mean = y_mean / np.linalg.norm(y_mean)

    Pyy = np.zeros((y_mean.size, y_mean.size))
    Pxy = np.zeros((x.size, y_mean.size))
    for i in range(sigma.shape[0]):
        dy = y_sigma[i] - y_mean
        dx = sigma[i] - x
        Pyy += w_c[i] * np.outer(dy, dy)
        Pxy += w_c[i] * np.outer(dx, dy)
    Pyy += R

    K = Pxy @ np.linalg.inv(Pyy)
    x_new = x + K @ (z - y_mean)
    P_new = P - K @ Pyy @ K.T
    return x_new, P_new


def _sample_state_ensemble(
    x0: np.ndarray,
    P0: np.ndarray,
    members: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if members < 2:
        raise ValueError("members must be >= 2 for EnKF")
    P_sym = 0.5 * (P0 + P0.T)
    jitter = 0.0
    for _ in range(6):
        try:
            L = np.linalg.cholesky(P_sym + jitter * np.eye(P_sym.shape[0]))
            break
        except np.linalg.LinAlgError:
            diag_max = float(np.max(np.diag(P_sym))) if P_sym.size else 1.0
            jitter = (1e-12 if jitter == 0.0 else jitter * 10.0) * max(1.0, diag_max)
    else:
        evals, evecs = np.linalg.eigh(P_sym)
        eps = 1e-12 * max(1.0, float(np.max(evals)))
        evals = np.maximum(evals, eps)
        L = evecs @ np.diag(np.sqrt(evals))
    X = x0 + rng.standard_normal((members, x0.size)) @ L.T
    q_ref = x0[6:10].copy()
    if np.linalg.norm(q_ref) == 0.0:
        q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    for i in range(members):
        q = X[i, 6:10]
        n = np.linalg.norm(q)
        if n == 0.0:
            q = q_ref.copy()
        else:
            q = q / n
        if np.dot(q, q_ref) < 0.0:
            q = -q
        X[i, 6:10] = q
    return X


def _ensemble_mean_cov(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(X, axis=0)
    q_ref = X[0, 6:10]
    q = X[:, 6:10].copy()
    for i in range(q.shape[0]):
        if np.dot(q[i], q_ref) < 0.0:
            q[i] = -q[i]
    q_mean = np.mean(q, axis=0)
    n = np.linalg.norm(q_mean)
    if n == 0.0:
        q_mean = np.array([1.0, 0.0, 0.0, 0.0])
    else:
        q_mean = q_mean / n
    mean[6:10] = q_mean
    cov = np.cov(X, rowvar=False) if X.shape[0] > 1 else np.zeros((X.shape[1], X.shape[1]), dtype=float)
    cov = 0.5 * (cov + cov.T)
    return mean, cov


def _enkf_scalar_update(
    X: np.ndarray,
    y_pred: np.ndarray,
    z: float,
    sigma: float,
    rng: np.random.Generator,
    inflation: float,
) -> np.ndarray:
    if sigma <= 0.0:
        raise ValueError("measurement sigma must be positive")
    if inflation > 0.0 and inflation != 1.0:
        mean_x = np.mean(X, axis=0)
        X = mean_x + inflation * (X - mean_x)
    xa = X - np.mean(X, axis=0)
    ya = y_pred - float(np.mean(y_pred))
    pyy = float((ya @ ya) / max(X.shape[0] - 1, 1) + sigma * sigma)
    if pyy <= 0.0 or not np.isfinite(pyy):
        return X
    pxy = (xa.T @ ya) / max(X.shape[0] - 1, 1)
    K = pxy / pyy
    noise = rng.normal(0.0, sigma, size=X.shape[0])
    innov = float(z) + noise - y_pred
    X = X + np.outer(innov, K)
    for i in range(X.shape[0]):
        X[i, 6:10] = quat_normalize(X[i, 6:10])
    return X


def run_mekf_fullstate(
    prop_stm: StmPropagator,
    x0: np.ndarray,
    P0: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence[EnvInputs],
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    star_sigma: float = 1e-4,
    gyro_meas: Optional[np.ndarray] = None,
    gyro_sigma: float = 1e-3,
    magnetometer_meas: Optional[np.ndarray] = None,
    magnetometer_sigma: float = 1e-7,
) -> FullStateFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if t_grid.size < 2:
        raise ValueError("t_grid must have at least 2 entries")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")
    if P0.shape != (prop_stm.state_size, prop_stm.state_size):
        raise ValueError("P0 must match state size")
    if gyro_meas is not None and gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    if magnetometer_meas is not None and magnetometer_meas.shape != (t_grid.size, 3):
        raise ValueError("magnetometer_meas must be (Nt, 3)")

    x = np.zeros((t_grid.size, prop_stm.state_size))
    P = np.zeros((t_grid.size, prop_stm.state_size, prop_stm.state_size))
    x[0] = np.array(x0, dtype=float)
    P[0] = np.array(P0, dtype=float)

    star_map = {int(idx): i for i, idx in enumerate(star_indices)}

    for i in range(1, t_grid.size):
        t_pair = np.array([t_grid[i - 1], t_grid[i]], dtype=float)
        env_pair = [env[i - 1], env[i]]
        mean, cov = prop_stm.propagate(x[i - 1], P[i - 1], t_pair, env_pair)
        x_pred = mean[-1].copy()
        P_pred = cov[-1].copy()

        if gyro_meas is not None:
            H = np.zeros((3, prop_stm.state_size))
            H[:, 10:13] = np.eye(3)
            R = np.eye(3) * (gyro_sigma ** 2)
            z = gyro_meas[i]
            x_pred, P_pred = _kalman_update(x_pred, P_pred, z, H, R)

        if magnetometer_meas is not None:
            z = np.array(magnetometer_meas[i], dtype=float)
            z_pred = _predict_magnetometer_body(x_pred[6:10], env[i])
            H = np.zeros((3, prop_stm.state_size))
            eps = 1e-7
            for j in range(4):
                x_pert = x_pred.copy()
                x_pert[6 + j] += eps
                x_pert[6:10] = quat_normalize(x_pert[6:10])
                z_pert = _predict_magnetometer_body(x_pert[6:10], env[i])
                H[:, 6 + j] = (z_pert - z_pred) / eps
            R = np.eye(3) * (magnetometer_sigma ** 2)
            x_pred, P_pred = _kalman_update(x_pred, P_pred, z, H, R)
            x_pred[6:10] = quat_normalize(x_pred[6:10])

        if i in star_map:
            z = np.array(star_meas[star_map[i]], dtype=float)
            q_pred = x_pred[6:10]
            if np.dot(z, q_pred) < 0.0:
                z = -z
            H = np.zeros((4, prop_stm.state_size))
            H[:, 6:10] = np.eye(4)
            R = np.eye(4) * (star_sigma ** 2)
            x_pred, P_pred = _kalman_update(x_pred, P_pred, z, H, R)
            x_pred[6:10] = quat_normalize(x_pred[6:10])

        x[i] = x_pred
        P[i] = 0.5 * (P_pred + P_pred.T)

    return FullStateFilterResult(states=x, cov=P)


def run_ukf_fullstate(
    prop_ut: SigmaPointPropagatorUT,
    x0: np.ndarray,
    P0: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence[EnvInputs],
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    star_sigma: float = 1e-4,
    gyro_meas: Optional[np.ndarray] = None,
    gyro_sigma: float = 1e-3,
    magnetometer_meas: Optional[np.ndarray] = None,
    magnetometer_sigma: float = 1e-7,
    alpha: float = 1e-3,
    beta: float = 2.0,
    kappa: float = 0.0,
) -> FullStateFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if t_grid.size < 2:
        raise ValueError("t_grid must have at least 2 entries")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")
    if gyro_meas is not None and gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    if magnetometer_meas is not None and magnetometer_meas.shape != (t_grid.size, 3):
        raise ValueError("magnetometer_meas must be (Nt, 3)")

    x = np.zeros((t_grid.size, prop_ut.state_size))
    P = np.zeros((t_grid.size, prop_ut.state_size, prop_ut.state_size))
    x[0] = np.array(x0, dtype=float)
    P[0] = np.array(P0, dtype=float)

    star_map = {int(idx): i for i, idx in enumerate(star_indices)}

    for i in range(1, t_grid.size):
        t_pair = np.array([t_grid[i - 1], t_grid[i]], dtype=float)
        env_pair = [env[i - 1], env[i]]
        mean, cov = prop_ut.propagate(x[i - 1], P[i - 1], t_pair, env_pair)
        x_pred = mean[-1].copy()
        P_pred = cov[-1].copy()

        if gyro_meas is not None:
            z = gyro_meas[i]
            R = np.eye(3) * (gyro_sigma ** 2)
            def h_gyro(xi: np.ndarray) -> np.ndarray:
                return xi[10:13]
            x_pred, P_pred = _ukf_measurement_update(
                x_pred, P_pred, z, h_gyro, R, alpha, beta, kappa
            )

        if magnetometer_meas is not None:
            z = np.array(magnetometer_meas[i], dtype=float)
            R = np.eye(3) * (magnetometer_sigma ** 2)

            def h_mag(xi: np.ndarray) -> np.ndarray:
                return _predict_magnetometer_body(xi[6:10], env[i])

            x_pred, P_pred = _ukf_measurement_update(
                x_pred, P_pred, z, h_mag, R, alpha, beta, kappa
            )
            x_pred[6:10] = quat_normalize(x_pred[6:10])

        if i in star_map:
            z = np.array(star_meas[star_map[i]], dtype=float)
            q_pred = x_pred[6:10]
            if np.dot(z, q_pred) < 0.0:
                z = -z
            R = np.eye(4) * (star_sigma ** 2)
            def h_star(xi: np.ndarray) -> np.ndarray:
                q = xi[6:10]
                return q / np.linalg.norm(q)
            x_pred, P_pred = _ukf_measurement_update(
                x_pred, P_pred, z, h_star, R, alpha, beta, kappa
            )
            x_pred[6:10] = quat_normalize(x_pred[6:10])

        x[i] = x_pred
        P[i] = 0.5 * (P_pred + P_pred.T)

    return FullStateFilterResult(states=x, cov=P)


def run_enkf_fullstate(
    prop_mc: EnsemblePropagatorMC,
    x0: np.ndarray,
    P0: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence[EnvInputs],
    star_indices: np.ndarray,
    star_meas: np.ndarray,
    star_sigma: float = 1e-4,
    gyro_meas: Optional[np.ndarray] = None,
    gyro_sigma: float = 1e-3,
    magnetometer_meas: Optional[np.ndarray] = None,
    magnetometer_sigma: float = 1e-7,
    members: int = 64,
    seed: int = 0,
    inflation: float = 1.0,
) -> FullStateFilterResult:
    t_grid = np.array(t_grid, dtype=float)
    if t_grid.size < 2:
        raise ValueError("t_grid must have at least 2 entries")
    if len(env) != t_grid.size:
        raise ValueError("env must align with t_grid")
    if gyro_meas is not None and gyro_meas.shape != (t_grid.size, 3):
        raise ValueError("gyro_meas must be (Nt, 3)")
    if magnetometer_meas is not None and magnetometer_meas.shape != (t_grid.size, 3):
        raise ValueError("magnetometer_meas must be (Nt, 3)")

    rng = np.random.default_rng(seed)
    X = _sample_state_ensemble(np.array(x0, dtype=float), np.array(P0, dtype=float), members, rng)

    x_mean = np.zeros((t_grid.size, prop_mc.state_size))
    x_cov = np.zeros((t_grid.size, prop_mc.state_size, prop_mc.state_size))
    mean0, cov0 = _ensemble_mean_cov(X)
    x_mean[0] = mean0
    x_cov[0] = cov0

    star_map = {int(idx): i for i, idx in enumerate(star_indices)}

    for i in range(1, t_grid.size):
        t_pair = np.array([t_grid[i - 1], t_grid[i]], dtype=float)
        env_pair = [env[i - 1], env[i]]
        X = prop_mc.propagate(X, t_pair, env_pair)[-1]
        for j in range(X.shape[0]):
            X[j, 6:10] = quat_normalize(X[j, 6:10])

        if gyro_meas is not None:
            z = np.array(gyro_meas[i], dtype=float)
            for axis in range(3):
                y_pred = X[:, 10 + axis]
                X = _enkf_scalar_update(
                    X,
                    y_pred=y_pred,
                    z=float(z[axis]),
                    sigma=gyro_sigma,
                    rng=rng,
                    inflation=inflation,
                )

        if magnetometer_meas is not None:
            z = np.array(magnetometer_meas[i], dtype=float)
            y_pred_all = np.zeros((X.shape[0], 3), dtype=float)
            for j in range(X.shape[0]):
                y_pred_all[j] = _predict_magnetometer_body(X[j, 6:10], env[i])
            for axis in range(3):
                X = _enkf_scalar_update(
                    X,
                    y_pred=y_pred_all[:, axis],
                    z=float(z[axis]),
                    sigma=magnetometer_sigma,
                    rng=rng,
                    inflation=inflation,
                )
            for j in range(X.shape[0]):
                X[j, 6:10] = quat_normalize(X[j, 6:10])

        if i in star_map:
            z = np.array(star_meas[star_map[i]], dtype=float)
            z = z / np.linalg.norm(z)
            q_pred = X[:, 6:10].copy()
            for j in range(q_pred.shape[0]):
                if np.dot(q_pred[j], z) < 0.0:
                    q_pred[j] *= -1.0
            for comp in range(4):
                y_pred = q_pred[:, comp]
                X = _enkf_scalar_update(
                    X,
                    y_pred=y_pred,
                    z=float(z[comp]),
                    sigma=star_sigma,
                    rng=rng,
                    inflation=inflation,
                )

        mean_i, cov_i = _ensemble_mean_cov(X)
        x_mean[i] = mean_i
        x_cov[i] = cov_i

    return FullStateFilterResult(states=x_mean, cov=x_cov)
