#!/usr/bin/env python3
"""Attitude dynamics + MEKF example."""

from __future__ import annotations

import argparse
import numpy as np

from vleo_uq import (
    run_mekf,
    simulate_attitude_truth,
    simulate_gyro_measurements,
    simulate_star_tracker_measurements,
)


def quat_angle_error(q_est: np.ndarray, q_true: np.ndarray) -> float:
    dot = np.abs(np.dot(q_est, q_true))
    dot = np.clip(dot, -1.0, 1.0)
    return 2.0 * np.arccos(dot)


def main() -> None:
    parser = argparse.ArgumentParser(description="Attitude dynamics + MEKF example.")
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=1.0)
    parser.add_argument("--gyro_sigma", type=float, default=1e-4)
    parser.add_argument("--star_sigma", type=float, default=5e-4)
    parser.add_argument("--star_cadence", type=float, default=5.0)
    args = parser.parse_args()

    t_grid = np.arange(0.0, args.duration_s + args.dt_s, args.dt_s)
    inertia = np.diag([1.0, 1.2, 0.8])
    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    w0 = np.array([0.0, 0.0, 0.01])

    truth = simulate_attitude_truth(
        t_grid,
        q0,
        w0,
        inertia,
        bias_tau_s=200.0,
        bias_sigma_rad_s=1e-5,
    )
    gyro = simulate_gyro_measurements(truth, gyro_sigma_rad_s=args.gyro_sigma)
    star_idx, star_q = simulate_star_tracker_measurements(
        truth, sigma_rad=args.star_sigma, cadence_s=args.star_cadence
    )

    P0 = np.diag([1e-6] * 6)
    result = run_mekf(
        t_grid,
        gyro,
        star_idx,
        star_q,
        q0,
        np.zeros(3),
        P0,
        gyro_noise_std=args.gyro_sigma,
        bias_rw_std=1e-6,
        meas_noise_std=args.star_sigma,
    )

    err = quat_angle_error(result.q_wxyz[-1], truth.q_wxyz[-1])
    print(f"Final attitude error (rad): {err:.6e}")


if __name__ == "__main__":
    main()
