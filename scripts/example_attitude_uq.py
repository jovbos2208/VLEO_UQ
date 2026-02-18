#!/usr/bin/env python3
"""Example attitude UQ: OU errors mapped to ground geolocation and smear."""

from __future__ import annotations

import argparse
import numpy as np

from vleo_uq import (
    dcm_to_quat_wxyz,
    boresight_eci_from_states,
    boresight_eci_with_errors,
    ground_points_from_los,
    geolocation_error_m,
    simulate_attitude_ou,
    smear_over_exposure_m,
)


def build_states_circular(r0_m: float, v0_m_s: float, t_grid: np.ndarray) -> np.ndarray:
    states = np.zeros((t_grid.size, 18))
    for i, t in enumerate(t_grid):
        theta = v0_m_s * t / r0_m
        r = np.array([r0_m * np.cos(theta), r0_m * np.sin(theta), 0.0])
        v = np.array([-v0_m_s * np.sin(theta), v0_m_s * np.cos(theta), 0.0])

        r_hat = r / np.linalg.norm(r)
        h_hat = np.cross(r, v)
        h_hat /= np.linalg.norm(h_hat)
        z_b = -r_hat
        y_b = -h_hat
        x_b = np.cross(y_b, z_b)
        R_BI = np.vstack((x_b, y_b, z_b))
        q = dcm_to_quat_wxyz(R_BI)

        states[i, 0:3] = r
        states[i, 3:6] = v
        states[i, 6:10] = q
    return states


def main() -> None:
    parser = argparse.ArgumentParser(description="Attitude UQ example.")
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=1.0)
    parser.add_argument("--tau_s", type=float, default=200.0)
    parser.add_argument("--sigma_rad", type=float, default=1e-4)
    parser.add_argument("--exposure_s", type=float, default=1.0)
    args = parser.parse_args()

    R = 6378137.0
    alt = 500e3
    r0 = R + alt
    v0 = np.sqrt(3.986004418e14 / r0)

    t_grid = np.arange(0.0, args.duration_s + args.dt_s, args.dt_s)
    states = build_states_circular(r0, v0, t_grid)

    errors = simulate_attitude_ou(t_grid, tau_s=args.tau_s, sigma_rad=args.sigma_rad)
    boresight_nom = boresight_eci_from_states(states)
    boresight_err = boresight_eci_with_errors(states, errors.errors_rad)

    ground_nom = ground_points_from_los(states[:, 0:3], boresight_nom)
    ground_err = ground_points_from_los(states[:, 0:3], boresight_err)
    geo_err = geolocation_error_m(ground_nom, ground_err)
    smear = smear_over_exposure_m(t_grid, ground_err, exposure_s=args.exposure_s)

    print(f"Geolocation RMS (m): {np.sqrt(np.mean(geo_err ** 2)):.3f}")
    print(f"Smear median (m): {np.median(smear):.3f}")


if __name__ == "__main__":
    main()
