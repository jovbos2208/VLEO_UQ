import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    EnvInputs,
    VehicleParams,
    default_geometry,
    run_mekf_aero,
    run_ukf_aero,
    simulate_attitude_truth_aero,
    simulate_gyro_measurements,
    simulate_star_tracker_measurements,
)


def quat_angle_error(q_est: np.ndarray, q_true: np.ndarray) -> float:
    dot = np.abs(np.dot(q_est, q_true))
    dot = np.clip(dot, -1.0, 1.0)
    return 2.0 * np.arccos(dot)


class TestAttitudeFilterAero(unittest.TestCase):
    def test_filters_run_with_aero(self):
        t_grid = np.linspace(0.0, 60.0, 61)
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 1e-12
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        v_I = np.zeros((t_grid.size, 3))
        v_I[:, 1] = 7500.0

        q0 = np.array([1.0, 0.0, 0.0, 0.0])
        w0 = np.array([0.0, 0.0, 0.01])

        truth = simulate_attitude_truth_aero(
            t_grid,
            q0,
            w0,
            v_I,
            env,
            vehicle,
            aero,
            bias_tau_s=200.0,
            bias_sigma_rad_s=1e-5,
        )
        gyro = simulate_gyro_measurements(truth, gyro_sigma_rad_s=1e-4)
        star_idx, star_q = simulate_star_tracker_measurements(truth, sigma_rad=5e-4, cadence_s=5.0)

        P0 = np.diag([1e-6] * 9)
        mekf = run_mekf_aero(
            t_grid,
            v_I,
            env,
            aero,
            vehicle,
            gyro,
            star_idx,
            star_q,
            q0,
            w0,
            np.zeros(3),
            P0,
            gyro_noise_std=1e-4,
            bias_rw_std=1e-6,
            meas_noise_std=5e-4,
        )
        ukf = run_ukf_aero(
            t_grid,
            v_I,
            env,
            aero,
            vehicle,
            gyro,
            star_idx,
            star_q,
            q0,
            w0,
            np.zeros(3),
            P0,
            gyro_noise_std=1e-4,
            bias_rw_std=1e-6,
            meas_noise_std=5e-4,
        )

        err_mekf = quat_angle_error(mekf.q_wxyz[-1], truth.q_wxyz[-1])
        err_ukf = quat_angle_error(ukf.q_wxyz[-1], truth.q_wxyz[-1])
        self.assertTrue(np.isfinite(err_mekf))
        self.assertTrue(np.isfinite(err_ukf))


if __name__ == "__main__":
    unittest.main()
