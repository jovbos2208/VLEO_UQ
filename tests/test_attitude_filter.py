import unittest

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


class TestAttitudeFilter(unittest.TestCase):
    def test_mekf_reduces_error(self):
        t_grid = np.linspace(0.0, 200.0, 201)
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
        gyro = simulate_gyro_measurements(truth, gyro_sigma_rad_s=1e-4)
        star_idx, star_q = simulate_star_tracker_measurements(truth, sigma_rad=5e-4, cadence_s=5.0)

        P0 = np.diag([1e-6] * 3 + [1e-6] * 3)
        result = run_mekf(
            t_grid,
            gyro,
            star_idx,
            star_q,
            q0,
            np.zeros(3),
            P0,
            gyro_noise_std=1e-4,
            bias_rw_std=1e-6,
            meas_noise_std=5e-4,
        )

        no_update = run_mekf(
            t_grid,
            gyro,
            np.array([], dtype=int),
            np.zeros((0, 4)),
            q0,
            np.zeros(3),
            P0,
            gyro_noise_std=1e-4,
            bias_rw_std=1e-6,
            meas_noise_std=5e-4,
        )

        err_filter = quat_angle_error(result.q_wxyz[-1], truth.q_wxyz[-1])
        err_open = quat_angle_error(no_update.q_wxyz[-1], truth.q_wxyz[-1])
        self.assertLess(err_filter, err_open)


if __name__ == "__main__":
    unittest.main()
