import unittest

import numpy as np

from vleo_uq import (
    AttitudeErrorSeries,
    boresight_eci_from_states,
    boresight_eci_with_errors,
    geolocation_error_m,
    ground_intersect_spherical,
    ground_points_from_los,
    range_bias_m,
    simulate_attitude_ou,
    smear_over_exposure_m,
)


class TestAttitudeUq(unittest.TestCase):
    def test_simulate_attitude_ou_shape(self):
        t_grid = np.linspace(0.0, 10.0, 6)
        series = simulate_attitude_ou(t_grid, tau_s=20.0, sigma_rad=1e-3)
        self.assertIsInstance(series, AttitudeErrorSeries)
        self.assertEqual(series.errors_rad.shape, (len(t_grid), 3))

    def test_ground_intersect_nadir(self):
        R = 6378137.0
        r_sc = np.array([R + 500e3, 0.0, 0.0])
        los = np.array([-1.0, 0.0, 0.0])
        ground = ground_intersect_spherical(r_sc, los, earth_radius_m=R)
        self.assertTrue(np.allclose(ground, np.array([R, 0.0, 0.0]), atol=1e-6))

    def test_geolocation_error_zero(self):
        ground = np.array([[6378137.0, 0.0, 0.0], [0.0, 6378137.0, 0.0]])
        err = geolocation_error_m(ground, ground)
        self.assertTrue(np.allclose(err, 0.0))

    def test_smear_positive(self):
        t_grid = np.array([0.0, 1.0, 2.0])
        points = np.array(
            [
                [6378137.0, 0.0, 0.0],
                [6378137.0, 10.0, 0.0],
                [6378137.0, 20.0, 0.0],
            ]
        )
        smear = smear_over_exposure_m(t_grid, points, exposure_s=2.0)
        self.assertGreater(smear[1], 0.0)

    def test_boresight_and_range_bias(self):
        states = np.zeros((2, 18))
        states[:, 6] = 1.0
        states[:, 0] = 6378137.0 + 500e3
        states[:, 3] = 7500.0
        boresight = boresight_eci_from_states(states)
        self.assertTrue(np.allclose(boresight[0], np.array([0.0, 0.0, 1.0])))

        errors = np.zeros((2, 3))
        errors[1, 0] = 1e-3
        boresight_err = boresight_eci_with_errors(states, errors)
        self.assertFalse(np.allclose(boresight_err[1], boresight[1]))

        r_sc = states[:, 0:3]
        los = boresight_err
        bias = range_bias_m(r_sc, los)
        self.assertEqual(bias.shape[0], 2)


if __name__ == "__main__":
    unittest.main()
