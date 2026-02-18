import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnsemblePropagatorMC,
    EnvInputs,
    PropagatorConfig,
    SigmaPointPropagatorUT,
    StmPropagator,
    VehicleParams,
    default_geometry,
    run_enkf_fullstate,
    run_mekf_fullstate,
    run_ukf_fullstate,
)


def quat_angle_error(q_est: np.ndarray, q_true: np.ndarray) -> float:
    dot = np.abs(np.dot(q_est, q_true))
    dot = np.clip(dot, -1.0, 1.0)
    return 2.0 * np.arccos(dot)


class TestFullStateFilter(unittest.TestCase):
    def test_fullstate_filters_track_truth(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()
        config = PropagatorConfig()

        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_mc = EnsemblePropagatorMC(aero, vehicle, config)
        prop_stm = StmPropagator(aero, vehicle, config)
        prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop_det.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 60.0, 61)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth = prop_det.propagate(x0, t_grid, env)
        star_indices = np.arange(0, t_grid.size, 5, dtype=int)
        star_meas = truth[star_indices, 6:10]
        gyro_meas = truth[:, 10:13]

        P0 = np.eye(prop_det.state_size) * 1e-12
        mekf = run_mekf_fullstate(
            prop_stm,
            x0,
            P0,
            t_grid,
            env,
            star_indices,
            star_meas,
            star_sigma=1e-6,
            gyro_meas=gyro_meas,
            gyro_sigma=1e-6,
        )
        ukf = run_ukf_fullstate(
            prop_ut,
            x0,
            P0,
            t_grid,
            env,
            star_indices,
            star_meas,
            star_sigma=1e-6,
            gyro_meas=gyro_meas,
            gyro_sigma=1e-6,
        )
        enkf = run_enkf_fullstate(
            prop_mc,
            x0,
            P0,
            t_grid,
            env,
            star_indices,
            star_meas,
            star_sigma=1e-6,
            gyro_meas=gyro_meas,
            gyro_sigma=1e-6,
            members=32,
            seed=13,
            inflation=1.0,
        )

        err_mekf = quat_angle_error(mekf.states[-1, 6:10], truth[-1, 6:10])
        err_ukf = quat_angle_error(ukf.states[-1, 6:10], truth[-1, 6:10])
        err_enkf = quat_angle_error(enkf.states[-1, 6:10], truth[-1, 6:10])
        self.assertLess(err_mekf, 1e-4)
        self.assertLess(err_ukf, 1e-4)
        self.assertLess(err_enkf, 1e-4)


if __name__ == "__main__":
    unittest.main()
