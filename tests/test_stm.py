import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    EnvInputs,
    PropagatorConfig,
    SigmaPointPropagatorUT,
    StmPropagator,
    VehicleParams,
    default_geometry,
)


class TestStmPropagator(unittest.TestCase):
    def test_stm_shapes_zero_noise(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop = StmPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        P0 = np.zeros((prop.state_size, prop.state_size))

        t_grid = np.array([0.0, 10.0, 20.0])
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        mean, cov, stm = prop.propagate(x0, P0, t_grid, env, return_stm=True)

        self.assertEqual(mean.shape, (len(t_grid), prop.state_size))
        self.assertEqual(cov.shape, (len(t_grid), prop.state_size, prop.state_size))
        self.assertEqual(stm.shape, (len(t_grid), prop.state_size, prop.state_size))
        self.assertLess(np.max(np.abs(cov)), 1e-12)
        self.assertTrue(np.allclose(stm[0], np.eye(prop.state_size)))

    def test_stm_matches_ut_small_uncertainty(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop_stm = StmPropagator(aero, vehicle, config)
        prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop_stm.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        diag = np.array([
            1.0, 1.0, 1.0,
            1e-4, 1e-4, 1e-4,
            1e-8, 1e-8, 1e-8, 1e-8,
            1e-6, 1e-6, 1e-6,
            1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
        ])
        P0 = np.diag(diag)

        t_grid = np.array([0.0, 30.0])
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        mean_stm, cov_stm, stm = prop_stm.propagate(x0, P0, t_grid, env, return_stm=True)
        mean_ut, cov_ut = prop_ut.propagate(x0, P0, t_grid, env)

        pos_err = np.linalg.norm(mean_stm[-1, 0:3] - mean_ut[-1, 0:3])
        pos_rel = pos_err / np.linalg.norm(mean_ut[-1, 0:3])
        self.assertLess(pos_rel, 1e-8)

        cov_block_stm = cov_stm[-1, 0:6, 0:6]
        cov_block_ut = cov_ut[-1, 0:6, 0:6]
        rel_cov = np.linalg.norm(cov_block_stm - cov_block_ut) / np.linalg.norm(cov_block_ut)
        self.assertLess(rel_cov, 0.05)

        P_from_phi = stm[-1] @ P0 @ stm[-1].T
        rel_phi = np.linalg.norm(P_from_phi[0:6, 0:6] - cov_block_stm) / np.linalg.norm(cov_block_stm)
        self.assertLess(rel_phi, 1e-6)


if __name__ == "__main__":
    unittest.main()
