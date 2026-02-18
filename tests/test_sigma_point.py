import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagatorConfig,
    SigmaPointPropagatorUT,
    VehicleParams,
    default_geometry,
)


class TestSigmaPointPropagatorUT(unittest.TestCase):
    def test_ut_mean_matches_deterministic(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop_ut.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0  # unit quaternion w

        period = 2.0 * math.pi * math.sqrt(np.linalg.norm(r0) ** 3 / mu)
        t_grid = np.linspace(0.0, 0.25 * period, 81)

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        X_det = prop_det.propagate(x0, t_grid, env)

        P0 = np.eye(prop_ut.state_size) * 1e-6
        mean, cov = prop_ut.propagate(x0, P0, t_grid, env)

        r_det = X_det[-1, 0:3]
        r_ut = mean[-1, 0:3]
        rel_pos_err = np.linalg.norm(r_ut - r_det) / np.linalg.norm(r_det)

        self.assertLess(rel_pos_err, 1e-6)
        self.assertEqual(mean.shape, (len(t_grid), prop_ut.state_size))
        self.assertEqual(cov.shape, (len(t_grid), prop_ut.state_size, prop_ut.state_size))

    def test_ut_process_noise_inflates_latent_cov(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_tau_s = 1200.0
        config.rho_fast_sigma = 0.2
        config.wind_tau_s = 800.0
        config.wind_sigma = 0.5

        prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0 = np.zeros(prop_ut.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0
        P0 = np.zeros((prop_ut.state_size, prop_ut.state_size))

        t_grid = np.linspace(0.0, 600.0, 7)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        mean, cov = prop_ut.propagate(x0, P0, t_grid, env, augment_process_noise=True)

        self.assertGreater(cov[-1, 13, 13], 0.0)
        self.assertGreater(cov[-1, 15, 15], 0.0)


if __name__ == "__main__":
    unittest.main()
