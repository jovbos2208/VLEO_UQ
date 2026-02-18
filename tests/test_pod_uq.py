import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagatorConfig,
    StmPropagator,
    VehicleParams,
    default_geometry,
    run_pod_uq,
)


class TestPodUqHarness(unittest.TestCase):
    def test_pod_od_reduces_position_error(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_stm = StmPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0_truth = np.zeros(prop_stm.state_size)
        x0_truth[0:3] = r0
        x0_truth[3:6] = v0
        x0_truth[6] = 1.0

        x0_guess = x0_truth.copy()
        x0_guess[0:3] += np.array([20.0, -10.0, 5.0])

        diag = np.array([
            100.0, 100.0, 100.0,
            1e-2, 1e-2, 1e-2,
            1e-6, 1e-6, 1e-6, 1e-6,
            1e-6, 1e-6, 1e-6,
            1e-4, 1e-4, 1e-4, 1e-4, 1e-4,
        ])
        P0_guess = np.diag(diag)

        t_grid = np.linspace(0.0, 60.0, 7)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        result = run_pod_uq(
            prop_det,
            prop_stm,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            meas_sigma_m=0.5,
            meas_cadence_s=10.0,
            max_iter=5,
        )

        err0 = np.linalg.norm(x0_guess[0:3] - x0_truth[0:3])
        err1 = np.linalg.norm(result.x0_est[0:3] - x0_truth[0:3])
        self.assertLess(err1, 0.3 * err0)
        self.assertIn("1sigma", result.coverage_sigma)
        self.assertIn("2sigma", result.coverage_sigma)
        self.assertIn("3sigma", result.coverage_sigma)


if __name__ == "__main__":
    unittest.main()
