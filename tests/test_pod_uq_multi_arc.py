import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    GroundStation,
    PropagatorConfig,
    StmPropagator,
    VehicleParams,
    default_geometry,
    run_pod_uq_multi_arc,
    simulate_slr_measurements,
    slice_slr_measurements,
)


class TestPodMultiArc(unittest.TestCase):
    def test_multi_arc_coupled_runs(self):
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
        x0_guess[0:3] += np.array([15.0, -5.0, 2.0])

        P0_guess = np.eye(prop_stm.state_size) * 1e-3
        t_grid = np.linspace(0.0, 120.0, 13)

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth_states = prop_det.propagate(x0_truth, t_grid, env)
        stations = [GroundStation("S1", 0.0, 0.0, 0.0)]
        availability = np.ones((len(t_grid), 1), dtype=bool)
        slr_full = simulate_slr_measurements(
            truth_states,
            t_grid,
            stations,
            sigma_m=0.05,
            rng=np.random.default_rng(7),
            cadence_s=10.0,
            availability_mask=availability,
        )

        arcs = [(0, 7), (6, 13)]
        slr_sets = [slice_slr_measurements(slr_full, s, e) for s, e in arcs]
        results = run_pod_uq_multi_arc(
            prop_det,
            prop_stm,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            arcs,
            slr_sets=slr_sets,
            max_iter=4,
            couple_arcs=True,
            carry_covariance=True,
        )

        self.assertEqual(len(results), 2)
        self.assertIn("1sigma", results[0].coverage_sigma)
        self.assertIn("2sigma", results[0].coverage_sigma)
        self.assertIn("3sigma", results[0].coverage_sigma)

    def test_multi_arc_enkf_runs(self):
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
        x0_guess[0:3] += np.array([20.0, -8.0, 3.0])
        P0_guess = np.eye(prop_stm.state_size) * 1e-3
        t_grid = np.linspace(0.0, 120.0, 13)

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth_states = prop_det.propagate(x0_truth, t_grid, env)
        stations = [GroundStation("S1", 0.0, 0.0, 0.0)]
        availability = np.ones((len(t_grid), 1), dtype=bool)
        slr_full = simulate_slr_measurements(
            truth_states,
            t_grid,
            stations,
            sigma_m=0.05,
            rng=np.random.default_rng(9),
            cadence_s=10.0,
            availability_mask=availability,
        )

        arcs = [(0, 7), (6, 13)]
        slr_sets = [slice_slr_measurements(slr_full, s, e) for s, e in arcs]
        results = run_pod_uq_multi_arc(
            prop_det,
            prop_stm,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            arcs,
            slr_sets=slr_sets,
            couple_arcs=True,
            carry_covariance=True,
            estimator="enkf",
            enkf_members=24,
            enkf_seed=17,
            enkf_inflation=1.0,
        )

        self.assertEqual(len(results), 2)
        self.assertTrue(np.all(np.isfinite(results[0].x0_est)))
        self.assertTrue(np.all(np.isfinite(results[1].x0_est)))


if __name__ == "__main__":
    unittest.main()
