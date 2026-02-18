import math
import os
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
)


def _build_env(t_grid: np.ndarray) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = 0.0
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def _make_x0(state_size: int, mu_earth: float) -> np.ndarray:
    r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
    v0 = np.array([0.0, math.sqrt(mu_earth / np.linalg.norm(r0)), 0.0], dtype=float)
    x0 = np.zeros(state_size, dtype=float)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    return x0


class TestMethodParity(unittest.TestCase):
    def _build_stack(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        cfg = PropagatorConfig()
        cfg.rho_fast_sigma = 0.0
        cfg.rho_bias_sigma = 0.0
        cfg.wind_sigma = 0.0
        cfg.rtol = 1e-7
        cfg.atol = 1e-9
        cfg.rng_seed = 1234
        return aero, vehicle, cfg

    def test_mc_shape_and_mean_matches_det_when_particles_identical(self):
        aero, vehicle, cfg = self._build_stack()
        prop_det = DeterministicPropagator(aero, vehicle, cfg)
        prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
        x0 = _make_x0(prop_det.state_size, cfg.mu_earth_m3_s2)

        t_grid = np.linspace(0.0, 600.0, 13)
        env = _build_env(t_grid)

        n = 24
        X0 = np.repeat(x0[None, :], n, axis=0)
        det = prop_det.propagate(x0, t_grid, env)
        mc = prop_mc.propagate(X0, t_grid, env)

        self.assertEqual(mc.shape, (len(t_grid), n, prop_det.state_size))
        self.assertLess(float(np.max(np.abs(det - mc.mean(axis=1)))), 1e-8)

    def test_mc_reproducible_across_thread_counts(self):
        aero, vehicle, cfg = self._build_stack()
        prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
        x0 = _make_x0(prop_mc.state_size, cfg.mu_earth_m3_s2)
        t_grid = np.linspace(0.0, 300.0, 7)
        env = _build_env(t_grid)

        rng = np.random.default_rng(2026)
        n = 32
        X0 = np.repeat(x0[None, :], n, axis=0)
        X0[:, 0:3] += rng.normal(0.0, 2.0, size=(n, 3))
        X0[:, 3:6] += rng.normal(0.0, 0.02, size=(n, 3))
        X0[:, 6] = 1.0

        prev = os.environ.get("OMP_NUM_THREADS")
        try:
            os.environ["OMP_NUM_THREADS"] = "1"
            mc_1 = prop_mc.propagate(X0, t_grid, env)
            os.environ["OMP_NUM_THREADS"] = "2"
            mc_2 = prop_mc.propagate(X0, t_grid, env)
        finally:
            if prev is None:
                os.environ.pop("OMP_NUM_THREADS", None)
            else:
                os.environ["OMP_NUM_THREADS"] = prev

        self.assertTrue(np.array_equal(mc_1, mc_2))

    def test_ut_stm_covariance_agree_with_mc_short_arc(self):
        aero, vehicle, cfg = self._build_stack()
        prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
        prop_ut = SigmaPointPropagatorUT(aero, vehicle, cfg)
        prop_stm = StmPropagator(aero, vehicle, cfg)
        x0 = _make_x0(prop_mc.state_size, cfg.mu_earth_m3_s2)

        t_grid = np.linspace(0.0, 300.0, 11)
        env = _build_env(t_grid)

        P0 = np.diag(
            [
                10.0,
                10.0,
                10.0,
                1e-2,
                1e-2,
                1e-2,
                1e-8,
                1e-8,
                1e-8,
                1e-8,
                1e-6,
                1e-6,
                1e-6,
                1e-4,
                1e-4,
                1e-4,
                1e-4,
                1e-4,
            ]
        )
        rng = np.random.default_rng(99)
        n = 160
        X0 = rng.multivariate_normal(mean=x0, cov=P0, size=n)
        X0[:, 6:10] = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)

        mc = prop_mc.propagate(X0, t_grid, env)
        mc_mean = mc.mean(axis=1)
        mc_cov = np.cov(mc[-1], rowvar=False)

        ut_mean, ut_cov = prop_ut.propagate(x0, P0, t_grid, env)
        stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)

        ut_mean_err = float(np.linalg.norm(ut_mean[-1, 0:6] - mc_mean[-1, 0:6]))
        stm_mean_err = float(np.linalg.norm(stm_mean[-1, 0:6] - mc_mean[-1, 0:6]))
        ut_cov_err = float(np.linalg.norm(ut_cov[-1, 0:6, 0:6] - mc_cov[0:6, 0:6]) / np.linalg.norm(mc_cov[0:6, 0:6]))
        stm_cov_err = float(np.linalg.norm(stm_cov[-1, 0:6, 0:6] - mc_cov[0:6, 0:6]) / np.linalg.norm(mc_cov[0:6, 0:6]))

        self.assertLess(ut_mean_err, 25.0)
        self.assertLess(stm_mean_err, 25.0)
        self.assertLess(ut_cov_err, 0.35)
        self.assertLess(stm_cov_err, 0.35)


if __name__ == "__main__":
    unittest.main()
