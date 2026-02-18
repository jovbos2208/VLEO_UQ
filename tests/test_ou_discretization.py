import math
import os
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    EnsemblePropagatorMC,
    EnvInputs,
    PropagatorConfig,
    VehicleParams,
    default_geometry,
)

os.environ.setdefault("OMP_NUM_THREADS", "1")


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


def _build_mc_setup(n_particles: int, cfg: PropagatorConfig):
    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()
    prop = EnsemblePropagatorMC(aero, vehicle, cfg)

    mu = cfg.mu_earth_m3_s2
    r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
    if mu > 0.0:
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0], dtype=float)
    else:
        v0 = np.array([7500.0, 0.0, 0.0], dtype=float)
    X0 = np.zeros((n_particles, prop.state_size), dtype=float)
    X0[:, 0:3] = r0
    X0[:, 3:6] = v0
    X0[:, 6] = 1.0
    return prop, X0


class TestOuDiscretization(unittest.TestCase):
    def test_ou_latent_statistics_match_theory(self):
        cfg = PropagatorConfig()
        cfg.freeze_attitude = True
        cfg.rtol = 1e-4
        cfg.atol = 1e-6
        cfg.mu_earth_m3_s2 = 0.0
        cfg.rho_fast_tau_s = 600.0
        cfg.rho_fast_sigma = 0.2
        cfg.rho_bias_sigma = 0.0
        cfg.wind_sigma = 0.0
        cfg.rng_seed = 12345

        n_particles = 128
        prop, X0 = _build_mc_setup(n_particles, cfg)
        dt = 60.0
        t_grid = np.arange(0.0, 2400.0 + dt, dt)
        env = _build_env(t_grid)

        X = prop.propagate(X0, t_grid, env)
        x = X[:, :, 13]

        t_final = float(t_grid[-1] - t_grid[0])
        var_theory = cfg.rho_fast_sigma**2 * (1.0 - math.exp(-2.0 * t_final / cfg.rho_fast_tau_s))

        mean_emp = float(np.mean(x[-1]))
        var_emp = float(np.var(x[-1], ddof=1))
        self.assertLess(abs(mean_emp), 0.05)
        self.assertLess(abs(var_emp - var_theory) / max(var_theory, 1e-12), 0.35)

        burn_idx = int(np.searchsorted(t_grid, 2.0 * cfg.rho_fast_tau_s))
        burn_idx = min(burn_idx, max(0, len(t_grid) - 3))
        x_now = x[burn_idx:-1].reshape(-1)
        x_next = x[burn_idx + 1 :].reshape(-1)
        corr = float(np.corrcoef(x_now, x_next)[0, 1])
        phi = math.exp(-dt / cfg.rho_fast_tau_s)
        self.assertLess(abs(corr - phi), 0.14)

    def test_random_walk_fallback_variance(self):
        cfg = PropagatorConfig()
        cfg.freeze_attitude = True
        cfg.rtol = 1e-4
        cfg.atol = 1e-6
        cfg.mu_earth_m3_s2 = 0.0
        cfg.rho_fast_sigma = 0.0
        cfg.wind_sigma = 0.0
        cfg.rho_bias_tau_s = 0.0
        cfg.rho_bias_sigma = 0.02
        cfg.rng_seed = 2026

        n_particles = 128
        prop, X0 = _build_mc_setup(n_particles, cfg)
        dt = 10.0
        t_grid = np.arange(0.0, 1200.0 + dt, dt)
        env = _build_env(t_grid)
        X = prop.propagate(X0, t_grid, env)

        x_bias = X[-1, :, 14]
        t_final = float(t_grid[-1] - t_grid[0])
        var_theory = cfg.rho_bias_sigma**2 * t_final
        var_emp = float(np.var(x_bias, ddof=1))
        self.assertLess(abs(var_emp - var_theory) / max(var_theory, 1e-12), 0.35)

    def test_mc_seed_reproducibility_for_latent_states(self):
        cfg = PropagatorConfig()
        cfg.freeze_attitude = True
        cfg.rtol = 1e-4
        cfg.atol = 1e-6
        cfg.mu_earth_m3_s2 = 0.0
        cfg.rho_fast_tau_s = 900.0
        cfg.rho_fast_sigma = 0.1
        cfg.rng_seed = 777

        n_particles = 8
        prop_a, X0_a = _build_mc_setup(n_particles, cfg)
        prop_b, X0_b = _build_mc_setup(n_particles, cfg)

        t_grid = np.arange(0.0, 120.0 + 30.0, 30.0)
        env = _build_env(t_grid)

        Xa = prop_a.propagate(X0_a, t_grid, env)
        Xb = prop_b.propagate(X0_b, t_grid, env)
        self.assertTrue(np.array_equal(Xa, Xb))


if __name__ == "__main__":
    unittest.main()
