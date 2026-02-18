import math
import os
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagatorConfig,
    VehicleParams,
    default_geometry,
    default_ilrs_stations,
    simulate_gnss_measurements,
    simulate_slr_measurements,
)


def _build_env(t_grid: np.ndarray, density: float) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = float(density)
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def _make_sat_positions(t_grid: np.ndarray, n_sats: int = 8) -> np.ndarray:
    radius_m = 26560e3
    mu = 3.986004418e14
    n = math.sqrt(mu / radius_m**3)
    sat_positions = np.zeros((t_grid.size, n_sats, 3), dtype=float)
    phases = np.linspace(0.0, 2.0 * math.pi, n_sats, endpoint=False)
    for s in range(n_sats):
        for i, t in enumerate(t_grid):
            ang = phases[s] + n * float(t)
            sat_positions[i, s, 0] = radius_m * math.cos(ang)
            sat_positions[i, s, 1] = radius_m * math.sin(ang)
    return sat_positions


class TestPropagatorRegression(unittest.TestCase):
    def _make_propagator(self, *, rtol: float = 1e-9, atol: float = 1e-10):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()
        cfg = PropagatorConfig()
        cfg.rtol = rtol
        cfg.atol = atol
        return aero, DeterministicPropagator(aero, vehicle, cfg), cfg

    def test_two_body_matches_circular_reference(self):
        _aero, prop, cfg = self._make_propagator()
        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 1200.0, 121)
        env = _build_env(t_grid, density=0.0)
        states = prop.propagate(x0, t_grid, env)

        n = math.sqrt(cfg.mu_earth_m3_s2 / np.linalg.norm(r0) ** 3)
        ref = np.zeros((t_grid.size, 6), dtype=float)
        for i, t in enumerate(t_grid):
            th = n * float(t)
            c = math.cos(th)
            s = math.sin(th)
            ref[i, 0:3] = np.array([r0[0] * c, r0[0] * s, 0.0])
            ref[i, 3:6] = np.array([-np.linalg.norm(v0) * s, np.linalg.norm(v0) * c, 0.0])

        pos_err = np.linalg.norm(states[:, 0:3] - ref[:, 0:3], axis=1)
        vel_err = np.linalg.norm(states[:, 3:6] - ref[:, 3:6], axis=1)

        self.assertLess(float(np.max(pos_err)), 50.0)
        self.assertLess(float(np.max(vel_err)), 0.08)

    def test_aero_shadowing_toggle_changes_force_and_torque(self):
        aero, _prop, _cfg = self._make_propagator(rtol=1e-6, atol=1e-8)

        q_wxyz = np.array([0.965925826, 0.0, 0.258819045, 0.0], dtype=float)
        w_b = np.array([0.02, -0.03, 0.01], dtype=float)
        v_i = np.array([7600.0, 120.0, -80.0], dtype=float)
        wind_i = np.zeros(3, dtype=float)
        density = 2e-11
        temp = 950.0
        pmass = 28.0 * 1.6605390689252e-27
        eta1 = math.radians(22.0)
        eta2 = math.radians(-17.0)

        prev_flag = os.environ.pop("FM_DISABLE_SHADOWING", None)
        try:
            f_on, tau_on = aero.compute_ft(
                q_wxyz,
                w_b,
                v_i,
                wind_i,
                density,
                temp,
                pmass,
                eta1,
                eta2,
                1,
            )
            os.environ["FM_DISABLE_SHADOWING"] = "1"
            f_off, tau_off = aero.compute_ft(
                q_wxyz,
                w_b,
                v_i,
                wind_i,
                density,
                temp,
                pmass,
                eta1,
                eta2,
                1,
            )
        finally:
            if prev_flag is None:
                os.environ.pop("FM_DISABLE_SHADOWING", None)
            else:
                os.environ["FM_DISABLE_SHADOWING"] = prev_flag

        fn_on = float(np.linalg.norm(f_on))
        fn_off = float(np.linalg.norm(f_off))
        tn_on = float(np.linalg.norm(tau_on))
        tn_off = float(np.linalg.norm(tau_off))

        self.assertTrue(np.isfinite(fn_on) and np.isfinite(fn_off))
        self.assertTrue(np.isfinite(tn_on) and np.isfinite(tn_off))
        self.assertGreater(fn_on, 1e-9)
        self.assertLess(fn_on, 0.1)
        self.assertLess(tn_on, 0.1)

        # Shadowing OFF should change integrated force/torque in this asymmetric attitude.
        self.assertGreater(np.linalg.norm(f_off - f_on), 1e-10)
        self.assertGreater(np.linalg.norm(tau_off - tau_on), 1e-11)

    def test_gnss_slr_measurement_sanity_fixed_seed(self):
        _aero, prop, cfg = self._make_propagator(rtol=1e-8, atol=1e-10)
        r0 = np.array([7050e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.arange(0.0, 601.0, 60.0, dtype=float)
        env = _build_env(t_grid, density=0.0)
        truth = prop.propagate(x0, t_grid, env)

        sat_positions = _make_sat_positions(t_grid, n_sats=8)
        gnss = simulate_gnss_measurements(
            truth,
            t_grid,
            sat_positions,
            rng=np.random.default_rng(123),
            cadence_s=None,
            require_los=False,
            fov_half_angle_deg=None,
            dropout_prob=0.0,
            cycle_slip_prob_per_min=0.0,
        )
        self.assertGreater(int(gnss.values.size), 0)
        self.assertEqual(gnss.values.shape, gnss.sigmas.shape)
        self.assertTrue(np.all(np.isfinite(gnss.values)))
        self.assertTrue(np.all(gnss.sigmas > 0.0))

        stations = default_ilrs_stations()[:3]
        avail = np.ones((t_grid.size, len(stations)), dtype=bool)
        slr = simulate_slr_measurements(
            truth,
            t_grid,
            stations,
            rng=np.random.default_rng(321),
            cadence_s=None,
            availability_mask=avail,
            sigma_m=0.02,
        )
        self.assertGreater(int(slr.values.size), 0)
        self.assertEqual(slr.values.shape, slr.sigmas.shape)
        self.assertTrue(np.all(np.isfinite(slr.values)))
        self.assertTrue(np.all(slr.values > 0.0))


if __name__ == "__main__":
    unittest.main()
