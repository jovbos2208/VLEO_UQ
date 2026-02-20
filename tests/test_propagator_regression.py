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


def _build_env(
    t_grid: np.ndarray,
    density: float,
    sun_position_m: np.ndarray | None = None,
    moon_position_m: np.ndarray | None = None,
) -> list:
    env = []
    if sun_position_m is None:
        sun_position_m = np.zeros(3, dtype=float)
    if moon_position_m is None:
        moon_position_m = np.zeros(3, dtype=float)
    for _ in t_grid:
        e = EnvInputs()
        e.density = float(density)
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        e.sun_position_I_m = sun_position_m
        e.moon_position_I_m = moon_position_m
        e.albedo_ir_scale = 1.0
        e.tide_loading_accel_I_m_s2 = np.zeros(3, dtype=float)
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

    def test_j2_toggle_changes_trajectory_for_inclined_orbit(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_off = PropagatorConfig()
        cfg_off.rtol = 1e-9
        cfg_off.atol = 1e-10
        self.assertFalse(cfg_off.use_j2_perturbation)

        cfg_on = PropagatorConfig()
        cfg_on.rtol = 1e-9
        cfg_on.atol = 1e-10
        cfg_on.use_j2_perturbation = True

        prop_off = DeterministicPropagator(aero, vehicle, cfg_off)
        prop_on = DeterministicPropagator(aero, vehicle, cfg_on)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v_circ = math.sqrt(cfg_off.mu_earth_m3_s2 / np.linalg.norm(r0))
        inc = math.radians(63.4)
        v0 = np.array([0.0, v_circ * math.cos(inc), v_circ * math.sin(inc)], dtype=float)
        x0 = np.zeros(prop_off.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 3.0 * 3600.0, 181)
        env = _build_env(t_grid, density=0.0)

        states_off = prop_off.propagate(x0, t_grid, env)
        states_on = prop_on.propagate(x0, t_grid, env)

        final_pos_delta = float(np.linalg.norm(states_on[-1, 0:3] - states_off[-1, 0:3]))
        final_vel_delta = float(np.linalg.norm(states_on[-1, 3:6] - states_off[-1, 3:6]))
        self.assertGreater(final_pos_delta, 1.0)
        self.assertGreater(final_vel_delta, 1e-3)

    def test_third_body_toggles_change_trajectory(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_none = PropagatorConfig()
        cfg_none.rtol = 1e-9
        cfg_none.atol = 1e-10

        cfg_sun = PropagatorConfig()
        cfg_sun.rtol = 1e-9
        cfg_sun.atol = 1e-10
        cfg_sun.use_sun_third_body = True

        cfg_sun_scaled = PropagatorConfig()
        cfg_sun_scaled.rtol = 1e-9
        cfg_sun_scaled.atol = 1e-10
        cfg_sun_scaled.use_sun_third_body = True
        cfg_sun_scaled.sun_ephemeris_scale = 1.01

        cfg_moon = PropagatorConfig()
        cfg_moon.rtol = 1e-9
        cfg_moon.atol = 1e-10
        cfg_moon.use_moon_third_body = True

        cfg_both = PropagatorConfig()
        cfg_both.rtol = 1e-9
        cfg_both.atol = 1e-10
        cfg_both.use_sun_third_body = True
        cfg_both.use_moon_third_body = True

        prop_none = DeterministicPropagator(aero, vehicle, cfg_none)
        prop_sun = DeterministicPropagator(aero, vehicle, cfg_sun)
        prop_sun_scaled = DeterministicPropagator(aero, vehicle, cfg_sun_scaled)
        prop_moon = DeterministicPropagator(aero, vehicle, cfg_moon)
        prop_both = DeterministicPropagator(aero, vehicle, cfg_both)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v_circ = math.sqrt(cfg_none.mu_earth_m3_s2 / np.linalg.norm(r0))
        inc = math.radians(30.0)
        v0 = np.array([0.0, v_circ * math.cos(inc), v_circ * math.sin(inc)], dtype=float)
        x0 = np.zeros(prop_none.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 6.0 * 3600.0, 241)
        sun_position_m = np.array([149597870700.0, 1.0e9, 0.0], dtype=float)
        moon_position_m = np.array([3.844e8, 2.0e7, -1.0e7], dtype=float)
        env = _build_env(
            t_grid,
            density=0.0,
            sun_position_m=sun_position_m,
            moon_position_m=moon_position_m,
        )

        states_none = prop_none.propagate(x0, t_grid, env)
        states_sun = prop_sun.propagate(x0, t_grid, env)
        states_sun_scaled = prop_sun_scaled.propagate(x0, t_grid, env)
        states_moon = prop_moon.propagate(x0, t_grid, env)
        states_both = prop_both.propagate(x0, t_grid, env)

        delta_sun_pos = float(np.linalg.norm(states_sun[-1, 0:3] - states_none[-1, 0:3]))
        delta_moon_pos = float(np.linalg.norm(states_moon[-1, 0:3] - states_none[-1, 0:3]))
        delta_both_pos = float(np.linalg.norm(states_both[-1, 0:3] - states_none[-1, 0:3]))
        delta_sun_scale_pos = float(np.linalg.norm(states_sun_scaled[-1, 0:3] - states_sun[-1, 0:3]))

        self.assertGreater(delta_sun_pos, 5.0)
        self.assertGreater(delta_moon_pos, 10.0)
        self.assertGreater(delta_both_pos, delta_sun_pos)
        self.assertGreater(delta_both_pos, delta_moon_pos)
        self.assertGreater(delta_sun_scale_pos, 0.01)

    def test_j3_j4_toggle_changes_trajectory(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_base = PropagatorConfig()
        cfg_base.rtol = 1e-9
        cfg_base.atol = 1e-10

        cfg_zonal = PropagatorConfig()
        cfg_zonal.rtol = 1e-9
        cfg_zonal.atol = 1e-10
        cfg_zonal.use_j2_perturbation = True
        cfg_zonal.use_j3_perturbation = True
        cfg_zonal.use_j4_perturbation = True

        prop_base = DeterministicPropagator(aero, vehicle, cfg_base)
        prop_zonal = DeterministicPropagator(aero, vehicle, cfg_zonal)

        r0 = np.array([7100e3, 0.0, 0.0], dtype=float)
        v_circ = math.sqrt(cfg_base.mu_earth_m3_s2 / np.linalg.norm(r0))
        inc = math.radians(74.0)
        v0 = np.array([0.0, v_circ * math.cos(inc), v_circ * math.sin(inc)], dtype=float)
        x0 = np.zeros(prop_base.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 5.0 * 3600.0, 201)
        env = _build_env(t_grid, density=0.0)
        states_base = prop_base.propagate(x0, t_grid, env)
        states_zonal = prop_zonal.propagate(x0, t_grid, env)
        delta_pos = float(np.linalg.norm(states_zonal[-1, 0:3] - states_base[-1, 0:3]))
        delta_vel = float(np.linalg.norm(states_zonal[-1, 3:6] - states_base[-1, 3:6]))
        self.assertGreater(delta_pos, 1.0)
        self.assertGreater(delta_vel, 1e-4)

    def test_srp_toggle_changes_trajectory(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_off = PropagatorConfig()
        cfg_off.rtol = 1e-9
        cfg_off.atol = 1e-10
        cfg_off.use_srp_acceleration = False

        cfg_on = PropagatorConfig()
        cfg_on.rtol = 1e-9
        cfg_on.atol = 1e-10
        cfg_on.use_srp_acceleration = True
        cfg_on.srp_area_m2 = 0.35
        cfg_on.srp_cr = 1.5

        prop_off = DeterministicPropagator(aero, vehicle, cfg_off)
        prop_on = DeterministicPropagator(aero, vehicle, cfg_on)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg_on.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop_off.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0
        t_grid = np.linspace(0.0, 8.0 * 3600.0, 321)
        env = _build_env(t_grid, density=0.0, sun_position_m=np.array([149597870700.0, 0.0, 0.0], dtype=float))

        states_off = prop_off.propagate(x0, t_grid, env)
        states_on = prop_on.propagate(x0, t_grid, env)
        final_pos_delta = float(np.linalg.norm(states_on[-1, 0:3] - states_off[-1, 0:3]))
        self.assertGreater(final_pos_delta, 0.1)

    def test_albedo_ir_toggle_changes_trajectory(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_off = PropagatorConfig()
        cfg_off.rtol = 1e-9
        cfg_off.atol = 1e-10
        cfg_off.use_albedo_ir_acceleration = False

        cfg_on = PropagatorConfig()
        cfg_on.rtol = 1e-9
        cfg_on.atol = 1e-10
        cfg_on.use_albedo_ir_acceleration = True
        cfg_on.albedo_ir_area_m2 = 0.35
        cfg_on.albedo_ir_cr = 1.3

        prop_off = DeterministicPropagator(aero, vehicle, cfg_off)
        prop_on = DeterministicPropagator(aero, vehicle, cfg_on)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg_on.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop_off.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0
        t_grid = np.linspace(0.0, 8.0 * 3600.0, 321)
        env = _build_env(t_grid, density=0.0)
        for e in env:
            e.albedo_ir_scale = 1.0

        states_off = prop_off.propagate(x0, t_grid, env)
        states_on = prop_on.propagate(x0, t_grid, env)
        final_pos_delta = float(np.linalg.norm(states_on[-1, 0:3] - states_off[-1, 0:3]))
        self.assertGreater(final_pos_delta, 0.1)

    def test_tide_loading_toggle_changes_trajectory(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_off = PropagatorConfig()
        cfg_off.rtol = 1e-9
        cfg_off.atol = 1e-10
        cfg_off.use_tide_loading_acceleration = False

        cfg_on = PropagatorConfig()
        cfg_on.rtol = 1e-9
        cfg_on.atol = 1e-10
        cfg_on.use_tide_loading_acceleration = True
        cfg_on.tide_loading_scale = 1.0

        prop_off = DeterministicPropagator(aero, vehicle, cfg_off)
        prop_on = DeterministicPropagator(aero, vehicle, cfg_on)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg_on.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop_off.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0
        t_grid = np.linspace(0.0, 4.0 * 3600.0, 241)
        env = _build_env(t_grid, density=0.0)
        for i, e in enumerate(env):
            phase = 2.0 * math.pi * (i / max(1, len(env) - 1))
            e.tide_loading_accel_I_m_s2 = np.array([1e-8 * math.sin(phase), 8e-9 * math.cos(phase), 0.0], dtype=float)

        states_off = prop_off.propagate(x0, t_grid, env)
        states_on = prop_on.propagate(x0, t_grid, env)
        final_pos_delta = float(np.linalg.norm(states_on[-1, 0:3] - states_off[-1, 0:3]))
        self.assertGreater(final_pos_delta, 1e-2)

    def test_magnetic_torque_changes_angular_rate(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()

        cfg_off = PropagatorConfig()
        cfg_off.rtol = 1e-9
        cfg_off.atol = 1e-10
        cfg_off.use_magnetic_torque = False

        cfg_on = PropagatorConfig()
        cfg_on.rtol = 1e-9
        cfg_on.atol = 1e-10
        cfg_on.use_magnetic_torque = True
        cfg_on.residual_dipole_B_A_m2 = np.array([0.15, -0.08, 0.04], dtype=float)

        prop_off = DeterministicPropagator(aero, vehicle, cfg_off)
        prop_on = DeterministicPropagator(aero, vehicle, cfg_on)

        r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
        v0 = np.array([0.0, math.sqrt(cfg_on.mu_earth_m3_s2 / np.linalg.norm(r0)), 0.0], dtype=float)
        x0 = np.zeros(prop_off.state_size, dtype=float)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0
        t_grid = np.linspace(0.0, 1200.0, 121)
        env = _build_env(t_grid, density=0.0)
        for e in env:
            e.magnetic_field_I_T = np.array([2.2e-5, -1.1e-5, 3.4e-5], dtype=float)

        states_off = prop_off.propagate(x0, t_grid, env)
        states_on = prop_on.propagate(x0, t_grid, env)
        delta_w = float(np.linalg.norm(states_on[-1, 10:13] - states_off[-1, 10:13]))
        self.assertGreater(delta_w, 1e-6)

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
