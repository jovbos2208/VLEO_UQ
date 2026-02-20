import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagationEvent,
    PropagatorConfig,
    VehicleParams,
    default_geometry,
    propagate_with_events,
    sample_event_realizations,
)


class TestPropagationEvents(unittest.TestCase):
    def _make_case(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)
        vehicle = VehicleParams()
        config = PropagatorConfig()
        prop = DeterministicPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])
        x0 = np.zeros(prop.state_size)
        x0[0:3] = r0
        x0[3:6] = v0
        x0[6] = 1.0

        t_grid = np.linspace(0.0, 200.0, 21)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)
        return prop, x0, t_grid, env

    def test_segmented_no_events_matches_direct(self):
        prop, x0, t_grid, env = self._make_case()
        x_direct = prop.propagate(x0, t_grid, env)
        x_events = propagate_with_events(prop, x0, t_grid, env, [])
        self.assertLess(np.max(np.abs(x_direct - x_events)), 1e-6)

    def test_delta_v_event_applies_velocity_jump(self):
        prop, x0, t_grid, env = self._make_case()
        x_base = prop.propagate(x0, t_grid, env)

        dv = np.array([0.0, 0.5, 0.0])
        event = PropagationEvent(
            t_s=float(t_grid[10]),
            kind="delta_v",
            params={"dv_eci_m_s": dv},
        )
        x_events = propagate_with_events(prop, x0, t_grid, env, [event])

        delta = x_events[10, 3:6] - x_base[10, 3:6]
        self.assertLess(np.linalg.norm(delta - dv), 1e-9)

    def test_event_uncertainty_sampling_is_reproducible(self):
        events = [
            PropagationEvent(
                t_s=100.0,
                kind="delta_v",
                params={
                    "dv_eci_m_s": [0.0, 0.2, 0.0],
                    "sigma_time_s": 0.5,
                    "sigma_scale": 0.1,
                },
            )
        ]
        a = sample_event_realizations(events, np.random.default_rng(123))
        b = sample_event_realizations(events, np.random.default_rng(123))
        self.assertAlmostEqual(float(a[0].t_s), float(b[0].t_s), places=12)
        np.testing.assert_allclose(np.asarray(a[0].params["dv_eci_m_s"]), np.asarray(b[0].params["dv_eci_m_s"]))

    def test_finite_burn_event_accumulates_velocity(self):
        prop, x0, t_grid, env = self._make_case()
        x_base = prop.propagate(x0, t_grid, env)
        accel = np.array([0.0, 2.0e-3, 0.0], dtype=float)
        duration_s = 60.0
        ev = PropagationEvent(
            t_s=40.0,
            kind="finite_burn",
            params={
                "accel_eci_m_s2": accel,
                "duration_s": duration_s,
            },
        )
        x_burn = propagate_with_events(prop, x0, t_grid, env, [ev])
        dv_expected = accel * duration_s
        burn_end_idx = int(np.argmin(np.abs(t_grid - (40.0 + duration_s))))
        dv_actual = x_burn[burn_end_idx, 3:6] - x_base[burn_end_idx, 3:6]
        self.assertLess(float(np.linalg.norm(dv_actual - dv_expected)), 5e-4)

    def test_finite_burn_sampling_reproducible(self):
        events = [
            PropagationEvent(
                t_s=80.0,
                kind="finite_burn",
                params={
                    "accel_eci_m_s2": [0.0, 3e-4, 0.0],
                    "duration_s": 45.0,
                    "sigma_time_s": 0.5,
                    "sigma_scale": 0.15,
                    "sigma_duration_s": 2.0,
                },
            )
        ]
        a = sample_event_realizations(events, np.random.default_rng(9))
        b = sample_event_realizations(events, np.random.default_rng(9))
        self.assertAlmostEqual(float(a[0].t_s), float(b[0].t_s), places=12)
        self.assertAlmostEqual(float(a[0].params["duration_s"]), float(b[0].params["duration_s"]), places=12)
        np.testing.assert_allclose(
            np.asarray(a[0].params["accel_eci_m_s2"]),
            np.asarray(b[0].params["accel_eci_m_s2"]),
        )


if __name__ == "__main__":
    unittest.main()
