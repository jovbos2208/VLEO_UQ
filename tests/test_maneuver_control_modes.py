import unittest

import numpy as np

from scripts.run_mc_job_array import (
    _extract_maneuver_runtime_config,
    _maneuver_control_mode,
    _propagate_member_with_axis_runtime,
)


class _DummyEnv:
    eta1_rad: float = 0.0
    eta2_rad: float = 0.0


class _DummyDetPropagator:
    def propagate(self, x, t_grid, env):
        x = np.asarray(x, dtype=float).reshape(-1)
        nt = int(np.asarray(t_grid).size)
        out = np.zeros((nt, x.size), dtype=float)
        out[:] = x
        return out


class TestManeuverControlModes(unittest.TestCase):
    def test_control_mode_from_toggle(self):
        self.assertEqual(_maneuver_control_mode({"catalog_control_toggles": ["CTL_ATT_AERO_RATE"]}), "rate")
        self.assertEqual(_maneuver_control_mode({"catalog_control_toggles": ["CTL_ATT_AERO_POINT"]}), "point")
        self.assertEqual(_maneuver_control_mode({}), "open_loop")

    def test_extract_runtime_cfg_infers_mode(self):
        cfg = _extract_maneuver_runtime_config(
            {
                "maneuver_axis": "roll",
                "wing_constant_deg": 8.0,
                "catalog_control_toggles": ["CTL_ATT_AERO_RATE"],
            }
        )
        self.assertIsNotNone(cfg)
        self.assertEqual(cfg["control_mode"], "rate")
        self.assertGreater(cfg["rate_kp_cmd_per_deg_s"], 0.0)

    def test_point_mode_target_command(self):
        x0 = np.zeros(13, dtype=float)
        x0[6] = 1.0  # identity quaternion [w, x, y, z]
        t_grid = np.array([0.0, 1.0, 2.0], dtype=float)
        env = [_DummyEnv(), _DummyEnv(), _DummyEnv()]
        out, diag = _propagate_member_with_axis_runtime(
            _DummyDetPropagator(),
            x0,
            t_grid,
            env,
            axis="roll",
            control_mode="point",
            wing_constant_deg=10.0,
            half_turn_deg=180.0,
            full_turn_deg=360.0,
            eta_max_deg=110.0,
            eta_rate_max_deg_s=100.0,
            eta_accel_max_deg_s2=1000.0,
            point_kp_cmd_per_deg=0.05,
            point_kd_cmd_per_deg_s=0.0,
            rate_target_deg_s=1.0,
            rate_kp_cmd_per_deg_s=12.0,
            inertia_B=np.eye(3),
        )
        self.assertEqual(out.shape, (3, 13))
        self.assertEqual(diag["control_mode"], "point")
        self.assertAlmostEqual(float(diag["eta_target_deg"][0]), 9.0, places=6)

    def test_rate_mode_target_command(self):
        x0 = np.zeros(13, dtype=float)
        x0[6] = 1.0  # identity quaternion [w, x, y, z]
        t_grid = np.array([0.0, 1.0, 2.0], dtype=float)
        env = [_DummyEnv(), _DummyEnv(), _DummyEnv()]
        out, diag = _propagate_member_with_axis_runtime(
            _DummyDetPropagator(),
            x0,
            t_grid,
            env,
            axis="roll",
            control_mode="rate",
            wing_constant_deg=10.0,
            half_turn_deg=180.0,
            full_turn_deg=360.0,
            eta_max_deg=110.0,
            eta_rate_max_deg_s=100.0,
            eta_accel_max_deg_s2=1000.0,
            point_kp_cmd_per_deg=0.5,
            point_kd_cmd_per_deg_s=2.0,
            rate_target_deg_s=0.2,
            rate_kp_cmd_per_deg_s=20.0,
            inertia_B=np.eye(3),
        )
        self.assertEqual(out.shape, (3, 13))
        self.assertEqual(diag["control_mode"], "rate")
        self.assertAlmostEqual(float(diag["eta_target_deg"][0]), 4.0, places=6)


if __name__ == "__main__":
    unittest.main()
