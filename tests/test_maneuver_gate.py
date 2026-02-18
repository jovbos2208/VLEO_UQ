import json
import tempfile
import unittest
from pathlib import Path

from scripts.validate_maneuver_realism_gate import (
    _authority_eval,
    _evaluate_post_stats,
    _run_post_gate_for_scenario,
)


class TestManeuverGate(unittest.TestCase):
    def test_authority_eval_pass_and_fail(self):
        ok = _authority_eval(
            max_abs_torque_axis_nm=1e-5,
            inertia_axis_kgm2=0.15,
            full_turn_deg=360.0,
            duration_s=1000.0,
            min_authority_ratio=0.5,
        )
        self.assertTrue(ok["pass"])
        bad = _authority_eval(
            max_abs_torque_axis_nm=1e-9,
            inertia_axis_kgm2=0.15,
            full_turn_deg=360.0,
            duration_s=1000.0,
            min_authority_ratio=0.5,
        )
        self.assertFalse(bad["pass"])

    def test_post_stats_eval(self):
        stats_ok = {
            "saturation_time_fraction": 0.2,
            "command_rate_max_abs_deg_s": 4.5,
            "eta_rate_max_deg_s": 5.0,
            "reached_fraction_180": 0.9,
            "full_turn_reached": True,
        }
        passed, reasons = _evaluate_post_stats(
            stats=stats_ok,
            max_saturation_fraction=0.85,
            max_rate_limit_ratio=1.05,
            min_reached_fraction_180=0.5,
            require_full_turn=True,
        )
        self.assertTrue(passed)
        self.assertEqual(reasons, [])

        stats_bad = {
            "saturation_time_fraction": 0.95,
            "command_rate_max_abs_deg_s": 7.0,
            "eta_rate_max_deg_s": 5.0,
            "reached_fraction_180": 0.2,
            "full_turn_reached": False,
        }
        passed, reasons = _evaluate_post_stats(
            stats=stats_bad,
            max_saturation_fraction=0.85,
            max_rate_limit_ratio=1.05,
            min_reached_fraction_180=0.5,
            require_full_turn=True,
        )
        self.assertFalse(passed)
        self.assertGreaterEqual(len(reasons), 3)

    def test_run_post_gate_for_scenario(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            case_dir = root / "att_a8_roll_rotation_baseline"
            case_dir.mkdir(parents=True, exist_ok=True)
            summary = {
                "maneuver_actuator_stats": {
                    "saturation_time_fraction": 0.3,
                    "command_rate_max_abs_deg_s": 4.0,
                    "eta_rate_max_deg_s": 5.0,
                    "reached_fraction_180": 0.8,
                    "full_turn_reached": True,
                }
            }
            (case_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
            scenario = {
                "name": "att_a8_roll_rotation_baseline",
                "type": "attitude",
                "maneuver_axis": "roll",
                "wing_constant_deg": 8.0,
            }
            out = _run_post_gate_for_scenario(
                scenario,
                results_dir=root,
                pre_result_by_scenario={"att_a8_roll_rotation_baseline": {"pass": True}},
                max_saturation_fraction=0.85,
                max_rate_limit_ratio=1.05,
                min_reached_fraction_180=0.5,
                require_full_turn=True,
                require_pre_authority_pass=True,
            )
            self.assertTrue(out["pass"])


if __name__ == "__main__":
    unittest.main()
