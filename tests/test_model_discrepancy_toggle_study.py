import json
import tempfile
import unittest
from pathlib import Path

from scripts.model_discrepancy_toggle_study import run_study


class TestModelDiscrepancyToggleStudy(unittest.TestCase):
    def test_pairs_case_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            off = root / "off" / "att_a1"
            on = root / "on" / "att_a1"
            off.mkdir(parents=True, exist_ok=True)
            on.mkdir(parents=True, exist_ok=True)
            (off / "summary.json").write_text(
                json.dumps(
                    {
                        "case": "att_a1",
                        "ut_mean_err": 1.0,
                        "payload_impact": {"final_p95_geolocation_error_m": 10.0},
                    }
                ),
                encoding="utf-8",
            )
            (on / "summary.json").write_text(
                json.dumps(
                    {
                        "case": "att_a1",
                        "ut_mean_err": 1.5,
                        "payload_impact": {"final_p95_geolocation_error_m": 11.0},
                    }
                ),
                encoding="utf-8",
            )
            rep = run_study(root / "off", root / "on")
            self.assertTrue(bool(rep["pass"]))
            self.assertEqual(rep["paired_case_count"], 1)
            self.assertAlmostEqual(rep["cases"][0]["delta"]["ut_mean_err"], 0.5)


if __name__ == "__main__":
    unittest.main()
