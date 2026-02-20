import json
import tempfile
import unittest
from pathlib import Path

from scripts.aggregate_run_summary import build_standardized_run_summary


class TestAggregateRunSummary(unittest.TestCase):
    def test_metadata_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            case_dir = run_dir / "mission"
            case_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "case": "mission",
                "metadata": {
                    "git_hash": "abc1234",
                    "scenario_id": "MISSION",
                    "seed": 42,
                    "toggles": {},
                },
            }
            (case_dir / "summary.json").write_text(json.dumps(payload), encoding="utf-8")
            report = build_standardized_run_summary(run_dir)
            self.assertEqual(report["case_count"], 1)
            self.assertTrue(bool(report["pass"]))


if __name__ == "__main__":
    unittest.main()
