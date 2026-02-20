import json
import tempfile
import unittest
from pathlib import Path

from scripts.generate_maneuver_realism_report import generate_report


class TestManeuverRealismReport(unittest.TestCase):
    def test_report_reads_gate_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "maneuver_gate_summary.json").write_text(json.dumps({"pass": True}), encoding="utf-8")
            (d / "maneuver_gate_pre.json").write_text(json.dumps({"results": [], "all_pass": True}), encoding="utf-8")
            (d / "maneuver_gate_post.json").write_text(
                json.dumps({"results": [{"scenario": "att_a1", "pass": True}], "all_pass": True}),
                encoding="utf-8",
            )
            out = generate_report(d)
            self.assertTrue(bool(out["pass"]))
            self.assertEqual(out["post_count"], 1)


if __name__ == "__main__":
    unittest.main()
