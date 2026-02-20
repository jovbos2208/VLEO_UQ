import json
import tempfile
import unittest
from pathlib import Path

from scripts.generate_campaign_report import build_campaign_report


class TestCampaignReport(unittest.TestCase):
    def test_phase_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            for case in (
                "att_a1_detumble_aero_only",
                "att_a2_nadir_pointing",
                "od_c1_gnss_continuous_bestcase",
                "od_c2_gnss_realistic_outages",
                "od_c3_gnss_plus_accel",
            ):
                d = run_dir / case
                d.mkdir(parents=True, exist_ok=True)
                (d / "summary.json").write_text(json.dumps({"case": case}), encoding="utf-8")
            report = build_campaign_report(run_dir)
            self.assertTrue(bool(report["phases"]["phase_a"]["pass"]))
            self.assertFalse(bool(report["pass"]))  # B/C still missing


if __name__ == "__main__":
    unittest.main()
