import json
import tempfile
import unittest
from pathlib import Path

from scripts.run_phase_campaign import PHASE_SCENARIO_IDS, _apply_config_overrides, _selected_phases


class TestPhaseCampaign(unittest.TestCase):
    def test_phase_mapping(self):
        self.assertEqual(_selected_phases("all"), ["phase_a", "phase_b", "phase_c"])
        self.assertEqual(_selected_phases("phase_b"), ["phase_b"])

    def test_required_case_counts(self):
        self.assertEqual(len(PHASE_SCENARIO_IDS["phase_a"]), 5)
        self.assertEqual(len(PHASE_SCENARIO_IDS["phase_b"]), 4)
        self.assertEqual(len(PHASE_SCENARIO_IDS["phase_c"]), 4)

    def test_config_overrides(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = Path(tmp) / "cfg.json"
            payload = {
                "scenarios": [
                    {"name": "a", "particles": 256, "pod_uq": True},
                    {"name": "b", "particles": 12, "pod_uq": False},
                    {"name": "c", "pod_uq": True},
                ]
            }
            cfg.write_text(json.dumps(payload), encoding="utf-8")
            _apply_config_overrides(cfg, max_particles=32, disable_pod=True)
            out = json.loads(cfg.read_text(encoding="utf-8"))
            self.assertEqual(out["scenarios"][0]["particles"], 32)
            self.assertEqual(out["scenarios"][1]["particles"], 12)
            self.assertEqual(out["scenarios"][2]["particles"], 32)
            self.assertFalse(bool(out["scenarios"][0]["pod_uq"]))
            self.assertFalse(bool(out["scenarios"][2]["pod_uq"]))


if __name__ == "__main__":
    unittest.main()
