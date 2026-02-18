import json
import tempfile
import unittest
from pathlib import Path

from scripts.validate_measurement_realism_gate import _validate_case
from scripts.validate_pod_uq_harness import _validate_case_pod_summary


class TestPodValidationGates(unittest.TestCase):
    def test_mp4_harness_case_validation(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            pod_dir = root / "case_a" / "pod"
            pod_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "pod_summaries": [
                    {
                        "radial_rms_m": 1.2,
                        "coverage_sigma": {"1sigma": {}, "2sigma": {}, "3sigma": {}},
                        "coverage_3sigma": {},
                    }
                ]
            }
            path = pod_dir / "summary.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            out = _validate_case_pod_summary(path)
            self.assertTrue(out["pass"])

    def test_vg5_case_validation(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            pod_dir = root / "case_b" / "pod"
            pod_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "pod_summaries": [
                    {
                        "measurement_residuals": {
                            "gnss_code": {
                                "count": 120,
                                "mean_norm": 0.01,
                                "std_norm": 1.05,
                                "rms_norm": 1.1,
                            }
                        }
                    }
                ]
            }
            path = pod_dir / "summary.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            out = _validate_case(
                path,
                min_count=20,
                max_abs_norm_mean=0.35,
                min_norm_std=0.5,
                max_norm_std=1.8,
                max_norm_rms=2.0,
            )
            self.assertTrue(out["pass"])


if __name__ == "__main__":
    unittest.main()
