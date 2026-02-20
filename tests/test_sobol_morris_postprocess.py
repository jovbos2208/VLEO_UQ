import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.postprocess_sobol_morris import (
    analyze_records,
    build_records,
    discover_summary_paths,
)


class TestSobolMorrisPostprocess(unittest.TestCase):
    def test_dominant_parameter_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(123)
            case_dir = root / "run_a" / "case_x"
            case_dir.mkdir(parents=True, exist_ok=True)
            for i in range(64):
                x1 = float(rng.normal())
                x2 = float(rng.normal())
                y = 4.0 * x1 + 0.1 * x2 + 0.05 * float(rng.normal())
                summary = {
                    "case": "case_x",
                    "uq_parameter_draw": {
                        "seed": i,
                        "p1": {"value": x1},
                        "p2": {"value": x2},
                    },
                    "ut_mean_err": y,
                }
                (case_dir / f"summary_{i}.json").write_text(json.dumps(summary), encoding="utf-8")
            # Also include canonical summary.json paths for discovery.
            for i in range(64):
                p = root / f"run_{i:03d}" / "case_x"
                p.mkdir(parents=True, exist_ok=True)
                x1 = float(rng.normal())
                x2 = float(rng.normal())
                y = 4.0 * x1 + 0.1 * x2 + 0.05 * float(rng.normal())
                summary = {
                    "case": "case_x",
                    "uq_parameter_draw": {
                        "seed": i,
                        "p1": {"value": x1},
                        "p2": {"value": x2},
                    },
                    "ut_mean_err": y,
                }
                (p / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

            paths = discover_summary_paths([root])
            records = build_records(paths)
            report = analyze_records(records, qois=["ut_mean_err"], min_samples=20, seed=7)
            rows = report["cases"]["case_x"]["ut_mean_err"]["sobol_proxy"]["parameters"]
            self.assertGreaterEqual(len(rows), 2)
            self.assertEqual(rows[0]["name"], "p1")


if __name__ == "__main__":
    unittest.main()
