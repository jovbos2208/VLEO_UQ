import unittest

from scripts.validate_coverage_short_arc import run_validation


class TestCoverageValidation(unittest.TestCase):
    def test_ut_stm_coverage_against_mc(self):
        result = run_validation(
            mc_particles=128,
            duration_s=300.0,
            dt_s=30.0,
            seed=21,
        )
        self.assertTrue(result["pass"])
        self.assertTrue(result["ut_pass"])
        self.assertTrue(result["stm_pass"])

        for key in ("ut_coverage", "stm_coverage"):
            cov = result[key]
            self.assertIn("1sigma", cov)
            self.assertIn("2sigma", cov)
            self.assertIn("3sigma", cov)
            self.assertLessEqual(cov["1sigma"], cov["2sigma"])
            self.assertLessEqual(cov["2sigma"], cov["3sigma"])


if __name__ == "__main__":
    unittest.main()
