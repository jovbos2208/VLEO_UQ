import unittest

from scripts.validate_tolerance_sweep import run_sweep


class TestToleranceSweep(unittest.TestCase):
    def test_tolerance_sweep_converges_for_mission_and_attitude(self):
        result = run_sweep(
            duration_s=300.0,
            dt_s=20.0,
            rtol_values=[1e-5, 1e-6, 1e-7],
            ref_rtol=1e-9,
            ref_atol=1e-11,
        )
        self.assertTrue(result["pass"])
        self.assertIn("mission", result["cases"])
        self.assertIn("attitude", result["cases"])
        for name in ("mission", "attitude"):
            self.assertTrue(result["cases"][name]["pass"])
            metrics = result["cases"][name]["metrics"]
            self.assertGreaterEqual(len(metrics), 2)
            coarse = metrics[0]
            fine = metrics[-1]
            self.assertLessEqual(fine["pos_final_m"], coarse["pos_final_m"] + 1e-12)
            self.assertLessEqual(fine["vel_final_mps"], coarse["vel_final_mps"] + 1e-15)
            if "att_final_deg" in fine:
                self.assertLessEqual(fine["att_final_deg"], coarse["att_final_deg"] + 1e-12)


if __name__ == "__main__":
    unittest.main()
