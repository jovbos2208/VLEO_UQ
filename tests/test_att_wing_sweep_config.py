import unittest

from scripts.build_att_wing_sweep_config import enforce_recommended_range, parse_angles


class TestAttWingSweepConfig(unittest.TestCase):
    def test_parse_angles(self):
        self.assertEqual(parse_angles("0, 5,10"), [0.0, 5.0, 10.0])

    def test_enforce_recommended_range_blocks_out_of_range(self):
        with self.assertRaises(ValueError):
            enforce_recommended_range(
                angles_deg=[0.0, 60.0, 65.0],
                recommended_abs_max_deg=60.0,
                allow_out_of_recommended_range=False,
            )

    def test_enforce_recommended_range_allows_override(self):
        vals = enforce_recommended_range(
            angles_deg=[0.0, 60.0, 65.0],
            recommended_abs_max_deg=60.0,
            allow_out_of_recommended_range=True,
        )
        self.assertEqual(vals, [0.0, 60.0, 65.0])


if __name__ == "__main__":
    unittest.main()
