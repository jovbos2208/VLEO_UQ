import datetime as dt
import tempfile
import unittest

import numpy as np

from vleo_uq import parse_sp3


class TestSp3Parser(unittest.TestCase):
    def test_parse_sp3_minimal(self):
        sp3 = "\n".join(
            [
                "#cP2020  1  1  0  0  0.00000000      96 ORBIT IGS14 HLM 0 0 0",
                "*  2020  1  1  0  0  0.00000000",
                "PG01  10000.000000 20000.000000 30000.000000 0.000000",
                "PG02  11000.000000 21000.000000 31000.000000 0.000000",
                "*  2020  1  1  0 15  0.00000000",
                "PG01  10010.000000 20010.000000 30010.000000 0.000000",
                "PG02  11010.000000 21010.000000 31010.000000 0.000000",
            ]
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
            tmp.write(sp3)
            path = tmp.name

        eph = parse_sp3(path, constellation_prefix=("G",))
        self.assertEqual(len(eph.epochs), 2)
        self.assertEqual(eph.positions_ecef_m.shape, (2, 2, 3))
        self.assertIsNotNone(eph.clocks_s)

        t_grid = [dt.datetime(2020, 1, 1, 0, 7, 30)]
        sample = eph.sample(t_grid, frame="ecef")
        self.assertEqual(sample.shape, (1, 2, 3))
        self.assertTrue(np.all(np.isfinite(sample)))

        clocks = eph.sample_clock(t_grid)
        self.assertEqual(clocks.shape, (1, 2))


if __name__ == "__main__":
    unittest.main()
