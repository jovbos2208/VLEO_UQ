import unittest

import numpy as np

from scripts.attitude_uq_modes import summarize_attitude_uq_modes


class TestAttitudeUqModes(unittest.TestCase):
    def test_fast_and_full_modes_available(self):
        nt = 10
        npart = 5
        t = np.linspace(0.0, 90.0, nt)
        det = np.zeros((nt, 18), dtype=float)
        det[:, 0] = 7000e3
        det[:, 6] = 1.0
        mc = np.repeat(det[:, None, :], npart, axis=1)
        mc[:, :, 7] += np.linspace(0.0, 1e-3, npart)[None, :]
        out = summarize_attitude_uq_modes(
            t_grid=t,
            det_states=det,
            mc_states=mc,
            seed=123,
            fast_tau_s=150.0,
            fast_sigma_rad=2e-4,
        )
        self.assertTrue(bool(out.get("available")))
        self.assertTrue(bool(out["fast_mode"].get("available")))
        self.assertTrue(bool(out["full_mode"].get("available")))


if __name__ == "__main__":
    unittest.main()
