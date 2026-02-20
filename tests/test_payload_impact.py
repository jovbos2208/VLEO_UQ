import unittest

import numpy as np

from scripts.payload_impact import (
    compute_payload_impact_metrics,
    compute_payload_impact_proxy_from_stats,
)


class TestPayloadImpact(unittest.TestCase):
    def test_proxy_metrics_present(self):
        mc_mean = np.zeros((5, 18), dtype=float)
        mc_std = np.zeros((5, 18), dtype=float)
        mc_std[-1, 0:3] = np.array([3.0, 4.0, 0.0])
        mc_std[-1, 7:10] = np.array([1e-6, 2e-6, 0.0])
        out = compute_payload_impact_proxy_from_stats(mc_mean=mc_mean, mc_std=mc_std)
        self.assertTrue(bool(out.get("available")))
        self.assertGreater(float(out.get("final_p95_geolocation_error_m", 0.0)), 0.0)

    def test_ensemble_metrics_structure(self):
        nt = 6
        npart = 4
        t = np.linspace(0.0, 50.0, nt)
        det = np.zeros((nt, 18), dtype=float)
        det[:, 0] = 7000e3
        det[:, 6] = 1.0
        mc = np.repeat(det[:, None, :], npart, axis=1)
        mc[:, :, 1] += np.linspace(0.0, 5.0, npart)[None, :]
        out = compute_payload_impact_metrics(det_states=det, mc_states=mc, t_grid=t, max_samples=4)
        self.assertIn("available", out)
        self.assertIn("samples", out)
        self.assertGreaterEqual(len(out["samples"]), 1)


if __name__ == "__main__":
    unittest.main()
