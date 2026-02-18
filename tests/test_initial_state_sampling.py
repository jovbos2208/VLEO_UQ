import unittest

import numpy as np

from scripts.run_case_studies import sample_mc_initial_states


class TestInitialStateSampling(unittest.TestCase):
    def test_sampling_is_reproducible_for_same_seed(self):
        x0 = np.zeros(18, dtype=float)
        x0[6] = 1.0
        P0 = np.diag([4.0, 4.0, 4.0] + [0.01] * 3 + [1e-6] * 4 + [1e-4] * 3 + [0.01] * 2 + [0.25] * 3)
        rng_a = np.random.default_rng(123)
        rng_b = np.random.default_rng(123)
        Xa = sample_mc_initial_states(x0, P0, particles=32, rng=rng_a, freeze_attitude=False)
        Xb = sample_mc_initial_states(x0, P0, particles=32, rng=rng_b, freeze_attitude=False)
        self.assertTrue(np.allclose(Xa, Xb))

    def test_freeze_attitude_keeps_nominal_attitude_and_rates(self):
        x0 = np.zeros(18, dtype=float)
        x0[6] = 1.0
        x0[10:13] = np.array([0.1, -0.2, 0.3])
        P0 = np.eye(18) * 0.1
        X = sample_mc_initial_states(
            x0,
            P0,
            particles=16,
            rng=np.random.default_rng(7),
            freeze_attitude=True,
        )
        self.assertTrue(np.allclose(X[:, 6:10], x0[6:10]))
        self.assertTrue(np.allclose(X[:, 10:13], x0[10:13]))


if __name__ == "__main__":
    unittest.main()
