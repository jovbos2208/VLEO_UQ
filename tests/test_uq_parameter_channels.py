import unittest

from scripts.uq_parameter_channels import apply_propagator_overrides, apply_uq_parameter_channels


class _Cfg:
    rho_bias_tau_s = 0.0
    rho_bias_sigma = 0.0


class TestUqParameterChannels(unittest.TestCase):
    def test_reproducible_parameter_draws(self):
        scenario = {
            "name": "od_c5_gnss_slr_hybrid_anchor",
            "spacecraft_mass_kg": 12.0,
            "spacecraft_inertia_kgm2": [0.15, 0.12, 0.2],
            "density_scale": 1.0,
            "propagator_overrides": {"rho_bias_tau_s": 7200.0, "rho_bias_sigma": 0.08},
            "uq_mass_rel_sigma": 0.1,
            "uq_inertia_rel_sigma": 0.15,
            "uq_density_scale_rel_sigma": 0.2,
            "uq_rho_bias_tau_rel_sigma": 0.25,
            "uq_rho_bias_sigma_rel_sigma": 0.3,
        }
        a, draw_a = apply_uq_parameter_channels(scenario, seed=42)
        b, draw_b = apply_uq_parameter_channels(scenario, seed=42)
        self.assertEqual(draw_a, draw_b)
        self.assertEqual(a["uq_parameter_draw"], b["uq_parameter_draw"])
        self.assertGreater(a["spacecraft_mass_kg"], 0.0)
        self.assertGreater(a["density_scale"], 0.0)
        self.assertGreater(a["propagator_overrides"]["rho_bias_tau_s"], 0.0)
        self.assertGreater(a["propagator_overrides"]["rho_bias_sigma"], 0.0)

    def test_no_sigma_no_draw(self):
        scenario = {"name": "att_a2_nadir_pointing", "spacecraft_mass_kg": 12.0}
        out, draw = apply_uq_parameter_channels(scenario, seed=7)
        self.assertEqual(draw, {})
        self.assertNotIn("uq_parameter_draw", out)

    def test_apply_propagator_overrides(self):
        cfg = _Cfg()
        apply_propagator_overrides(cfg, {"rho_bias_tau_s": 8000.0, "rho_bias_sigma": 0.1, "unknown": 2.0})
        self.assertAlmostEqual(cfg.rho_bias_tau_s, 8000.0)
        self.assertAlmostEqual(cfg.rho_bias_sigma, 0.1)


if __name__ == "__main__":
    unittest.main()
