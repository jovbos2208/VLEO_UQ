import unittest

import numpy as np

from scripts.model_discrepancy import apply_density_model_discrepancy


class _Env:
    def __init__(self, density: float, particles_mass_kg: float = 28.0 * 1.6605390689252e-27):
        self.density = float(density)
        self.particles_mass_kg = float(particles_mass_kg)
        self.eta1_rad = 0.95
        self.eta2_rad = 0.95
        self.srp_scale = 1.0
        self.magnetic_field_I_T = np.array([2.1e-5, -1.2e-5, 3.3e-5], dtype=float)
        self.temperature_ratio_method = 1


class TestModelDiscrepancy(unittest.TestCase):
    def test_disabled_no_change(self):
        env = [_Env(1.0e-12) for _ in range(5)]
        t_grid = np.linspace(0.0, 40.0, 5)
        before = np.array([e.density for e in env], dtype=float)
        meta = apply_density_model_discrepancy(
            env,
            t_grid,
            scenario={"name": "att_a1", "model_discrepancy_on": False, "model_discrepancy_sigma_rel": 0.2},
            seed=42,
        )
        after = np.array([e.density for e in env], dtype=float)
        np.testing.assert_allclose(after, before)
        self.assertFalse(bool(meta.get("applied")))

    def test_reproducible_when_enabled(self):
        t_grid = np.linspace(0.0, 120.0, 13)
        scenario = {
            "name": "att_a1",
            "model_discrepancy_on": True,
            "model_discrepancy_sigma_rel": 0.25,
            "model_discrepancy_tau_s": 60.0,
            "model_discrepancy_df": 4.0,
        }
        env_a = [_Env(1.0e-12) for _ in range(t_grid.size)]
        env_b = [_Env(1.0e-12) for _ in range(t_grid.size)]
        meta_a = apply_density_model_discrepancy(env_a, t_grid, scenario=scenario, seed=77)
        meta_b = apply_density_model_discrepancy(env_b, t_grid, scenario=scenario, seed=77)
        dens_a = np.array([e.density for e in env_a], dtype=float)
        dens_b = np.array([e.density for e in env_b], dtype=float)
        np.testing.assert_allclose(dens_a, dens_b)
        self.assertTrue(bool(meta_a.get("applied")))
        self.assertEqual(int(meta_a["seed"]), int(meta_b["seed"]))

    def test_composition_channel_reproducible(self):
        t_grid = np.linspace(0.0, 80.0, 9)
        scenario = {
            "name": "od_c1",
            "model_discrepancy_on": False,
            "composition_discrepancy_on": True,
            "composition_sigma_rel": 0.2,
            "composition_tau_s": 50.0,
            "composition_df": 4.0,
        }
        env_a = [_Env(1.0e-12) for _ in range(t_grid.size)]
        env_b = [_Env(1.0e-12) for _ in range(t_grid.size)]
        meta_a = apply_density_model_discrepancy(env_a, t_grid, scenario=scenario, seed=19)
        meta_b = apply_density_model_discrepancy(env_b, t_grid, scenario=scenario, seed=19)
        comp_a = np.array([e.particles_mass_kg for e in env_a], dtype=float)
        comp_b = np.array([e.particles_mass_kg for e in env_b], dtype=float)
        np.testing.assert_allclose(comp_a, comp_b)
        self.assertTrue(bool(meta_a.get("applied_composition")))
        self.assertTrue(bool(meta_a.get("applied")))

    def test_storm_jump_channel_applies_piecewise_density(self):
        t_grid = np.linspace(0.0, 6.0 * 3600.0, 61)
        scenario = {
            "name": "od_c2",
            "model_discrepancy_on": False,
            "storm_jump_on": True,
            "storm_jump_rate_per_day": 24.0,
            "storm_jump_duration_s": 1800.0,
            "storm_jump_sigma_rel": 0.4,
            "storm_jump_mean_rel": 0.15,
            "storm_jump_df": 4.0,
        }
        env = [_Env(1.0e-12) for _ in range(t_grid.size)]
        meta = apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=7)
        dens = np.array([e.density for e in env], dtype=float)
        self.assertTrue(bool(meta.get("applied")))
        self.assertIn("storm_jump", meta)
        self.assertGreaterEqual(int(meta["storm_jump"].get("event_count", 0)), 0)
        self.assertGreater(float(np.std(dens)), 0.0)

    def test_ao_plasma_structural_surrogates_apply(self):
        t_grid = np.linspace(0.0, 3600.0, 13)
        scenario = {
            "name": "att_b3",
            "model_discrepancy_on": False,
            "ao_erosion_on": True,
            "ao_erosion_sigma_rel": 0.1,
            "ao_erosion_rate_rel_per_day": 0.2,
            "plasma_charging_on": True,
            "plasma_charging_sigma_rel": 0.2,
            "structural_flex_on": True,
            "structural_flex_sigma_rel": 0.05,
        }
        env = [_Env(1.0e-12) for _ in range(t_grid.size)]
        eta_before = np.array([[e.eta1_rad, e.eta2_rad] for e in env], dtype=float)
        srp_before = np.array([e.srp_scale for e in env], dtype=float)
        b_before = np.array([e.magnetic_field_I_T for e in env], dtype=float)

        meta = apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=23)
        eta_after = np.array([[e.eta1_rad, e.eta2_rad] for e in env], dtype=float)
        srp_after = np.array([e.srp_scale for e in env], dtype=float)
        b_after = np.array([e.magnetic_field_I_T for e in env], dtype=float)

        self.assertTrue(bool(meta.get("applied_ao_erosion")))
        self.assertTrue(bool(meta.get("applied_plasma_charging")))
        self.assertTrue(bool(meta.get("applied_structural_flex")))
        self.assertGreater(float(np.linalg.norm(eta_after - eta_before)), 0.0)
        self.assertGreater(float(np.linalg.norm(srp_after - srp_before)), 0.0)
        self.assertGreater(float(np.linalg.norm(b_after - b_before)), 0.0)

    def test_rarefied_regime_switch_applies(self):
        t_grid = np.linspace(0.0, 2.0 * 3600.0, 25)
        scenario = {
            "name": "att_c3",
            "model_discrepancy_on": False,
            "rarefied_regime_switch_on": True,
            "rarefied_switch_rate_per_hour": 1.0e4,
            "rarefied_switch_mean_duration_s": 300.0,
            "rarefied_switch_alt_temperature_ratio_method": 2,
        }
        env = [_Env(1.0e-12) for _ in range(t_grid.size)]
        meta = apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=77)
        modes = np.array([e.temperature_ratio_method for e in env], dtype=int)
        self.assertTrue(bool(meta.get("rarefied_regime_switch", {}).get("enabled")))
        self.assertGreater(int(meta.get("rarefied_regime_switch", {}).get("switched_count", 0)), 0)
        self.assertTrue(np.any(modes == 2))


if __name__ == "__main__":
    unittest.main()
