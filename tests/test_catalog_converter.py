import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


class TestCatalogConverter(unittest.TestCase):
    def test_convert_block_rephase_event_uses_scaled_duration(self):
        from scripts.catalog_to_case_config import convert_block

        block = {
            "scenario_id": "FORM_B4_REPHASE_DIFFERENTIAL_DRAG",
            "area": "FORM",
            "duration_s": 1000.0,
            "dt_control_s": 5.0,
            "environment": {"env_toggle": ["ENV_S2_STORM_PULSE"]},
            "control": {"control_toggle": "CTL_DD_REPHASE"},
            "uq": {"method": "MC", "n_mc": 10},
            "sensors": {"sensor_toggle": ["SENS_GNSS_RAW_DUAL"]},
        }
        out = convert_block(
            block,
            default_particles=16,
            duration_scale=0.1,
            start_utc=None,
            use_env_sources=False,
            dt_floor=1.0,
        )
        self.assertEqual(out["type"], "formation")
        self.assertAlmostEqual(out["duration_s"], 100.0)
        self.assertAlmostEqual(out["events"][0]["t_s"], 50.0)
        self.assertEqual(out["catalog_scenario_id"], "FORM_B4_REPHASE_DIFFERENTIAL_DRAG")
        self.assertEqual(out["catalog_env_toggles"], ["ENV_S2_STORM_PULSE"])
        self.assertEqual(out["catalog_sensor_toggles"], ["SENS_GNSS_RAW_DUAL"])
        self.assertEqual(out["catalog_est_toggles"], [])
        self.assertEqual(out["catalog_control_toggles"], ["CTL_DD_REPHASE"])

    def test_summary_metadata_has_required_keys(self):
        from scripts.summary_metadata import build_run_metadata

        scenario = {
            "name": "att_a1_detumble_aero_only",
            "catalog_scenario_id": "ATT_A1_DETUMBLE_AERO_ONLY",
            "catalog_env_toggles": ["ENV_S1_OU_LOGRHO"],
            "catalog_gsi_toggles": ["GSI_M1_MAXWELL"],
            "catalog_sensor_toggles": ["SENS_ATT_GYRO"],
            "catalog_ground_toggles": ["GRD_REG_EU"],
            "catalog_est_toggles": ["EST_EKF"],
            "catalog_control_toggles": ["CTL_ATT_AERO_RATE"],
        }
        meta = build_run_metadata(scenario, seed=42)
        self.assertIn("git_hash", meta)
        self.assertIn("scenario_id", meta)
        self.assertIn("seed", meta)
        self.assertIn("toggles", meta)
        self.assertEqual(meta["scenario_id"], "ATT_A1_DETUMBLE_AERO_ONLY")
        self.assertEqual(meta["seed"], 42)
        self.assertEqual(meta["toggles"]["env"], ["ENV_S1_OU_LOGRHO"])
        self.assertEqual(meta["toggles"]["ground"], ["GRD_REG_EU"])

    def test_convert_block_maps_c8_c9_measurement_and_weather_controls(self):
        from scripts.catalog_to_case_config import convert_block

        c8 = {
            "scenario_id": "OD_C8_STATION_WEATHER_AVAILABILITY",
            "area": "OD",
            "duration_s": 1000.0,
            "dt_meas_s": 5.0,
            "environment": {"env_toggle": ["ENV_S1_OU_LOGRHO"]},
            "sensors": {"sensor_toggle": ["SENS_SLR"]},
            "ground": {"ground_toggle": ["GRD_GLB_20", "GRD_WEATHER_ON"]},
            "uq": {"method": "MC", "n_mc": 10},
        }
        c9 = {
            "scenario_id": "OD_C9_MEAS_ERROR_STRESS_TESTS",
            "area": "OD",
            "duration_s": 1000.0,
            "dt_meas_s": 1.0,
            "environment": {"env_toggle": ["ENV_S1_OU_LOGRHO"]},
            "sensors": {"sensor_toggle": ["SENS_GNSS_RAW_DUAL", "SENS_SLR"]},
            "ground": {"ground_toggle": ["GRD_REG_EU", "GRD_WEATHER_ON"]},
            "uq": {"method": "MC", "n_mc": 10},
        }
        out_c8 = convert_block(
            c8,
            default_particles=16,
            duration_scale=1.0,
            start_utc=None,
            use_env_sources=False,
            dt_floor=1.0,
        )
        out_c9 = convert_block(
            c9,
            default_particles=16,
            duration_scale=1.0,
            start_utc=None,
            use_env_sources=False,
            dt_floor=1.0,
        )
        self.assertTrue(out_c8["pod_slr_weather_on"])
        self.assertIn("GRD_WEATHER_ON", out_c8["catalog_ground_toggles"])
        self.assertGreater(out_c9["pod_gnss_code_outlier_prob"], 0.0)
        self.assertGreater(out_c9["pod_gnss_cycle_slip_prob_per_min"], 0.0)


if __name__ == "__main__":
    unittest.main()
