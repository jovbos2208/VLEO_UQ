import math
import unittest

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    apply_measurement_latency,
    GroundStation,
    PropagatorConfig,
    StmPropagator,
    VehicleParams,
    default_geometry,
    run_pod_uq_measurements,
    run_pod_uq_measurements_enkf,
    simulate_gnss_measurements,
    simulate_slr_measurements,
    slr_availability_mask,
    summarize_arc_metrics,
)


class TestPodMeasurements(unittest.TestCase):
    def test_measurement_latency_index_mapping(self):
        t_grid = np.arange(0.0, 110.0, 10.0, dtype=float)
        idx = np.array([1, 3, 5], dtype=int)
        eff = apply_measurement_latency(t_grid, idx, latency_s=15.0, jitter_s=0.0, rng=np.random.default_rng(1))
        np.testing.assert_array_equal(eff, np.array([0, 1, 3], dtype=int))

    def test_gnss_outlier_and_cycle_slip_flags(self):
        t_grid = np.arange(0.0, 60.0, 10.0, dtype=float)
        truth_states = np.zeros((len(t_grid), 16), dtype=float)
        truth_states[:, 0] = 7000e3

        sat_positions_eci = np.zeros((len(t_grid), 1, 3), dtype=float)
        sat_positions_eci[:, 0, 0] = 26560e3

        gnss = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            frequencies_hz=(1575.42e6,),
            code_sigma_m=0.5,
            carrier_sigma_m=0.01,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            dropout_prob=0.0,
            cadence_s=None,
            cycle_slip_gap_s=0.0,
            cycle_slip_prob_per_min=1.0,
            code_outlier_prob=1.0,
            carrier_outlier_prob=1.0,
            rng=np.random.default_rng(7),
        )

        self.assertEqual(gnss.values.shape, gnss.is_outlier.shape)
        self.assertEqual(gnss.values.shape, gnss.is_cycle_slip.shape)
        self.assertTrue(np.all(gnss.is_outlier))

        carrier_mask = gnss.ambiguity_ids >= 0
        self.assertTrue(np.any(carrier_mask))
        self.assertTrue(np.all(~gnss.is_cycle_slip[~carrier_mask]))
        self.assertTrue(np.any(gnss.is_cycle_slip[carrier_mask]))
        self.assertGreater(int(gnss.ambiguities_cycles.shape[0]), 1)

    def test_ops_outage_reduces_gnss_measurements(self):
        t_grid = np.arange(0.0, 600.0, 10.0, dtype=float)
        truth_states = np.zeros((len(t_grid), 16), dtype=float)
        truth_states[:, 0] = 7000e3
        sat_positions_eci = np.zeros((len(t_grid), 2, 3), dtype=float)
        sat_positions_eci[:, 0, 0] = 26560e3
        sat_positions_eci[:, 1, 1] = 26560e3

        gnss_nom = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            frequencies_hz=(1575.42e6,),
            code_sigma_m=0.2,
            carrier_sigma_m=0.01,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            dropout_prob=0.0,
            cadence_s=10.0,
            include_carrier=False,
            rng=np.random.default_rng(2),
        )
        gnss_out = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            frequencies_hz=(1575.42e6,),
            code_sigma_m=0.2,
            carrier_sigma_m=0.01,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            dropout_prob=0.0,
            cadence_s=10.0,
            include_carrier=False,
            ops_outage_on=True,
            ops_outage_rate_per_hour=8.0,
            ops_outage_mean_duration_s=120.0,
            rng=np.random.default_rng(2),
        )
        self.assertLess(int(gnss_out.values.shape[0]), int(gnss_nom.values.shape[0]))

    def test_slr_weather_mask_reproducible(self):
        t_grid = np.arange(0.0, 600.0, 20.0, dtype=float)
        truth_states = np.zeros((len(t_grid), 16), dtype=float)
        stations = [
            GroundStation("S1", 0.0, 0.0, 0.0),
            GroundStation("S2", np.deg2rad(20.0), np.deg2rad(30.0), 0.0),
        ]

        base = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=-90.0,
            require_night=False,
            weather_enabled=False,
        )
        mask_a = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=-90.0,
            require_night=False,
            weather_enabled=True,
            weather_clear_prob=0.5,
            weather_p_stay_clear=0.9,
            weather_p_stay_blocked=0.8,
            weather_seed=1234,
        )
        mask_b = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=-90.0,
            require_night=False,
            weather_enabled=True,
            weather_clear_prob=0.5,
            weather_p_stay_clear=0.9,
            weather_p_stay_blocked=0.8,
            weather_seed=1234,
        )
        mask_c = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=-90.0,
            require_night=False,
            weather_enabled=True,
            weather_clear_prob=0.5,
            weather_p_stay_clear=0.9,
            weather_p_stay_blocked=0.8,
            weather_seed=999,
        )

        self.assertTrue(np.array_equal(mask_a, mask_b))
        self.assertFalse(np.array_equal(base, mask_a))
        self.assertFalse(np.array_equal(mask_a, mask_c))
        self.assertTrue(np.any(mask_a))
        self.assertTrue(np.any(~mask_a))

    def test_gnss_slr_od_reduces_error(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_stm = StmPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0_truth = np.zeros(prop_stm.state_size)
        x0_truth[0:3] = r0
        x0_truth[3:6] = v0
        x0_truth[6] = 1.0

        x0_guess = x0_truth.copy()
        x0_guess[0:3] += np.array([50.0, -20.0, 10.0])

        P0_guess = np.diag(np.ones(prop_stm.state_size))

        t_grid = np.linspace(0.0, 60.0, 7)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth_states = prop_det.propagate(x0_truth, t_grid, env)

        sat_alt = 26560e3
        sats = np.array(
            [
                [sat_alt, 0.0, 0.0],
                [0.0, sat_alt, 0.0],
                [0.0, 0.0, sat_alt],
                [-sat_alt, 0.0, 0.0],
            ],
            dtype=float,
        )
        sat_positions_eci = np.tile(sats[None, :, :], (len(t_grid), 1, 1))

        gnss = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            code_sigma_m=0.5,
            carrier_sigma_m=0.02,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            frequencies_hz=(1575.42e6,),
            dropout_prob=0.0,
            cycle_slip_gap_s=None,
            cycle_slip_prob_per_min=0.0,
            include_carrier=False,
        )

        stations = [GroundStation("STA", 0.0, 0.0, 0.0)]
        slr = simulate_slr_measurements(
            truth_states,
            t_grid,
            stations,
            sigma_m=0.05,
            min_elevation_deg=0.0,
        )

        result = run_pod_uq_measurements(
            prop_det,
            prop_stm,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            gnss=gnss,
            slr=slr,
            max_iter=5,
            estimate_clock_bias=False,
            estimate_ambiguity=False,
        )

        err0 = np.linalg.norm(x0_guess[0:3] - x0_truth[0:3])
        err1 = np.linalg.norm(result.x0_est[0:3] - x0_truth[0:3])
        self.assertLess(err1, 0.3 * err0)

        summary = summarize_arc_metrics(result, horizons_s=(10.0, 30.0))
        self.assertGreater(len(summary["horizons_s"]), 0)
        self.assertIn("measurement_residuals", summary)
        self.assertIn("gnss_all", summary["measurement_residuals"])
        self.assertGreater(summary["measurement_residuals"]["gnss_all"]["count"], 0)

    def test_gnss_slr_od_enkf_reduces_error(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0

        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_stm = StmPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])

        x0_truth = np.zeros(prop_stm.state_size)
        x0_truth[0:3] = r0
        x0_truth[3:6] = v0
        x0_truth[6] = 1.0

        x0_guess = x0_truth.copy()
        x0_guess[0:3] += np.array([60.0, -25.0, 12.0])

        P0_guess = np.diag(np.ones(prop_stm.state_size))

        t_grid = np.linspace(0.0, 60.0, 7)
        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth_states = prop_det.propagate(x0_truth, t_grid, env)

        sat_alt = 26560e3
        sats = np.array(
            [
                [sat_alt, 0.0, 0.0],
                [0.0, sat_alt, 0.0],
                [0.0, 0.0, sat_alt],
                [-sat_alt, 0.0, 0.0],
            ],
            dtype=float,
        )
        sat_positions_eci = np.tile(sats[None, :, :], (len(t_grid), 1, 1))

        gnss = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            code_sigma_m=0.5,
            carrier_sigma_m=0.02,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            frequencies_hz=(1575.42e6,),
            dropout_prob=0.0,
            cycle_slip_gap_s=None,
            cycle_slip_prob_per_min=0.0,
            include_carrier=False,
        )

        stations = [GroundStation("STA", 0.0, 0.0, 0.0)]
        slr = simulate_slr_measurements(
            truth_states,
            t_grid,
            stations,
            sigma_m=0.05,
            min_elevation_deg=0.0,
        )

        result = run_pod_uq_measurements_enkf(
            prop_det,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            gnss=gnss,
            slr=slr,
            members=48,
            seed=11,
            estimate_clock_bias=False,
            estimate_clock_drift=False,
            estimate_tropo=False,
            estimate_ambiguity=False,
            estimate_slr_bias=False,
            use_carrier=False,
            inflation=1.0,
        )

        err0 = np.linalg.norm(x0_guess[0:3] - x0_truth[0:3])
        err1 = np.linalg.norm(result.x0_est[0:3] - x0_truth[0:3])
        self.assertLess(err1, 0.8 * err0)
        self.assertTrue(np.all(np.isfinite(result.P0_post)))

    def test_batch_od_accepts_latency_controls(self):
        geom = default_geometry()
        aero = AeroAdapter()
        aero.init(geom)

        vehicle = VehicleParams()
        config = PropagatorConfig()
        config.rho_fast_sigma = 0.0
        config.rho_bias_sigma = 0.0
        config.wind_sigma = 0.0
        prop_det = DeterministicPropagator(aero, vehicle, config)
        prop_stm = StmPropagator(aero, vehicle, config)

        mu = config.mu_earth_m3_s2
        r0 = np.array([7000e3, 0.0, 0.0])
        v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0])
        x0_truth = np.zeros(prop_stm.state_size)
        x0_truth[0:3] = r0
        x0_truth[3:6] = v0
        x0_truth[6] = 1.0
        x0_guess = x0_truth.copy()
        x0_guess[0:3] += np.array([30.0, -15.0, 8.0])
        P0_guess = np.diag(np.ones(prop_stm.state_size))
        t_grid = np.linspace(0.0, 60.0, 7)

        env = []
        for _ in t_grid:
            e = EnvInputs()
            e.density = 0.0
            e.temperature_K = 1000.0
            e.particles_mass_kg = 28.0 * 1.6605390689252e-27
            e.wind_I = np.zeros(3)
            env.append(e)

        truth_states = prop_det.propagate(x0_truth, t_grid, env)
        sat_positions_eci = np.zeros((len(t_grid), 1, 3), dtype=float)
        sat_positions_eci[:, 0, 0] = 26560e3
        gnss = simulate_gnss_measurements(
            truth_states,
            t_grid,
            sat_positions_eci,
            frequencies_hz=(1575.42e6,),
            code_sigma_m=0.3,
            carrier_sigma_m=0.01,
            clock_bias_sigma_rw_m=0.0,
            clock_drift_sigma_rw_m_s=0.0,
            tropo_sigma_rw_m=0.0,
            require_los=False,
            dropout_prob=0.0,
            include_carrier=False,
            cadence_s=10.0,
            rng=np.random.default_rng(41),
        )
        result = run_pod_uq_measurements(
            prop_det,
            prop_stm,
            x0_truth,
            x0_guess,
            P0_guess,
            t_grid,
            env,
            gnss=gnss,
            slr=None,
            max_iter=4,
            estimate_clock_bias=False,
            estimate_ambiguity=False,
            measurement_latency_s=5.0,
            measurement_latency_jitter_s=1.0,
            measurement_latency_seed=9,
        )
        self.assertIn("measurement_latency_s", result.nuisance)
        self.assertIn("gnss_time_effective", result.nuisance)
        self.assertTrue(np.all(np.asarray(result.nuisance["gnss_time_effective"]) >= 0))


if __name__ == "__main__":
    unittest.main()
