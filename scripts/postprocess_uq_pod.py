from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_case_studies import plot_case, save_json
from scripts.model_discrepancy import apply_density_model_discrepancy
from scripts.payload_impact import compute_payload_impact_proxy_from_stats
from scripts.summary_metadata import build_run_metadata
import gc


def load_mc_stats_from_stats(path: Path):
    with np.load(path) as data:
        return data["t_grid"], None, data["mc_mean"], data["mc_std"], data["mc_cov_final"]


def load_mc_stats_from_shards(indir: Path):
    paths = sorted(indir.glob("ensemble_rank*.npz"))
    if not paths:
        raise RuntimeError(f"no shard files found in {indir}")
    t_grid = None
    sum_x = None
    sum_x2 = None
    n_total = 0
    sum_final = None
    sum_xx_final = None
    for path in paths:
        with np.load(path) as data:
            X = data["X"]
            t_local = data["t_grid"]
        if t_grid is None:
            t_grid = t_local
            nt, _, nx = X.shape
            sum_x = np.zeros((nt, nx), dtype=float)
            sum_x2 = np.zeros((nt, nx), dtype=float)
            sum_final = np.zeros(nx, dtype=float)
            sum_xx_final = np.zeros((nx, nx), dtype=float)
        sum_x += X.sum(axis=1)
        sum_x2 += (X ** 2).sum(axis=1)
        n_local = X.shape[1]
        n_total += n_local
        Xf = X[-1]
        sum_final += Xf.sum(axis=0)
        sum_xx_final += Xf.T @ Xf
    if n_total <= 1:
        raise RuntimeError("not enough MC samples to compute statistics")
    mc_mean = sum_x / n_total
    mc_var = (sum_x2 - (sum_x ** 2) / n_total) / (n_total - 1)
    mc_std = np.sqrt(np.maximum(mc_var, 0.0))
    mc_cov_final = (sum_xx_final - np.outer(sum_final, sum_final) / n_total) / (n_total - 1)
    return t_grid, None, mc_mean, mc_std, mc_cov_final


def load_ut_stats(path: Path):
    with np.load(path) as data:
        return data["t_grid"], data["ut_mean"], data["ut_cov"]


def load_det_stm(path: Path):
    with np.load(path) as data:
        return data["t_grid"], data["det_states"], data["stm_mean"], data["stm_cov"]


def run_pod(
    det_states,
    t_grid,
    pod_sp3,
    pod_eop,
    pod_use_sat_clock,
    pod_arc_s,
    pod_overlap_s,
    pod_skip_slr,
    seed,
    t0_utc=None,
    use_env_sources=False,
    omni_path="data/space_weather/omni/omni2_all_years.dat",
    hwm14_lib="data/hwm14/libhwm14.so",
    hwm14_data="data/hwm14",
    env_interpolation="nearest",
    freeze_attitude=False,
    pod_estimator: str = "batch",
    pod_enkf_members: int = 64,
    pod_enkf_seed: int | None = None,
    pod_enkf_inflation: float = 1.0,
    pod_enkf_use_carrier: bool = False,
    pod_gnss_cycle_slip_gap_s: float | None = None,
    pod_gnss_cycle_slip_prob_per_min: float = 0.0,
    pod_gnss_code_outlier_prob: float = 0.0,
    pod_gnss_code_outlier_sigma_scale: float = 25.0,
    pod_gnss_carrier_outlier_prob: float = 0.0,
    pod_gnss_carrier_outlier_sigma_scale: float = 25.0,
    pod_slr_weather_on: bool = False,
    pod_slr_weather_clear_prob: float = 0.8,
    pod_slr_weather_p_stay_clear: float = 0.985,
    pod_slr_weather_p_stay_blocked: float = 0.93,
    pod_slr_weather_seed: int | None = None,
    spacecraft_mass_kg: float | None = None,
    spacecraft_inertia_kgm2: list[float] | None = None,
    scenario: dict | None = None,
):
    from vleo_uq import (
        default_ilrs_stations,
        parse_eop_all,
        parse_sp3,
        run_pod_uq_measurements,
        run_pod_uq_measurements_enkf,
        simulate_gnss_measurements,
        simulate_slr_measurements,
        slr_availability_mask,
        slice_gnss_measurements,
        slice_slr_measurements,
        summarize_arc_metrics,
    )
    from scripts.run_case_studies import synthetic_gnss_positions
    import datetime as dt
    rng = np.random.default_rng(seed)
    scenario_dict = scenario if isinstance(scenario, dict) else {}
    latency_s = float(scenario_dict.get("pod_measurement_latency_s", 0.0))
    latency_jitter_s = float(scenario_dict.get("pod_measurement_latency_jitter_s", 0.0))
    latency_seed = int(scenario_dict.get("pod_measurement_latency_seed", seed))
    ops_outage_on = bool(scenario_dict.get("pod_ops_outage_on", False))
    ops_outage_rate_per_hour = float(scenario_dict.get("pod_ops_outage_rate_per_hour", 0.0))
    ops_outage_mean_duration_s = float(scenario_dict.get("pod_ops_outage_mean_duration_s", 0.0))
    gnss_positions = None
    sat_clock_bias_m = None
    if pod_sp3:
        sp3_path = Path(pod_sp3)
        if sp3_path.exists():
            eph = parse_sp3(str(sp3_path), constellation_prefix=("G",))
            t0 = eph.epochs[0]
            t_grid_dt = [t0 + dt.timedelta(seconds=float(sec)) for sec in t_grid]
            eop = None
            if pod_eop:
                eop_path = Path(pod_eop)
                if eop_path.exists():
                    eop = parse_eop_all(str(eop_path))
            gnss_positions = eph.sample(t_grid_dt, frame="eci", eop=eop)
            if pod_use_sat_clock:
                clocks = eph.sample_clock(t_grid_dt)
                if clocks is not None:
                    sat_clock_bias_m = clocks * 299792458.0
    if gnss_positions is None:
        gnss_positions = synthetic_gnss_positions(t_grid, 3.986004418e14)
    dt_grid = float(t_grid[1] - t_grid[0]) if len(t_grid) > 1 else 1.0
    gnss_cadence_s = max(10.0, dt_grid)
    gnss = simulate_gnss_measurements(
        det_states,
        t_grid,
        gnss_positions,
        rng=rng,
        cadence_s=gnss_cadence_s,
        frequencies_hz=(1575.42e6, 1227.60e6),
        code_sigma_m=0.5,
        carrier_sigma_m=0.003,
        require_los=True,
        earth_margin_m=50e3,
        fov_half_angle_deg=80.0,
        dropout_prob=0.01,
        cycle_slip_gap_s=pod_gnss_cycle_slip_gap_s,
        cycle_slip_prob_per_min=float(pod_gnss_cycle_slip_prob_per_min),
        code_outlier_prob=float(pod_gnss_code_outlier_prob),
        code_outlier_sigma_scale=float(pod_gnss_code_outlier_sigma_scale),
        carrier_outlier_prob=float(pod_gnss_carrier_outlier_prob),
        carrier_outlier_sigma_scale=float(pod_gnss_carrier_outlier_sigma_scale),
        sat_clock_bias_m=sat_clock_bias_m,
        ops_outage_on=ops_outage_on,
        ops_outage_rate_per_hour=ops_outage_rate_per_hour,
        ops_outage_mean_duration_s=ops_outage_mean_duration_s,
    )
    print(
        f"[pod] gnss generated cadence_s={gnss_cadence_s:.1f} meas={int(gnss.values.shape[0])} "
        f"carrier={int(np.sum(gnss.types == 1))} amb_total={int(gnss.ambiguities_cycles.shape[0])}",
        flush=True,
    )
    slr = None
    if not pod_skip_slr:
        stations = default_ilrs_stations()
        weather_seed = pod_slr_weather_seed
        if weather_seed is None:
            weather_seed = int(seed) + 7919
        slr_mask = slr_availability_mask(
            det_states,
            t_grid,
            stations,
            t0=t0_utc,
            weather_enabled=bool(pod_slr_weather_on),
            weather_clear_prob=float(pod_slr_weather_clear_prob),
            weather_p_stay_clear=float(pod_slr_weather_p_stay_clear),
            weather_p_stay_blocked=float(pod_slr_weather_p_stay_blocked),
            weather_seed=int(weather_seed),
        )
        slr = simulate_slr_measurements(
            det_states,
            t_grid,
            stations,
            rng=rng,
            sigma_m=0.01,
            availability_mask=slr_mask,
            ops_outage_on=ops_outage_on,
            ops_outage_rate_per_hour=ops_outage_rate_per_hour,
            ops_outage_mean_duration_s=ops_outage_mean_duration_s,
        )
    arcs = []
    step = max(1, int(round((pod_arc_s - pod_overlap_s) / (t_grid[1] - t_grid[0]))))
    arc_len = max(2, int(round(pod_arc_s / (t_grid[1] - t_grid[0]))))
    n = len(t_grid)
    start = 0
    while start < n - 1:
        end = min(n, start + arc_len)
        arcs.append((start, end))
        start += step
        if end == n:
            break
    # build propagators/env for POD
    from vleo_uq import (
        AeroAdapter,
        DeterministicPropagator,
        EnvInputs,
        PropagatorConfig,
        StmPropagator,
        VehicleParams,
        default_geometry,
    )
    from scripts.run_case_studies import build_env, build_env_from_sources

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()
    if spacecraft_mass_kg is not None:
        try:
            mass = float(spacecraft_mass_kg)
            if np.isfinite(mass) and mass > 0.0:
                vehicle.mass_kg = mass
        except Exception:
            pass
    if isinstance(spacecraft_inertia_kgm2, list):
        try:
            arr = np.asarray(spacecraft_inertia_kgm2, dtype=float)
        except Exception:
            arr = np.array([], dtype=float)
        if arr.size == 3:
            vehicle.inertia_B = np.diag(arr)
        elif arr.size == 9:
            vehicle.inertia_B = arr.reshape(3, 3)

    cfg = PropagatorConfig()
    cfg.rtol = 1e-6
    cfg.atol = 1e-6
    cfg.max_step_s = 120.0
    cfg.rng_seed = seed
    cfg.freeze_attitude = freeze_attitude
    prop_det = DeterministicPropagator(aero, vehicle, cfg)
    prop_stm = StmPropagator(aero, vehicle, cfg)

    x0_truth = det_states[0].copy()
    x0_guess = x0_truth.copy()
    x0_guess[0:3] += np.array([50.0, -20.0, 10.0])
    P0_guess = np.diag(np.ones(prop_stm.state_size))

    if use_env_sources:
        if t0_utc is None:
            raise ValueError("start_utc required for env sources in POD")
        env_series = build_env_from_sources(
            t_grid_s=t_grid,
            t0_utc=t0_utc,
            x0=x0_truth,
            prop_det=prop_det,
            EnvInputs=EnvInputs,
            eta1_rad=None,
            eta2_rad=None,
            space_weather_source="omni",
            omni_path=omni_path,
            hwm14_lib=hwm14_lib,
            hwm14_data=hwm14_data,
            interpolation=env_interpolation,
        )
        env = env_series.env_inputs
    else:
        env = build_env(
            t_grid,
            EnvInputs,
            density=1e-12,
            temperature_K=1000.0,
            particle_mass_kg=28.0 * 1.6605390689252e-27,
        )
    apply_density_model_discrepancy(env, t_grid, scenario=scenario, seed=int(seed))

    # Streaming arc processing: summarize and release each arc immediately
    # to keep memory bounded on long OD scenarios.
    guess_full = prop_det.propagate(x0_guess, t_grid, env)
    arc_summaries = []
    estimator = str(pod_estimator).strip().lower()
    if estimator not in {"batch", "enkf"}:
        raise ValueError("pod_estimator must be 'batch' or 'enkf'")

    prev_start = None
    prev_x0_est = None
    prev_P0_post = None
    for i, (s_idx, e_idx) in enumerate(arcs):
        gnss_slice = None
        slr_slice = None
        try:
            gnss_slice = slice_gnss_measurements(gnss, s_idx, e_idx)
        except ValueError:
            gnss_slice = None
        if slr is not None:
            try:
                slr_slice = slice_slr_measurements(slr, s_idx, e_idx)
            except ValueError:
                slr_slice = None
        if gnss_slice is None and slr_slice is None:
            arc_summaries.append(
                {
                    "arc_index": i,
                    "start_idx": int(s_idx),
                    "end_idx": int(e_idx),
                    "status": "skipped_no_measurements",
                }
            )
            continue

        t_arc = t_grid[s_idx:e_idx]
        env_arc = env[s_idx:e_idx]
        if i % 10 == 0:
            n_gnss = 0 if gnss_slice is None else int(gnss_slice.values.shape[0])
            n_slr = 0 if slr_slice is None else int(slr_slice.values.shape[0])
            n_amb = 0 if gnss_slice is None else int(gnss_slice.ambiguities_cycles.shape[0])
            print(
                f"[pod] arc={i+1}/{len(arcs)} idx=({s_idx},{e_idx}) "
                f"gnss={n_gnss} slr={n_slr} amb={n_amb}",
                flush=True,
            )
        x0_truth_arc = det_states[s_idx].copy()
        x0_guess_arc = guess_full[s_idx].copy()
        P0_arc = P0_guess.copy()
        if prev_x0_est is not None and prev_P0_post is not None and prev_start is not None:
            if s_idx == prev_start:
                x0_guess_arc = prev_x0_est.copy()
                P0_arc = prev_P0_post.copy()
            elif s_idx > prev_start:
                t_link = t_grid[prev_start : s_idx + 1]
                env_link = env[prev_start : s_idx + 1]
                mean_link, cov_link = prop_stm.propagate(prev_x0_est, prev_P0_post, t_link, env_link)
                x0_guess_arc = mean_link[-1].copy()
                P0_arc = cov_link[-1].copy()

        if estimator == "enkf":
            arc_res = run_pod_uq_measurements_enkf(
                prop_det=prop_det,
                x0_truth=x0_truth_arc,
                x0_guess=x0_guess_arc,
                P0_guess=P0_arc,
                t_grid=t_arc,
                env=env_arc,
                gnss=gnss_slice,
                slr=slr_slice,
                members=int(pod_enkf_members),
                seed=int((seed if pod_enkf_seed is None else pod_enkf_seed) + i),
                estimate_clock_bias=True,
                estimate_clock_drift=False,
                estimate_tropo=False,
                estimate_ambiguity=bool(pod_enkf_use_carrier),
                estimate_slr_bias=False,
                use_carrier=bool(pod_enkf_use_carrier),
                inflation=float(pod_enkf_inflation),
                measurement_latency_s=latency_s,
                measurement_latency_jitter_s=latency_jitter_s,
                measurement_latency_seed=latency_seed + i,
            )
        else:
            arc_res = run_pod_uq_measurements(
                prop_det,
                prop_stm,
                x0_truth_arc,
                x0_guess_arc,
                P0_arc,
                t_arc,
                env_arc,
                gnss=gnss_slice,
                slr=slr_slice,
                max_iter=6,
                measurement_latency_s=latency_s,
                measurement_latency_jitter_s=latency_jitter_s,
                measurement_latency_seed=latency_seed + i,
            )
        summary = summarize_arc_metrics(arc_res, t_grid=t_arc)
        summary["arc_index"] = i
        summary["start_idx"] = int(s_idx)
        summary["end_idx"] = int(e_idx)
        summary["estimator"] = estimator
        arc_summaries.append(summary)

        if estimator == "enkf":
            # Match batch coupling semantics (estimate at arc start) for overlap-safe chaining.
            prev_start = s_idx
            prev_x0_est = arc_res.x0_est.copy()
            prev_P0_post = arc_res.P0_post.copy()
        else:
            prev_start = s_idx
            prev_x0_est = arc_res.x0_est.copy()
            prev_P0_post = arc_res.P0_post.copy()
        del arc_res
        if (i + 1) % 10 == 0:
            gc.collect()

    return [], arc_summaries


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 2: postprocess UQ + POD from stored runs.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--pod_uq", action="store_true")
    parser.add_argument("--pod_arc_s", type=float, default=1800.0)
    parser.add_argument("--pod_overlap_s", type=float, default=600.0)
    parser.add_argument("--pod_skip_slr", action="store_true")
    parser.add_argument("--pod_sp3", default=None)
    parser.add_argument("--pod_eop", default=None)
    parser.add_argument("--pod_use_sat_clock", action="store_true")
    parser.add_argument("--pod_estimator", choices=("batch", "enkf"), default="batch")
    parser.add_argument("--pod_enkf_members", type=int, default=64)
    parser.add_argument("--pod_enkf_seed", type=int, default=None)
    parser.add_argument("--pod_enkf_inflation", type=float, default=1.0)
    parser.add_argument("--pod_enkf_use_carrier", action="store_true")
    parser.add_argument("--pod_gnss_cycle_slip_gap_s", type=float, default=None)
    parser.add_argument("--pod_gnss_cycle_slip_prob_per_min", type=float, default=0.0)
    parser.add_argument("--pod_gnss_code_outlier_prob", type=float, default=0.0)
    parser.add_argument("--pod_gnss_code_outlier_sigma_scale", type=float, default=25.0)
    parser.add_argument("--pod_gnss_carrier_outlier_prob", type=float, default=0.0)
    parser.add_argument("--pod_gnss_carrier_outlier_sigma_scale", type=float, default=25.0)
    parser.add_argument("--pod_slr_weather_on", action="store_true")
    parser.add_argument("--pod_slr_weather_clear_prob", type=float, default=0.8)
    parser.add_argument("--pod_slr_weather_p_stay_clear", type=float, default=0.985)
    parser.add_argument("--pod_slr_weather_p_stay_blocked", type=float, default=0.93)
    parser.add_argument("--pod_slr_weather_seed", type=int, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    args = parser.parse_args()
    t0_utc = None
    if args.start_utc:
        import datetime as dt

        t0_utc = dt.datetime.fromisoformat(args.start_utc)

    raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenarios = raw.get("scenarios", [])
    if not scenarios:
        raise ValueError("config must include scenarios")

    outdir = Path(args.outdir)
    for scenario in scenarios:
        name = scenario["name"]
        typ = scenario["type"]
        if typ == "formation":
            offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
            results = []
            for idx, _ in enumerate(offsets):
                obj = f"{name}/obj_{idx+1:02d}"
                det_path = outdir / obj / "det_stm.npz"
                mc_path = outdir / f"{obj}_mc_merged.npz"
                mc_stats = outdir / f"{obj}_mc_stats.npz"
                mc_shards = outdir / "mc_shards" / f"formation_obj_{idx+1:02d}"
                ut_path = outdir / f"{obj}_ut_merged.npz"
                t_grid, det_states, stm_mean, stm_cov = load_det_stm(det_path)
                if mc_stats.exists():
                    _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_stats(mc_stats)
                elif mc_shards.exists():
                    _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_shards(mc_shards)
                else:
                    raise RuntimeError(f"MC stats not found for {obj}; expected stats or shard files")
                _, ut_mean, ut_cov = load_ut_stats(ut_path)
                case_dir = outdir / obj
                payload = dict(
                    t_grid=t_grid,
                    det_states=det_states,
                    mc_mean=mc_mean,
                    mc_std=mc_std,
                    mc_cov_final=mc_cov_final,
                    ut_mean=ut_mean,
                    ut_cov=ut_cov,
                    stm_mean=stm_mean,
                    stm_cov=stm_cov,
                )
                np.savez(case_dir / "mc_ut_stm.npz", **payload)
                summary = {
                    "case": obj,
                    "particles": int(scenario.get("particles", 0)),
                    "ut_mean_err": float(np.linalg.norm(mc_mean[-1] - ut_mean[-1])),
                    "stm_mean_err": float(np.linalg.norm(mc_mean[-1] - stm_mean[-1])),
                    "ut_cov_rel_err": float(np.linalg.norm(mc_cov_final - ut_cov[-1]) / max(np.linalg.norm(mc_cov_final), 1e-12)),
                    "stm_cov_rel_err": float(np.linalg.norm(mc_cov_final - stm_cov[-1]) / max(np.linalg.norm(mc_cov_final), 1e-12)),
                    "uq_parameter_draw": dict(scenario.get("uq_parameter_draw", {}))
                    if isinstance(scenario.get("uq_parameter_draw"), dict)
                    else {},
                    "payload_impact": compute_payload_impact_proxy_from_stats(mc_mean=mc_mean, mc_std=mc_std),
                    "metadata": build_run_metadata(scenario, args.seed + idx),
                }
                save_json(case_dir / "summary.json", summary)
                if args.plot:
                    plot_case(
                        outdir,
                        obj,
                        t_grid,
                        det_states,
                        mc_mean,
                        mc_std,
                        ut_mean,
                        stm_mean,
                        False,
                    )
                if args.pod_uq:
                    arc_results, pod_summary = run_pod(
                        det_states,
                        t_grid,
                        scenario.get("pod_sp3", args.pod_sp3),
                        scenario.get("pod_eop", args.pod_eop),
                        bool(scenario.get("pod_use_sat_clock", args.pod_use_sat_clock)),
                        float(scenario.get("pod_arc_s", args.pod_arc_s)),
                        float(scenario.get("pod_overlap_s", args.pod_overlap_s)),
                        bool(scenario.get("pod_skip_slr", args.pod_skip_slr)),
                        args.seed + idx,
                        t0_utc=t0_utc,
                        use_env_sources=args.use_env_sources,
                        omni_path=args.omni_path,
                        hwm14_lib=args.hwm14_lib,
                        hwm14_data=args.hwm14_data,
                        env_interpolation=args.env_interpolation,
                        freeze_attitude=True,
                        pod_estimator=str(scenario.get("pod_estimator", args.pod_estimator)),
                        pod_enkf_members=int(scenario.get("pod_enkf_members", args.pod_enkf_members)),
                        pod_enkf_seed=scenario.get("pod_enkf_seed", args.pod_enkf_seed),
                        pod_enkf_inflation=float(scenario.get("pod_enkf_inflation", args.pod_enkf_inflation)),
                        pod_enkf_use_carrier=bool(
                            scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)
                        ),
                        pod_gnss_cycle_slip_gap_s=scenario.get(
                            "pod_gnss_cycle_slip_gap_s",
                            args.pod_gnss_cycle_slip_gap_s,
                        ),
                        pod_gnss_cycle_slip_prob_per_min=float(
                            scenario.get(
                                "pod_gnss_cycle_slip_prob_per_min",
                                args.pod_gnss_cycle_slip_prob_per_min,
                            )
                        ),
                        pod_gnss_code_outlier_prob=float(
                            scenario.get("pod_gnss_code_outlier_prob", args.pod_gnss_code_outlier_prob)
                        ),
                        pod_gnss_code_outlier_sigma_scale=float(
                            scenario.get(
                                "pod_gnss_code_outlier_sigma_scale",
                                args.pod_gnss_code_outlier_sigma_scale,
                            )
                        ),
                        pod_gnss_carrier_outlier_prob=float(
                            scenario.get("pod_gnss_carrier_outlier_prob", args.pod_gnss_carrier_outlier_prob)
                        ),
                        pod_gnss_carrier_outlier_sigma_scale=float(
                            scenario.get(
                                "pod_gnss_carrier_outlier_sigma_scale",
                                args.pod_gnss_carrier_outlier_sigma_scale,
                            )
                        ),
                        pod_slr_weather_on=bool(scenario.get("pod_slr_weather_on", args.pod_slr_weather_on)),
                        pod_slr_weather_clear_prob=float(
                            scenario.get("pod_slr_weather_clear_prob", args.pod_slr_weather_clear_prob)
                        ),
                        pod_slr_weather_p_stay_clear=float(
                            scenario.get("pod_slr_weather_p_stay_clear", args.pod_slr_weather_p_stay_clear)
                        ),
                        pod_slr_weather_p_stay_blocked=float(
                            scenario.get("pod_slr_weather_p_stay_blocked", args.pod_slr_weather_p_stay_blocked)
                        ),
                        pod_slr_weather_seed=scenario.get("pod_slr_weather_seed", args.pod_slr_weather_seed),
                        scenario=scenario,
                    )
                    pod_dir = case_dir / "pod"
                    pod_dir.mkdir(parents=True, exist_ok=True)
                    save_json(pod_dir / "summary.json", pod_summary)
                results.append((t_grid, det_states, mc_mean))
            if len(results) > 1:
                t_grid = results[0][0]
                chief_det = results[0][1]
                chief_mc = results[0][2]
                det_sep = []
                mc_sep = []
                for (_, det, mc) in results[1:]:
                    det_sep.append(np.linalg.norm(det[:, 0:3] - chief_det[:, 0:3], axis=1))
                    mc_sep.append(np.linalg.norm(mc[:, 0:3] - chief_mc[:, 0:3], axis=1))
                det_sep = np.vstack(det_sep)
                mc_sep = np.vstack(mc_sep)
                np.savez(
                    outdir / name / "formation_metrics.npz",
                    t_grid=t_grid,
                    det_separation_m=det_sep,
                    mc_mean_separation_m=mc_sep,
                )
        else:
            det_path = outdir / name / "det_stm.npz"
            mc_path = outdir / f"{name}_mc_merged.npz"
            mc_stats = outdir / f"{name}_mc_stats.npz"
            mc_shards = outdir / "mc_shards" / name
            ut_path = outdir / f"{name}_ut_merged.npz"
            t_grid, det_states, stm_mean, stm_cov = load_det_stm(det_path)
            if mc_stats.exists():
                _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_stats(mc_stats)
            elif mc_shards.exists():
                _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_shards(mc_shards)
            else:
                raise RuntimeError(f"MC stats not found for {name}; expected stats or shard files")
            _, ut_mean, ut_cov = load_ut_stats(ut_path)
            case_dir = outdir / name
            payload = dict(
                t_grid=t_grid,
                det_states=det_states,
                mc_mean=mc_mean,
                mc_std=mc_std,
                mc_cov_final=mc_cov_final,
                ut_mean=ut_mean,
                ut_cov=ut_cov,
                stm_mean=stm_mean,
                stm_cov=stm_cov,
            )
            np.savez(case_dir / "mc_ut_stm.npz", **payload)
            summary = {
                "case": name,
                "particles": int(scenario.get("particles", 0)),
                "ut_mean_err": float(np.linalg.norm(mc_mean[-1] - ut_mean[-1])),
                "stm_mean_err": float(np.linalg.norm(mc_mean[-1] - stm_mean[-1])),
                "ut_cov_rel_err": float(np.linalg.norm(mc_cov_final - ut_cov[-1]) / max(np.linalg.norm(mc_cov_final), 1e-12)),
                "stm_cov_rel_err": float(np.linalg.norm(mc_cov_final - stm_cov[-1]) / max(np.linalg.norm(mc_cov_final), 1e-12)),
                "uq_parameter_draw": dict(scenario.get("uq_parameter_draw", {}))
                if isinstance(scenario.get("uq_parameter_draw"), dict)
                else {},
                "payload_impact": compute_payload_impact_proxy_from_stats(mc_mean=mc_mean, mc_std=mc_std),
                "metadata": build_run_metadata(scenario, args.seed),
            }
            save_json(case_dir / "summary.json", summary)
            if args.plot:
                plot_case(
                    outdir,
                    name,
                    t_grid,
                    det_states,
                    mc_mean,
                    mc_std,
                    ut_mean,
                    stm_mean,
                    False,
                )
            if args.pod_uq:
                arc_results, pod_summary = run_pod(
                    det_states,
                    t_grid,
                    scenario.get("pod_sp3", args.pod_sp3),
                    scenario.get("pod_eop", args.pod_eop),
                    bool(scenario.get("pod_use_sat_clock", args.pod_use_sat_clock)),
                    float(scenario.get("pod_arc_s", args.pod_arc_s)),
                    float(scenario.get("pod_overlap_s", args.pod_overlap_s)),
                    bool(scenario.get("pod_skip_slr", args.pod_skip_slr)),
                    args.seed,
                    t0_utc=t0_utc,
                    use_env_sources=args.use_env_sources,
                    omni_path=args.omni_path,
                    hwm14_lib=args.hwm14_lib,
                    hwm14_data=args.hwm14_data,
                    env_interpolation=args.env_interpolation,
                    freeze_attitude=(typ != "attitude"),
                    pod_estimator=str(scenario.get("pod_estimator", args.pod_estimator)),
                    pod_enkf_members=int(scenario.get("pod_enkf_members", args.pod_enkf_members)),
                    pod_enkf_seed=scenario.get("pod_enkf_seed", args.pod_enkf_seed),
                    pod_enkf_inflation=float(scenario.get("pod_enkf_inflation", args.pod_enkf_inflation)),
                    pod_enkf_use_carrier=bool(scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)),
                    pod_gnss_cycle_slip_gap_s=scenario.get(
                        "pod_gnss_cycle_slip_gap_s",
                        args.pod_gnss_cycle_slip_gap_s,
                    ),
                    pod_gnss_cycle_slip_prob_per_min=float(
                        scenario.get(
                            "pod_gnss_cycle_slip_prob_per_min",
                            args.pod_gnss_cycle_slip_prob_per_min,
                        )
                    ),
                    pod_gnss_code_outlier_prob=float(
                        scenario.get("pod_gnss_code_outlier_prob", args.pod_gnss_code_outlier_prob)
                    ),
                    pod_gnss_code_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_code_outlier_sigma_scale",
                            args.pod_gnss_code_outlier_sigma_scale,
                        )
                    ),
                    pod_gnss_carrier_outlier_prob=float(
                        scenario.get("pod_gnss_carrier_outlier_prob", args.pod_gnss_carrier_outlier_prob)
                    ),
                    pod_gnss_carrier_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_carrier_outlier_sigma_scale",
                            args.pod_gnss_carrier_outlier_sigma_scale,
                        )
                    ),
                    pod_slr_weather_on=bool(scenario.get("pod_slr_weather_on", args.pod_slr_weather_on)),
                    pod_slr_weather_clear_prob=float(
                        scenario.get("pod_slr_weather_clear_prob", args.pod_slr_weather_clear_prob)
                    ),
                    pod_slr_weather_p_stay_clear=float(
                        scenario.get("pod_slr_weather_p_stay_clear", args.pod_slr_weather_p_stay_clear)
                    ),
                    pod_slr_weather_p_stay_blocked=float(
                        scenario.get("pod_slr_weather_p_stay_blocked", args.pod_slr_weather_p_stay_blocked)
                    ),
                    pod_slr_weather_seed=scenario.get("pod_slr_weather_seed", args.pod_slr_weather_seed),
                    scenario=scenario,
                )
                pod_dir = case_dir / "pod"
                pod_dir.mkdir(parents=True, exist_ok=True)
                save_json(pod_dir / "summary.json", pod_summary)
        gc.collect()


if __name__ == "__main__":
    main()
