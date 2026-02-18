from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.postprocess_uq_pod import (
    load_det_stm,
    load_mc_stats_from_shards,
    load_mc_stats_from_stats,
    load_ut_stats,
    run_pod,
)
from scripts.run_case_studies import plot_case, save_json


def main() -> None:
    parser = argparse.ArgumentParser(description="Postprocess one scenario (UQ + POD).")
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--plot_all_states", action="store_true")
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
    parser.add_argument(
        "--pod_force",
        action="store_true",
        help="Force POD even if scenario config has pod_uq=false.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    args = parser.parse_args()

    raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenarios = raw.get("scenarios", [])
    if not scenarios:
        raise ValueError("config must include scenarios")

    t0_utc = None
    if args.start_utc:
        import datetime as dt

        t0_utc = dt.datetime.fromisoformat(args.start_utc)

    outdir = Path(args.outdir)
    scenario = next((s for s in scenarios if s.get("name") == args.scenario), None)
    if scenario is None:
        raise ValueError(f"scenario '{args.scenario}' not found")
    do_pod = bool(args.pod_uq) and (bool(scenario.get("pod_uq", False)) or bool(args.pod_force))
    print(
        f"[post] scenario={args.scenario} do_pod={do_pod} "
        f"(arg_pod_uq={args.pod_uq}, scenario_pod_uq={bool(scenario.get('pod_uq', False))}, pod_force={args.pod_force})",
        flush=True,
    )

    name = scenario["name"]
    typ = scenario["type"]
    if typ == "formation":
        offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
        results = []
        for idx, _ in enumerate(offsets):
            obj = f"{name}/obj_{idx+1:02d}"
            det_path = outdir / obj / "det_stm.npz"
            mc_stats = outdir / f"{obj}_mc_stats.npz"
            mc_shards = outdir / "mc_shards" / f"formation_obj_{idx+1:02d}"
            ut_path = outdir / f"{obj}_ut_merged.npz"
            t_grid, det_states, stm_mean, stm_cov = load_det_stm(det_path)
            if mc_stats.exists():
                _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_stats(mc_stats)
            elif mc_shards.exists():
                _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_shards(mc_shards)
            else:
                raise RuntimeError(f"MC stats not found for {obj}")
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
            save_json(case_dir / "summary.json", {"case": obj})
            if args.plot:
                plot_case(outdir, obj, t_grid, det_states, mc_mean, mc_std, ut_mean, stm_mean, args.plot_all_states)
            if do_pod:
                arc_results, pod_summary = run_pod(
                    det_states,
                    t_grid,
                    args.pod_sp3,
                    args.pod_eop,
                    args.pod_use_sat_clock,
                    args.pod_arc_s,
                    args.pod_overlap_s,
                    args.pod_skip_slr,
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
                    spacecraft_mass_kg=scenario.get("spacecraft_mass_kg"),
                    spacecraft_inertia_kgm2=scenario.get("spacecraft_inertia_kgm2"),
                )
                pod_dir = case_dir / "pod"
                pod_dir.mkdir(parents=True, exist_ok=True)
                save_json(pod_dir / "summary.json", {"pod_summaries": pod_summary})
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
        mc_stats = outdir / f"{name}_mc_stats.npz"
        mc_shards = outdir / "mc_shards" / name
        ut_path = outdir / f"{name}_ut_merged.npz"
        t_grid, det_states, stm_mean, stm_cov = load_det_stm(det_path)
        if mc_stats.exists():
            _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_stats(mc_stats)
        elif mc_shards.exists():
            _, mc_states, mc_mean, mc_std, mc_cov_final = load_mc_stats_from_shards(mc_shards)
        else:
            raise RuntimeError(f"MC stats not found for {name}")
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
        save_json(case_dir / "summary.json", {"case": name})
        if args.plot:
            plot_case(outdir, name, t_grid, det_states, mc_mean, mc_std, ut_mean, stm_mean, args.plot_all_states)
        if do_pod:
            arc_results, pod_summary = run_pod(
                det_states,
                t_grid,
                args.pod_sp3,
                args.pod_eop,
                args.pod_use_sat_clock,
                args.pod_arc_s,
                args.pod_overlap_s,
                args.pod_skip_slr,
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
                pod_enkf_use_carrier=bool(
                    scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)
                ),
                spacecraft_mass_kg=scenario.get("spacecraft_mass_kg"),
                spacecraft_inertia_kgm2=scenario.get("spacecraft_inertia_kgm2"),
            )
            pod_dir = case_dir / "pod"
            pod_dir.mkdir(parents=True, exist_ok=True)
            save_json(pod_dir / "summary.json", {"pod_summaries": pod_summary})


if __name__ == "__main__":
    main()
