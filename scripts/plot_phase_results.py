#!/usr/bin/env python3
"""Generate phase-level paper plots from run outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


PHASE_CASES = {
    "phase_a": [
        "att_a1_detumble_aero_only",
        "att_a2_nadir_pointing",
        "od_c1_gnss_continuous_bestcase",
        "od_c2_gnss_realistic_outages",
        "od_c3_gnss_plus_accel",
    ],
    "phase_b": [
        "att_a7_eclipse_thermal_transition",
        "od_c9_meas_error_stress_tests",
        "od_c8_station_weather_availability",
        "form_b6_robust_formation_storm",
    ],
    "phase_c": [
        "form_b2_alongtrack_2sat_dd_keeping",
        "form_b3_lvlh_box_constraint_3sat",
        "form_b5_close_approach_risk",
        "form_b7_cluster_sparse_ground_contact",
    ],
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _as_finite_float(value: Any) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    if not np.isfinite(out):
        return None
    return out


def _load_summaries(run_dir: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for path in sorted(run_dir.rglob("summary.json")):
        if path.parent.name == "pod":
            continue
        rel = str(path.parent.relative_to(run_dir))
        try:
            out[rel] = _load_json(path)
        except Exception:
            continue
    return out


def _case_children(summaries: dict[str, dict[str, Any]], case: str) -> list[tuple[str, dict[str, Any]]]:
    pref = case + "/obj_"
    rows = [(k, v) for k, v in summaries.items() if k.startswith(pref)]
    rows.sort(key=lambda x: x[0])
    return rows


def _case_metric(
    summaries: dict[str, dict[str, Any]],
    case: str,
    key: str,
) -> dict[str, Any]:
    base = summaries.get(case, {})
    base_val = _as_finite_float(base.get(key))
    if base_val is not None:
        return {"available": True, "source": "base", "value": float(base_val)}

    vals = []
    for _, payload in _case_children(summaries, case):
        v = _as_finite_float(payload.get(key))
        if v is not None:
            vals.append(v)
    if not vals:
        return {"available": False}
    arr = np.asarray(vals, dtype=float)
    return {
        "available": True,
        "source": "objects",
        "value": float(np.mean(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "count": int(arr.size),
    }


def _case_payload_valid_fraction(summaries: dict[str, dict[str, Any]], case: str) -> dict[str, Any]:
    base = summaries.get(case, {})
    payload = base.get("payload_impact")
    if isinstance(payload, dict):
        v = _as_finite_float(payload.get("valid_sample_fraction"))
        if v is not None:
            return {"available": True, "source": "base", "value": float(v)}

    vals = []
    for _, entry in _case_children(summaries, case):
        p = entry.get("payload_impact")
        if not isinstance(p, dict):
            continue
        v = _as_finite_float(p.get("valid_sample_fraction"))
        if v is not None:
            vals.append(v)
    if not vals:
        return {"available": False}
    arr = np.asarray(vals, dtype=float)
    return {
        "available": True,
        "source": "objects",
        "value": float(np.mean(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "count": int(arr.size),
    }


def _case_attitude_full_p95(summaries: dict[str, dict[str, Any]], case: str) -> dict[str, Any]:
    base = summaries.get(case, {})
    modes = base.get("attitude_uq_modes")
    if not isinstance(modes, dict):
        return {"available": False}
    full_mode = modes.get("full_mode")
    if not isinstance(full_mode, dict):
        return {"available": False}
    v = _as_finite_float(full_mode.get("final_p95_error_deg"))
    if v is None:
        return {"available": False}
    return {"available": True, "value": float(v)}


def _phase_case_list(phase: str) -> list[str]:
    if phase == "all":
        out: list[str] = []
        for key in ("phase_a", "phase_b", "phase_c"):
            out.extend(PHASE_CASES[key])
        return out
    return list(PHASE_CASES[phase])


def _safe_import_pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    return plt


def _plot_method_errors(
    outdir: Path,
    *,
    labels: list[str],
    ut_values: list[float],
    stm_values: list[float],
    title: str,
    filename: str,
) -> str | None:
    plt = _safe_import_pyplot()
    if plt is None or not labels:
        return None
    x = np.arange(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.8))
    ax.bar(x - w / 2.0, ut_values, width=w, label="UT", color="#1f77b4")
    ax.bar(x + w / 2.0, stm_values, width=w, label="STM", color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("error metric")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = outdir / filename
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _plot_payload_fraction(
    outdir: Path,
    *,
    labels: list[str],
    values: list[float],
) -> str | None:
    plt = _safe_import_pyplot()
    if plt is None or not labels:
        return None
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.8))
    ax.barh(y, values, color="#2ca02c")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("payload valid sample fraction")
    ax.set_xlim(0.0, 1.0)
    ax.set_title("Payload Impact Geometry Validity")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    out = outdir / "payload_valid_fraction.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _plot_attitude_p95(
    outdir: Path,
    *,
    labels: list[str],
    values: list[float],
) -> str | None:
    plt = _safe_import_pyplot()
    if plt is None or not labels:
        return None
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.8))
    ax.bar(x, values, color="#9467bd")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("final p95 attitude error [deg]")
    ax.set_title("Attitude UQ Full-Mode Final P95")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    out = outdir / "attitude_uq_fullmode_p95.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _plot_formation_separation(
    run_dir: Path,
    outdir: Path,
    *,
    cases: list[str],
) -> str | None:
    plt = _safe_import_pyplot()
    if plt is None:
        return None
    rows: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []
    for case in cases:
        path = run_dir / case / "formation_metrics.npz"
        if not path.exists():
            continue
        try:
            with np.load(path) as data:
                t_grid = np.asarray(data["t_grid"], dtype=float)
                det = np.asarray(data["det_separation_m"], dtype=float)
                mc = np.asarray(data["mc_mean_separation_m"], dtype=float)
            if det.ndim != 2 or mc.ndim != 2:
                continue
            rows.append((case, t_grid, np.mean(det, axis=0), np.mean(mc, axis=0)))
        except Exception:
            continue
    if not rows:
        return None

    n = len(rows)
    fig, axes = plt.subplots(n, 1, figsize=(9, max(3.2, 2.8 * n)), sharex=False)
    if n == 1:
        axes = [axes]
    for ax, (name, t_grid, det_m, mc_m) in zip(axes, rows):
        ax.plot(t_grid, det_m, label="DET separation mean", color="#1f77b4")
        ax.plot(t_grid, mc_m, label="MC mean separation", color="#d62728")
        ax.set_ylabel("separation [m]")
        ax.set_title(name)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
    axes[-1].set_xlabel("time [s]")
    fig.tight_layout()
    out = outdir / "formation_separation_overview.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def build_phase_plot_bundle(run_dir: Path, *, phase: str, outdir: Path) -> dict[str, Any]:
    summaries = _load_summaries(run_dir)
    cases = _phase_case_list(phase)
    outdir.mkdir(parents=True, exist_ok=True)

    method_rows = []
    for case in cases:
        ut = _case_metric(summaries, case, "ut_mean_err")
        stm = _case_metric(summaries, case, "stm_mean_err")
        ut_cov = _case_metric(summaries, case, "ut_cov_rel_err")
        stm_cov = _case_metric(summaries, case, "stm_cov_rel_err")
        method_rows.append(
            {
                "case": case,
                "ut_mean_err": ut,
                "stm_mean_err": stm,
                "ut_cov_rel_err": ut_cov,
                "stm_cov_rel_err": stm_cov,
            }
        )

    labels_mean = []
    ut_vals_mean = []
    stm_vals_mean = []
    labels_cov = []
    ut_vals_cov = []
    stm_vals_cov = []
    for row in method_rows:
        ut_m = row["ut_mean_err"]
        stm_m = row["stm_mean_err"]
        if ut_m.get("available") and stm_m.get("available"):
            labels_mean.append(row["case"])
            ut_vals_mean.append(float(ut_m["value"]))
            stm_vals_mean.append(float(stm_m["value"]))
        ut_c = row["ut_cov_rel_err"]
        stm_c = row["stm_cov_rel_err"]
        if ut_c.get("available") and stm_c.get("available"):
            labels_cov.append(row["case"])
            ut_vals_cov.append(float(ut_c["value"]))
            stm_vals_cov.append(float(stm_c["value"]))

    payload_rows = []
    payload_labels = []
    payload_vals = []
    for case in cases:
        row = _case_payload_valid_fraction(summaries, case)
        payload_rows.append({"case": case, "payload_valid_fraction": row})
        if row.get("available"):
            payload_labels.append(case)
            payload_vals.append(float(row["value"]))

    att_rows = []
    att_labels = []
    att_vals = []
    for case in cases:
        row = _case_attitude_full_p95(summaries, case)
        att_rows.append({"case": case, "attitude_fullmode_final_p95_deg": row})
        if row.get("available"):
            att_labels.append(case)
            att_vals.append(float(row["value"]))

    plots = {}
    plots["method_mean_error"] = _plot_method_errors(
        outdir,
        labels=labels_mean,
        ut_values=ut_vals_mean,
        stm_values=stm_vals_mean,
        title=f"{phase}: UT/STM Mean Error vs MC",
        filename="method_mean_error.png",
    )
    plots["method_cov_error"] = _plot_method_errors(
        outdir,
        labels=labels_cov,
        ut_values=ut_vals_cov,
        stm_values=stm_vals_cov,
        title=f"{phase}: UT/STM Covariance Relative Error vs MC",
        filename="method_cov_error.png",
    )
    plots["payload_valid_fraction"] = _plot_payload_fraction(outdir, labels=payload_labels, values=payload_vals)
    plots["attitude_fullmode_p95"] = _plot_attitude_p95(outdir, labels=att_labels, values=att_vals)
    plots["formation_separation"] = _plot_formation_separation(run_dir, outdir, cases=cases)

    report = {
        "run_dir": str(run_dir),
        "phase": phase,
        "cases": cases,
        "method_metrics": method_rows,
        "payload_metrics": payload_rows,
        "attitude_metrics": att_rows,
        "plots": plots,
    }
    report_path = outdir / "phase_plot_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate phase-level summary plots.")
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--phase", choices=("phase_a", "phase_b", "phase_c", "all"), default="all")
    parser.add_argument("--outdir", default=None)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    outdir = Path(args.outdir) if args.outdir else (run_dir / f"plots_{args.phase}")
    report = build_phase_plot_bundle(run_dir, phase=args.phase, outdir=outdir)
    print(f"[phase-plots] wrote {outdir / 'phase_plot_report.json'}")
    has_any = any(v for v in report.get("plots", {}).values())
    raise SystemExit(0 if has_any else 1)


if __name__ == "__main__":
    main()
