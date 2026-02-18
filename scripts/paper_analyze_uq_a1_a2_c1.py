#!/usr/bin/env python3
"""Aggregate paper-ready metrics/plots for A1/A2/C1 runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


def quat_angle_deg(q_a: np.ndarray, q_b: np.ndarray) -> np.ndarray:
    qa = np.array(q_a, dtype=float, copy=True)
    qb = np.array(q_b, dtype=float, copy=True)
    qa /= np.maximum(np.linalg.norm(qa, axis=1, keepdims=True), 1e-15)
    qb /= np.maximum(np.linalg.norm(qb, axis=1, keepdims=True), 1e-15)
    dots = np.clip(np.abs(np.sum(qa * qb, axis=1)), -1.0, 1.0)
    return np.degrees(2.0 * np.arccos(dots))


def rms(values: np.ndarray) -> float:
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(arr * arr)))


def p95(values: np.ndarray) -> float:
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return float("nan")
    return float(np.percentile(arr, 95.0))


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_pod_summaries(path: Path) -> list[dict[str, Any]]:
    payload = load_json(path)
    if isinstance(payload, dict) and isinstance(payload.get("pod_summaries"), list):
        return payload["pod_summaries"]
    if isinstance(payload, list):
        return payload
    return []


def scalar_or_nan(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_table(path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "scenario",
        "type",
        "mc_pos_rms_m",
        "ut_vs_mc_final_pos_m",
        "stm_vs_mc_final_pos_m",
        "attitude_err_p95_deg",
        "pod_mean_radial_rms_m",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        vals = []
        for h in headers:
            v = row.get(h, "")
            if isinstance(v, float):
                if np.isfinite(v):
                    vals.append(f"{v:.4g}")
                else:
                    vals.append("nan")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze A1/A2/C1 outputs for a paper workflow.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--outdir", default=None)
    args = parser.parse_args()

    config = load_json(Path(args.config))
    scenarios = config.get("scenarios", [])
    if not scenarios:
        raise ValueError("config has no scenarios")

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir) if args.outdir else (results_dir / "paper_analysis")
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    att_plot_data: list[dict[str, Any]] = []
    c1_pod_arcs: list[dict[str, Any]] = []

    for scenario in scenarios:
        name = str(scenario["name"])
        typ = str(scenario.get("type", ""))
        case_dir = results_dir / name
        npz_path = case_dir / "mc_ut_stm.npz"
        if not npz_path.exists():
            raise FileNotFoundError(f"missing {npz_path}")
        with np.load(npz_path) as data:
            t_grid = np.array(data["t_grid"], dtype=float)
            det = np.array(data["det_states"], dtype=float)
            mc_mean = np.array(data["mc_mean"], dtype=float)
            ut_mean = np.array(data["ut_mean"], dtype=float)
            stm_mean = np.array(data["stm_mean"], dtype=float)

        det_pos = det[:, 0:3]
        mc_pos = mc_mean[:, 0:3]
        ut_pos = ut_mean[:, 0:3]
        stm_pos = stm_mean[:, 0:3]

        mc_pos_err = np.linalg.norm(mc_pos - det_pos, axis=1)
        ut_vs_mc_err = np.linalg.norm(ut_pos - mc_pos, axis=1)
        stm_vs_mc_err = np.linalg.norm(stm_pos - mc_pos, axis=1)

        row: dict[str, Any] = {
            "scenario": name,
            "type": typ,
            "duration_s": scalar_or_nan(scenario.get("duration_s")),
            "dt_s": scalar_or_nan(scenario.get("dt_s")),
            "particles": int(scenario.get("particles", 0)),
            "mc_pos_rms_m": rms(mc_pos_err),
            "mc_final_pos_err_m": scalar_or_nan(mc_pos_err[-1]),
            "ut_vs_mc_pos_rms_m": rms(ut_vs_mc_err),
            "ut_vs_mc_final_pos_m": scalar_or_nan(ut_vs_mc_err[-1]),
            "stm_vs_mc_pos_rms_m": rms(stm_vs_mc_err),
            "stm_vs_mc_final_pos_m": scalar_or_nan(stm_vs_mc_err[-1]),
            "attitude_err_rms_deg": float("nan"),
            "attitude_err_p95_deg": float("nan"),
            "rate_err_rms_deg_s": float("nan"),
            "pod_estimator": scenario.get("pod_estimator", ""),
            "pod_n_arcs": 0,
            "pod_mean_radial_rms_m": float("nan"),
            "pod_mean_cov3_radial": float("nan"),
        }

        if typ == "attitude":
            att_err_deg = quat_angle_deg(mc_mean[:, 6:10], det[:, 6:10])
            rate_err_deg_s = np.degrees(np.linalg.norm(mc_mean[:, 10:13] - det[:, 10:13], axis=1))
            row["attitude_err_rms_deg"] = rms(att_err_deg)
            row["attitude_err_p95_deg"] = p95(att_err_deg)
            row["rate_err_rms_deg_s"] = rms(rate_err_deg_s)
            att_plot_data.append(
                {
                    "scenario": name,
                    "t_grid": t_grid,
                    "att_err_deg": att_err_deg,
                    "rate_err_deg_s": rate_err_deg_s,
                }
            )

        pod_dir = case_dir / "pod"
        pod_summary_path = pod_dir / "summary.json"
        if pod_summary_path.exists():
            pod_summaries = load_pod_summaries(pod_summary_path)
            row["pod_n_arcs"] = len(pod_summaries)
            radial_rms_vals = []
            cov3_vals = []
            for summary in pod_summaries:
                cov3 = summary.get("coverage_3sigma", {})
                if isinstance(cov3, dict):
                    cov3_vals.append(scalar_or_nan(cov3.get("radial")))
            for pod_arc_path in sorted(pod_dir.glob("pod_arc_*.npz")):
                with np.load(pod_arc_path) as pod_data:
                    radial_rms_vals.append(scalar_or_nan(pod_data["radial_rms_m"]))
                    if "rtn_error" in pod_data and "rtn_sigma" in pod_data:
                        c1_pod_arcs.append(
                            {
                                "scenario": name,
                                "arc": pod_arc_path.stem,
                                "t": np.array(pod_data["t_grid"], dtype=float),
                                "radial_err": np.array(pod_data["rtn_error"], dtype=float)[:, 0],
                                "radial_sig3": 3.0 * np.array(pod_data["rtn_sigma"], dtype=float)[:, 0],
                            }
                        )
            if radial_rms_vals:
                row["pod_mean_radial_rms_m"] = float(np.nanmean(radial_rms_vals))
            finite_cov3 = [x for x in cov3_vals if np.isfinite(x)]
            if finite_cov3:
                row["pod_mean_cov3_radial"] = float(np.mean(finite_cov3))

        rows.append(row)

    write_csv(outdir / "paper_metrics.csv", rows)
    (outdir / "paper_metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    write_markdown_table(outdir / "paper_metrics.md", rows)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = [row["scenario"] for row in rows]
        ut_final = [scalar_or_nan(row["ut_vs_mc_final_pos_m"]) for row in rows]
        stm_final = [scalar_or_nan(row["stm_vs_mc_final_pos_m"]) for row in rows]
        x = np.arange(len(names))
        width = 0.38
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.bar(x - width / 2.0, ut_final, width=width, label="UT vs MC")
        ax.bar(x + width / 2.0, stm_final, width=width, label="STM vs MC")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=20, ha="right")
        ax.set_ylabel("Final position error [m]")
        ax.set_title("Surrogate final position error by scenario")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "surrogate_final_position_error.png", dpi=180)
        plt.close(fig)

        if att_plot_data:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
            for rec in att_plot_data:
                t_h = rec["t_grid"] / 3600.0
                ax1.plot(t_h, rec["att_err_deg"], label=rec["scenario"])
                ax2.plot(t_h, rec["rate_err_deg_s"], label=rec["scenario"])
            ax1.set_ylabel("Attitude error [deg]")
            ax2.set_ylabel("Rate error [deg/s]")
            ax2.set_xlabel("Time [h]")
            ax1.set_title("A1/A2 attitude error traces")
            ax1.grid(True, alpha=0.3)
            ax2.grid(True, alpha=0.3)
            ax1.legend()
            fig.tight_layout()
            fig.savefig(outdir / "attitude_error_timeseries.png", dpi=180)
            plt.close(fig)

        if c1_pod_arcs:
            fig, ax = plt.subplots(figsize=(10, 4.5))
            for rec in c1_pod_arcs:
                t_h = rec["t"] / 3600.0
                ax.plot(t_h, rec["radial_err"], label=f"{rec['arc']} R err")
                ax.plot(t_h, rec["radial_sig3"], "--", alpha=0.7, label=f"{rec['arc']} 3sigma")
            ax.set_xlabel("Time [h]")
            ax.set_ylabel("Radial [m]")
            ax.set_title("C1 POD radial error and 3sigma envelopes")
            ax.grid(True, alpha=0.3)
            ax.legend(ncol=2, fontsize=8)
            fig.tight_layout()
            fig.savefig(outdir / "c1_pod_radial_errors.png", dpi=180)
            plt.close(fig)
    except Exception as exc:
        print(f"[paper-analysis] plotting skipped: {exc}")

    print(f"[paper-analysis] wrote metrics -> {outdir}")


if __name__ == "__main__":
    main()
