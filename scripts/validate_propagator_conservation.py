#!/usr/bin/env python3
"""Validate propagator behavior with conservation diagnostics and plots."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from pathlib import Path

import numpy as np


def parse_dt_values(raw: str) -> list[float]:
    vals = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        vals.append(float(token))
    vals = sorted(set(vals))
    if not vals:
        raise ValueError("dt_values is empty")
    return vals


def build_env(t_grid: np.ndarray, EnvInputs, density: float) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = float(density)
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def build_t_grid(duration_s: float, dt_s: float) -> np.ndarray:
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    n = int(np.floor(duration_s / dt_s))
    t = np.arange(n + 1, dtype=float) * dt_s
    if t[-1] < duration_s:
        t = np.concatenate([t, np.array([duration_s], dtype=float)])
    else:
        t[-1] = duration_s
    return t


def specific_energy(mu: float, r: np.ndarray, v: np.ndarray) -> np.ndarray:
    return 0.5 * np.sum(v * v, axis=1) - mu / np.linalg.norm(r, axis=1)


def angular_momentum(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    return np.cross(r, v)


def angle_between_series(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    an = np.linalg.norm(a, axis=1, keepdims=True)
    bn = np.linalg.norm(b, axis=1, keepdims=True)
    an[an == 0.0] = 1.0
    bn[bn == 0.0] = 1.0
    aa = a / an
    bb = b / bn
    dots = np.clip(np.sum(aa * bb, axis=1), -1.0, 1.0)
    return np.rad2deg(np.arccos(dots))


def orbital_elements(mu: float, r: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = angular_momentum(r, v)
    h_norm = np.linalg.norm(h, axis=1)
    eps = specific_energy(mu, r, v)
    a = -mu / np.maximum(2.0 * eps, -1e-30)
    e_vec = np.cross(v, h) / mu - r / np.linalg.norm(r, axis=1, keepdims=True)
    e = np.linalg.norm(e_vec, axis=1)
    return a, e


def summarize_run(mu: float, t_grid: np.ndarray, X: np.ndarray, x0: np.ndarray) -> tuple[dict, dict]:
    r = X[:, 0:3]
    v = X[:, 3:6]
    q = X[:, 6:10]
    eps = specific_energy(mu, r, v)
    h = angular_momentum(r, v)
    h_norm = np.linalg.norm(h, axis=1)
    a, e = orbital_elements(mu, r, v)

    eps0 = float(eps[0])
    h0 = float(h_norm[0])
    a0 = float(a[0])
    e0 = float(e[0])
    h0_vec = h[0:1].repeat(len(h), axis=0)

    rel_eps = (eps - eps0) / max(abs(eps0), 1e-30)
    rel_h = (h_norm - h0) / max(abs(h0), 1e-30)
    h_dir_deg = angle_between_series(h, h0_vec)
    rel_a = (a - a0) / max(abs(a0), 1e-30)
    abs_e = np.abs(e - e0)
    q_norm_err = np.abs(np.linalg.norm(q, axis=1) - 1.0)

    metrics = {
        "energy_rel_rms": float(np.sqrt(np.mean(rel_eps * rel_eps))),
        "energy_rel_final": float(rel_eps[-1]),
        "energy_rel_max_abs": float(np.max(np.abs(rel_eps))),
        "h_rel_rms": float(np.sqrt(np.mean(rel_h * rel_h))),
        "h_rel_final": float(rel_h[-1]),
        "h_rel_max_abs": float(np.max(np.abs(rel_h))),
        "h_dir_rms_deg": float(np.sqrt(np.mean(h_dir_deg * h_dir_deg))),
        "h_dir_max_deg": float(np.max(np.abs(h_dir_deg))),
        "a_rel_rms": float(np.sqrt(np.mean(rel_a * rel_a))),
        "e_abs_rms": float(np.sqrt(np.mean(abs_e * abs_e))),
        "closure_pos_m": float(np.linalg.norm(X[-1, 0:3] - x0[0:3])),
        "closure_vel_mps": float(np.linalg.norm(X[-1, 3:6] - x0[3:6])),
        "quat_norm_max_err": float(np.max(q_norm_err)),
    }
    series = {
        "t_grid": t_grid,
        "energy_rel": rel_eps,
        "h_rel": rel_h,
        "h_dir_deg": h_dir_deg,
        "a_rel": rel_a,
        "e_abs": abs_e,
        "energy": eps,
        "h_norm": h_norm,
    }
    return metrics, series


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = sorted({k for row in rows for k in row.keys()})
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def setup_plot_style():
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10.5,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_conservation_vs_dt(rows: list[dict], outdir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    dts = np.array([r["dt_s"] for r in rows], dtype=float)
    e_rms = np.array([r["energy_rel_rms"] for r in rows], dtype=float)
    h_rms = np.array([r["h_rel_rms"] for r in rows], dtype=float)
    cpos = np.array([r["closure_pos_m"] for r in rows], dtype=float)
    cvel = np.array([r["closure_vel_mps"] for r in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2), constrained_layout=True)
    axes[0].plot(dts, e_rms, "o-", label=r"$\mathrm{RMS}(\Delta \epsilon/\epsilon_0)$")
    axes[0].plot(dts, h_rms, "s-", label=r"$\mathrm{RMS}(\Delta h/h_0)$")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$\Delta t$ [s]")
    axes[0].set_ylabel("Relative error")
    axes[0].set_title("Conservation Error vs Timestep")
    axes[0].legend(frameon=False)

    axes[1].plot(dts, cpos, "o-", label="orbit closure position [m]")
    axes[1].plot(dts, cvel, "s-", label="orbit closure velocity [m/s]")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$\Delta t$ [s]")
    axes[1].set_title("Closure Error vs Timestep")
    axes[1].legend(frameon=False)

    p_png = outdir / "fig01_conservation_vs_dt.png"
    p_pdf = outdir / "fig01_conservation_vs_dt.pdf"
    fig.savefig(p_png, dpi=300, bbox_inches="tight")
    fig.savefig(p_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [p_png, p_pdf]


def plot_timeseries(series_by_dt: dict[float, dict], outdir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    dt_min = min(series_by_dt.keys())
    dt_max = max(series_by_dt.keys())
    s_min = series_by_dt[dt_min]
    s_max = series_by_dt[dt_max]

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.0), constrained_layout=True)
    for s, label, color in [(s_min, f"dt={dt_min:g}s", "#1f77b4"), (s_max, f"dt={dt_max:g}s", "#d62728")]:
        t_orb = s["t_grid"] / s["t_grid"][-1]
        axes[0, 0].plot(t_orb, s["energy_rel"], color=color, label=label)
        axes[0, 1].plot(t_orb, s["h_rel"], color=color, label=label)
        axes[1, 0].plot(t_orb, s["h_dir_deg"], color=color, label=label)
        axes[1, 1].plot(t_orb, s["a_rel"], color=color, label=label)

    axes[0, 0].set_title(r"Relative specific energy drift $\Delta \epsilon/\epsilon_0$")
    axes[0, 1].set_title(r"Relative angular momentum norm drift $\Delta h/h_0$")
    axes[1, 0].set_title("Angular momentum direction error [deg]")
    axes[1, 1].set_title(r"Semi-major axis drift $\Delta a/a_0$")
    for ax in axes.flat:
        ax.set_xlabel("Normalized time [0, 1]")
    axes[0, 0].legend(frameon=False)

    p_png = outdir / "fig02_conservation_timeseries.png"
    p_pdf = outdir / "fig02_conservation_timeseries.pdf"
    fig.savefig(p_png, dpi=300, bbox_inches="tight")
    fig.savefig(p_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [p_png, p_pdf]


def plot_drag_comparison(conservative: dict, drag: dict, outdir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    t0 = conservative["t_grid"] / 3600.0
    t1 = drag["t_grid"] / 3600.0

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2), constrained_layout=True)
    axes[0].plot(t0, conservative["energy"], label="Vacuum (rho=0)")
    axes[0].plot(t1, drag["energy"], label="With drag (rho>0)")
    axes[0].set_xlabel("Time [h]")
    axes[0].set_ylabel("Specific energy [J/kg]")
    axes[0].set_title("Energy trend")
    axes[0].legend(frameon=False)

    axes[1].plot(t0, conservative["h_norm"], label="Vacuum (rho=0)")
    axes[1].plot(t1, drag["h_norm"], label="With drag (rho>0)")
    axes[1].set_xlabel("Time [h]")
    axes[1].set_ylabel("Specific angular momentum [m$^2$/s]")
    axes[1].set_title("Angular momentum trend")
    axes[1].legend(frameon=False)

    p_png = outdir / "fig03_drag_vs_conservative.png"
    p_pdf = outdir / "fig03_drag_vs_conservative.pdf"
    fig.savefig(p_png, dpi=300, bbox_inches="tight")
    fig.savefig(p_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [p_png, p_pdf]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run propagator conservation validation and generate plots.")
    parser.add_argument("--dt_values", default="1,2,5,10,30,60")
    parser.add_argument("--orbits", type=float, default=3.0)
    parser.add_argument("--radius_m", type=float, default=7000e3)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--drag_density", type=float, default=1e-12)
    parser.add_argument("--skip_drag_case", action="store_true")
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--max_step_s", type=float, default=120.0)
    args = parser.parse_args()

    from vleo_uq import AeroAdapter, DeterministicPropagator, EnvInputs, PropagatorConfig, VehicleParams, default_geometry

    dt_values = parse_dt_values(args.dt_values)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else (Path("results") / f"propagator_validation_{stamp}")
    outdir.mkdir(parents=True, exist_ok=True)

    setup_plot_style()

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()

    cfg = PropagatorConfig()
    cfg.rtol = float(args.rtol)
    cfg.atol = float(args.atol)
    cfg.max_step_s = float(args.max_step_s)
    cfg.rho_fast_sigma = 0.0
    cfg.rho_bias_sigma = 0.0
    cfg.wind_sigma = 0.0
    cfg.rng_seed = int(args.seed)
    cfg.freeze_attitude = True
    prop = DeterministicPropagator(aero, vehicle, cfg)

    mu = float(cfg.mu_earth_m3_s2)
    r0 = np.array([float(args.radius_m), 0.0, 0.0], dtype=float)
    v0 = np.array([0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0], dtype=float)
    x0 = np.zeros(prop.state_size, dtype=float)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0

    period_s = float(2.0 * np.pi * np.sqrt(np.linalg.norm(r0) ** 3 / mu))
    duration_s = float(args.orbits) * period_s

    rows: list[dict] = []
    series_by_dt: dict[float, dict] = {}
    for dt_s in dt_values:
        print(f"[prop-val] running dt={dt_s:g}s ...", flush=True)
        t_grid = build_t_grid(duration_s, dt_s)
        env = build_env(t_grid, EnvInputs, density=0.0)
        X = prop.propagate(x0, t_grid, env)
        metrics, series = summarize_run(mu, t_grid, X, x0)
        metrics["dt_s"] = float(dt_s)
        metrics["num_steps"] = int(t_grid.size)
        rows.append(metrics)
        series_by_dt[float(dt_s)] = series

    rows = sorted(rows, key=lambda r: r["dt_s"])
    (outdir / "conservation_metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    write_csv(outdir / "conservation_metrics.csv", rows)

    generated = []
    generated += [str(p) for p in plot_conservation_vs_dt(rows, outdir)]
    generated += [str(p) for p in plot_timeseries(series_by_dt, outdir)]

    drag_metrics = None
    if not args.skip_drag_case:
        print("[prop-val] running drag comparison case ...", flush=True)
        dt_ref = min(dt_values)
        t_grid = build_t_grid(duration_s, dt_ref)
        env_vac = build_env(t_grid, EnvInputs, density=0.0)
        env_drag = build_env(t_grid, EnvInputs, density=float(args.drag_density))
        X_vac = prop.propagate(x0, t_grid, env_vac)
        X_drag = prop.propagate(x0, t_grid, env_drag)
        _, s_vac = summarize_run(mu, t_grid, X_vac, x0)
        m_drag, s_drag = summarize_run(mu, t_grid, X_drag, x0)
        drag_metrics = m_drag
        generated += [str(p) for p in plot_drag_comparison(s_vac, s_drag, outdir)]
        np.savez(
            outdir / "drag_comparison_series.npz",
            t_grid=t_grid,
            energy_vac=s_vac["energy"],
            energy_drag=s_drag["energy"],
            h_vac=s_vac["h_norm"],
            h_drag=s_drag["h_norm"],
        )

    summary = {
        "mu_earth_m3_s2": mu,
        "period_s": period_s,
        "duration_s": duration_s,
        "dt_values_s": dt_values,
        "rows": rows,
        "drag_case": drag_metrics,
        "figures": generated,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[prop-val] wrote validation outputs -> {outdir}")


if __name__ == "__main__":
    main()
