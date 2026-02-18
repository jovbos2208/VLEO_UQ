#!/usr/bin/env python3
"""Generate thesis-ready plots from convergence study outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


THESIS_TEXT_WIDTH_MM = 160.0  # A4 with 25 mm margins on both sides.


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def finite_or_nan(values: np.ndarray) -> np.ndarray:
    arr = np.array(values, dtype=float)
    arr[~np.isfinite(arr)] = np.nan
    return arr


def metric_series(summary_rows: list[dict], key: str) -> tuple[np.ndarray, np.ndarray]:
    rows = sorted(summary_rows, key=lambda r: float(r["dt_floor_s"]))
    x = np.array([float(r["dt_floor_s"]) for r in rows], dtype=float)
    y = finite_or_nan(np.array([float(r.get(key, np.nan)) for r in rows], dtype=float))
    return x, y


def save_figure(fig, outbase: Path, formats: list[str], dpi: int) -> list[str]:
    paths = []
    for fmt in formats:
        out = outbase.with_suffix(f".{fmt}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        paths.append(str(out))
    return paths


def _mm_to_in(value_mm: float) -> float:
    return float(value_mm) / 25.4


def setup_plot_style(base_font_size: float = 11.0):
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.figsize": (_mm_to_in(THESIS_TEXT_WIDTH_MM), _mm_to_in(95.0)),
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": base_font_size,
            "axes.labelsize": base_font_size,
            "axes.titlesize": base_font_size + 1.0,
            "legend.fontsize": base_font_size - 0.5,
            "xtick.labelsize": base_font_size - 0.5,
            "ytick.labelsize": base_font_size - 0.5,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": "-",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.linewidth": 0.9,
            "lines.linewidth": 2.0,
            "lines.markersize": 5.0,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_summary_comparison(
    att_summary: list[dict],
    mission_summary: list[dict],
    outdir: Path,
    formats: list[str],
    dpi: int,
) -> list[str]:
    import matplotlib.pyplot as plt

    files = []
    x_att, y_att_pos = metric_series(att_summary, "mean_pos_rms_m")
    _, y_att_vel = metric_series(att_summary, "mean_vel_rms_mps")
    _, y_att_energy = metric_series(att_summary, "mean_energy_drift_rms")
    x_mis, y_mis_pos = metric_series(mission_summary, "mean_pos_rms_m")
    _, y_mis_vel = metric_series(mission_summary, "mean_vel_rms_mps")
    _, y_mis_energy = metric_series(mission_summary, "mean_energy_drift_rms")

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(_mm_to_in(THESIS_TEXT_WIDTH_MM), _mm_to_in(62.0)),
        constrained_layout=True,
    )
    panels = [
        ("Mean Position RMS [m]", y_att_pos, y_mis_pos),
        ("Mean Velocity RMS [m/s]", y_att_vel, y_mis_vel),
        ("Mean Energy Drift RMS", y_att_energy, y_mis_energy),
    ]
    for ax, (title, y_att, y_mis) in zip(axes, panels):
        ax.plot(x_att, y_att, marker="o", lw=2.0, color="#1f77b4", label="ATT_A*")
        ax.plot(x_mis, y_mis, marker="s", lw=2.0, color="#d62728", label="OD_C*")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Timestep floor $\\Delta t$ [s]")
        ax.set_title(title)
        ax.legend(frameon=False)
    fig.suptitle("Convergence Summary vs Timestep", fontsize=12.5, fontweight="bold")
    files += save_figure(fig, outdir / "fig01_convergence_summary", formats, dpi)
    plt.close(fig)

    return files


def plot_attitude_only(
    att_summary: list[dict],
    outdir: Path,
    formats: list[str],
    dpi: int,
) -> list[str]:
    import matplotlib.pyplot as plt

    files = []
    x, y_point = metric_series(att_summary, "mean_boresight_rms_deg")
    _, y_point_p95 = metric_series(att_summary, "mean_boresight_p95_deg")
    _, y_rate = metric_series(att_summary, "mean_rate_rms_deg_s")
    _, y_stm = metric_series(att_summary, "mean_stm_cov_final_rel_frob")

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(_mm_to_in(THESIS_TEXT_WIDTH_MM), _mm_to_in(62.0)),
        constrained_layout=True,
    )
    axes[0].plot(x, y_point, marker="o", lw=2.0, color="#0B6E4F", label="RMS")
    axes[0].plot(x, y_point_p95, marker="s", lw=1.8, color="#1F9D55", label="p95")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("$\\Delta t$ [s]")
    axes[0].set_ylabel("deg")
    axes[0].set_title("Mean Boresight Error")
    axes[0].legend(frameon=False)

    axes[1].plot(x, y_rate, marker="o", lw=2.0, color="#BA3B46")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("$\\Delta t$ [s]")
    axes[1].set_ylabel("deg/s")
    axes[1].set_title("Mean Rate RMS")

    y_stm_plot = np.where(np.isfinite(y_stm) & (y_stm > 0.0), y_stm, np.nan)
    axes[2].plot(x, y_stm_plot, marker="o", lw=2.0, color="#3A506B")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("$\\Delta t$ [s]")
    axes[2].set_ylabel("relative Frobenius")
    axes[2].set_title("Mean STM Covariance Error")

    fig.suptitle("Attitude Convergence Metrics", fontsize=13, fontweight="bold")
    files += save_figure(fig, outdir / "fig02_attitude_metrics", formats, dpi)
    plt.close(fig)
    return files


def _scenario_dt_matrix(
    metrics_rows: list[dict],
    metric_key: str,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    scenarios = sorted({str(r["scenario"]) for r in metrics_rows})
    dts = sorted({float(r["dt_floor_s"]) for r in metrics_rows})
    mat = np.full((len(scenarios), len(dts)), np.nan, dtype=float)
    idx_s = {name: i for i, name in enumerate(scenarios)}
    idx_d = {dt: j for j, dt in enumerate(dts)}
    for row in metrics_rows:
        s = str(row["scenario"])
        dt = float(row["dt_floor_s"])
        val = float(row.get(metric_key, np.nan))
        if np.isfinite(val):
            mat[idx_s[s], idx_d[dt]] = val
    return scenarios, np.array(dts, dtype=float), mat


def _shorten_scenario(name: str) -> str:
    tokens = name.split("_")
    if len(tokens) >= 2:
        return "_".join(tokens[:2]).upper()
    return name


def plot_heatmaps(
    metrics_rows: list[dict],
    outdir: Path,
    formats: list[str],
    dpi: int,
    prefix: str,
    include_attitude_metric: bool,
) -> list[str]:
    import matplotlib.pyplot as plt

    files = []
    targets = [("pos_rms_m", "Position RMS [m]")]
    if include_attitude_metric:
        targets.append(("boresight_rms_deg", "Boresight RMS [deg]"))

    for key, title in targets:
        scenarios, dts, mat = _scenario_dt_matrix(metrics_rows, key)
        if np.all(~np.isfinite(mat)):
            continue
        with np.errstate(invalid="ignore"):
            z = np.log10(np.where(mat > 0.0, mat, np.nan))
        fig_h_mm = max(95.0, 11.0 * len(scenarios))
        fig, ax = plt.subplots(
            figsize=(_mm_to_in(THESIS_TEXT_WIDTH_MM), _mm_to_in(fig_h_mm)),
            constrained_layout=True,
        )
        im = ax.imshow(z, aspect="auto", cmap="viridis")
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label(f"log10({title})")
        ax.set_title(f"{prefix}: {title} by Scenario and dt")
        ax.set_xlabel("$\\Delta t$ [s]")
        ax.set_ylabel("Scenario")
        ax.set_xticks(np.arange(len(dts)))
        ax.set_xticklabels([f"{x:g}" for x in dts])
        ax.set_yticks(np.arange(len(scenarios)))
        ax.set_yticklabels([_shorten_scenario(s) for s in scenarios])
        files += save_figure(fig, outdir / f"{prefix.lower()}_{key}_heatmap", formats, dpi)
        plt.close(fig)
    return files


def write_latex_table(
    summary_rows: list[dict],
    type_name: str,
    outpath: Path,
) -> None:
    def fmt(x: float) -> str:
        if not np.isfinite(x):
            return "--"
        if x == 0.0:
            return "0"
        ax = abs(x)
        if ax >= 1e4 or ax < 1e-2:
            return f"{x:.3e}"
        return f"{x:.4g}"

    rows = sorted(summary_rows, key=lambda r: float(r["dt_floor_s"]))
    is_att = any("mean_boresight_rms_deg" in r for r in rows)
    if is_att:
        header = (
            "\\begin{tabular}{rcccccc}\n"
            "\\toprule\n"
            "$\\Delta t$ [s] & Mean pos RMS [m] & Mean vel RMS [m/s] & Mean boresight RMS [deg] & "
            "Mean boresight p95 [deg] & Mean rate RMS [deg/s] & Mean STM cov err \\\\\n"
            "\\midrule\n"
        )
    else:
        header = (
            "\\begin{tabular}{rcccc}\n"
            "\\toprule\n"
            "$\\Delta t$ [s] & Mean pos RMS [m] & Mean vel RMS [m/s] & Mean energy drift & Mean STM cov err \\\\\n"
            "\\midrule\n"
        )
    body = ""
    for r in rows:
        dt = float(r["dt_floor_s"])
        p = float(r.get("mean_pos_rms_m", np.nan))
        v = float(r.get("mean_vel_rms_mps", np.nan))
        s = float(r.get("mean_stm_cov_final_rel_frob", np.nan))
        if is_att:
            b = float(r.get("mean_boresight_rms_deg", np.nan))
            bp95 = float(r.get("mean_boresight_p95_deg", np.nan))
            wr = float(r.get("mean_rate_rms_deg_s", np.nan))
            body += f"{dt:g} & {fmt(p)} & {fmt(v)} & {fmt(b)} & {fmt(bp95)} & {fmt(wr)} & {fmt(s)} \\\\\n"
        else:
            e = float(r.get("mean_energy_drift_rms", np.nan))
            body += f"{dt:g} & {fmt(p)} & {fmt(v)} & {fmt(e)} & {fmt(s)} \\\\\n"
    footer = "\\bottomrule\n\\end{tabular}\n"
    caption = f"% Convergence summary table for {type_name}\n"
    outpath.write_text(caption + header + body + footer, encoding="utf-8")


def write_latex_figure_snippets(manifest: dict, outdir: Path) -> None:
    lines = []
    lines.append("% Auto-generated include snippets for convergence thesis figures")
    lines.append("% Requires: \\usepackage{graphicx}")
    lines.append("")
    for fp in manifest["files"]:
        p = Path(fp)
        if p.suffix.lower() != ".pdf":
            continue
        rel = p.relative_to(outdir.parent).as_posix()
        label = p.stem.replace("_", "-")
        lines.append("\\begin{figure}[H]")
        lines.append("  \\centering")
        lines.append(f"  \\includegraphics[width=\\textwidth]{{{rel}}}")
        lines.append(f"  \\caption{{{p.stem.replace('_', ' ').title()}}}")
        lines.append(f"  \\label{{fig:{label}}}")
        lines.append("\\end{figure}")
        lines.append("")
    (outdir / "figures_include.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate thesis plots from convergence results.")
    parser.add_argument("--results_dir", required=True, help="Path to results/convergence_* directory.")
    parser.add_argument("--outdir", default=None, help="Default: <results_dir>/thesis_plots")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated formats, e.g. png,pdf")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font_size", type=float, default=11.0)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        raise FileNotFoundError(results_dir)
    outdir = Path(args.outdir) if args.outdir else (results_dir / "thesis_plots")
    ensure_dir(outdir)
    formats = [x.strip().lower() for x in args.formats.split(",") if x.strip()]
    if not formats:
        raise ValueError("no output formats specified")

    setup_plot_style(base_font_size=float(args.font_size))

    att_metrics = load_json(results_dir / "attitude" / "metrics.json")
    mission_metrics = load_json(results_dir / "mission_od" / "metrics.json")
    att_summary = load_json(results_dir / "attitude" / "summary.json")
    mission_summary = load_json(results_dir / "mission_od" / "summary.json")

    generated: list[str] = []
    generated += plot_summary_comparison(att_summary, mission_summary, outdir, formats, args.dpi)
    generated += plot_attitude_only(att_summary, outdir, formats, args.dpi)
    generated += plot_heatmaps(
        mission_metrics,
        outdir=outdir,
        formats=formats,
        dpi=args.dpi,
        prefix="MissionOD",
        include_attitude_metric=False,
    )
    generated += plot_heatmaps(
        att_metrics,
        outdir=outdir,
        formats=formats,
        dpi=args.dpi,
        prefix="Attitude",
        include_attitude_metric=True,
    )

    write_latex_table(mission_summary, "mission_od", outdir / "table_convergence_mission_od.tex")
    write_latex_table(att_summary, "attitude", outdir / "table_convergence_attitude.tex")

    manifest = {
        "results_dir": str(results_dir),
        "outdir": str(outdir),
        "formats": formats,
        "files": generated,
        "tables": [
            str(outdir / "table_convergence_mission_od.tex"),
            str(outdir / "table_convergence_attitude.tex"),
        ],
    }
    write_latex_figure_snippets(manifest, outdir)
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[thesis-plots] wrote {len(generated)} figures + 2 tables -> {outdir}")


if __name__ == "__main__":
    main()
