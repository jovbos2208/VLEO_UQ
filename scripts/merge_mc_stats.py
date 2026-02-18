from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute MC stats from shard files without merging.")
    parser.add_argument("--indir", required=True, help="input directory with ensemble_rank*.npz")
    parser.add_argument("--out", required=True, help="output npz with mc_mean/mc_std/mc_cov_final")
    args = parser.parse_args()

    indir = Path(args.indir)
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

    if n_total == 0:
        # All members failed in upstream propagation for this object/shard set.
        # Keep pipeline alive and emit NaN stats so downstream can flag degraded quality.
        nt = sum_x.shape[0]
        nx = sum_x.shape[1]
        mc_mean = np.full((nt, nx), np.nan, dtype=float)
        mc_std = np.full((nt, nx), np.nan, dtype=float)
        mc_cov_final = np.full((nx, nx), np.nan, dtype=float)
        print(f"[warn] {indir}: no valid MC samples; writing NaN statistics")
    elif n_total == 1:
        # Variance is undefined for one sample; use deterministic fallback.
        mc_mean = sum_x.copy()
        mc_std = np.zeros_like(sum_x, dtype=float)
        mc_cov_final = np.zeros_like(sum_xx_final, dtype=float)
        print(f"[warn] {indir}: only one valid MC sample; covariance set to zeros")
    else:
        mc_mean = sum_x / n_total
        mc_var = (sum_x2 - (sum_x ** 2) / n_total) / (n_total - 1)
        mc_std = np.sqrt(np.maximum(mc_var, 0.0))
        mc_cov_final = (sum_xx_final - np.outer(sum_final, sum_final) / n_total) / (n_total - 1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_path,
        t_grid=t_grid,
        mc_mean=mc_mean,
        mc_std=mc_std,
        mc_cov_final=mc_cov_final,
        n_total=n_total,
    )
    print(f"Wrote MC stats to {out_path}")


if __name__ == "__main__":
    main()
