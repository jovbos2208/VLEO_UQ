from __future__ import annotations

import argparse
import importlib
import os
import time
from pathlib import Path
from typing import Callable, Dict, Tuple

import numpy as np


def parse_scenario(spec: str) -> Callable[[], Dict[str, object]]:
    if ":" not in spec:
        raise ValueError("scenario must be in module:function form")
    module_name, func_name = spec.split(":", 1)
    module = importlib.import_module(module_name)
    func = getattr(module, func_name)
    if not callable(func):
        raise ValueError("scenario is not callable")
    return func


def shard_from_env() -> Tuple[int, int]:
    if "SLURM_ARRAY_TASK_ID" in os.environ and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        return int(os.environ["SLURM_ARRAY_TASK_ID"]), int(os.environ["SLURM_ARRAY_TASK_COUNT"])
    if "SLURM_PROCID" in os.environ and "SLURM_NTASKS" in os.environ:
        return int(os.environ["SLURM_PROCID"]), int(os.environ["SLURM_NTASKS"])
    if "OMPI_COMM_WORLD_RANK" in os.environ and "OMPI_COMM_WORLD_SIZE" in os.environ:
        return int(os.environ["OMPI_COMM_WORLD_RANK"]), int(os.environ["OMPI_COMM_WORLD_SIZE"])
    return -1, -1


def generate_sigma_points(x0: np.ndarray, P0: np.ndarray, alpha: float, beta: float, kappa: float):
    n = x0.size
    lam = alpha ** 2 * (n + kappa) - n
    scale = n + lam
    if scale <= 0.0:
        raise ValueError("invalid UT scaling; n+lambda must be > 0")
    P = 0.5 * (P0 + P0.T)
    jitter = 0.0
    for _ in range(6):
        try:
            L = np.linalg.cholesky(scale * (P + jitter * np.eye(n)))
            break
        except np.linalg.LinAlgError:
            diag_max = float(np.max(np.diag(P))) if P.size else 1.0
            jitter = (1e-12 if jitter == 0.0 else jitter * 10.0) * max(1.0, diag_max)
    else:
        evals, evecs = np.linalg.eigh(P)
        eps = 1e-12 * max(1.0, float(np.max(evals)))
        evals = np.maximum(evals, eps)
        L = np.sqrt(scale) * (evecs @ np.diag(np.sqrt(evals)))
    sigmas = np.zeros((2 * n + 1, n), dtype=float)
    sigmas[0] = x0
    for i in range(n):
        col = L[:, i]
        sigmas[1 + i] = x0 + col
        sigmas[1 + n + i] = x0 - col
    w0m = lam / scale
    w0c = w0m + (1.0 - alpha ** 2 + beta)
    wi = 1.0 / (2.0 * scale)
    return sigmas, w0m, w0c, wi


def normalize_quaternions(sigmas: np.ndarray, ref: np.ndarray) -> None:
    q_ref = ref[6:10].copy()
    if np.linalg.norm(q_ref) == 0.0:
        q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    for i in range(sigmas.shape[0]):
        q = sigmas[i, 6:10]
        norm = np.linalg.norm(q)
        if norm == 0.0:
            q = q_ref.copy()
        else:
            q = q / norm
        if float(np.dot(q, q_ref)) < 0.0:
            q = -q
        sigmas[i, 6:10] = q


def propagate_chunked(
    propagator,
    x0: np.ndarray,
    t_grid: np.ndarray,
    env,
    chunk_steps: int,
    tag: str,
) -> np.ndarray:
    nt = len(t_grid)
    if nt == 0:
        raise ValueError("t_grid must have at least 1 entry")
    if nt == 1:
        out = np.zeros((1, x0.shape[0], x0.shape[1]), dtype=float)
        out[0] = x0
        return out

    out = np.zeros((nt, x0.shape[0], x0.shape[1]), dtype=float)
    x_curr = x0
    i0 = 0
    while i0 < nt - 1:
        i1 = min(nt - 1, i0 + chunk_steps)
        t_chunk = t_grid[i0 : i1 + 1] - t_grid[i0]
        env_chunk = env[i0 : i1 + 1]
        t_chunk_start = time.time()
        out_chunk = propagator.propagate(x_curr, t_chunk, env_chunk)
        out[i0 : i1 + 1] = out_chunk
        x_curr = out_chunk[-1]
        elapsed = time.time() - t_chunk_start
        print(f"{tag} step {i1+1}/{nt} chunk_s={elapsed:.1f}")
        i0 = i1
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run UT sigma-point shard for job arrays.")
    parser.add_argument("--scenario", required=True, help="module:function returning scenario dict")
    parser.add_argument("--outdir", default="results/ut_shards", help="output directory")
    parser.add_argument("--shard", type=int, default=None, help="shard index (0-based)")
    parser.add_argument("--shards", type=int, default=None, help="total shards")
    parser.add_argument("--ut_alpha", type=float, default=1.0)
    parser.add_argument("--ut_beta", type=float, default=2.0)
    parser.add_argument("--ut_kappa", type=float, default=0.0)
    parser.add_argument("--chunk_steps", type=int, default=0, help="chunk size in steps for progress")
    args = parser.parse_args()

    shard, shards = shard_from_env()
    if args.shard is not None:
        shard = args.shard
    if args.shards is not None:
        shards = args.shards
    if shard < 0 or shards <= 0:
        raise RuntimeError("provide shard/shards or set SLURM_ARRAY_TASK_ID/COUNT")

    scenario = parse_scenario(args.scenario)()
    for key in ("propagator", "t_grid", "env", "x0", "P0"):
        if key not in scenario:
            raise RuntimeError(f"scenario missing key '{key}'")

    propagator = scenario["propagator"]
    t_grid = np.asarray(scenario["t_grid"], dtype=float)
    env = scenario["env"]
    x0 = np.asarray(scenario["x0"], dtype=float)
    P0 = np.asarray(scenario["P0"], dtype=float)
    freeze_attitude = bool(scenario.get("freeze_attitude", False))

    sigmas, w0m, w0c, wi = generate_sigma_points(x0, P0, args.ut_alpha, args.ut_beta, args.ut_kappa)
    if freeze_attitude:
        sigmas[:, 6:10] = x0[6:10]
        sigmas[:, 10:13] = x0[10:13]
    else:
        normalize_quaternions(sigmas, x0)

    n_sig = sigmas.shape[0]
    start = (n_sig * shard) // shards
    end = (n_sig * (shard + 1)) // shards
    if start >= end:
        print(f"[ut] rank {shard}/{shards} empty slice (sigmas={n_sig}), skipping")
        return

    sigma_slice = sigmas[start:end]
    t0 = time.time()
    print(
        f"[ut] rank {shard}/{shards} start={start} end={end} "
        f"sigmas={end-start} steps={len(t_grid)}"
    )

    chunk_steps = int(args.chunk_steps)
    if chunk_steps > 0:
        out = propagate_chunked(
            propagator,
            sigma_slice,
            t_grid,
            env,
            chunk_steps,
            f"[ut] rank {shard}/{shards}",
        )
    else:
        out = propagator.propagate(sigma_slice, t_grid, env)

    weights_m = np.full(end - start, wi, dtype=float)
    weights_c = np.full(end - start, wi, dtype=float)
    if start == 0:
        weights_m[0] = w0m
        weights_c[0] = w0c

    sum_wm_x = np.tensordot(weights_m, out, axes=(0, 1))
    sum_wc_x = np.tensordot(weights_c, out, axes=(0, 1))
    sum_wc_xxt = np.einsum("i,tin,tim->tnm", weights_c, out, out)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"ut_rank{shard:04d}.npz"
    np.savez_compressed(
        out_path,
        sum_wm_x=sum_wm_x,
        sum_wc_x=sum_wc_x,
        sum_wc_xxt=sum_wc_xxt,
        t_grid=t_grid,
        shard=shard,
        shards=shards,
        w0m=w0m,
        w0c=w0c,
        wi=wi,
    )
    elapsed = time.time() - t0
    print(f"[ut] rank {shard}/{shards} done in {elapsed:.1f}s")
    print(f"Saved UT shard {shard}/{shards} to {out_path}")


if __name__ == "__main__":
    main()
