from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge UT sigma-point shards.")
    parser.add_argument("--indir", default="results/ut_shards", help="input directory")
    parser.add_argument("--out", default="results/ut_merged.npz", help="output npz path")
    args = parser.parse_args()

    indir = Path(args.indir)
    paths = sorted(indir.glob("ut_rank*.npz"))
    if not paths:
        raise RuntimeError(f"no shard files found in {indir}")

    shards = []
    for path in paths:
        with np.load(path) as data:
            shards.append(
                {
                    "path": path,
                    "sum_wm_x": data["sum_wm_x"],
                    "sum_wc_x": data["sum_wc_x"],
                    "sum_wc_xxt": data["sum_wc_xxt"],
                    "t_grid": data["t_grid"],
                    "w0m": float(data["w0m"]),
                    "w0c": float(data["w0c"]),
                    "wi": float(data["wi"]),
                    "shards": int(data["shards"]),
                }
            )

    t_grid = shards[0]["t_grid"]
    sum_wm_x = np.zeros_like(shards[0]["sum_wm_x"])
    sum_wc_x = np.zeros_like(shards[0]["sum_wc_x"])
    sum_wc_xxt = np.zeros_like(shards[0]["sum_wc_xxt"])
    for s in shards:
        sum_wm_x += s["sum_wm_x"]
        sum_wc_x += s["sum_wc_x"]
        sum_wc_xxt += s["sum_wc_xxt"]

    n = sum_wm_x.shape[1]
    w0m = shards[0]["w0m"]
    w0c = shards[0]["w0c"]
    wi = shards[0]["wi"]
    sum_wc = w0c + 2.0 * n * wi

    mean = sum_wm_x
    cov = sum_wc_xxt - np.einsum("ti,tj->tij", mean, sum_wc_x) - np.einsum(
        "ti,tj->tij", sum_wc_x, mean
    ) + sum_wc * np.einsum("ti,tj->tij", mean, mean)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, ut_mean=mean, ut_cov=cov, t_grid=t_grid)
    print(f"Merged {len(shards)} UT shards into {out_path}")


if __name__ == "__main__":
    main()
