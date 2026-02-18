from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge MC ensemble shards.")
    parser.add_argument("--indir", default="results/mc_shards", help="input directory")
    parser.add_argument("--out", default="results/ensemble_merged.npz", help="output npz path")
    args = parser.parse_args()

    indir = Path(args.indir)
    paths = sorted(indir.glob("ensemble_rank*.npz"))
    if not paths:
        raise RuntimeError(f"no shard files found in {indir}")

    shards = []
    for path in paths:
        with np.load(path) as data:
            shards.append(
                {
                    "path": path,
                    "X": data["X"],
                    "start": int(data["start"]),
                    "end": int(data["end"]),
                    "shard": int(data["shard"]),
                    "shards": int(data["shards"]),
                    "t_grid": data["t_grid"],
                }
            )

    shards.sort(key=lambda s: s["start"])
    t_grid = shards[0]["t_grid"]
    nx = shards[0]["X"].shape[-1]
    n_total = max(s["end"] for s in shards)

    X_merged = np.zeros((len(t_grid), n_total, nx), dtype=float)
    for s in shards:
        X_merged[:, s["start"] : s["end"], :] = s["X"]

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, X=X_merged, t_grid=t_grid)
    print(f"Merged {len(shards)} shards into {out_path}")


if __name__ == "__main__":
    main()
