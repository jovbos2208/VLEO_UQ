#!/usr/bin/env python3
"""Proxy Sobol/Morris sensitivity postprocess from per-run summary files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _as_finite_float(value: Any) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    if not np.isfinite(out):
        return None
    return float(out)


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _flatten_uq_draw(draw: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, val in draw.items():
        if str(key) == "seed":
            continue
        if isinstance(val, dict):
            v = _as_finite_float(val.get("value"))
            if v is not None:
                out[str(key)] = v
                continue
        v_scalar = _as_finite_float(val)
        if v_scalar is not None:
            out[str(key)] = v_scalar
    return out


def _extract_payload_p95(summary: dict[str, Any]) -> float | None:
    payload = summary.get("payload_impact")
    if not isinstance(payload, dict):
        return None
    return _as_finite_float(payload.get("final_p95_geolocation_error_m"))


def extract_summary_record(path: Path) -> dict[str, Any] | None:
    summary = _load_json(path)
    case = str(summary.get("case", "")).strip()
    if not case:
        return None
    draw = summary.get("uq_parameter_draw")
    if not isinstance(draw, dict):
        return None
    params = _flatten_uq_draw(draw)
    if not params:
        return None

    qois = {}
    for key in ("ut_mean_err", "stm_mean_err", "ut_cov_rel_err", "stm_cov_rel_err"):
        v = _as_finite_float(summary.get(key))
        if v is not None:
            qois[key] = v
    payload_p95 = _extract_payload_p95(summary)
    if payload_p95 is not None:
        qois["payload_p95_geolocation_error_m"] = payload_p95
    if not qois:
        return None
    return {
        "path": str(path),
        "case": case,
        "params": params,
        "qois": qois,
    }


def _corr_safe(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if x.size != y.size or x.size < 3:
        return 0.0
    sx = float(np.std(x))
    sy = float(np.std(y))
    if sx <= 0.0 or sy <= 0.0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def compute_morris_proxy(X: np.ndarray, y: np.ndarray, names: list[str]) -> dict[str, Any]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    std_y = float(np.std(y))
    rows = []
    for j, name in enumerate(names):
        x = X[:, j]
        order = np.argsort(x)
        xs = x[order]
        ys = y[order]
        dx = np.diff(xs)
        dy = np.diff(ys)
        mask = np.abs(dx) > 1e-12
        if not np.any(mask):
            rows.append({"name": name, "mu_star": 0.0, "sigma": 0.0, "mu_star_norm": 0.0})
            continue
        ee = dy[mask] / dx[mask]
        mu_star = float(np.mean(np.abs(ee)))
        sigma = float(np.std(ee))
        std_x = float(np.std(x))
        mu_star_norm = float(mu_star * std_x / max(std_y, 1e-12)) if std_y > 0.0 else 0.0
        rows.append({"name": name, "mu_star": mu_star, "sigma": sigma, "mu_star_norm": mu_star_norm})
    rows.sort(key=lambda r: float(r["mu_star_norm"]), reverse=True)
    return {"parameters": rows}


def compute_sobol_proxy(X: np.ndarray, y: np.ndarray, names: list[str], seed: int = 0) -> dict[str, Any]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    n, d = X.shape
    xm = np.mean(X, axis=0)
    xs = np.std(X, axis=0)
    xs = np.where(xs > 0.0, xs, 1.0)
    ym = float(np.mean(y))
    ys = float(np.std(y))
    if ys <= 0.0:
        ys = 1.0
    Xz = (X - xm) / xs
    yz = (y - ym) / ys

    corr = np.array([_corr_safe(Xz[:, j], yz) for j in range(d)], dtype=float)
    s1 = np.square(corr)
    s1_sum = float(np.sum(s1))
    if s1_sum > 1.0:
        s1 = s1 / s1_sum

    beta, *_ = np.linalg.lstsq(Xz, yz, rcond=None)
    pred = Xz @ beta
    full_mse = float(np.mean((yz - pred) ** 2))
    rng = np.random.default_rng(int(seed))
    st = np.zeros(d, dtype=float)
    for j in range(d):
        Xp = np.array(Xz, copy=True)
        Xp[:, j] = Xp[rng.permutation(n), j]
        pred_p = Xp @ beta
        mse_p = float(np.mean((yz - pred_p) ** 2))
        st[j] = max(0.0, (mse_p - full_mse) / max(float(np.var(yz)), 1e-12))
    st = np.clip(st, 0.0, 1.0)

    rows = []
    for j, name in enumerate(names):
        rows.append({"name": name, "S1_proxy": float(s1[j]), "ST_proxy": float(st[j])})
    rows.sort(key=lambda r: float(r["ST_proxy"]), reverse=True)
    return {"parameters": rows}


def _group_name(param: str) -> str:
    p = str(param).lower()
    if "mass" in p:
        return "mass"
    if "inertia" in p:
        return "inertia"
    if ("rho" in p) or ("density" in p):
        return "density"
    if "wind" in p:
        return "wind"
    if ("gsi" in p) or ("alpha" in p) or ("accommodation" in p):
        return "gsi"
    return "other"


def grouped_attribution(sobol_proxy: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in sobol_proxy.get("parameters", []):
        name = str(row.get("name", ""))
        s1 = _as_finite_float(row.get("S1_proxy")) or 0.0
        grp = _group_name(name)
        out[grp] = out.get(grp, 0.0) + float(s1)
    total = float(sum(out.values()))
    if total > 0.0:
        for k in list(out.keys()):
            out[k] = float(out[k] / total)
    return dict(sorted(out.items(), key=lambda kv: kv[1], reverse=True))


def analyze_case_records(
    records: list[dict[str, Any]],
    *,
    qoi: str,
    min_samples: int = 12,
    seed: int = 0,
) -> dict[str, Any] | None:
    if len(records) < min_samples:
        return None
    all_params = sorted({k for r in records for k in r["params"].keys()})
    if not all_params:
        return None
    X_rows = []
    y_rows = []
    for rec in records:
        y = _as_finite_float(rec["qois"].get(qoi))
        if y is None:
            continue
        row = []
        ok = True
        for p in all_params:
            v = _as_finite_float(rec["params"].get(p))
            if v is None:
                ok = False
                break
            row.append(v)
        if not ok:
            continue
        X_rows.append(row)
        y_rows.append(y)
    if len(X_rows) < min_samples:
        return None
    X = np.asarray(X_rows, dtype=float)
    y = np.asarray(y_rows, dtype=float)
    keep = np.std(X, axis=0) > 0.0
    X = X[:, keep]
    names = [p for p, k in zip(all_params, keep) if bool(k)]
    if X.shape[1] == 0:
        return None

    morris = compute_morris_proxy(X, y, names)
    sobol = compute_sobol_proxy(X, y, names, seed=seed)
    grouped = grouped_attribution(sobol)
    return {
        "n_samples": int(X.shape[0]),
        "qoi": qoi,
        "morris_proxy": morris,
        "sobol_proxy": sobol,
        "grouped_attribution": grouped,
    }


def discover_summary_paths(results_dirs: list[Path]) -> list[Path]:
    out: list[Path] = []
    for base in results_dirs:
        if not base.exists():
            continue
        out.extend(sorted(base.rglob("summary.json")))
    return out


def build_records(paths: list[Path]) -> list[dict[str, Any]]:
    rows = []
    for path in paths:
        try:
            rec = extract_summary_record(path)
        except Exception:
            rec = None
        if rec is not None:
            rows.append(rec)
    return rows


def analyze_records(
    records: list[dict[str, Any]],
    *,
    qois: list[str],
    min_samples: int = 12,
    seed: int = 0,
) -> dict[str, Any]:
    by_case: dict[str, list[dict[str, Any]]] = {}
    for rec in records:
        by_case.setdefault(rec["case"], []).append(rec)

    cases_out = {}
    for case, rows in sorted(by_case.items()):
        qoi_out = {}
        for qoi in qois:
            result = analyze_case_records(rows, qoi=qoi, min_samples=min_samples, seed=seed)
            if result is not None:
                qoi_out[qoi] = result
        if qoi_out:
            cases_out[case] = qoi_out

    return {
        "cases": cases_out,
        "n_records": int(len(records)),
        "qois": [str(q) for q in qois],
        "min_samples": int(min_samples),
        "notes": [
            "Sobol/Morris values are proxy estimators from random-seed campaign summaries (not Saltelli design exact Sobol).",
            "Grouped attribution is normalized sum of first-order proxy contributions by uncertainty channel.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Proxy Sobol/Morris sensitivity postprocess.")
    parser.add_argument(
        "--results_dir",
        action="append",
        default=["results"],
        help="Results root directory (can be repeated).",
    )
    parser.add_argument(
        "--qoi",
        action="append",
        default=["ut_mean_err", "stm_mean_err", "ut_cov_rel_err", "stm_cov_rel_err", "payload_p95_geolocation_error_m"],
        help="QoI key to analyze (can be repeated).",
    )
    parser.add_argument("--min_samples", type=int, default=12)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", default="results/sensitivity/sobol_morris_summary.json")
    args = parser.parse_args()

    roots = [Path(p) for p in args.results_dir]
    paths = discover_summary_paths(roots)
    records = build_records(paths)
    report = analyze_records(
        records,
        qois=[str(q) for q in args.qoi],
        min_samples=max(4, int(args.min_samples)),
        seed=int(args.seed),
    )
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[sobol-morris] summaries={len(paths)} usable={len(records)} wrote {out_path}")


if __name__ == "__main__":
    main()
