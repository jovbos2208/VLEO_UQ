from __future__ import annotations

import functools
import subprocess
from pathlib import Path
from typing import Any, Mapping


def _as_list(val: Any) -> list[str]:
    if val is None:
        return []
    if isinstance(val, (list, tuple, set)):
        return [str(x) for x in val]
    return [str(val)]


@functools.lru_cache(maxsize=1)
def detect_git_hash() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo_root,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return out or "unknown"
    except Exception:
        return "unknown"


def build_run_metadata(scenario: Mapping[str, Any] | None, seed: int | None) -> dict[str, Any]:
    s = scenario or {}
    scenario_name = str(s.get("name", "unknown"))
    scenario_id = str(s.get("catalog_scenario_id", scenario_name))
    toggles = {
        "env": _as_list(s.get("catalog_env_toggles")),
        "gsi": _as_list(s.get("catalog_gsi_toggles")),
        "sensor": _as_list(s.get("catalog_sensor_toggles")),
        "ground": _as_list(s.get("catalog_ground_toggles")),
        "est": _as_list(s.get("catalog_est_toggles")),
        "control": _as_list(s.get("catalog_control_toggles")),
    }
    return {
        "git_hash": detect_git_hash(),
        "scenario_id": scenario_id,
        "seed": None if seed is None else int(seed),
        "toggles": toggles,
    }
