#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-/tmp/vleo_ci_smoke}"

python scripts/run_case_studies.py \
  --case mission \
  --mission_duration_s 60 \
  --mission_dt_s 30 \
  --mission_particles 4 \
  --outdir "$OUTDIR"

export VLEO_CI_SMOKE_OUTDIR="$OUTDIR"
python - <<'PY'
import json
import os
from pathlib import Path

summary_path = Path(os.environ["VLEO_CI_SMOKE_OUTDIR"]) / "mission" / "summary.json"
if not summary_path.exists():
    raise SystemExit(f"missing summary: {summary_path}")
summary = json.loads(summary_path.read_text(encoding="utf-8"))
meta = summary.get("metadata", {})
required_meta = ["git_hash", "scenario_id", "seed", "toggles"]
missing = [k for k in required_meta if k not in meta]
if missing:
    raise SystemExit(f"missing metadata keys: {missing}")
if not isinstance(meta.get("toggles"), dict):
    raise SystemExit("metadata.toggles must be an object")
print("[ci-smoke] summary metadata keys present")
PY

python -m unittest tests/test_catalog_converter.py tests/test_propagator_regression.py -v

echo "[ci-smoke] deterministic smoke passed"
