#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_scenarios.sh <ensemble_size> [threads] [config_in] [config_out]

Defaults:
  ensemble_size: 1000
  threads: 4
  config_in: configs/case_studies.json
  config_out: configs/case_studies_<ensemble_size>.json
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

ENSEMBLE_SIZE="${1:-1000}"
THREADS="${2:-4}"
CONFIG_IN="${3:-configs/case_studies.json}"
CONFIG_OUT="${4:-configs/case_studies_${ENSEMBLE_SIZE}.json}"

if [[ ! -f "$CONFIG_IN" ]]; then
  echo "Config not found: $CONFIG_IN" >&2
  exit 1
fi

python3 - <<PY
import json
from pathlib import Path

cfg_in = Path("${CONFIG_IN}")
cfg_out = Path("${CONFIG_OUT}")
data = json.loads(cfg_in.read_text())
for scenario in data.get("scenarios", []):
    scenario["particles"] = int("${ENSEMBLE_SIZE}")
cfg_out.write_text(json.dumps(data, indent=2))
print(f"[run_scenarios] wrote {cfg_out}")
PY

if [[ -f "/home/jovan/venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source /home/jovan/venv/bin/activate
fi

python scripts/run_case_studies.py --config "$CONFIG_OUT" --threads "$THREADS" --plot --pod_uq
