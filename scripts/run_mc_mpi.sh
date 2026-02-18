#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_mc_mpi.sh <scenario_func> <num_tasks> <outdir>

Example:
  export VLEO_CASE_CONFIG=configs/case_studies_1000.json
  export VLEO_CASE_NAME=mission
  export VLEO_USE_ENV_SOURCES=1
  export VLEO_START_UTC=2024-05-10T00:00:00
  ./scripts/run_mc_mpi.sh scripts.mc_scenarios:mission_scenario 8 results/mc_shards/mission
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 3 ]]; then
  usage
  exit 0
fi

SCENARIO_FUNC="$1"
TASKS="$2"
OUTDIR="$3"

if command -v mpirun >/dev/null 2>&1; then
  mpirun -n "$TASKS" python scripts/run_mc_job_array.py --scenario "$SCENARIO_FUNC" --outdir "$OUTDIR"
elif command -v srun >/dev/null 2>&1; then
  srun -n "$TASKS" python scripts/run_mc_job_array.py --scenario "$SCENARIO_FUNC" --outdir "$OUTDIR"
else
  echo "mpirun or srun not found; please launch via your scheduler." >&2
  exit 1
fi
