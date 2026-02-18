#!/usr/bin/env bash
set -euo pipefail

CATALOG="${1:-scenario_catalog.md}"
CONFIG_OUT="${2:-configs/case_studies_catalog.json}"
IDS="${3:-}"
DURATION_SCALE="${DURATION_SCALE:-1.0}"
DEFAULT_PARTICLES="${DEFAULT_PARTICLES:-256}"
THREADS="${THREADS:-4}"
DT_FLOOR="${DT_FLOOR:-1.0}"
START_UTC="${START_UTC:-}"
NO_ENV_SOURCES="${NO_ENV_SOURCES:-0}"

CMD=(python scripts/catalog_to_case_config.py
  --catalog "$CATALOG"
  --out "$CONFIG_OUT"
  --duration_scale "$DURATION_SCALE"
  --default_particles "$DEFAULT_PARTICLES"
  --dt_floor "$DT_FLOOR")

if [[ -n "$IDS" ]]; then
  CMD+=(--ids "$IDS")
fi
if [[ -n "$START_UTC" ]]; then
  CMD+=(--start_utc "$START_UTC")
fi
if [[ "$NO_ENV_SOURCES" == "1" ]]; then
  CMD+=(--no_env_sources)
fi

"${CMD[@]}"

python scripts/run_case_studies.py \
  --config "$CONFIG_OUT" \
  --threads "$THREADS"
