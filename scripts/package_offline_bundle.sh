#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-vleo_uq_offline_bundle.tar.gz}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

tar --exclude="results" \
  --exclude="__pycache__" \
  --exclude="*.pyc" \
  -czf "$OUT" \
  -C "$ROOT" \
  CMakeLists.txt pyproject.toml OFFLINE.md Plan.md EnvInputs.md EnvIInputs.md \
  cpp python scripts configs data tests wheelhouse

echo "Wrote offline bundle to $OUT"
