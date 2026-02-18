#!/usr/bin/env bash
set -euo pipefail

WHEELHOUSE="${1:-wheelhouse}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if [ ! -d "$WHEELHOUSE" ]; then
  echo "Wheelhouse directory not found: $WHEELHOUSE" >&2
  exit 1
fi

"$PYTHON_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" \
  pip setuptools wheel

"$PYTHON_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" \
  numpy scipy pydantic xarray h5py SALib pyyaml pymsis matplotlib scikit-build-core pybind11 cmake ninja

"$PYTHON_BIN" -m pip install --no-index --find-links "$WHEELHOUSE" \
  -e . --no-build-isolation

echo "Offline install complete"
