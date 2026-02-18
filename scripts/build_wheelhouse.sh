#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-wheelhouse}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

REQS=(
  numpy
  scipy
  pydantic
  xarray
  h5py
  SALib
  pyyaml
  pymsis
  matplotlib
  scikit-build-core
  pybind11
  setuptools
  wheel
  cmake
  ninja
  pip
)

mkdir -p "$OUT_DIR"
"$PYTHON_BIN" -m pip download -d "$OUT_DIR" "${REQS[@]}"

echo "Wheelhouse ready in $OUT_DIR"
