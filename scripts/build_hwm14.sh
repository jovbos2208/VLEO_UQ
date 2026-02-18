#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 /path/to/HWM14 [output_dir]" >&2
  exit 1
fi

SRC_DIR="$1"
OUT_DIR="${2:-data/hwm14}"

if [ ! -d "$SRC_DIR/src" ]; then
  echo "HWM14 source directory must contain src/ with Fortran files" >&2
  exit 1
fi

if ! command -v gfortran >/dev/null 2>&1; then
  echo "gfortran is required to build HWM14" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

SRC_ORDER=(
  "hwm.f90"
  "alf.f90"
  "gd2qdc.f90"
  "qwm.f90"
  "dwm.f90"
  "hwm14_utilities.f90"
  "hwm14.f90"
)

for src in "${SRC_ORDER[@]}"; do
  gfortran -O2 -fPIC -J "$BUILD_DIR" -c "$SRC_DIR/src/$src" -o "$BUILD_DIR/${src%.f90}.o"
done

gfortran -shared -o "$OUT_DIR/libhwm14.so" "$BUILD_DIR"/*.o

cp "$SRC_DIR"/dwm07b104i.dat "$OUT_DIR"/
cp "$SRC_DIR"/gd2qd.dat "$OUT_DIR"/
cp "$SRC_DIR"/hwm123114.bin "$OUT_DIR"/

echo "HWM14 library and data staged in $OUT_DIR"
