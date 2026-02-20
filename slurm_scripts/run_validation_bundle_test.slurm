#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=vleo_validation_test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=vleo_validation_test_%j.out
#SBATCH --error=vleo_validation_test_%j.err

set -euo pipefail

echo "[vleo_validation_test] $(date -u +%Y-%m-%dT%H:%M:%SZ) start job=${SLURM_JOB_ID:-na}"
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "${SLURM_SUBMIT_DIR}"
fi
echo "[vleo_validation_test] cwd=$(pwd)"

module purge
module load GCCcore/13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Ninja/1.12.1-GCCcore-13.3.0

echo "[vleo_validation_test] modules loaded"
source venv/bin/activate

echo "[vleo_validation_test] venv activated: $(python -c 'import sys; print(sys.executable)')"
python -c 'import vleo_uq; print("[vleo_validation_test] vleo_uq OK")'

export PYTHONNOUSERSITE=1
export PYTHONPATH="$(pwd):$(pwd)/python"
export MPLBACKEND=Agg

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

export PYTHON_BIN="${PYTHON_BIN:-$(pwd)/venv/bin/python3}"
OUTDIR_DEFAULT="results/reports/validation_bundle_test_$(date -u +%Y%m%d_%H%M%S)"
export OUTDIR="${OUTDIR:-${OUTDIR_DEFAULT}}"

CMD=("${PYTHON_BIN}" scripts/run_validation_bundle.py --outdir "${OUTDIR}")

echo "[vleo_validation_test] cmd: ${CMD[*]}"
"${CMD[@]}"

echo "[vleo_validation_test] $(date -u +%Y-%m-%dT%H:%M:%SZ) done"
