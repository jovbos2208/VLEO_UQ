#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=vleo_att_full
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=vleo_att_full_%j.out
#SBATCH --error=vleo_att_full_%j.err

set -euo pipefail

echo "[vleo_att_full] $(date -u +%Y-%m-%dT%H:%M:%SZ) start job=${SLURM_JOB_ID:-na}"
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "${SLURM_SUBMIT_DIR}"
fi
echo "[vleo_att_full] cwd=$(pwd)"

module purge
module load GCCcore/13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Ninja/1.12.1-GCCcore-13.3.0

echo "[vleo_att_full] modules loaded"
source venv/bin/activate

echo "[vleo_att_full] venv activated: $(python -c 'import sys; print(sys.executable)')"
python -c 'import vleo_uq; print("[vleo_att_full] vleo_uq OK")'

export PYTHONNOUSERSITE=1
export PYTHONPATH="$(pwd):$(pwd)/python"
export MPLBACKEND=Agg

export THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${THREADS}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

export PYTHON_BIN="${PYTHON_BIN:-$(pwd)/venv/bin/python3}"
export CATALOG="${CATALOG:-scenario_catalog.md}"
export PHASE="${PHASE:-att}"
export SEED="${SEED:-42}"
export DEFAULT_PARTICLES="${DEFAULT_PARTICLES:-256}"
export MAX_PARTICLES="${MAX_PARTICLES:-}"
export DURATION_SCALE="${DURATION_SCALE:-1.0}"
export DT_FLOOR="${DT_FLOOR:-1.0}"
export DISABLE_POD="${DISABLE_POD:-0}"
export USE_ENV_SOURCES="${USE_ENV_SOURCES:-1}"
export START_UTC="${START_UTC:-2025-01-24T00:00:00}"
export STRICT_CAMPAIGN="${STRICT_CAMPAIGN:-1}"
export IDS="${IDS:-}"

OUTDIR_DEFAULT="results/run_att_full_$(date -u +%Y%m%d_%H%M%S)"
export OUTDIR="${OUTDIR:-${OUTDIR_DEFAULT}}"

CMD=(
  "${PYTHON_BIN}" scripts/run_phase_campaign.py
  --catalog "${CATALOG}"
  --phase "${PHASE}"
  --outdir "${OUTDIR}"
  --threads "${THREADS}"
  --seed "${SEED}"
  --default_particles "${DEFAULT_PARTICLES}"
  --duration_scale "${DURATION_SCALE}"
  --dt_floor "${DT_FLOOR}"
)

if [ -n "${MAX_PARTICLES}" ]; then
  CMD+=(--max_particles "${MAX_PARTICLES}")
fi
if [ "${DISABLE_POD}" = "1" ]; then
  CMD+=(--disable_pod)
fi
if [ "${USE_ENV_SOURCES}" = "1" ]; then
  CMD+=(--use_env_sources --start_utc "${START_UTC}")
fi
if [ "${STRICT_CAMPAIGN}" = "1" ]; then
  CMD+=(--strict_campaign)
fi
if [ -n "${IDS}" ]; then
  CMD+=(--ids "${IDS}")
fi

echo "[vleo_att_full] cmd: ${CMD[*]}"
"${CMD[@]}"

echo "[vleo_att_full] $(date -u +%Y-%m-%dT%H:%M:%SZ) done"
