#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=vleo_ATT_baseline
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --output=ATT_baseline_%j.out
#SBATCH --error=ATT_baseline_%j.err

set -euo pipefail

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) start job=${SLURM_JOB_ID:-na}"

# Use srun as task launcher only; avoid cluster-specific MPI plugin mismatches.
export SLURM_MPI_TYPE="${VLEO_SLURM_MPI_TYPE:-none}"
echo "[job] SLURM_MPI_TYPE=${SLURM_MPI_TYPE}"

module purge
module load GCCcore/13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Ninja/1.12.1-GCCcore-13.3.0

source venv/bin/activate

export PYTHONNOUSERSITE=1
export PYTHONPATH="$(pwd)"
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
VENV_PY="$(pwd)/venv/bin/python3"

echo "[job] python=$(which python3)"
echo "[job] venv_py=${VENV_PY}"
"${VENV_PY}" -c "import sys; print(sys.executable)"
"${VENV_PY}" -c "import yaml; print('PyYAML OK', yaml.__version__)"
"${VENV_PY}" -c "import vleo_uq; print('vleo_uq OK')"

# Offline environment sources
export VLEO_USE_ENV_SOURCES=1
export VLEO_START_UTC="${VLEO_START_UTC:-2025-01-24T00:00:00}"
export VLEO_OMNI_PATH="${VLEO_OMNI_PATH:-data/space_weather/omni/omni2_all_years.dat}"
export VLEO_HWM14_LIB="${VLEO_HWM14_LIB:-data/hwm14/libhwm14.so}"
export VLEO_HWM14_DATA="${VLEO_HWM14_DATA:-data/hwm14}"
export VLEO_ENV_INTERP="${VLEO_ENV_INTERP:-nearest}"

# Baseline maneuver set (pitch removed: A8 + A10).
ATT_BASE_IDS="${ATT_BASE_IDS:-ATT_A8_ROLL_ROTATION_BASELINE,ATT_A10_YAW_ROTATION_BASELINE}"
ATT_SWEEP_WING_ANGLES_DEG="${ATT_SWEEP_WING_ANGLES_DEG:-0,5,10,15,20,25,30,35,40,45,50,55,60}"
ATT_RECOMMENDED_ABS_MAX_DEG="${ATT_RECOMMENDED_ABS_MAX_DEG:-60}"
ATT_ALLOW_OUT_OF_RECOMMENDED_RANGE="${ATT_ALLOW_OUT_OF_RECOMMENDED_RANGE:-0}"
ATT_SWEEP_DURATION_S="${ATT_SWEEP_DURATION_S:-120}"
ATT_SWEEP_PARTICLES="${ATT_SWEEP_PARTICLES:-1000}"
ATT_SWEEP_DT_S="${ATT_SWEEP_DT_S:-1.0}"
ETA_MAX_DEG="${ETA_MAX_DEG:-110}"
ETA_RATE_MAX_DEG_S="${ETA_RATE_MAX_DEG_S:-5}"
ETA_ACCEL_MAX_DEG_S2="${ETA_ACCEL_MAX_DEG_S2:-1}"

ENABLE_PLOT="${ENABLE_PLOT:-1}"
PLOT_ALL_STATES="${PLOT_ALL_STATES:-0}"
MC_CHUNK_STEPS="${MC_CHUNK_STEPS:-50}"
ATT_MC_SAFE_MODE="${ATT_MC_SAFE_MODE:-1}"
CATALOG_DT_FLOOR="${CATALOG_DT_FLOOR:-1.0}"
export VLEO_MC_SAFE_MODE="${ATT_MC_SAFE_MODE}"

RUN_ID="${RUN_ID:-run_att_baseline_$(date -u +%Y%m%d_%H%M%S)}"
RESULTS_DIR="results/${RUN_ID}"
BASE_CFG="configs/generated/catalog_ATT_baseline_base.json"
CFG="configs/generated/catalog_ATT_baseline_sweep.json"
mkdir -p "${RESULTS_DIR}" configs/generated

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) build base config ids=${ATT_BASE_IDS}"
"${VENV_PY}" scripts/catalog_to_case_config.py \
  --catalog scenario_catalog.md \
  --out "${BASE_CFG}" \
  --ids "${ATT_BASE_IDS}" \
  --start_utc "${VLEO_START_UTC}" \
  --dt_floor "${CATALOG_DT_FLOOR}"

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) build wing sweep angles=${ATT_SWEEP_WING_ANGLES_DEG} duration_s=${ATT_SWEEP_DURATION_S} particles=${ATT_SWEEP_PARTICLES}"
ATT_SWEEP_CMD=(
  "${VENV_PY}" scripts/build_att_wing_sweep_config.py
  --in_config "${BASE_CFG}"
  --out "${CFG}"
  --angles_deg "${ATT_SWEEP_WING_ANGLES_DEG}"
  --recommended_abs_max_deg "${ATT_RECOMMENDED_ABS_MAX_DEG}"
  --duration_s "${ATT_SWEEP_DURATION_S}"
  --particles "${ATT_SWEEP_PARTICLES}"
  --dt_s "${ATT_SWEEP_DT_S}"
  --eta_max_deg "${ETA_MAX_DEG}"
  --eta_rate_max_deg_s "${ETA_RATE_MAX_DEG_S}"
  --eta_accel_max_deg_s2 "${ETA_ACCEL_MAX_DEG_S2}"
)
if [ "${ATT_ALLOW_OUT_OF_RECOMMENDED_RANGE}" = "1" ]; then
  ATT_SWEEP_CMD+=(--allow_out_of_recommended_range)
fi
"${ATT_SWEEP_CMD[@]}"

mapfile -t ATT_SCENARIOS < <("${VENV_PY}" - <<'PY'
import json
cfg=json.load(open("configs/generated/catalog_ATT_baseline_sweep.json","r",encoding="utf-8"))
for s in cfg["scenarios"]:
    print(s["name"])
PY
)

N_CASES="${#ATT_SCENARIOS[@]}"
if [ "${N_CASES}" -lt 1 ]; then
  echo "[job] no ATT baseline scenarios resolved from config"
  exit 1
fi

TOTAL_TASKS="${SLURM_NTASKS:-$(( ${SLURM_JOB_NUM_NODES:-1} * 128 ))}"
MC_TASKS_PER_CASE_DEFAULT=$(( TOTAL_TASKS / N_CASES ))
if [ "${MC_TASKS_PER_CASE_DEFAULT}" -lt 1 ]; then MC_TASKS_PER_CASE_DEFAULT=1; fi
MC_TASKS_PER_CASE="${MC_TASKS_PER_CASE:-$MC_TASKS_PER_CASE_DEFAULT}"
if [ "${MC_TASKS_PER_CASE}" -lt 1 ]; then MC_TASKS_PER_CASE=1; fi

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) layout n_cases=${N_CASES} total_tasks=${TOTAL_TASKS} mc_tasks_per_case=${MC_TASKS_PER_CASE}"

export VLEO_CASE_CONFIG="${CFG}"

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) phase1 MC shards (scenario-parallel)"
MC_PARALLEL=$(( TOTAL_TASKS / MC_TASKS_PER_CASE ))
if [ "${MC_PARALLEL}" -lt 1 ]; then MC_PARALLEL=1; fi
PIDS=()
RUNNING=0
for CASE in "${ATT_SCENARIOS[@]}"; do
  (
    PARTICLES=$("${VENV_PY}" - <<PY
import json
cfg=json.load(open("${CFG}","r",encoding="utf-8"))
s=next(x for x in cfg["scenarios"] if x["name"]=="${CASE}")
print(int(s.get("particles", 2)))
PY
)
    TASKS="${MC_TASKS_PER_CASE}"
    if [ "${TASKS}" -gt "${PARTICLES}" ]; then TASKS="${PARTICLES}"; fi
    if [ "${TASKS}" -lt 1 ]; then TASKS=1; fi
    export VLEO_CASE_NAME="${CASE}"
    echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) phase1 ${CASE} MC tasks=${TASKS} particles=${PARTICLES}"
    srun --exclusive -n "${TASKS}" --kill-on-bad-exit=1 \
      --output="${RESULTS_DIR}/mc_${CASE}_%j_%t.log" \
      --error="${RESULTS_DIR}/mc_${CASE}_%j_%t.err" \
      "${VENV_PY}" scripts/run_mc_job_array.py \
      --scenario scripts.mc_scenarios:attitude_scenario \
      --outdir "${RESULTS_DIR}/mc_shards/${CASE}" \
      --chunk_steps "${MC_CHUNK_STEPS}"

    "${VENV_PY}" scripts/merge_mc_stats.py \
      --indir "${RESULTS_DIR}/mc_shards/${CASE}" \
      --out "${RESULTS_DIR}/${CASE}_mc_stats.npz"
  ) &
  PIDS+=($!)
  RUNNING=$((RUNNING + 1))
  if [ "${RUNNING}" -ge "${MC_PARALLEL}" ]; then
    wait -n
    RUNNING=$((RUNNING - 1))
  fi
done
FAIL=0
for PID in "${PIDS[@]}"; do
  wait "${PID}" || FAIL=1
done
if [ "${FAIL}" -ne 0 ]; then
  echo "[job] phase1 MC failed"
  exit 1
fi

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) phase2 MC-only postprocess"
for CASE in "${ATT_SCENARIOS[@]}"; do
  echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) phase2 ${CASE}"
  CMD=(
    "${VENV_PY}" scripts/postprocess_mc_only.py
    --config "${CFG}"
    --outdir "${RESULTS_DIR}"
    --scenario "${CASE}"
  )
  if [ "${ENABLE_PLOT}" = "1" ]; then
    CMD+=(--plot)
  fi
  if [ "${PLOT_ALL_STATES}" = "1" ]; then
    CMD+=(--plot_all_states)
  fi
  "${CMD[@]}"
done

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) phase3 sweep table"
"${VENV_PY}" scripts/summarize_att_wing_sweep.py \
  --config "${CFG}" \
  --results_dir "${RESULTS_DIR}" \
  --out_csv "${RESULTS_DIR}/att_wing_sweep_times.csv" \
  --out_md "${RESULTS_DIR}/att_wing_sweep_times.md"

echo "[job] $(date -u +%Y-%m-%dT%H:%M:%SZ) done"
echo "[job] results: ${RESULTS_DIR}"
