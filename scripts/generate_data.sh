#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: generate_data.sh <config> <outdir>

Environment:
  VLEO_START_UTC   (required if VLEO_USE_ENV_SOURCES=1)
  VLEO_USE_ENV_SOURCES=1|0
  VLEO_OMNI_PATH, VLEO_HWM14_LIB, VLEO_HWM14_DATA, VLEO_ENV_INTERP
  MC_TASKS (default 256), MC_TPN (default 64)
EOF
}

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

CONFIG="$1"
OUTDIR="$2"
NODES="${SLURM_JOB_NUM_NODES:-4}"
TPN="${SLURM_NTASKS_PER_NODE:-128}"
MC_TASKS="${MC_TASKS:-$((NODES * TPN / 2))}"
MC_TPN="${MC_TPN:-$((MC_TASKS / NODES))}"

mkdir -p "$OUTDIR"

echo "[gen] phase1: DET+STM"
srun -N 1 -n 1 python scripts/run_det_stm_phase.py \
  --config "$CONFIG" \
  --outdir "$OUTDIR" \
  ${VLEO_USE_ENV_SOURCES:+--use_env_sources} \
  --start_utc "${VLEO_START_UTC:-}" \
  --omni_path "${VLEO_OMNI_PATH:-data/space_weather/omni/omni2_all_years.dat}" \
  --hwm14_lib "${VLEO_HWM14_LIB:-data/hwm14/libhwm14.so}" \
  --hwm14_data "${VLEO_HWM14_DATA:-data/hwm14}" \
  --env_interpolation "${VLEO_ENV_INTERP:-nearest}"

echo "[gen] phase1: MC+UT shards"
export VLEO_CASE_CONFIG="$CONFIG"

# Mission MC/UT
export VLEO_CASE_NAME=mission
srun -N ${NODES} -n ${MC_TASKS} --ntasks-per-node=${MC_TPN} --kill-on-bad-exit=0 --wait=0 \
  python scripts/run_mc_job_array.py \
  --scenario scripts.mc_scenarios:mission_scenario \
  --outdir "${OUTDIR}/mc_shards/mission" \
  --chunk_steps 50
python scripts/merge_mc_stats.py --indir "${OUTDIR}/mc_shards/mission" --out "${OUTDIR}/mission_mc_stats.npz"
srun -N 1 -n 37 --kill-on-bad-exit=0 --wait=0 \
  python scripts/run_ut_job_array.py \
  --scenario scripts.ut_scenarios:mission_ut_scenario \
  --outdir "${OUTDIR}/ut_shards/mission" \
  --chunk_steps 50
python scripts/merge_ut_shards.py --indir "${OUTDIR}/ut_shards/mission" --out "${OUTDIR}/mission_ut_merged.npz"

# Formation MC/UT
export VLEO_CASE_NAME=formation
for IDX in 0 1 2; do
  export VLEO_FORMATION_INDEX=$IDX
  OBJ=$(printf "obj_%02d" $((IDX+1)))
  srun -N ${NODES} -n ${MC_TASKS} --ntasks-per-node=${MC_TPN} --kill-on-bad-exit=0 --wait=0 \
    python scripts/run_mc_job_array.py \
    --scenario scripts.mc_scenarios:formation_scenario \
    --outdir "${OUTDIR}/mc_shards/formation_${OBJ}" \
    --chunk_steps 50
  python scripts/merge_mc_stats.py --indir "${OUTDIR}/mc_shards/formation_${OBJ}" --out "${OUTDIR}/formation/obj_$(printf "%02d" $((IDX+1)))_mc_stats.npz"
  srun -N 1 -n 37 --kill-on-bad-exit=0 --wait=0 \
    python scripts/run_ut_job_array.py \
    --scenario scripts.ut_scenarios:formation_ut_scenario \
    --outdir "${OUTDIR}/ut_shards/formation_${OBJ}" \
    --chunk_steps 50
  python scripts/merge_ut_shards.py --indir "${OUTDIR}/ut_shards/formation_${OBJ}" --out "${OUTDIR}/formation/obj_$(printf "%02d" $((IDX+1)))_ut_merged.npz"
done

# Attitude MC/UT
export VLEO_CASE_NAME=attitude
srun -N ${NODES} -n ${MC_TASKS} --ntasks-per-node=${MC_TPN} --kill-on-bad-exit=0 --wait=0 \
  python scripts/run_mc_job_array.py \
  --scenario scripts.mc_scenarios:attitude_scenario \
  --outdir "${OUTDIR}/mc_shards/attitude" \
  --chunk_steps 50
python scripts/merge_mc_stats.py --indir "${OUTDIR}/mc_shards/attitude" --out "${OUTDIR}/attitude_mc_stats.npz"
srun -N 1 -n 37 --kill-on-bad-exit=0 --wait=0 \
  python scripts/run_ut_job_array.py \
  --scenario scripts.ut_scenarios:attitude_ut_scenario \
  --outdir "${OUTDIR}/ut_shards/attitude" \
  --chunk_steps 50
python scripts/merge_ut_shards.py --indir "${OUTDIR}/ut_shards/attitude" --out "${OUTDIR}/attitude_ut_merged.npz"

echo "[gen] done"
