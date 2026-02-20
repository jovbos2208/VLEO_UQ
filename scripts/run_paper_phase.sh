#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "usage: $0 <phase_a|phase_b|phase_c|att|mis|orb|aero>"
  exit 2
fi

PHASE="$1"
case "${PHASE}" in
  phase_a|phase_b|phase_c|att|mis|orb|aero) ;;
  *)
    echo "unsupported phase: ${PHASE}"
    exit 2
    ;;
esac

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-8}"
SEED="${SEED:-42}"
START_UTC="${START_UTC:-2025-01-24T00:00:00}"
OMNI_PATH="${OMNI_PATH:-data/space_weather/omni/omni2_all_years.dat}"
HWM14_LIB="${HWM14_LIB:-data/hwm14/libhwm14.so}"
HWM14_DATA="${HWM14_DATA:-data/hwm14}"
ENV_INTERPOLATION="${ENV_INTERPOLATION:-nearest}"
POD_SP3="${POD_SP3:-data/gnss/quiet/IGS0OPSFIN_20250240000_01D_15M_ORB.SP3}"
POD_EOP="${POD_EOP:-data/eop/EOP-All.txt}"
POD_ESTIMATOR="${POD_ESTIMATOR:-enkf}"
POD_ENKF_MEMBERS="${POD_ENKF_MEMBERS:-128}"
POD_ENKF_INFLATION="${POD_ENKF_INFLATION:-1.02}"
DEFAULT_PARTICLES="${DEFAULT_PARTICLES:-256}"
MAX_PARTICLES="${MAX_PARTICLES:-}"
DURATION_SCALE="${DURATION_SCALE:-1.0}"
DT_FLOOR="${DT_FLOOR:-1.0}"
USE_ENV_SOURCES="${USE_ENV_SOURCES:-1}"
DISABLE_POD="${DISABLE_POD:-0}"
POD_USE_SAT_CLOCK="${POD_USE_SAT_CLOCK:-1}"
POD_ENKF_USE_CARRIER="${POD_ENKF_USE_CARRIER:-1}"
PLOT_CASE_TRACES="${PLOT_CASE_TRACES:-0}"
PLOT_ALL_STATES="${PLOT_ALL_STATES:-0}"
MAKE_PHASE_PLOTS="${MAKE_PHASE_PLOTS:-1}"

STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUTDIR="${OUTDIR:-results/paper_${PHASE}_${STAMP}}"
CONFIG_OUT="${CONFIG_OUT:-configs/generated/paper_${PHASE}_${STAMP}.json}"

CMD=(
  "${PYTHON_BIN}" scripts/run_paper_phase.py
  --phase "${PHASE}"
  --threads "${THREADS}"
  --seed "${SEED}"
  --outdir "${OUTDIR}"
  --config_out "${CONFIG_OUT}"
  --default_particles "${DEFAULT_PARTICLES}"
  --duration_scale "${DURATION_SCALE}"
  --dt_floor "${DT_FLOOR}"
  --pod_estimator "${POD_ESTIMATOR}"
  --pod_enkf_members "${POD_ENKF_MEMBERS}"
  --pod_enkf_inflation "${POD_ENKF_INFLATION}"
  --pod_eop "${POD_EOP}"
)

if [ -n "${MAX_PARTICLES}" ]; then
  CMD+=(--max_particles "${MAX_PARTICLES}")
fi
if [ "${USE_ENV_SOURCES}" = "1" ]; then
  CMD+=(
    --use_env_sources
    --start_utc "${START_UTC}"
    --omni_path "${OMNI_PATH}"
    --hwm14_lib "${HWM14_LIB}"
    --hwm14_data "${HWM14_DATA}"
    --env_interpolation "${ENV_INTERPOLATION}"
  )
fi
if [ "${DISABLE_POD}" = "1" ]; then
  CMD+=(--disable_pod)
fi
if [ -n "${POD_SP3}" ] && [ -f "${POD_SP3}" ]; then
  CMD+=(--pod_sp3 "${POD_SP3}")
fi
if [ "${POD_USE_SAT_CLOCK}" = "1" ]; then
  CMD+=(--pod_use_sat_clock)
fi
if [ "${POD_ENKF_USE_CARRIER}" = "1" ]; then
  CMD+=(--pod_enkf_use_carrier)
fi
if [ "${PLOT_CASE_TRACES}" = "1" ]; then
  CMD+=(--plot)
fi
if [ "${PLOT_ALL_STATES}" = "1" ]; then
  CMD+=(--plot_all_states)
fi
if [ "${MAKE_PHASE_PLOTS}" = "1" ]; then
  CMD+=(--make_phase_plots)
fi

echo "[paper-phase] phase=${PHASE} outdir=${OUTDIR}"
echo "[paper-phase] cmd: ${CMD[*]}"
"${CMD[@]}"
