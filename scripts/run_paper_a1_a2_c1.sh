#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-8}"
SEED="${SEED:-42}"
START_UTC="${START_UTC:-2025-01-24T00:00:00}"
POD_SP3="${POD_SP3:-data/gnss/quiet/IGS0OPSFIN_20250240000_01D_15M_ORB.SP3}"
POD_EOP="${POD_EOP:-data/eop/EOP-All.txt}"
POD_ESTIMATOR="${POD_ESTIMATOR:-enkf}"
POD_ENKF_MEMBERS="${POD_ENKF_MEMBERS:-128}"
POD_ENKF_INFLATION="${POD_ENKF_INFLATION:-1.02}"
DURATION_SCALE="${DURATION_SCALE:-1.0}"

STAMP="$(date -u +%Y%m%d_%H%M%S)"
RESULTS_DIR="${RESULTS_DIR:-results/paper_a1_a2_c1_${STAMP}}"
CONFIG_PATH="${CONFIG_PATH:-configs/generated/paper_a1_a2_c1_${STAMP}.json}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${RESULTS_DIR}/paper_analysis}"

echo "[paper] building config -> ${CONFIG_PATH}"
"${PYTHON_BIN}" scripts/paper_build_config_a1_a2_c1.py \
  --out "${CONFIG_PATH}" \
  --start_utc "${START_UTC}" \
  --duration_scale "${DURATION_SCALE}" \
  --pod_estimator "${POD_ESTIMATOR}" \
  --pod_enkf_members "${POD_ENKF_MEMBERS}" \
  --pod_enkf_inflation "${POD_ENKF_INFLATION}"

echo "[paper] running case studies -> ${RESULTS_DIR}"
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/run_case_studies.py \
  --config "${CONFIG_PATH}" \
  --outdir "${RESULTS_DIR}" \
  --threads "${THREADS}" \
  --seed "${SEED}" \
  --pod_sp3 "${POD_SP3}" \
  --pod_eop "${POD_EOP}" \
  --pod_use_sat_clock \
  --start_utc "${START_UTC}" \
  --use_env_sources \
  --omni_path "data/space_weather/omni/omni2_all_years.dat" \
  --hwm14_lib "data/hwm14/libhwm14.so" \
  --hwm14_data "data/hwm14" \
  --env_interpolation "nearest"

echo "[paper] analyzing outputs -> ${ANALYSIS_DIR}"
MPLBACKEND=Agg "${PYTHON_BIN}" scripts/paper_analyze_uq_a1_a2_c1.py \
  --config "${CONFIG_PATH}" \
  --results_dir "${RESULTS_DIR}" \
  --outdir "${ANALYSIS_DIR}"

echo "[paper] done"
echo "  config:  ${CONFIG_PATH}"
echo "  results: ${RESULTS_DIR}"
echo "  analysis:${ANALYSIS_DIR}"
