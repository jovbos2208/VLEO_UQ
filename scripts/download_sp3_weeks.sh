#!/usr/bin/env bash
set -euo pipefail

# Downloads IGS Rapid (igr) and Final (igs) SP3 files for a quiet week (2025-01-24..30)
# and a storm week (2024-05-10..16) from ftp://igs.ign.fr/pub/igs/products/.

ROOT_URL="ftp://igs.ign.fr/pub/igs/products"
OUTDIR="data/gnss"

quiet_dir="${OUTDIR}/quiet_2025w2350_2351"
storm_dir="${OUTDIR}/storm_2024w2313_2314"

mkdir -p "${quiet_dir}" "${storm_dir}"

download_file() {
  local week="$1"
  local dir="$2"
  local file="$3"
  if ! wget -c "${ROOT_URL}/${week}/${file}" -P "${dir}"; then
    echo "[warn] missing: ${week}/${file}" >&2
    return 1
  fi
}

download_week() {
  local week="$1"
  shift
  local dir="$1"
  shift
  local failed=0
  for f in "$@"; do
    if ! download_file "${week}" "${dir}" "${f}"; then
      failed=1
    fi
  done
  if [[ "${failed}" -ne 0 ]]; then
    echo "[warn] some files missing in week ${week}" >&2
  fi
}

download_long() {
  local week="$1"
  local year="$2"
  local doy="$3"
  local dir="$4"
  local prefix="$5"  # IGS0OPSRAP or IGS0OPSFIN
  local fname="${prefix}_${year}${doy}0000_01D_15M_ORB.SP3.gz"
  if ! wget -c "${ROOT_URL}/${week}/${fname}" -P "${dir}"; then
    echo "[warn] missing long file: ${week}/${fname}" >&2
    return 1
  fi
}

# Quiet week: 2025-01-24..30
download_week 2350 "${quiet_dir}" \
  igr23505.sp3.Z igr23506.sp3.Z \
  igs23505.sp3.Z igs23506.sp3.Z

download_week 2351 "${quiet_dir}" \
  igr23510.sp3.Z igr23511.sp3.Z igr23512.sp3.Z igr23513.sp3.Z igr23514.sp3.Z \
  igs23510.sp3.Z igs23511.sp3.Z igs23512.sp3.Z igs23513.sp3.Z igs23514.sp3.Z

# Try long-name FIN/RAP files (quiet week)
download_long 2350 2025 024 "${quiet_dir}" IGS0OPSRAP
download_long 2350 2025 024 "${quiet_dir}" IGS0OPSFIN
download_long 2350 2025 025 "${quiet_dir}" IGS0OPSRAP
download_long 2350 2025 025 "${quiet_dir}" IGS0OPSFIN
download_long 2351 2025 026 "${quiet_dir}" IGS0OPSRAP
download_long 2351 2025 026 "${quiet_dir}" IGS0OPSFIN
download_long 2351 2025 027 "${quiet_dir}" IGS0OPSRAP
download_long 2351 2025 027 "${quiet_dir}" IGS0OPSFIN
download_long 2351 2025 028 "${quiet_dir}" IGS0OPSRAP
download_long 2351 2025 028 "${quiet_dir}" IGS0OPSFIN
download_long 2351 2025 029 "${quiet_dir}" IGS0OPSRAP
download_long 2351 2025 029 "${quiet_dir}" IGS0OPSFIN
download_long 2351 2025 030 "${quiet_dir}" IGS0OPSRAP
download_long 2351 2025 030 "${quiet_dir}" IGS0OPSFIN

# Storm week: 2024-05-10..16
download_week 2313 "${storm_dir}" \
  igr23135.sp3.Z igr23136.sp3.Z \
  igs23135.sp3.Z igs23136.sp3.Z

download_week 2314 "${storm_dir}" \
  igr23140.sp3.Z igr23141.sp3.Z igr23142.sp3.Z igr23143.sp3.Z igr23144.sp3.Z \
  igs23140.sp3.Z igs23141.sp3.Z igs23142.sp3.Z igs23143.sp3.Z igs23144.sp3.Z

# Try long-name FIN/RAP files (storm week)
download_long 2313 2024 131 "${storm_dir}" IGS0OPSRAP
download_long 2313 2024 131 "${storm_dir}" IGS0OPSFIN
download_long 2313 2024 132 "${storm_dir}" IGS0OPSRAP
download_long 2313 2024 132 "${storm_dir}" IGS0OPSFIN
download_long 2314 2024 133 "${storm_dir}" IGS0OPSRAP
download_long 2314 2024 133 "${storm_dir}" IGS0OPSFIN
download_long 2314 2024 134 "${storm_dir}" IGS0OPSRAP
download_long 2314 2024 134 "${storm_dir}" IGS0OPSFIN
download_long 2314 2024 135 "${storm_dir}" IGS0OPSRAP
download_long 2314 2024 135 "${storm_dir}" IGS0OPSFIN
download_long 2314 2024 136 "${storm_dir}" IGS0OPSRAP
download_long 2314 2024 136 "${storm_dir}" IGS0OPSFIN
download_long 2314 2024 137 "${storm_dir}" IGS0OPSRAP
download_long 2314 2024 137 "${storm_dir}" IGS0OPSFIN

echo "Done. Files saved under ${OUTDIR}/quiet_2025w2350_2351 and ${OUTDIR}/storm_2024w2313_2314"
