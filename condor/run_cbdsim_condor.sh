#!/bin/bash
set -eo pipefail

# Usage:
#   run_cbdsim_condor.sh <macro_path_or_name> <seed> <prefix>
# Example:
#   run_cbdsim_condor.sh /u/user/rmsvlf000/Trigger/build/CBDsim/run_ele_60gev.mac 0 v2_60GeV_e-_noLG_3000

REPO_ROOT="/u/user/rmsvlf000/Trigger"
BUILD_DIR="${REPO_ROOT}/build/CBDsim"
ENVSET="${REPO_ROOT}/envset.sh"

MACRO="${1:-${BUILD_DIR}/run_ele_60gev.mac}"
SEED="${2:-0}"
PREFIX="${3:-v2_60GeV_e-_noLG_3000}"

# Allow both absolute path and short macro name.
if [[ "${MACRO}" != /* ]]; then
  MACRO="${BUILD_DIR}/${MACRO}"
fi

echo "[INFO] host=$(hostname)"
echo "[INFO] start=$(date -Is)"
echo "[INFO] macro=${MACRO} seed=${SEED} prefix=${PREFIX}"

if [[ ! -f "${MACRO}" ]]; then
  echo "[ERROR] macro file not found: ${MACRO}" >&2
  exit 4
fi

if [[ ! -f "${ENVSET}" ]]; then
  echo "[ERROR] envset.sh not found at ${ENVSET}" >&2
  exit 2
fi

# Load the same LCG stack used by interactive runs.
# shellcheck disable=SC1090
source "${ENVSET}"

LIBSTDCXX="$(g++ -print-file-name=libstdc++.so.6 || true)"
if [[ -z "${LIBSTDCXX}" || "${LIBSTDCXX}" == "libstdc++.so.6" || ! -f "${LIBSTDCXX}" ]]; then
  echo "[ERROR] failed to resolve LCG libstdc++.so.6 via g++" >&2
  exit 3
fi

LIBSTDCXX_DIR="$(dirname "${LIBSTDCXX}")"
export LD_LIBRARY_PATH="${LIBSTDCXX_DIR}:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="${LIBSTDCXX}${LD_PRELOAD:+:${LD_PRELOAD}}"

echo "[INFO] g++=$(command -v g++)"
echo "[INFO] libstdcxx=${LIBSTDCXX}"
echo "[INFO] ld_library_path_head=$(echo "${LD_LIBRARY_PATH}" | awk -F: '{print $1}')"

cd "${BUILD_DIR}"
LDD_STDCPP="$(ldd ./CBDsim | rg 'libstdc\\+\\+' || true)"
if [[ -n "${LDD_STDCPP}" ]]; then
  echo "[INFO] ldd(libstdc++)=${LDD_STDCPP}"
else
  echo "[WARN] ldd(libstdc++) not found in filtered output; ldd head follows:"
  ldd ./CBDsim | sed -n '1,20p'
fi

if [[ "${CBDSIM_CONDOR_DIAG_ONLY:-0}" == "1" ]]; then
  echo "[INFO] CBDSIM_CONDOR_DIAG_ONLY=1, skip CBDsim run"
  exit 0
fi

./CBDsim "${MACRO}" "${SEED}" "${PREFIX}"

echo "[INFO] done=$(date -Is)"
