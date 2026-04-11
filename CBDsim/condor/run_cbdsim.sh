#!/usr/bin/env bash
# CBDsim 배치 1 job = 시드 $(Process) 로 한 번 실행.
# 사용: condor_submit 에서 arguments = $(Process) 로 넘김.
#
# 환경변수(선택):
#   MACRO          기본 run_ele.mac  (build/CBDsim 에 복사된 매크로 이름)
#   OUT_SUBDIR     stats 아래 출력 하위 경로, 기본 condor
#   EXTRA_ARGS     CBDsim 에 넘길 추가 인자 없음 — 필요 시 스크립트 수정
#
set -euo pipefail

JOB_ID="${1:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# .../Trigger/CBDsim/condor -> 저장소 루트 Trigger
TRIGGER_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_CBD="${TRIGGER_ROOT}/build/CBDsim"
MACRO="${MACRO:-run_ele.mac}"
OUT_SUBDIR="${OUT_SUBDIR:-condor}"

if [[ ! -x "${BUILD_CBD}/CBDsim" ]]; then
  echo "ERROR: ${BUILD_CBD}/CBDsim 없거나 실행 불가. cmake 빌드 후 제출하세요." >&2
  exit 1
fi

# LCG/ROOT (envset.sh 가 있으면)
if [[ -f "${TRIGGER_ROOT}/envset.sh" ]]; then
  # shellcheck source=/dev/null
  source "${TRIGGER_ROOT}/envset.sh"
fi
export LD_LIBRARY_PATH="${TRIGGER_ROOT}/build/rootIO:${LD_LIBRARY_PATH:-}"

mkdir -p "${BUILD_CBD}/stats/${OUT_SUBDIR}"

cd "${BUILD_CBD}"
# RunAction: 파일명 = prefix + "_" + seed + ".root"
OUT_PREFIX="stats/${OUT_SUBDIR}/job_${JOB_ID}"
exec ./CBDsim "${MACRO}" "${JOB_ID}" "${OUT_PREFIX}"
