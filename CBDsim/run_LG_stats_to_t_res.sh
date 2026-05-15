#!/usr/bin/env bash
# LG 지오메트리 빌드된 CBDsim 으로, 에너지별 3000이벤트씩 시뮬 후
# analysis/t_res/data/LG_stats_{에너지}GeV_0.root 생성
# (NoLG_stats_*GeV_0.root 과 이름 규칙만 LG_ 로 맞춤)
#
# 사용 전:
#   source envset.sh   # Geant4
#   cd build && cmake .. && make -j
#
# 실행 (저장소 루트에서):
#   bash CBDsim/run_LG_stats_to_t_res.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT="${REPO}/analysis/t_res/data"
EXE="${REPO}/build/CBDsim/CBDsim"
MAC="${REPO}/CBDsim"

# 필요 시 수정 (예: 90 100 110)
ENERGIES=(80 100 120)
EVENTS=3000
SEED=0

mkdir -p "${OUT}"

if [[ ! -x "${EXE}" ]]; then
  echo "CBDsim 실행 파일 없음: ${EXE}  →  build 먼저 하세요." >&2
  exit 1
fi

export LD_LIBRARY_PATH="${REPO}/build/rootIO:${LD_LIBRARY_PATH:-}"

for E in "${ENERGIES[@]}"; do
  MACFILE="${MAC}/run_ele_${E}gev_3000.mac"
  if [[ ! -f "${MACFILE}" ]]; then
    echo "매크로 없음: ${MACFILE}" >&2
    exit 1
  fi
  STEM="${OUT}/LG_stats_${E}GeV"
  OUTROOT="${STEM}_${SEED}.root"
  echo "========================================"
  echo " ${E} GeV  /  ${EVENTS} events  →  ${OUTROOT}"
  echo "========================================"
  "${EXE}" "${MACFILE}" "${SEED}" "${STEM}"
done

echo "완료. 출력: ${OUT}/LG_stats_*GeV_${SEED}.root"
