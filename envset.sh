source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

# 저장소 루트 (이 스크립트가 Trigger/ 에 있을 때)
export TRIGGER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# CBDsim 통계 ROOT 모음 (cmake 후 build/CBDsim/stats)
export CBDsim_STATS="${TRIGGER_ROOT}/build/CBDsim/stats"
