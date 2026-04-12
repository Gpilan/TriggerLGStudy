source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

# 저장소 루트 (이 스크립트가 Trigger/ 에 있을 때)
export TRIGGER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# plot_trigger_timing.py 첫 번째 입력 디렉터리 (미설정 시 analysis/t_res/data 와 동일 역할)
# 스크립트는 추가로 build/CBDsim/stats, build/CBDsim 도 자동 검색
export CBDsim_STATS="${TRIGGER_ROOT}/build/CBDsim/stats"
