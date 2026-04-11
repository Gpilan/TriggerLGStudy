#!/usr/bin/env bash
# verify_root.C 를 절대 경로로 실행 (어느 디렉터리에서든 동작)
#
#   chmod +x CBDsim/verify_root_run.sh   # 최초 1회
#   cd ~/Trigger && source envset.sh
#   ./CBDsim/verify_root_run.sh build/CBDsim/stats/verify_smoke_0.root
#
# 또는 (stats 아래 파일은 basename 만으로도 검색):
#   ./CBDsim/verify_root_run.sh verify_smoke_0.root
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# 소스: .../CBDsim  |  빌드 복사본: .../build/CBDsim  → 저장소 루트는 각각 .. / ../..
if [[ "$SCRIPT_DIR" == *"/build/CBDsim" ]]; then
  TRIGGER_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
else
  TRIGGER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
MACRO="$SCRIPT_DIR/verify_root.C"
ROOTIO_SO="${ROOTIO_LIB:-$TRIGGER_ROOT/build/rootIO/librootIO.so}"

resolve_rootfile() {
  local f="${1:-verify_smoke_0.root}"
  if [[ -f "$f" ]]; then
    readlink -f "$f"
    return
  fi
  if [[ -f "$TRIGGER_ROOT/$f" ]]; then
    readlink -f "$TRIGGER_ROOT/$f"
    return
  fi
  if [[ -f "$TRIGGER_ROOT/build/CBDsim/stats/$f" ]]; then
    readlink -f "$TRIGGER_ROOT/build/CBDsim/stats/$f"
    return
  fi
  if [[ -f "$TRIGGER_ROOT/build/CBDsim/$f" ]]; then
    readlink -f "$TRIGGER_ROOT/build/CBDsim/$f"
    return
  fi
  echo ""
}

if [[ ! -f "$MACRO" ]]; then
  echo "ERROR: missing $MACRO" >&2
  exit 1
fi
if [[ ! -f "$ROOTIO_SO" ]]; then
  echo "ERROR: librootIO not found: $ROOTIO_SO" >&2
  echo "  export ROOTIO_LIB=/path/to/Trigger/build/rootIO/librootIO.so" >&2
  exit 1
fi

RF="$(resolve_rootfile "${1:-verify_smoke_0.root}")"
if [[ -z "$RF" ]] || [[ ! -f "$RF" ]]; then
  echo "ERROR: ROOT file not found: ${1:-verify_smoke_0.root}" >&2
  echo "  Pass basename (stats/ 또는 build/CBDsim) 또는 전체 경로." >&2
  exit 1
fi

export ROOTIO_LIB="$ROOTIO_SO"
echo "[verify_root_run] ROOTIO_LIB=$ROOTIO_LIB"
echo "[verify_root_run] FILE=$RF"

exec root -l -b -e ".x $MACRO(\"$RF\")" -q
