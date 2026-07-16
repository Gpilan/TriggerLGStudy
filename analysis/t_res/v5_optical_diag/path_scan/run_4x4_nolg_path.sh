#!/usr/bin/env bash
# 4x4 v5 no-LG optical path diag (100 ev). Requires CBDsim_PROTO_NO_LG=1 at runtime.
set -eo pipefail

TRIGGER="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="$TRIGGER/analysis/t_res/v5_optical_diag/path_scan"
BUILD="$TRIGGER/build"
CBD="$BUILD/CBDsim/CBDsim"
MACRO="run_diag_100ev.mac"
TAG="4x4_nolg"
OUTBASE="$OUT/path_${TAG}"

ts() { date '+%Y-%m-%d %H:%M:%S KST'; }

if [[ -s "${OUTBASE}_100ev.path_hist.txt" ]]; then
  echo "[$(ts)] skip — ${OUTBASE}_100ev.path_hist.txt exists"
  exit 0
fi

echo "[$(ts)] rebuild CBDsim (no-LG geometry in same binary)..."
(cd "$TRIGGER" && source envset.sh && cd build && cmake --build . -j 4 --target CBDsim)

export CBDsim_PROTO_NO_LG=1
export CBDsim_OPTICAL_DIAG=1
export CBDsim_OPTICAL_DIAG_HIST=1
export CBDsim_OPTICAL_DIAG_OUT="${OUTBASE}_100ev.txt"
export CBDsim_OPTICAL_DIAG_HIST_OUT="${OUTBASE}_100ev.path_hist.txt"

echo "=== [$(ts)] START 4x4 no-LG path diag (100 ev) ==="
t0=$SECONDS
(cd "$BUILD/CBDsim" && "$CBD" "$MACRO" 0 "$OUTBASE")
echo "=== [$(ts)] DONE 4x4 no-LG path diag ($((SECONDS - t0))s) ==="

unset CBDsim_PROTO_NO_LG

echo "[$(ts)] plot (includes 4x4 LG vs no-LG overlay if both hist files exist)"
(cd "$TRIGGER" && source envset.sh && python3 analysis/plot_optical_path_scan.py --scan-dir "$OUT")
