#!/usr/bin/env bash
# V5 geometry: trigger size 1x1..4x4 optical photon path-length / fate scan (100 ev each).
set -eo pipefail

TRIGGER="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="$TRIGGER/analysis/t_res/v5_optical_diag/path_scan"
PROTO="$TRIGGER/CBDsim/src/CBDsimDetectorConstructionProto.cc"
BUILD="$TRIGGER/build"
CBD="$BUILD/CBDsim/CBDsim"
MACRO="run_diag_100ev.mac"
LOG="$OUT/path_scan.log"
PLOT="$TRIGGER/analysis/plot_optical_path_scan.py"

mkdir -p "$OUT"
exec > >(tee -a "$LOG") 2>&1

ts() { date '+%Y-%m-%d %H:%M:%S KST'; }

# size_label hx_mm
SIZES=(
  "1x1 5"
  "1p5x1p5 7.5"
  "2x2 10"
  "2p5x2p5 12.5"
  "3x3 15"
  "3p5x3p5 17.5"
  "4x4 20"
)

patch_hx() {
  local hx="$1"
  sed -i -e "s/constexpr G4double kHxWide = [0-9.]* \* mm;/constexpr G4double kHxWide = ${hx} * mm;/" "$PROTO"
  echo "[$(ts)] patched kHxWide=${hx} mm"
}

rebuild() {
  echo "[$(ts)] rebuild CBDsim..."
  (cd "$TRIGGER" && source envset.sh && cd build && cmake --build . -j 4 --target CBDsim)
}

run_one() {
  local tag="$1"
  local outbase="$OUT/path_${tag}"
  export CBDsim_OPTICAL_DIAG=1
  export CBDsim_OPTICAL_DIAG_HIST=1
  export CBDsim_OPTICAL_DIAG_OUT="${outbase}_100ev.txt"
  export CBDsim_OPTICAL_DIAG_HIST_OUT="${outbase}_100ev.path_hist.txt"
  echo "=== [$(ts)] START path scan: $tag ==="
  local t0=$SECONDS
  (cd "$BUILD/CBDsim" && "$CBD" "$MACRO" 0 "$outbase")
  local elapsed=$((SECONDS - t0))
  echo "=== [$(ts)] DONE path scan: $tag (${elapsed}s) ==="
}

echo "=== path scan pipeline start $(ts) ==="
rebuild

for entry in "${SIZES[@]}"; do
  tag="${entry%% *}"
  hx="${entry##* }"
  if [[ -s "$OUT/path_${tag}_100ev.path_hist.txt" ]]; then
    echo "[$(ts)] skip $tag — hist exists"
    continue
  fi
  patch_hx "$hx"
  rebuild
  run_one "$tag"
done

# restore 4x4 default
patch_hx "20"
rebuild
echo "[$(ts)] restored kHxWide=20 mm (4x4)"

echo "=== plot summary $(ts) ==="
(cd "$TRIGGER" && source envset.sh && python3 "$PLOT" --scan-dir "$OUT")

echo "=== path scan pipeline finished $(ts) ==="
