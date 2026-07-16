#!/usr/bin/env bash
# V5 steps ②~⑤ overnight pipeline. Logs to pipeline.log and updates RUNLOG.md
set -eo pipefail

TRIGGER="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="$TRIGGER/analysis/t_res/v5_optical_diag"
PROTO="$TRIGGER/CBDsim/src/CBDsimDetectorConstructionProto.cc"
BUILD="$TRIGGER/build"
CBD="$BUILD/CBDsim/CBDsim"
MACRO="run_diag_100ev.mac"
LOG="$OUT/pipeline.log"
RUNLOG="$OUT/RUNLOG.md"
COMPARE="$TRIGGER/analysis/compare_v5_optical_diag.py"

exec > >(tee -a "$LOG") 2>&1

ts() { date '+%Y-%m-%d %H:%M:%S KST'; }
append_runlog() { echo -e "\n$1" >> "$RUNLOG"; }

run_diag() {
  local tag="$1"
  local outbase="$OUT/diag_${tag}"
  export CBDsim_OPTICAL_DIAG=1
  export CBDsim_OPTICAL_DIAG_OUT="${outbase}_100ev.txt"
  echo "=== [$(ts)] START diag: $tag ==="
  local t0=$SECONDS
  (cd "$BUILD/CBDsim" && "$CBD" "$MACRO" 0 "$outbase")
  local elapsed=$((SECONDS - t0))
  echo "=== [$(ts)] DONE diag: $tag (${elapsed}s) -> ${outbase}_100ev.txt ==="
  append_runlog "### run $tag\n- finished: $(ts), wall ${elapsed}s\n- file: \`${outbase}_100ev.txt\`\n\`\`\`\n$(cat "${outbase}_100ev.txt")\n\`\`\`"
}

patch_proto() {
  local nphi="$1" gap_val="$2" gap_unit="$3" hx="$4"
  sed -i \
    -e "s/constexpr G4int kNPhi = [0-9]*;/constexpr G4int kNPhi = ${nphi};/" \
    -e "s/constexpr G4double kLGZGap = [^;]*;/constexpr G4double kLGZGap = ${gap_val} * ${gap_unit};/" \
    -e "s/constexpr G4double kHxWide = [0-9.]* \* mm;/constexpr G4double kHxWide = ${hx} * mm;/" \
    "$PROTO"
  echo "[$(ts)] patched kNPhi=$nphi kLGZGap=${gap_val}${gap_unit} kHxWide=${hx}mm"
}

rebuild() {
  echo "[$(ts)] rebuild CBDsim..."
  (cd "$TRIGGER" && source envset.sh && cd build && cmake --build . -j 4 --target CBDsim)
}

compare_pair() {
  local a="$1" b="$2" label="$3"
  echo "=== COMPARE $label ==="
  python3 "$COMPARE" "$OUT/${a}" "$OUT/${b}" | tee -a "$LOG"
  append_runlog "\n**Compare $label** (\`$a\` vs \`$b\`):\n\`\`\`\n$(python3 "$COMPARE" "$OUT/${a}" "$OUT/${b}")\n\`\`\`"
}

source "$TRIGGER/envset.sh" || { echo "FATAL: envset.sh failed"; exit 1; }

echo "======== V5 pipeline start $(ts) ========"

# Wait for manual/baseline nphi256 if still running
while pgrep -f "CBDsim run_diag_100ev.mac.*nphi256" >/dev/null 2>&1; do
  echo "[$(ts)] waiting for existing nphi256 job..."
  sleep 60
done

# ② nphi256 (skip if output exists and non-empty)
if [[ ! -s "$OUT/diag_4x4_nphi256_100ev.txt" ]]; then
  patch_proto 256 0.0 mm 20.0
  rebuild
  run_diag "4x4_nphi256"
else
  echo "[$(ts)] skip ② — diag_4x4_nphi256_100ev.txt exists"
fi
compare_pair "diag_4x4_baseline_100ev.txt" "diag_4x4_nphi256_100ev.txt" "① vs ②"

# ③ gap scan (keep nphi=256, hx=20)
for gap in 5 10; do
  tag="4x4_nphi256_gap${gap}um"
  if [[ ! -s "$OUT/diag_${tag}_100ev.txt" ]]; then
    patch_proto 256 "$gap" um 20.0
    rebuild
    run_diag "$tag"
  else
    echo "[$(ts)] skip ③ gap${gap}um — exists"
  fi
  compare_pair "diag_4x4_nphi256_100ev.txt" "diag_${tag}_100ev.txt" "② vs ③ gap${gap}um"
done

# pick best gap by lowest scint_to_world (python one-liner)
BEST_GAP_VAL="0.0"
BEST_GAP_UNIT="mm"
BEST_FILE="diag_4x4_nphi256_100ev.txt"
for f in diag_4x4_nphi256_100ev.txt diag_4x4_nphi256_gap5um_100ev.txt diag_4x4_nphi256_gap10um_100ev.txt; do
  [[ -f "$OUT/$f" ]] || continue
  sw=$(grep '^scint_to_world ' "$OUT/$f" | awk '{print $2}')
  best_sw=$(grep '^scint_to_world ' "$OUT/$BEST_FILE" | awk '{print $2}')
  if awk "BEGIN{exit !($sw < $best_sw)}"; then
    BEST_FILE="$f"
    if [[ "$f" == *gap5um* ]]; then BEST_GAP_VAL=5; BEST_GAP_UNIT=um; fi
    if [[ "$f" == *gap10um* ]]; then BEST_GAP_VAL=10; BEST_GAP_UNIT=um; fi
  fi
done
echo "[$(ts)] best 4x4 config file: $BEST_FILE"
append_runlog "\n## best 4x4 after ③\n- file: \`$BEST_FILE\`"

# ④ 1x1 diag (kHxWide=5 mm)
if [[ ! -s "$OUT/diag_1x1_nphi256_100ev.txt" ]]; then
  patch_proto 256 "$BEST_GAP_VAL" "$BEST_GAP_UNIT" 5.0
  rebuild
  run_diag "1x1_nphi256"
else
  echo "[$(ts)] skip ④ 1x1 — exists"
fi
compare_pair "diag_1x1_nphi256_100ev.txt" "$BEST_FILE" "1x1 vs best 4x4"

# ④ MPV quick compare on 100 ev roots (if plot script exists)
MPV_SCRIPT="$TRIGGER/analysis/plot_event_optical_mpv_matrix_root.py"
if [[ -f "$MPV_SCRIPT" ]]; then
  mkdir -p "$OUT/mpv_100ev"
  for tag in 1x1_nphi256 4x4_nphi256; do
    r="$OUT/diag_${tag}_100_0.root"
    [[ -f "$r" ]] && python3 "$MPV_SCRIPT" "$r" --outdir "$OUT/mpv_100ev" --tag "$tag" 2>/dev/null || true
  done
fi

# ⑤ decision: scint_to_world vs scint_to_lg ratio
read_sw() { grep "^scint_to_world " "$OUT/$1" | awk '{print $2}'; }
read_sl() { grep "^scint_to_lg " "$OUT/$1" | awk '{print $2}'; }
SW=$(read_sw "$BEST_FILE")
SL=$(read_sl "$BEST_FILE")
RATIO=$(awk "BEGIN{printf \"%.1f\", $SW/$SL}")
echo "[$(ts)] ⑤ check: scint_to_world/scint_to_lg = $RATIO"
append_runlog "\n## ⑤ foil decision\n- scint_to_world/scint_to_lg = $RATIO on \`$BEST_FILE\`"

if awk "BEGIN{exit !($SW/$SL > 5)}"; then
  append_runlog "- **ACTION**: scint↔world still dominant → foil hollow shell deferred to manual PR (geometry subtraction non-trivial overnight). Document as open item."
  echo "[$(ts)] ⑤ foil shell NOT auto-implemented — ratio $RATIO > 5, needs daytime review"
else
  append_runlog "- **SKIP foil**: ratio <= 5 after ②③"
fi

# final summary csv
SUMMARY="$OUT/v5_compare_summary.csv"
echo "tag,scint_to_world,world_to_scint,scint_to_lg,lg_to_world,lg_to_sipm,sipm_detect,path_scint_mm_per_event,path_lg_mm_per_event" > "$SUMMARY"
for f in "$OUT"/diag_*_100ev.txt; do
  [[ -f "$f" ]] || continue
  tag=$(basename "$f" _100ev.txt)
  get() { grep "^$1 " "$f" | awk '{print $2}'; }
  echo "$(basename $tag),$(get scint_to_world),$(get world_to_scint),$(get scint_to_lg),$(get lg_to_world),$(get lg_to_sipm),$(get sipm_detect),$(get path_scint_mm_per_event),$(get path_lg_mm_per_event)" >> "$SUMMARY"
done

append_runlog "\n## pipeline finished $(ts)\n- summary: \`v5_compare_summary.csv\`\n- full log: \`pipeline.log\`"
echo "======== V5 pipeline done $(ts) ========"
