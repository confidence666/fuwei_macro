#!/usr/bin/env bash
set -euo pipefail
ARCHIVE=${ARCHIVE:-/workspace/archive/ASPDAC2026_Timing_Driven_MacroPlacer_mpl2_archive}
OPENROAD=${OPENROAD:-/workspace/OpenROAD-flow-scripts/tools/OpenROAD/build/src/openroad}
export ASPDAC_MPL2_ARCHIVE_DIR="$ARCHIVE"
MAX_JOBS=${MAX_JOBS:-3}
TIMEOUT_SECONDS=${TIMEOUT_SECONDS:-14400}
BENCHMARKS=(bp bp_be ariane81 swerv_wrapper swerv_wrapper43 mor1kx or1200 ariane133 bp_multi bp_multi57 isa_npu dft68 ethernet bp_fe ariane136)

has_fatal_error_log() {
  local log=$1
  grep -Eiv '\[ERROR GPL-0305\]' "$log" | grep -Eiq '(\[ERROR|^Error:|Signal [0-9]+|segmentation fault|command failed|timeout)'
}

run_global() {
  local design=$1
  local tcl="$ARCHIVE/03_benchmarks/$design/run_global.tcl"
  local outdir="$ARCHIVE/04_results/$design/global"
  local log="$outdir/run_global.log"
  local expected_def="$outdir/${design}_global.def"
  local expected_odb="$outdir/${design}_global.odb"
  mkdir -p "$outdir"
  rm -f "$outdir"/*
  echo "[RUN] $design global"
  local rc=0
  if command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SECONDS" "$OPENROAD" -no_init "$tcl" > "$log" 2>&1 || rc=$?
  else
    "$OPENROAD" -no_init "$tcl" > "$log" 2>&1 || rc=$?
  fi
  if [ "$rc" -ne 0 ]; then
    echo "FAILED global $design: openroad exit code $rc" > "$outdir/FAILED"
    return 0
  fi
  if has_fatal_error_log "$log"; then
    echo "FAILED global $design: fatal error pattern found in log" > "$outdir/FAILED"
    return 0
  fi
  if [ ! -s "$expected_def" ] || [ ! -s "$expected_odb" ]; then
    echo "FAILED global $design: missing expected DEF/ODB" > "$outdir/FAILED"
    return 0
  fi
  echo "[OK] $design global"
}

active=0
for design in "${BENCHMARKS[@]}"; do
  run_global "$design" &
  active=$((active + 1))
  if [ "$active" -ge "$MAX_JOBS" ]; then
    wait -n || true
    active=$((active - 1))
  fi
done
wait || true
python3 "$ARCHIVE/scripts/collect_metrics.py" "$ARCHIVE" || true
