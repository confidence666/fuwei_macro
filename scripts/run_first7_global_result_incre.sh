#!/usr/bin/env bash
set -euo pipefail
ARCHIVE=${ARCHIVE:-/workspace/archive/ASPDAC2026_Timing_Driven_MacroPlacer_mpl2_archive}
OPENROAD=${OPENROAD:-/workspace/OpenROAD-flow-scripts/tools/OpenROAD/build/src/openroad}
export ASPDAC_MPL2_ARCHIVE_DIR="$ARCHIVE"
MAX_JOBS=${MAX_JOBS:-3}
TIMEOUT_SECONDS=${TIMEOUT_SECONDS:-14400}
BENCHMARKS=(bp bp_be ariane81 swerv_wrapper swerv_wrapper43 mor1kx or1200)

has_fatal_error_log() {
  local log=$1
  grep -Eiv '\[ERROR GPL-0305\]' "$log" | grep -Eiq '(\[ERROR|^Error:|Signal [0-9]+|segmentation fault|command failed|timeout)'
}

run_incre() {
  local design=$1
  local tcl="$ARCHIVE/03_benchmarks/$design/run_global_result_incre.tcl"
  local outdir="$ARCHIVE/04_results/$design/global_incre"
  local log="$outdir/run_global_incre.log"
  local expected_def="$outdir/${design}_global_incre.def"
  local expected_odb="$outdir/${design}_global_incre.odb"
  mkdir -p "$outdir"
  rm -f "$outdir"/*
  echo "[RUN] $design global_result_incre"
  local rc=0
  if command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SECONDS" "$OPENROAD" -no_init "$tcl" > "$log" 2>&1 || rc=$?
  else
    "$OPENROAD" -no_init "$tcl" > "$log" 2>&1 || rc=$?
  fi
  if [ "$rc" -ne 0 ]; then
    echo "FAILED global_result_incre $design: openroad exit code $rc" > "$outdir/FAILED"
    return 0
  fi
  if has_fatal_error_log "$log"; then
    echo "FAILED global_result_incre $design: fatal error pattern found in log" > "$outdir/FAILED"
    return 0
  fi
  if [ ! -s "$expected_def" ] || [ ! -s "$expected_odb" ]; then
    echo "FAILED global_result_incre $design: missing expected DEF/ODB" > "$outdir/FAILED"
    return 0
  fi
  echo "[OK] $design global_result_incre"
}

active=0
for design in "${BENCHMARKS[@]}"; do
  run_incre "$design" &
  active=$((active + 1))
  if [ "$active" -ge "$MAX_JOBS" ]; then
    wait -n || true
    active=$((active - 1))
  fi
done
wait || true
