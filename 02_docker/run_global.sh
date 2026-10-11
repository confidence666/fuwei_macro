#!/usr/bin/env bash
set -euo pipefail
DESIGN=${1:-bp_be}
cd "$(dirname "$0")/.."
docker run --rm \
  -e ASPDAC_MPL2_ARCHIVE_DIR=/workspace/archive \
  -v "$PWD":/workspace/archive \
  aspdac2026-mpl2:archive \
  -no_init "/workspace/archive/03_benchmarks/$DESIGN/run_global.tcl"
