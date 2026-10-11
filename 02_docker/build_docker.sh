#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
docker build -f 02_docker/Dockerfile -t aspdac2026-mpl2:archive .
