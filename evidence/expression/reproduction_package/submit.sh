#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT}"
export PROJECT_ROOT
export PACKAGE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python3 -c 'import numpy,matplotlib'
cd "$PROJECT_ROOT"
sbatch --parsable --export=ALL "$PACKAGE_DIR/job.slurm"
