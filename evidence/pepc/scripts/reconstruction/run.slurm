#!/bin/bash
#SBATCH --job-name=pepc_cds_tree
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=pepc_cds_tree_%j.out
#SBATCH --error=pepc_cds_tree_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
cd "$PROJECT_ROOT"
python3 -u pepc_cds_rebuild_v1/run.py
