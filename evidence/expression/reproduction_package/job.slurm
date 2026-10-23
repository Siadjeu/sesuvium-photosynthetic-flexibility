#!/usr/bin/env bash
#SBATCH --job-name=family_style
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=family_style_%j.out
#SBATCH --error=family_style_%j.err
set -euo pipefail
python3 "${PACKAGE_DIR:?}/run.py"
