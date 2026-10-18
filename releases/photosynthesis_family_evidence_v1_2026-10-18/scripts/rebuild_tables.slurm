#!/usr/bin/env bash
#SBATCH --job-name=family_freeze_qc
#SBATCH --partition=all_cpu.p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=family_freeze_qc_%j.out
#SBATCH --error=family_freeze_qc_%j.err
set -euo pipefail
: "${FREEZE_ROOT:?Export FREEZE_ROOT as the absolute frozen release directory}"
: "${OUTPUT_DIR:?Export OUTPUT_DIR as a new output directory}"
python3 "$FREEZE_ROOT/scripts/build_family_tables.py" --release "$FREEZE_ROOT" --out "$OUTPUT_DIR"
(cd "$OUTPUT_DIR" && sha256sum -c SHA256SUMS.txt)
