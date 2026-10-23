#!/usr/bin/env bash
#SBATCH --job-name=ss_evidence_freeze
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=ss_evidence_freeze_%j.out
#SBATCH --error=ss_evidence_freeze_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
python3 "$PROJECT_ROOT/ss_transcript_metabolite_freeze_v1/freeze.py" --project-root "$PROJECT_ROOT"
