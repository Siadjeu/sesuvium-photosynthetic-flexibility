#!/usr/bin/env bash
#SBATCH --job-name=supported_integrate
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=supported_integrate_%j.out
#SBATCH --error=supported_integrate_%j.err
set -euo pipefail
: "${PROJECT_ROOT:?Export PROJECT_ROOT before submission}"
: "${INTEGRATION_SCRIPT_DIR:?Export INTEGRATION_SCRIPT_DIR before submission}"
SOURCE_RELEASE="${SOURCE_RELEASE:-$PROJECT_ROOT/08_focused_families/downstream_supported/full_sources/release_20559577_20261018T184854Z}"
OUT="$PROJECT_ROOT/08_focused_families/downstream_supported/integration/release_${SLURM_JOB_ID:-local}_$(date -u +%Y%m%dT%H%M%SZ)"
python3 "$INTEGRATION_SCRIPT_DIR/integrate_supported.py" --source-release "$SOURCE_RELEASE" --output "$OUT"
cp "$INTEGRATION_SCRIPT_DIR/integrate_supported.py" "$INTEGRATION_SCRIPT_DIR/03_integrate_supported.slurm" "$OUT/"
cp "$INTEGRATION_SCRIPT_DIR/METHODS_RESULTS.txt" "$OUT/"
(cd "$OUT"; sha256sum integrate_supported.py 03_integrate_supported.slurm METHODS_RESULTS.txt >> SHA256SUMS.txt; sha256sum -c SHA256SUMS.txt)
tar -czf "${OUT}_bundle.tar.gz" -C "$(dirname "$OUT")" "$(basename "$OUT")"
printf 'Output: %s\nBundle: %s_bundle.tar.gz\n' "$OUT" "$OUT"
