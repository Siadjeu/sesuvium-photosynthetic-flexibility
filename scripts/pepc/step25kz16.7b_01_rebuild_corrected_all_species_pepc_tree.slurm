#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_7b
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.7b/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.7b/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate focused_family_env

cd "$PROJECT_ROOT"

RELEASE="08_focused_families/validation/step25kz16.7b/release"
WORKDIR="$RELEASE/work"

mkdir -p "$RELEASE" "$WORKDIR"

for command_name in python3 mafft iqtree; do
    command -v "$command_name" >/dev/null || {
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    }
done

CANDIDATES="08_focused_families/validation/step25kz16.7/release/work/all_species_PEPC_primary_candidates.faa"

REFERENCE="08_focused_families/validation/step25kz16.7/release/reference_curation/05_corrected_PEPC_references_only.faa"

if [[ ! -s "$CANDIDATES" ]]; then
    echo "ERROR: primary candidate FASTA not found: $CANDIDATES" >&2
    exit 1
fi

if [[ ! -s "$REFERENCE" ]]; then
    echo "ERROR: corrected reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

export RELEASE WORKDIR CANDIDATES REFERENCE

python3 <<'PY'
from pathlib import Path
import os
import re

release = Path(os.environ["RELEASE"])
workdir = Path(os.environ["WORKDIR"])
candidates = Path(os.environ["CANDIDATES"])
reference = Path(os.environ["REFERENCE"])

def read_fasta(path):
    header = None
    sequence = []

    with path.open(
        encoding="utf-8",
        errors="ignore",
    ) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)

                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        yield header, "".join(sequence)

output = workdir / "corrected_all_species_PEPC_primary_with_refs.faa"

used = set()
candidate_count = 0
reference_count = 0

with output.open("w", encoding="utf-8") as handle:
    for header, sequence in read_fasta(reference):
        clean_id = "REF__" + re.sub(
            r"[^A-Za-z0-9_.-]+",
            "_",
            header.split()[0],
        )

        if clean_id in used:
            continue

        used.add(clean_id)
        reference_count += 1
        handle.write(f">{clean_id}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

    for header, sequence in read_fasta(candidates):
        clean_id = header.split()[0]

        if clean_id in used:
            continue

        used.add(clean_id)
        candidate_count += 1
        handle.write(f">{clean_id}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

if "SS__g24599" not in used:
    raise SystemExit(
        "ERROR: SS__g24599 is absent from corrected tree input."
    )

with (release / "00_corrected_tree_input_qc.tsv").open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"primary_candidate_count\t{candidate_count}\n")
    handle.write(f"corrected_reference_count\t{reference_count}\n")
    handle.write(f"total_tree_tips\t{len(used)}\n")
    handle.write("SS_g24599_present\tYES\n")
    handle.write(
        "reference_set\tcanonical_CAPP1_to_CAPP4_PPC1_to_PPC4\n"
    )

print(f"Primary candidates: {candidate_count}")
print(f"Corrected references: {reference_count}")
print(f"Total tree tips: {len(used)}")
print("Confirmed: SS__g24599 is present.")
PY

mafft \
    --auto \
    "$WORKDIR/corrected_all_species_PEPC_primary_with_refs.faa" \
    > "$RELEASE/01_corrected_all_species_PEPC_alignment.faa"

iqtree \
    -s "$RELEASE/01_corrected_all_species_PEPC_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    -pre "$RELEASE/PEPC_all_species_corrected_primary" \
    -quiet

python3 <<'PY'
from pathlib import Path
import os

release = Path(os.environ["RELEASE"])

tree = release / "PEPC_all_species_corrected_primary.treefile"
alignment = release / "01_corrected_all_species_PEPC_alignment.faa"

if not tree.exists() or tree.stat().st_size == 0:
    raise SystemExit("ERROR: corrected PEPC tree was not created.")

if "SS__g24599" not in tree.read_text(encoding="utf-8"):
    raise SystemExit(
        "ERROR: SS__g24599 is absent from corrected PEPC tree."
    )

with (release / "FORMAL_QC.tsv").open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\tstep25kz16.7b_corrected_all_species_pepc_tree\n"
    )
    handle.write("status\tCOMPLETE\n")
    handle.write("reference_class\tcanonical_PPC_only\n")
    handle.write("noncanonical_references_excluded\tYES\n")
    handle.write("SS_g24599_in_tree\tYES\n")
    handle.write("bootstrap_replicates\t1000\n")
    handle.write("SH_aLRT_replicates\t1000\n")

print((release / "FORMAL_QC.tsv").read_text(), end="")
PY

(
    cd "$RELEASE"

    sha256sum \
        00_corrected_tree_input_qc.tsv \
        01_corrected_all_species_PEPC_alignment.faa \
        PEPC_all_species_corrected_primary.treefile \
        FORMAL_QC.tsv \
        > SHA256SUMS.txt

    sha256sum -c SHA256SUMS.txt
)

{
    echo "checkpoint=step25kz16.7b_corrected_all_species_pepc_tree"
    echo "status=COMPLETE"
    echo "date=$(date --iso-8601=seconds)"
    echo "tree=PEPC_all_species_corrected_primary.treefile"
    echo "reference_set=canonical_PPC_only"
} > "$RELEASE/STEP25KZ16.7B_COMPLETE.txt"

echo "Step25kz16.7b completed successfully."
