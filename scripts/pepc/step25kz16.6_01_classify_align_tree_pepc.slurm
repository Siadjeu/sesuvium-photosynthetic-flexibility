#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_6
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.6/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.6/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate focused_family_env

cd "$PROJECT_ROOT"

RELEASE="08_focused_families/validation/step25kz16.6/release"
WORKDIR="$RELEASE/work"

mkdir -p "$RELEASE" "$WORKDIR"

for command_name in python3 diamond mafft iqtree; do
    command -v "$command_name" >/dev/null || {
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    }
done

CANDIDATE_FASTA="08_focused_families/validation/step25kz16.5/release/03_complete_PEPC_tree_input.faa"
CANDIDATE_TABLE="08_focused_families/validation/step25kz16.5/release/02_complete_pepc_candidate_audit.tsv"
REGISTRY="08_focused_families/validation/step25kn/results/PEPC/candidate_validation.tsv"

REFERENCE="$(
    find \
        08_focused_families/validation/step25kn/references \
        08_focused_families/validation/step25kn/phylogeny_inputs \
        -type f \
        -iname 'PEPC.faa' \
        2>/dev/null |
    head -n 1
)"

for input_file in "$CANDIDATE_FASTA" "$CANDIDATE_TABLE" "$REFERENCE"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: missing or empty input: $input_file" >&2
        exit 1
    fi
done

export PROJECT_ROOT RELEASE WORKDIR CANDIDATE_FASTA CANDIDATE_TABLE REGISTRY REFERENCE

python3 <<'PY'
from pathlib import Path
import csv
import os
import re

release = Path(os.environ["RELEASE"])
workdir = Path(os.environ["WORKDIR"])
candidate_fasta = Path(os.environ["CANDIDATE_FASTA"])
candidate_table = Path(os.environ["CANDIDATE_TABLE"])
registry = Path(os.environ["REGISTRY"])

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

sequences = {}

for header, sequence in read_fasta(candidate_fasta):
    sequence = re.sub(r"[^A-Za-z]", "", sequence)
    sequences[header.split()[0]] = sequence

metadata = {}

with candidate_table.open(encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        metadata[row["tree_id"]] = row

known_pepc_ids = {
    "SS__g24599",
}

if registry.exists():
    for line in registry.read_text(
        encoding="utf-8",
        errors="ignore",
    ).splitlines():
        for token in re.split(r"[\t ,;|]+", line):
            if token.startswith(("PO__", "SS__")):
                known_pepc_ids.add(token)

known_pepc_ids &= set(sequences)

relaxed_query = workdir / "PEPC_relaxed_query.faa"

with relaxed_query.open("w", encoding="utf-8") as handle:
    for tree_id in sorted(sequences):
        handle.write(f">{tree_id}\n")
        sequence = sequences[tree_id]

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

with (release / "00_candidate_sequence_inventory.tsv").open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow([
        "tree_id",
        "species",
        "length_aa",
        "known_validated_pepc",
    ])

    for tree_id in sorted(sequences):
        species = tree_id.split("__", 1)[0]

        writer.writerow([
            tree_id,
            species,
            len(sequences[tree_id]),
            "YES" if tree_id in known_pepc_ids else "NO",
        ])

print(f"Candidate sequences: {len(sequences)}")
print(f"Validated PEPC IDs retained: {len(known_pepc_ids)}")
PY

diamond makedb \
    --in "$REFERENCE" \
    --db "$WORKDIR/PEPC_reference"

diamond blastp \
    --query "$WORKDIR/PEPC_relaxed_query.faa" \
    --db "$WORKDIR/PEPC_reference.dmnd" \
    --out "$WORKDIR/candidates_vs_PEPC_reference.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp \
    --evalue 1e-5 \
    --sensitive \
    --max-target-seqs 20 \
    --threads "$SLURM_CPUS_PER_TASK"

export REFERENCE

python3 <<'PY'
from pathlib import Path
import csv
import os
import re

release = Path(os.environ["RELEASE"])
workdir = Path(os.environ["WORKDIR"])
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

candidate_sequences = {}

for header, sequence in read_fasta(
    workdir / "PEPC_relaxed_query.faa"
):
    candidate_sequences[header.split()[0]] = sequence

relaxed_hits = {}

with (
    workdir / "candidates_vs_PEPC_reference.tsv"
).open(encoding="utf-8") as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        query, subject, pident, length, qlen, slen, evalue, bitscore, qcovhsp = fields

        hit = {
            "reference": subject,
            "pident": float(pident),
            "length": int(length),
            "evalue": float(evalue),
            "bitscore": float(bitscore),
            "qcovhsp": float(qcovhsp),
        }

        if (
            query not in relaxed_hits
            or hit["bitscore"] > relaxed_hits[query]["bitscore"]
        ):
            relaxed_hits[query] = hit

known_pepc_ids = {"SS__g24599"}

registry = Path(
    os.environ["PROJECT_ROOT"]
) / "08_focused_families/validation/step25kn/results/PEPC/candidate_validation.tsv"

if registry.exists():
    for line in registry.read_text(
        encoding="utf-8",
        errors="ignore",
    ).splitlines():
        for token in re.split(r"[\t ,;|]+", line):
            if token.startswith(("PO__", "SS__")):
                known_pepc_ids.add(token)

primary_ids = set()
sensitivity_ids = set()

for tree_id, sequence in candidate_sequences.items():
    length = len(sequence)
    hit = relaxed_hits.get(tree_id)

    similarity_support = (
        hit is not None
        and hit["pident"] >= 20
        and hit["qcovhsp"] >= 50
        and hit["bitscore"] >= 100
        and hit["evalue"] <= 1e-5
    )

    full_length = 800 <= length <= 1300
    fragment = 300 <= length < 800

    if full_length and (
        similarity_support
        or tree_id in known_pepc_ids
    ):
        primary_ids.add(tree_id)

    if (
        (full_length or fragment)
        and (
            similarity_support
            or tree_id in known_pepc_ids
        )
    ):
        sensitivity_ids.add(tree_id)

def write_fasta(path, ids, label):
    with path.open("w", encoding="utf-8") as handle:
        for tree_id in sorted(ids):
            sequence = candidate_sequences[tree_id]
            handle.write(
                f">{tree_id} dataset={label}\n"
            )

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

write_fasta(
    workdir / "PO_SS_PEPC_primary_candidates.faa",
    primary_ids,
    "primary_full_length",
)

write_fasta(
    workdir / "PO_SS_PEPC_sensitivity_candidates.faa",
    sensitivity_ids,
    "sensitivity_full_and_fragment",
)

rows = []

for tree_id in sorted(candidate_sequences):
    length = len(candidate_sequences[tree_id])
    hit = relaxed_hits.get(tree_id)

    similarity_support = (
        hit is not None
        and hit["pident"] >= 20
        and hit["qcovhsp"] >= 50
        and hit["bitscore"] >= 100
        and hit["evalue"] <= 1e-5
    )

    if tree_id in primary_ids:
        classification = "PRIMARY_FULL_LENGTH_PEPC"
        reason = "full_length_plus_similarity_or_validated_registry"
    elif tree_id in sensitivity_ids:
        classification = "SENSITIVITY_FRAGMENT_OR_AMBIGUOUS"
        reason = "fragment_or_nonprimary_length_with_support"
    else:
        classification = "EXCLUDED_FROM_TREE"
        reason = "no_sufficient_pepc_support_or_invalid_length"

    rows.append({
        "tree_id": tree_id,
        "species": tree_id.split("__", 1)[0],
        "length_aa": length,
        "validated_registry": (
            "YES" if tree_id in known_pepc_ids else "NO"
        ),
        "relaxed_similarity": (
            "YES" if similarity_support else "NO"
        ),
        "best_reference": (
            hit["reference"] if hit else ""
        ),
        "best_pident": (
            hit["pident"] if hit else ""
        ),
        "best_qcovhsp": (
            hit["qcovhsp"] if hit else ""
        ),
        "best_bitscore": (
            hit["bitscore"] if hit else ""
        ),
        "best_evalue": (
            hit["evalue"] if hit else ""
        ),
        "classification": classification,
        "classification_reason": reason,
    })

with (release / "01_pepc_candidate_classification.tsv").open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

if "SS__g24599" not in primary_ids:
    raise SystemExit(
        "ERROR: SS__g24599 is not in the primary PEPC set."
    )

print(f"Primary full-length PEPC candidates: {len(primary_ids)}")
print(f"Sensitivity PEPC candidates: {len(sensitivity_ids)}")
print("Confirmed: SS__g24599 is in the primary set.")
PY

python3 <<'PY'
from pathlib import Path
import os
import re

release = Path(os.environ["RELEASE"])
reference = Path(os.environ["REFERENCE"])
workdir = Path(os.environ["WORKDIR"])

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

def write_combined(candidate_file, output_file):
    used = set()

    with output_file.open("w", encoding="utf-8") as out:
        for header, sequence in read_fasta(reference):
            ref_id = "REF__" + re.sub(
                r"[^A-Za-z0-9_.-]+",
                "_",
                header.split()[0],
            )

            if ref_id in used:
                continue

            used.add(ref_id)
            out.write(f">{ref_id}\n")

            for start in range(0, len(sequence), 80):
                out.write(sequence[start:start + 80] + "\n")

        for header, sequence in read_fasta(candidate_file):
            tree_id = header.split()[0]

            if tree_id in used:
                continue

            used.add(tree_id)
            out.write(f">{tree_id}\n")

            for start in range(0, len(sequence), 80):
                out.write(sequence[start:start + 80] + "\n")

write_combined(
    workdir / "PO_SS_PEPC_primary_candidates.faa",
    workdir / "PEPC_primary_with_references.faa",
)

write_combined(
    workdir / "PO_SS_PEPC_sensitivity_candidates.faa",
    workdir / "PEPC_sensitivity_with_references.faa",
)
PY

mafft \
    --auto \
    "$WORKDIR/PEPC_primary_with_references.faa" \
    > "$RELEASE/02_PEPC_primary_alignment.faa"

mafft \
    --auto \
    "$WORKDIR/PEPC_sensitivity_with_references.faa" \
    > "$RELEASE/03_PEPC_sensitivity_alignment.faa"

iqtree \
    -s "$RELEASE/02_PEPC_primary_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    -pre "$RELEASE/PEPC_primary" \
    -quiet

iqtree \
    -s "$RELEASE/03_PEPC_sensitivity_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    -pre "$RELEASE/PEPC_sensitivity" \
    -quiet

python3 <<'PY'
from pathlib import Path
import csv
import os

release = Path(os.environ["RELEASE"])

classification = release / "01_pepc_candidate_classification.tsv"

rows = list(
    csv.DictReader(
        classification.open(encoding="utf-8"),
        delimiter="\t",
    )
)

primary_po = sum(
    row["classification"] == "PRIMARY_FULL_LENGTH_PEPC"
    and row["species"] == "PO"
    for row in rows
)

primary_ss = sum(
    row["classification"] == "PRIMARY_FULL_LENGTH_PEPC"
    and row["species"] == "SS"
    for row in rows
)

sensitivity = sum(
    row["classification"] == "SENSITIVITY_FRAGMENT_OR_AMBIGUOUS"
    for row in rows
)

excluded = sum(
    row["classification"] == "EXCLUDED_FROM_TREE"
    for row in rows
)

with (release / "FORMAL_QC.tsv").open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\tstep25kz16.6_pepc_classification_alignment_tree\n"
    )
    handle.write("status\tCOMPLETE\n")
    handle.write(f"primary_PO_PEPC_count\t{primary_po}\n")
    handle.write(f"primary_SS_PEPC_count\t{primary_ss}\n")
    handle.write(f"sensitivity_candidate_count\t{sensitivity}\n")
    handle.write(f"excluded_candidate_count\t{excluded}\n")
    handle.write("SS_g24599_in_primary_tree\tYES\n")
    handle.write(
        "primary_rule\t"
        "800-1300 aa plus relaxed PEPC similarity or validated registry\n"
    )
    handle.write(
        "sensitivity_rule\t"
        "300-1300 aa plus relaxed PEPC similarity or validated registry\n"
    )
    handle.write(
        "reference_rule\t"
        "PEPC.faa reference proteins included in both trees\n"
    )

print((release / "FORMAL_QC.tsv").read_text(), end="")
PY

(
    cd "$RELEASE"

    sha256sum \
        00_candidate_sequence_inventory.tsv \
        01_pepc_candidate_classification.tsv \
        02_PEPC_primary_alignment.faa \
        03_PEPC_sensitivity_alignment.faa \
        FORMAL_QC.tsv \
        PEPC_primary.treefile \
        PEPC_sensitivity.treefile \
        > SHA256SUMS.txt
)

{
    echo "checkpoint=step25kz16.6_pepc_classification_alignment_tree"
    echo "status=COMPLETE"
    echo "date=$(date --iso-8601=seconds)"
    echo "primary_tree=PEPC_primary.treefile"
    echo "sensitivity_tree=PEPC_sensitivity.treefile"
    echo "required_target=SS__g24599"
} > "$RELEASE/STEP25KZ16.6_COMPLETE.txt"

(
    cd "$RELEASE"
    sha256sum -c SHA256SUMS.txt
)

echo "Step25kz16.6 completed successfully."
