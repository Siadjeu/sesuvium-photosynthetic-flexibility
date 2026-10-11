#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_9c
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.9c/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.9c/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

cd "$PROJECT_ROOT"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate focused_family_env

for command_name in mafft iqtree python3 sha256sum; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    }
done

RELEASE="08_focused_families/validation/step25kz16.9c/release"
WORKDIR="$RELEASE/work"

mkdir -p "$WORKDIR"

INPUT_OUR="08_focused_families/validation/step25kz16.7b/release/work/corrected_all_species_PEPC_primary_with_refs.faa"

INPUT_WANG="08_focused_families/validation/step25kz12/references/pepc_36_genes.fasta"

INPUT_MAIZE="08_focused_families/validation/step25kz12/work/reference_inputs/Zm00001d046170_PEPC1.faa"

for input_file in "$INPUT_OUR" "$INPUT_WANG" "$INPUT_MAIZE"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: missing or empty input: $input_file" >&2
        exit 1
    fi
done

export RELEASE WORKDIR INPUT_OUR INPUT_WANG INPUT_MAIZE

python3 <<'PY'
import csv
import hashlib
import os
import re
from pathlib import Path

release = Path(os.environ["RELEASE"])
workdir = Path(os.environ["WORKDIR"])

input_our = Path(os.environ["INPUT_OUR"])
input_wang = Path(os.environ["INPUT_WANG"])
input_maize = Path(os.environ["INPUT_MAIZE"])

combined = workdir / "all_PEPC_sequences_deduplicated.faa"
manifest_path = release / "00_sequence_deduplication_manifest.tsv"


def read_fasta(path):
    records = []
    header = None
    sequence = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))

                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(
                    re.sub(r"\s+", "", line).upper()
                )

    if header is not None:
        records.append((header, "".join(sequence)))

    return records


def clean_id(header):
    identifier = header.split()[0]
    identifier = re.sub(
        r"[^A-Za-z0-9_.|:-]+",
        "_",
        identifier,
    )
    return identifier


# Priority is deliberate:
# 1. maize anchor, so it is never lost;
# 2. our corrected project sequences;
# 3. Wang S7 sequences.
source_inputs = [
    ("MAIZE_COORDINATE_REFERENCE", "ZM", input_maize),
    ("OUR_CORRECTED_ALL_SPECIES", "OUR", input_our),
    ("WANG_2023_SUPPLEMENTARY_TABLE_S7", "WANG_S7", input_wang),
]

retained_records = []
manifest_rows = []
sequence_owner = {}
used_ids = set()

for source, prefix, path in source_inputs:
    for index, (header, sequence) in enumerate(
        read_fasta(path),
        start=1,
    ):
        if not sequence:
            continue

        sequence_hash = hashlib.sha1(
            sequence.encode("utf-8")
        ).hexdigest()

        original_id = clean_id(header)
        candidate_id = f"{prefix}__{original_id}"

        if candidate_id in used_ids:
            candidate_id = f"{candidate_id}__{index}"

        duplicate_of = sequence_owner.get(sequence_hash)

        if duplicate_of is None:
            tree_id = candidate_id
            used_ids.add(tree_id)
            sequence_owner[sequence_hash] = tree_id
            retained = "YES"
            retained_records.append((tree_id, sequence))
        else:
            tree_id = candidate_id
            retained = "NO"

        manifest_rows.append({
            "tree_id": tree_id,
            "source": source,
            "original_header": header,
            "original_id": original_id,
            "length_aa": len(sequence),
            "sequence_sha1": sequence_hash,
            "retained_in_tree": retained,
            "duplicate_of": duplicate_of or "",
        })

with combined.open("w", encoding="utf-8") as handle:
    for tree_id, sequence in retained_records:
        handle.write(f">{tree_id}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

with manifest_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    fields = [
        "tree_id",
        "source",
        "original_header",
        "original_id",
        "length_aa",
        "sequence_sha1",
        "retained_in_tree",
        "duplicate_of",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(manifest_rows)

source_totals = {}
source_retained = {}
source_duplicates = {}

for row in manifest_rows:
    source = row["source"]

    source_totals[source] = source_totals.get(source, 0) + 1

    if row["retained_in_tree"] == "YES":
        source_retained[source] = source_retained.get(source, 0) + 1
    else:
        source_duplicates[source] = source_duplicates.get(source, 0) + 1

qc_path = release / "FORMAL_QC.tsv"

with qc_path.open("w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\tstep25kz16.9c_all_pepc_deduplicated_tree\n"
    )
    handle.write("status\tPASS_PENDING_TREE\n")

    handle.write(
        f"all_input_records\t"
        f"{source_totals.get('OUR_CORRECTED_ALL_SPECIES', 0)}\n"
    )

    handle.write(
        f"wang_s7_input_records\t"
        f"{source_totals.get('WANG_2023_SUPPLEMENTARY_TABLE_S7', 0)}\n"
    )

    handle.write(
        f"maize_input_records\t"
        f"{source_totals.get('MAIZE_COORDINATE_REFERENCE', 0)}\n"
    )

    handle.write(
        f"deduplicated_tree_records\t{len(retained_records)}\n"
    )

    handle.write(
        f"duplicate_records_removed\t"
        f"{len(manifest_rows) - len(retained_records)}\n"
    )

    handle.write(
        f"wang_s7_records_retained\t"
        f"{source_retained.get('WANG_2023_SUPPLEMENTARY_TABLE_S7', 0)}\n"
    )

    handle.write(
        f"wang_s7_duplicates_removed\t"
        f"{source_duplicates.get('WANG_2023_SUPPLEMENTARY_TABLE_S7', 0)}\n"
    )

    handle.write(
        "deduplication_rule\t"
        "identical ungapped amino-acid sequences collapsed globally\n"
    )

    handle.write(
        "retention_priority\t"
        "maize anchor, then our corrected sequences, then Wang S7\n"
    )

    handle.write(
        "wang_inclusion_rule\t"
        "all records from pepc_36_genes.fasta retained unless identical\n"
    )

print(f"Input OUR records: {len(read_fasta(input_our))}")
print(f"Input Wang S7 records: {len(read_fasta(input_wang))}")
print(f"Input maize records: {len(read_fasta(input_maize))}")
print(f"Final deduplicated records: {len(retained_records)}")
print(
    "Duplicates removed: "
    f"{len(manifest_rows) - len(retained_records)}"
)
print(f"Created: {combined}")
print(f"Created: {manifest_path}")
print(f"Created: {qc_path}")
PY

echo "===== Running MAFFT ====="

mafft \
    --auto \
    --thread "$SLURM_CPUS_PER_TASK" \
    "$WORKDIR/all_PEPC_sequences_deduplicated.faa" \
    > "$RELEASE/01_all_PEPC_deduplicated_alignment.faa"

if [[ ! -s "$RELEASE/01_all_PEPC_deduplicated_alignment.faa" ]]; then
    echo "ERROR: alignment was not created" >&2
    exit 1
fi

echo "===== Running IQ-TREE ====="

iqtree \
    -s "$RELEASE/01_all_PEPC_deduplicated_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "$RELEASE/PEPC_all_deduplicated_WangS7_maize" \
    -redo

TREE="$RELEASE/PEPC_all_deduplicated_WangS7_maize.treefile"
ALIGNMENT="$RELEASE/01_all_PEPC_deduplicated_alignment.faa"

if [[ ! -s "$TREE" ]]; then
    echo "ERROR: tree was not created" >&2
    exit 1
fi

export TREE ALIGNMENT RELEASE

python3 <<'PY'
import csv
import os
from pathlib import Path

from Bio import Phylo

tree_path = Path(os.environ["TREE"])
alignment_path = Path(os.environ["ALIGNMENT"])
release = Path(os.environ["RELEASE"])


def read_fasta(path):
    records = {}
    header = None
    sequence = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence)

                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence)

    return records


alignment = read_fasta(alignment_path)

maize_ids = [
    sequence_id
    for sequence_id in alignment
    if sequence_id.startswith("ZM__")
]

if len(maize_ids) != 1:
    raise SystemExit(
        f"ERROR: expected exactly one retained maize anchor, found "
        f"{len(maize_ids)}"
    )

maize_id = maize_ids[0]
maize_alignment = alignment[maize_id]

coordinates = [519, 780, 890]
coordinate_columns = {}

for coordinate in coordinates:
    residue_count = 0
    found = False

    for column_index, residue in enumerate(maize_alignment):
        if residue != "-":
            residue_count += 1

        if residue_count == coordinate:
            coordinate_columns[coordinate] = column_index
            found = True
            break

    if not found:
        raise SystemExit(
            f"ERROR: maize coordinate {coordinate} is unavailable"
        )

tree = Phylo.read(str(tree_path), "newick")

tree_order = [
    terminal.name
    for terminal in tree.get_terminals()
    if terminal.name in alignment
]

metadata_path = release / "02_deduplicated_tree_tip_metadata.tsv"

with metadata_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

    writer.writerow([
        "tree_order",
        "tree_id",
        "aligned_length",
        "ungapped_length",
    ])

    for order, sequence_id in enumerate(tree_order, start=1):
        sequence = alignment[sequence_id]

        writer.writerow([
            order,
            sequence_id,
            len(sequence),
            len(sequence.replace("-", "")),
        ])

windows_path = release / "03_deduplicated_coordinate_windows_tree_order.tsv"

with windows_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

    writer.writerow([
        "tree_order",
        "tree_id",
        "source_coordinate",
        "maize_alignment_column_1based",
        "window_alignment_columns_1based",
        "maize_window",
        "sequence_window",
        "coverage_status",
    ])

    for order, sequence_id in enumerate(tree_order, start=1):
        sequence = alignment[sequence_id]

        for coordinate in coordinates:
            center = coordinate_columns[coordinate]
            start = center - 2
            end = center + 3

            window = sequence[start:end]
            maize_window = maize_alignment[start:end]

            if len(window) < 5 or set(window) == {"-"}:
                status = "NO_COVERAGE"
            elif "-" in window:
                status = "PARTIAL_COVERAGE"
            else:
                status = "FULL_COVERAGE"

            writer.writerow([
                order,
                sequence_id,
                coordinate,
                center + 1,
                f"{start + 1}-{end}",
                maize_window,
                window,
                status,
            ])

qc_path = release / "FORMAL_QC.tsv"

with qc_path.open("a", encoding="utf-8") as handle:
    handle.write(f"alignment_tip_count\t{len(alignment)}\n")
    handle.write(f"tree_terminal_count\t{len(tree_order)}\n")
    handle.write(f"maize_anchor\t{maize_id}\n")
    handle.write(
        f"maize_coordinate_519_alignment_column\t"
        f"{coordinate_columns[519] + 1}\n"
    )
    handle.write(
        f"maize_coordinate_780_alignment_column\t"
        f"{coordinate_columns[780] + 1}\n"
    )
    handle.write(
        f"maize_coordinate_890_alignment_column\t"
        f"{coordinate_columns[890] + 1}\n"
    )
    handle.write(
        "coordinate_window_definition\t"
        "two residues before and after each maize coordinate\n"
    )
    handle.write(
        "all_gap_windows\t"
        "reported as NO_COVERAGE and not biological substitutions\n"
    )

print(f"Tree terminals: {len(tree_order)}")
print(f"Maize anchor: {maize_id}")
print(f"Created: {metadata_path}")
print(f"Created: {windows_path}")
PY

cd "$RELEASE"

sha256sum \
    00_sequence_deduplication_manifest.tsv \
    01_all_PEPC_deduplicated_alignment.faa \
    02_deduplicated_tree_tip_metadata.tsv \
    03_deduplicated_coordinate_windows_tree_order.tsv \
    FORMAL_QC.tsv \
    PEPC_all_deduplicated_WangS7_maize.treefile \
    > SHA256SUMS.txt

sha256sum -c SHA256SUMS.txt

sed -i 's/status\tPASS_PENDING_TREE/status\tCOMPLETE/' FORMAL_QC.tsv

{
    echo "checkpoint=step25kz16.9c_all_pepc_deduplicated_tree"
    echo "status=COMPLETE"
    echo "date=$(date --iso-8601=seconds)"
} > STEP25KZ16.9c_COMPLETE.txt

column -t -s $'\t' FORMAL_QC.tsv
