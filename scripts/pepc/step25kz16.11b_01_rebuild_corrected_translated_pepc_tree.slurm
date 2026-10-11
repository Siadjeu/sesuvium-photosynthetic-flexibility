#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_11b
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.11b/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.11b/logs/%x_%j.err

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

RELEASE="08_focused_families/validation/step25kz16.11b/release"
WORKDIR="$RELEASE/work"

mkdir -p "$WORKDIR"

BASE_INPUT="08_focused_families/validation/step25kz16.7b/release/work/corrected_all_species_PEPC_primary_with_refs.faa"

WANG_TRANSLATED="08_focused_families/validation/step25kz16.11/release/work/Wang_S7_translated_proteins.faa"

MAIZE_INPUT="08_focused_families/validation/step25kz12/work/reference_inputs/Zm00001d046170_PEPC1.faa"

for input_file in \
    "$BASE_INPUT" \
    "$WANG_TRANSLATED" \
    "$MAIZE_INPUT"; do

    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: missing input: $input_file" >&2
        exit 1
    fi
done

export RELEASE WORKDIR BASE_INPUT WANG_TRANSLATED MAIZE_INPUT

python3 <<'PY'
import csv
import hashlib
import os
import re
from pathlib import Path

release = Path(os.environ["RELEASE"])
workdir = Path(os.environ["WORKDIR"])

base_input = Path(os.environ["BASE_INPUT"])
wang_input = Path(os.environ["WANG_TRANSLATED"])
maize_input = Path(os.environ["MAIZE_INPUT"])

excluded_gene = "SS__g24599"


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


def first_id(header):
    return header.split()[0]


def normalize_project_id(raw_id):
    raw_id = raw_id.replace("OUR__", "")

    if raw_id.startswith("REF__"):
        return "OUR__" + raw_id

    return "OUR__" + raw_id


def normalize_wang_id(raw_id):
    raw_id = raw_id.replace(
        "WANG_S7__WANG_S7__",
        "WANG_S7__",
    )

    if raw_id.startswith("WANG_S7__"):
        return raw_id

    return "WANG_S7__" + raw_id


all_records = []
manifest = []
sequence_owner = {}
used_ids = set()
excluded_records = []

# 1. Complete project primary input.
for header, sequence in read_fasta(base_input):
    raw_id = first_id(header)
    normalized_id = normalize_project_id(raw_id)
    gene_id = normalized_id.replace("OUR__", "", 1)

    if gene_id == excluded_gene:
        excluded_records.append({
            "original_id": raw_id,
            "source": "OUR_CORRECTED_ALL_SPECIES",
            "length_aa": len(sequence),
            "reason": "NO_INDEPENDENT_PEPC_HOMOLOGY",
        })
        continue

    sequence = sequence.replace("-", "")

    if not sequence:
        continue

    all_records.append({
        "tree_id": normalized_id,
        "source": "OUR_CORRECTED_ALL_SPECIES",
        "original_header": header,
        "original_id": raw_id,
        "sequence": sequence,
    })

# 2. Explicit maize anchor.
for header, sequence in read_fasta(maize_input):
    sequence = sequence.replace("-", "")

    if not sequence:
        continue

    all_records.append({
        "tree_id": "ZM__" + first_id(header),
        "source": "MAIZE_COORDINATE_REFERENCE",
        "original_header": header,
        "original_id": first_id(header),
        "sequence": sequence,
    })

# 3. All translated Wang S7 records.
for header, sequence in read_fasta(wang_input):
    sequence = sequence.replace("-", "")

    if not sequence:
        continue

    all_records.append({
        "tree_id": normalize_wang_id(first_id(header)),
        "source": "WANG_2023_SUPPLEMENTARY_TABLE_S7_TRANSLATED",
        "original_header": header,
        "original_id": first_id(header),
        "sequence": sequence,
    })

retained = []

for record in all_records:
    sequence = record["sequence"]

    sequence_hash = hashlib.sha1(
        sequence.encode("utf-8")
    ).hexdigest()

    duplicate_of = sequence_owner.get(sequence_hash)

    candidate_id = record["tree_id"]

    if duplicate_of is None:
        tree_id = candidate_id
        counter = 2

        while tree_id in used_ids:
            tree_id = f"{candidate_id}_{counter}"
            counter += 1

        used_ids.add(tree_id)
        sequence_owner[sequence_hash] = tree_id
        retained_status = "YES"
        retained.append((tree_id, sequence, record))
    else:
        tree_id = candidate_id
        retained_status = "NO"

    manifest.append({
        "tree_id": tree_id,
        "source": record["source"],
        "original_header": record["original_header"],
        "original_id": record["original_id"],
        "length_aa": len(sequence),
        "sequence_sha1": sequence_hash,
        "retained_in_tree": retained_status,
        "duplicate_of": duplicate_of or "",
    })

combined = workdir / "all_corrected_translated_PEPC_deduplicated.faa"

with combined.open("w", encoding="utf-8") as handle:
    for tree_id, sequence, record in retained:
        handle.write(f">{tree_id}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

manifest_path = release / "00_sequence_manifest.tsv"

with manifest_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "tree_id",
            "source",
            "original_header",
            "original_id",
            "length_aa",
            "sequence_sha1",
            "retained_in_tree",
            "duplicate_of",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(manifest)

excluded_path = release / "01_excluded_candidates.tsv"

with excluded_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "original_id",
            "source",
            "length_aa",
            "reason",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(excluded_records)

qc_path = release / "FORMAL_QC.tsv"

base_count = sum(
    record["source"] == "OUR_CORRECTED_ALL_SPECIES"
    for record in all_records
)

maize_count = sum(
    record["source"] == "MAIZE_COORDINATE_REFERENCE"
    for record in all_records
)

wang_count = sum(
    record["source"].startswith("WANG_")
    for record in all_records
)

with qc_path.open("w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\tstep25kz16.11b_corrected_translated_pepc_tree\n"
    )
    handle.write("status\tPASS_PENDING_TREE\n")
    handle.write(f"project_records_retained_before_dedup\t{base_count}\n")
    handle.write(f"maize_records_added\t{maize_count}\n")
    handle.write(f"wang_s7_records_added\t{wang_count}\n")
    handle.write(
        f"excluded_ss_g24599\t"
        f"{len(excluded_records)}\n"
    )
    handle.write(
        f"final_deduplicated_records\t{len(retained)}\n"
    )
    handle.write(
        f"duplicate_records_removed\t"
        f"{len(all_records) - len(retained)}\n"
    )
    handle.write(
        "id_normalization\t"
        "OUR, WANG_S7, and ZM prefixes standardized\n"
    )
    handle.write(
        "translation_status\t"
        "Wang S7 CDS records were translated before this step\n"
    )

print(f"Project records: {base_count}")
print(f"Maize records added: {maize_count}")
print(f"Wang S7 records added: {wang_count}")
print(f"SS__g24599 excluded: {len(excluded_records)}")
print(f"Final deduplicated records: {len(retained)}")
print(f"Duplicates removed: {len(all_records) - len(retained)}")
print(f"Created: {combined}")
print(f"Created: {manifest_path}")
print(f"Created: {excluded_path}")
PY

grep -q '^>ZM__' \
    "$WORKDIR/all_corrected_translated_PEPC_deduplicated.faa" || {
    echo "ERROR: maize anchor was not added" >&2
    exit 1
}

if grep -q '^>OUR__SS__g24599' \
    "$WORKDIR/all_corrected_translated_PEPC_deduplicated.faa"; then
    echo "ERROR: SS__g24599 was not excluded" >&2
    exit 1
fi

mafft \
    --auto \
    --thread "$SLURM_CPUS_PER_TASK" \
    "$WORKDIR/all_corrected_translated_PEPC_deduplicated.faa" \
    > "$RELEASE/02_corrected_translated_PEPC_alignment.faa"

iqtree \
    -s "$RELEASE/02_corrected_translated_PEPC_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "$RELEASE/PEPC_corrected_translated_S7" \
    -redo

TREE="$RELEASE/PEPC_corrected_translated_S7.treefile"
ALIGNMENT="$RELEASE/02_corrected_translated_PEPC_alignment.faa"

export TREE ALIGNMENT RELEASE

python3 <<'PY'
import csv
import os
from pathlib import Path

from Bio import Phylo

release = Path(os.environ["RELEASE"])
tree_path = Path(os.environ["TREE"])
alignment_path = Path(os.environ["ALIGNMENT"])


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

if "OUR__SS__g24599" in alignment:
    raise SystemExit(
        "ERROR: SS__g24599 is present in the alignment"
    )

maize_ids = [
    sequence_id
    for sequence_id in alignment
    if sequence_id.startswith("ZM__")
]

if len(maize_ids) != 1:
    raise SystemExit(
        f"ERROR: expected one maize anchor, found {len(maize_ids)}"
    )

maize_id = maize_ids[0]
maize = alignment[maize_id]

coordinates = [519, 780, 890]
coordinate_columns = {}

for coordinate in coordinates:
    target_positions = set(
        range(coordinate - 2, coordinate + 3)
    )

    residue_count = 0
    columns = []

    for column_index, residue in enumerate(maize):
        if residue != "-":
            residue_count += 1

            if residue_count in target_positions:
                columns.append(column_index)

    if len(columns) != 5:
        raise SystemExit(
            f"ERROR: unable to map maize residues around {coordinate}"
        )

    coordinate_columns[coordinate] = columns

tree = Phylo.read(str(tree_path), "newick")

tree_ids = [
    terminal.name
    for terminal in tree.get_terminals()
    if terminal.name in alignment
]

metadata_path = release / "03_tree_tip_metadata.tsv"

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

    for order, sequence_id in enumerate(tree_ids, start=1):
        sequence = alignment[sequence_id]

        writer.writerow([
            order,
            sequence_id,
            len(sequence),
            len(sequence.replace("-", "")),
        ])

coordinate_path = release / "04_maize_coordinate_windows_tree_order.tsv"

with coordinate_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

    writer.writerow([
        "tree_order",
        "tree_id",
        "maize_coordinate",
        "maize_residue_range",
        "maize_window",
        "sequence_window",
        "coverage_status",
    ])

    for order, sequence_id in enumerate(tree_ids, start=1):
        sequence = alignment[sequence_id]

        for coordinate in coordinates:
            columns = coordinate_columns[coordinate]

            maize_window = "".join(
                maize[column]
                for column in columns
            )

            sequence_window = "".join(
                sequence[column]
                if column < len(sequence)
                else "-"
                for column in columns
            )

            if set(sequence_window) == {"-"}:
                status = "NO_COVERAGE"
            elif "-" in sequence_window:
                status = "PARTIAL_COVERAGE"
            elif "X" in sequence_window:
                status = "AMBIGUOUS_COVERAGE"
            else:
                status = "FULL_COVERAGE"

            writer.writerow([
                order,
                sequence_id,
                coordinate,
                f"{coordinate - 2}-{coordinate + 2}",
                maize_window,
                sequence_window,
                status,
            ])

qc_path = release / "FORMAL_QC.tsv"

with qc_path.open("a", encoding="utf-8") as handle:
    handle.write(f"alignment_tip_count\t{len(alignment)}\n")
    handle.write(f"tree_terminal_count\t{len(tree_ids)}\n")
    handle.write(f"maize_anchor\t{maize_id}\n")
    handle.write("ss_g24599_in_tree\tNO\n")
    handle.write("coordinate_positions\t519,780,890\n")
    handle.write(
        "coordinate_rule\t"
        "five maize residues from coordinate minus two through plus two\n"
    )
    handle.write("tree_status\tCOMPLETE\n")

print(f"Tree terminals: {len(tree_ids)}")
print(f"Maize anchor: {maize_id}")
print(f"Created: {metadata_path}")
print(f"Created: {coordinate_path}")
PY

cd "$RELEASE"

sed -i 's/status\tPASS_PENDING_TREE/status\tCOMPLETE/' FORMAL_QC.tsv

sha256sum \
    00_sequence_manifest.tsv \
    01_excluded_candidates.tsv \
    02_corrected_translated_PEPC_alignment.faa \
    03_tree_tip_metadata.tsv \
    04_maize_coordinate_windows_tree_order.tsv \
    FORMAL_QC.tsv \
    PEPC_corrected_translated_S7.treefile \
    > SHA256SUMS.txt

sha256sum -c SHA256SUMS.txt

{
    echo "checkpoint=step25kz16.11b_corrected_translated_pepc_tree"
    echo "status=COMPLETE"
    echo "date=$(date --iso-8601=seconds)"
} > STEP25KZ16.11b_COMPLETE.txt

column -t -s $'\t' FORMAL_QC.tsv
