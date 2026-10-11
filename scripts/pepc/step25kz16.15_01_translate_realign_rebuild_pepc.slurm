#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_15
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.15/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.15/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

cd "$PROJECT_ROOT"

RELEASE="08_focused_families/validation/step25kz16.15/release"
INPUT="08_focused_families/validation/step25kz16.11b/release/work/all_corrected_translated_PEPC_deduplicated.faa"

mkdir -p "$RELEASE/work"

export RELEASE INPUT

python3 <<'PY'
import csv
import os
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


release = Path(os.environ["RELEASE"])
input_path = Path(os.environ["INPUT"])

records = list(SeqIO.parse(input_path, "fasta"))

if not records:
    raise SystemExit("No sequences found.")

translated_records = []
qc_rows = []

dna_pattern = re.compile(r"^[ACGTN]+$", re.IGNORECASE)

for record in records:
    sequence = str(record.seq).upper().replace("-", "")
    original_length = len(sequence)

    if (
        original_length >= 90
        and dna_pattern.fullmatch(sequence)
        and original_length % 3 == 0
    ):
        translated = str(
            Seq(sequence).translate(
                table=1,
                to_stop=False,
            )
        )

        terminal_stop = translated.endswith("*")

        if terminal_stop:
            translated = translated[:-1]

        internal_stop_count = translated.count("*")

        translated = translated.replace("*", "X")

        status = "TRANSLATED_FRAME_1"

        qc_rows.append(
            {
                "id": record.id,
                "original_length": original_length,
                "final_length": len(translated),
                "input_type": "DNA",
                "translation_status": status,
                "terminal_stop_removed": (
                    "YES" if terminal_stop else "NO"
                ),
                "internal_stops_replaced_X": str(
                    internal_stop_count
                ),
            }
        )
    else:
        translated = sequence

        qc_rows.append(
            {
                "id": record.id,
                "original_length": original_length,
                "final_length": len(translated),
                "input_type": "PROTEIN",
                "translation_status": "ALREADY_PROTEIN",
                "terminal_stop_removed": "NO",
                "internal_stops_replaced_X": "0",
            }
        )

    translated_records.append(
        (record.id, translated)
    )


# Remove exact duplicates after translation.
kept = []
seen = {}
manifest_rows = []

for record_id, sequence in translated_records:
    if sequence in seen:
        manifest_rows.append(
            {
                "tree_id": record_id,
                "length_aa": len(sequence),
                "retained": "NO",
                "duplicate_of": seen[sequence],
            }
        )
        continue

    seen[sequence] = record_id
    kept.append((record_id, sequence))

    manifest_rows.append(
        {
            "tree_id": record_id,
            "length_aa": len(sequence),
            "retained": "YES",
            "duplicate_of": "",
        }
    )


translated_input = (
    release / "work/"
    "all_corrected_translated_PEPC_deduplicated_after_translation.faa"
)

with translated_input.open("w", encoding="utf-8") as handle:
    for record_id, sequence in kept:
        handle.write(f">{record_id}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")


with (
    release / "01_translation_QC.tsv"
).open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "id",
            "original_length",
            "final_length",
            "input_type",
            "translation_status",
            "terminal_stop_removed",
            "internal_stops_replaced_X",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(qc_rows)


with (
    release / "02_translation_deduplication_manifest.tsv"
).open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "tree_id",
            "length_aa",
            "retained",
            "duplicate_of",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)


target = next(
    (
        sequence
        for record_id, sequence in kept
        if record_id == "OUR__SS__g8156"
    ),
    None,
)

if target is None:
    raise SystemExit(
        "ERROR: OUR__SS__g8156 was not found."
    )

if not re.search(r"[^ACGTN]", target):
    raise SystemExit(
        "ERROR: OUR__SS__g8156 still appears DNA-like."
    )

with (
    release / "FORMAL_QC_pre_alignment.tsv"
).open("w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "SS_g8156_original_length\t1125\n"
    )
    handle.write(
        f"SS_g8156_final_length\t{len(target)}\n"
    )
    handle.write(
        "SS_g8156_translation_status\tTRANSLATED\n"
    )
    handle.write(
        f"input_records\t{len(records)}\n"
    )
    handle.write(
        f"retained_records_after_translation_dedup\t"
        f"{len(kept)}\n"
    )

print(f"Input records: {len(records)}")
print(f"Retained records: {len(kept)}")
print(f"SS__g8156 translated length: {len(target)} aa")
PY

mafft \
  --auto \
  --thread "$SLURM_CPUS_PER_TASK" \
  "$RELEASE/work/all_corrected_translated_PEPC_deduplicated_after_translation.faa" \
  > "$RELEASE/03_corrected_translated_PEPC_alignment.faa"

iqtree \
  -s "$RELEASE/03_corrected_translated_PEPC_alignment.faa" \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T "$SLURM_CPUS_PER_TASK" \
  --prefix "$RELEASE/PEPC_corrected_translated"

export ALIGN="$RELEASE/03_corrected_translated_PEPC_alignment.faa"
export TREE="$RELEASE/PEPC_corrected_translated.treefile"

python3 <<'PY'
import csv
import os
from pathlib import Path

from Bio import AlignIO
from Bio import Phylo


release = Path(os.environ["RELEASE"])
alignment_path = Path(os.environ["ALIGN"])
tree_path = Path(os.environ["TREE"])

alignment = AlignIO.read(
    alignment_path,
    "fasta",
)

records = {
    record.id: str(record.seq)
    for record in alignment
}

maize_id = next(
    (
        record_id
        for record_id in records
        if record_id.startswith("ZM__")
    ),
    None,
)

if maize_id is None:
    raise SystemExit(
        "Maize reference was not found."
    )

maize = records[maize_id]

coordinate_columns = {}

ungapped_position = 0

for column, residue in enumerate(maize):
    if residue != "-":
        ungapped_position += 1

        if ungapped_position in {
            519,
            780,
            890,
        }:
            coordinate_columns[ungapped_position] = column

missing = {
    position
    for position in {519, 780, 890}
    if position not in coordinate_columns
}

if missing:
    raise SystemExit(
        "Missing maize coordinates: "
        + ", ".join(map(str, sorted(missing)))
    )


tree = Phylo.read(tree_path, "newick")

tree_ids = [
    terminal.name
    for terminal in tree.get_terminals()
    if terminal.name in records
]

rows = []

for order, record_id in enumerate(tree_ids, start=1):
    sequence = records[record_id]

    for position in [519, 780, 890]:
        center = coordinate_columns[position]

        columns = range(
            max(0, center - 2),
            min(len(sequence), center + 3),
        )

        window = "".join(
            sequence[column]
            for column in columns
        )

        if len(window) < 5:
            window = window.ljust(5, "-")

        nongap = sum(
            residue != "-"
            for residue in window
        )

        if nongap == 5:
            coverage = "FULL_COVERAGE"
        elif nongap == 0:
            coverage = "NO_COVERAGE"
        else:
            coverage = "PARTIAL_COVERAGE"

        rows.append(
            {
                "tree_order": order,
                "tree_id": record_id,
                "maize_position": position,
                "alignment_columns_1based": (
                    f"{center - 1}-"
                    f"{center + 3}"
                ),
                "window": window,
                "coverage": coverage,
            }
        )


coordinate_output = (
    release
    / "04_maize_coordinate_windows_tree_order.tsv"
)

with coordinate_output.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "tree_order",
            "tree_id",
            "maize_position",
            "alignment_columns_1based",
            "window",
            "coverage",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


with (release / "FORMAL_QC.tsv").open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\t"
        "step25kz16.15_translated_pepc_rebuild\n"
    )
    handle.write("status\tCOMPLETE\n")
    handle.write(
        f"alignment_tip_count\t{len(records)}\n"
    )
    handle.write(
        f"tree_terminal_count\t{len(tree_ids)}\n"
    )
    handle.write(
        f"SS_g8156_in_alignment\t"
        f"{'YES' if 'OUR__SS__g8156' in records else 'NO'}\n"
    )
    handle.write(
        f"SS_g8156_alignment_length\t"
        f"{len(records['OUR__SS__g8156'])}\n"
    )
    handle.write(
        "SS_g24599\tEXCLUDED\n"
    )
    handle.write(
        "coordinate_reference\tmaize PEPC\n"
    )
    handle.write(
        "coordinate_positions\t519,780,890\n"
    )
    handle.write(
        "translation_rule\tDNA-like sequences translated frame 1\n"
    )
    handle.write(
        "window_rule\tfive aligned amino acids centered on maize position\n"
    )

print("Coordinate table created.")
print("Corrected PEPC tree rebuild completed.")
PY

(
  cd "$RELEASE"

  sha256sum \
    01_translation_QC.tsv \
    02_translation_deduplication_manifest.tsv \
    03_corrected_translated_PEPC_alignment.faa \
    04_maize_coordinate_windows_tree_order.tsv \
    FORMAL_QC.tsv \
    PEPC_corrected_translated.treefile \
    > SHA256SUMS.txt

  sha256sum -c SHA256SUMS.txt
)
