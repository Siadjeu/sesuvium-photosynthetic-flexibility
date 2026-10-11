#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_9d
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=08_focused_families/validation/step25kz16.9d/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.9d/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

cd "$PROJECT_ROOT"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate focused_family_env

for command_name in diamond hmmbuild hmmsearch python3 sha256sum; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    }
done

RELEASE="08_focused_families/validation/step25kz16.9d/release"
WORKDIR="$RELEASE/work"

mkdir -p "$WORKDIR"

INPUT_FASTA="08_focused_families/validation/step25kz16.9c/release/work/all_PEPC_sequences_deduplicated.faa"

PEPC_PROFILE_ALIGNMENT="08_focused_families/validation/step25kz12/release/PEPC.comprehensive.min500.alignment.faa"

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "ERROR: missing deduplicated PEPC FASTA" >&2
    exit 1
fi

if [[ ! -s "$PEPC_PROFILE_ALIGNMENT" ]]; then
    echo "ERROR: missing PEPC profile alignment" >&2
    exit 1
fi

export RELEASE WORKDIR INPUT_FASTA PEPC_PROFILE_ALIGNMENT

python3 <<'PY'
import os
import re
from pathlib import Path

workdir = Path(os.environ["WORKDIR"])
input_fasta = Path(os.environ["INPUT_FASTA"])

query_id = "OUR__SS__g24599"


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
                sequence.append(
                    re.sub(r"\s+", "", line).upper()
                )

    if header is not None:
        records[header] = "".join(sequence)

    return records


records = read_fasta(input_fasta)

if query_id not in records:
    raise SystemExit(
        f"ERROR: {query_id} not found in deduplicated FASTA"
    )

query = records[query_id].replace("-", "")

targets = {
    sequence_id: sequence.replace("-", "")
    for sequence_id, sequence in records.items()
    if sequence_id != query_id
}

query_path = workdir / "SS_g24599_query.faa"
target_path = workdir / "all_other_PEPC_targets.faa"

with query_path.open("w", encoding="utf-8") as handle:
    handle.write(f">{query_id}\n")

    for start in range(0, len(query), 80):
        handle.write(query[start:start + 80] + "\n")

with target_path.open("w", encoding="utf-8") as handle:
    for sequence_id, sequence in targets.items():
        handle.write(f">{sequence_id}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")

print(f"Query length: {len(query)} aa")
print(f"Target sequences: {len(targets)}")
PY

diamond makedb \
    --in "$WORKDIR/all_other_PEPC_targets.faa" \
    --db "$WORKDIR/all_other_PEPC_targets"

diamond blastp \
    --query "$WORKDIR/SS_g24599_query.faa" \
    --db "$WORKDIR/all_other_PEPC_targets" \
    --out "$RELEASE/01_SS_g24599_DIAMOND_hits.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp \
    --evalue 1e-20 \
    --max-target-seqs 50 \
    --threads "$SLURM_CPUS_PER_TASK"

hmmbuild \
    "$WORKDIR/PEPC_profile.hmm" \
    "$PEPC_PROFILE_ALIGNMENT" \
    > "$RELEASE/02_hmmbuild.log" 2>&1

hmmsearch \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --domtblout "$RELEASE/03_SS_g24599_HMMER_domains.tbl" \
    "$WORKDIR/PEPC_profile.hmm" \
    "$WORKDIR/SS_g24599_query.faa" \
    > "$RELEASE/03_SS_g24599_HMMER_search.txt" 2>&1

python3 <<'PY'
import csv
import os
from pathlib import Path

release = Path(os.environ["RELEASE"])

diamond_path = release / "01_SS_g24599_DIAMOND_hits.tsv"
summary_path = release / "04_SS_g24599_homology_summary.tsv"

hits = []

if diamond_path.exists():
    with diamond_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) != 10:
                continue

            (
                query,
                subject,
                pident,
                length,
                qlen,
                slen,
                evalue,
                bitscore,
                qcovhsp,
            ) = fields[0:9]

            hits.append({
                "subject": subject,
                "pident": float(pident),
                "length": int(length),
                "qlen": int(qlen),
                "slen": int(slen),
                "evalue": evalue,
                "bitscore": float(bitscore),
                "qcovhsp": float(qcovhsp),
            })

hits.sort(
    key=lambda row: (
        row["bitscore"],
        row["pident"],
        row["qcovhsp"],
    ),
    reverse=True,
)

with summary_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

    writer.writerow([
        "metric",
        "value",
    ])

    writer.writerow([
        "query",
        "OUR__SS__g24599",
    ])

    writer.writerow([
        "diamond_hit_count",
        len(hits),
    ])

    if hits:
        top = hits[0]

        writer.writerow([
            "top_hit",
            top["subject"],
        ])

        writer.writerow([
            "top_hit_pident",
            top["pident"],
        ])

        writer.writerow([
            "top_hit_bitscore",
            top["bitscore"],
        ])

        writer.writerow([
            "top_hit_query_coverage",
            top["qcovhsp"],
        ])

        strong_hits = [
            row for row in hits
            if row["pident"] >= 25
            and row["qcovhsp"] >= 50
            and row["bitscore"] >= 150
            and float(row["evalue"]) <= 1e-20
        ]

        writer.writerow([
            "strong_diamond_hit_count",
            len(strong_hits),
        ])

        writer.writerow([
            "diamond_interpretation",
            "PEPC_LIKE" if strong_hits else "NO_STRONG_PEPC_HIT",
        ])
    else:
        writer.writerow([
            "diamond_interpretation",
            "NO_DIAMOND_HITS",
        ])

    writer.writerow([
        "hmmer_result",
        "inspect_03_SS_g24599_HMMER_search.txt_and_domains_tbl",
    ])

    writer.writerow([
        "coordinate_interpretation",
        "NO_COVERAGE_windows_are_not_substitutions",
    ])

qc_path = release / "FORMAL_QC.tsv"

with qc_path.open("w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\tstep25kz16.9d_ss_g24599_homology_audit\n"
    )
    handle.write("status\tCOMPLETE_REVIEW_REQUIRED\n")
    handle.write(
        "purpose\t"
        "independent PEPC homology assessment for SS__g24599\n"
    )
    handle.write(
        "methods\tDIAMOND plus HMMER PEPC profile search\n"
    )
    handle.write(
        "decision_rule\t"
        "retain in tree but do not interpret unsupported coordinate windows\n"
    )

print(f"Created: {summary_path}")
print(f"Created: {qc_path}")
PY

cd "$RELEASE"

sha256sum \
    01_SS_g24599_DIAMOND_hits.tsv \
    02_hmmbuild.log \
    03_SS_g24599_HMMER_domains.tbl \
    03_SS_g24599_HMMER_search.txt \
    04_SS_g24599_homology_summary.tsv \
    FORMAL_QC.tsv \
    > SHA256SUMS.txt

sha256sum -c SHA256SUMS.txt

column -t -s $'\t' FORMAL_QC.tsv
column -t -s $'\t' 04_SS_g24599_homology_summary.tsv
