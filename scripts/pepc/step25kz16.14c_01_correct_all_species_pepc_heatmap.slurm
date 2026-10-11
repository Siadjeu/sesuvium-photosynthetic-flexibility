#!/usr/bin/env bash
#SBATCH --job-name=step25kz16_14c
#SBATCH --partition=all_cpu.p
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=08_focused_families/validation/step25kz16.14c/logs/%x_%j.out
#SBATCH --error=08_focused_families/validation/step25kz16.14c/logs/%x_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?PROJECT_ROOT must be exported}"

cd "$PROJECT_ROOT"

export RELEASE="08_focused_families/validation/step25kz16.14c/release"
export CURRENT="08_focused_families/validation/step25kz16.14b/release/01_all_species_PEPC_evidence_heatmap_corrected.tsv"
export OLD_METADATA="08_focused_families/validation/step25kz16.8b/release/01_PEPC_tree_tip_metadata.tsv"

python3 <<'PY'
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle


release = Path(os.environ["RELEASE"])
current_path = Path(os.environ["CURRENT"])
old_path = Path(os.environ["OLD_METADATA"])

species_order = [
    "AT", "BV", "MC", "PO",
    "SLP_Hap1", "SP_Hex", "SP_Tet", "SS",
]

evidence_order = [
    "WGD/segmental",
    "CAFE5",
    "Positive selection",
    "Differential expression",
]

columns = {
    "WGD/segmental": "WGD_segmental",
    "CAFE5": "CAFE5",
    "Positive selection": "positive_selection",
    "Differential expression": "DE_PO_SS_only",
}

duplicate_map = {
    ("SS", "g29856.t1"): "g11076.t1",
    ("SS", "g7519.t1"): "g17009.t1",
}


def read_tsv(path):
    if not path.exists():
        raise SystemExit(f"Missing input: {path}")

    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames or []


def normalize(species, value):
    value = str(value or "").strip().split()[0]

    for prefix in [
        f"OUR__{species}__",
        f"{species}__",
        "OUR__",
    ]:
        if value.startswith(prefix):
            value = value[len(prefix):]

    return value


current_rows, current_fields = read_tsv(current_path)
old_rows, old_fields = read_tsv(old_path)

required_current = {
    "species",
    "gene_id",
    "duplicate_of",
    "DE_PO_SS_only",
}

missing = required_current - set(current_fields)

if missing:
    raise SystemExit(
        "Missing current-table columns: "
        + ", ".join(sorted(missing))
    )

required_old = {
    "tip_id",
    "species",
    "WGD_segmental",
    "CAFE5",
    "positive_selection",
}

missing = required_old - set(old_fields)

if missing:
    raise SystemExit(
        "Missing old-metadata columns: "
        + ", ".join(sorted(missing))
    )


old_index = {}

for row in old_rows:
    species = row["species"].strip()

    if species in species_order:
        old_index[
            (species, normalize(species, row["tip_id"]))
        ] = row


def is_yes(value):
    return str(value).strip().upper() in {
        "YES", "TRUE", "Y", "1"
    }


final_rows = []

for row in current_rows:
    species = row["species"].strip()
    gene_id = row["gene_id"].strip()

    representative = (
        row["duplicate_of"].strip()
        or gene_id
    )

    source = old_index.get(
        (species, normalize(species, representative))
    )

    if source is None:
        raise SystemExit(
            f"Missing validated metadata for "
            f"{species} {representative}"
        )

    final_rows.append(
        {
            "species": species,
            "gene_id": gene_id,
            "duplicate_of": row["duplicate_of"],
            "WGD_segmental": (
                "YES"
                if is_yes(source["WGD_segmental"])
                else "NO"
            ),
            "CAFE5": (
                "YES"
                if is_yes(source["CAFE5"])
                else "NO"
            ),
            "positive_selection": (
                "YES"
                if is_yes(source["positive_selection"])
                else "NO"
            ),
            "DE_PO_SS_only": row["DE_PO_SS_only"],
        }
    )


gene_output = (
    release / "01_all_species_PEPC_evidence_heatmap_final.tsv"
)

with gene_output.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species",
            "gene_id",
            "duplicate_of",
            "WGD_segmental",
            "CAFE5",
            "positive_selection",
            "DE_PO_SS_only",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(final_rows)


summary_rows = []

for species in species_order:
    rows = [
        row for row in final_rows
        if row["species"] == species
    ]

    total = len(rows)

    for evidence in evidence_order:
        if (
            evidence == "Differential expression"
            and species not in {"PO", "SS"}
        ):
            summary_rows.append(
                {
                    "species": species,
                    "evidence": evidence,
                    "n": "NA",
                    "nt": "NA",
                    "percentage": "NA",
                }
            )
            continue

        field = columns[evidence]

        n = sum(
            row[field] == "YES"
            for row in rows
        )

        percentage = (
            100.0 * n / total
            if total
            else 0.0
        )

        summary_rows.append(
            {
                "species": species,
                "evidence": evidence,
                "n": str(n),
                "nt": str(total),
                "percentage": f"{percentage:.1f}",
            }
        )


summary_output = (
    release / "02_PEPC_species_summary_n_over_nt_final.tsv"
)

with summary_output.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species",
            "evidence",
            "n",
            "nt",
            "percentage",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(summary_rows)


cmap = LinearSegmentedColormap.from_list(
    "blue_white_red",
    [
        "#2166ac",
        "#f7f7f7",
        "#b2182b",
    ],
)

fig, axis = plt.subplots(
    figsize=(13.5, 7.0),
    facecolor="white",
)

for y, species in enumerate(species_order):
    for x, evidence in enumerate(evidence_order):
        item = next(
            row for row in summary_rows
            if row["species"] == species
            and row["evidence"] == evidence
        )

        if item["percentage"] == "NA":
            facecolor = "#d9d9d9"
            text_color = "black"
            label = "NA"
        else:
            percentage = float(item["percentage"])

            facecolor = cmap(
                Normalize(0, 100)(percentage)
            )

            text_color = (
                "white"
                if percentage < 25
                or percentage > 75
                else "black"
            )

            label = (
                f"{item['n']}/{item['nt']}\n"
                f"{item['percentage']}%"
            )

        axis.add_patch(
            Rectangle(
                (x, y),
                1,
                1,
                facecolor=facecolor,
                edgecolor="black",
                linewidth=1.15,
            )
        )

        axis.text(
            x + 0.5,
            y + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color=text_color,
            linespacing=1.15,
        )

axis.set_xlim(0, len(evidence_order))
axis.set_ylim(0, len(species_order))
axis.invert_yaxis()

axis.set_xticks(
    [x + 0.5 for x in range(len(evidence_order))]
)

axis.set_xticklabels(
    [
        "WGD/segmental",
        "CAFE5",
        "Positive selection",
        "Differential expression",
    ],
    rotation=25,
    ha="right",
    fontsize=11,
    fontweight="bold",
)

axis.set_yticks(
    [y + 0.5 for y in range(len(species_order))]
)

axis.set_yticklabels(
    species_order,
    fontsize=11,
    fontweight="bold",
)

axis.tick_params(
    axis="both",
    length=0,
)

axis.set_xlabel(
    "Evidence category",
    fontsize=12,
    fontweight="bold",
)

axis.set_ylabel(
    "Species",
    fontsize=12,
    fontweight="bold",
)

axis.set_title(
    "PEPC evidence across eight species",
    fontsize=14,
    fontweight="bold",
    pad=16,
)

scalar = ScalarMappable(
    norm=Normalize(0, 100),
    cmap=cmap,
)

scalar.set_array([])

colorbar = fig.colorbar(
    scalar,
    ax=axis,
    fraction=0.045,
    pad=0.03,
)

colorbar.set_label(
    "Genes supported (%)",
    fontsize=11,
    fontweight="bold",
)

fig.subplots_adjust(
    left=0.12,
    right=0.88,
    bottom=0.22,
    top=0.90,
)

figure_prefix = (
    release / "03_PEPC_species_summary_n_over_nt_final"
)

for extension in ["pdf", "svg", "png"]:
    fig.savefig(
        figure_prefix.with_suffix(f".{extension}"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

plt.close(fig)


with (release / "FORMAL_QC.tsv").open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("metric\tvalue\n")
    handle.write(
        "checkpoint\t"
        "step25kz16.14c_final_pepc_evidence_heatmap\n"
    )
    handle.write("status\tCOMPLETE\n")
    handle.write("species_count\t8\n")
    handle.write(
        "cell_format\tn/nt first line; percentage second line\n"
    )
    handle.write(
        "color_scale\tblue-white-red based on percentage\n"
    )
    handle.write(
        "WGD_source\tStep25kz16.8b validated metadata\n"
    )
    handle.write(
        "DE_source\tStep25kz16.14b PO/SS values\n"
    )
    handle.write(
        "DE_assessed_species\tPO,SS\n"
    )
    handle.write(
        "DE_other_species\tNA\n"
    )
    handle.write(
        "SS_g24599\tEXCLUDED\n"
    )
    handle.write(
        "duplicate_gene_rule\t"
        "identical SS copies retained separately\n"
    )
    handle.write(
        "output_formats\tPDF, SVG, PNG, 300 dpi\n"
    )

print("Step25kz16.14c outputs created.")
PY

(
  cd "$RELEASE"

  sha256sum \
    01_all_species_PEPC_evidence_heatmap_final.tsv \
    02_PEPC_species_summary_n_over_nt_final.tsv \
    03_PEPC_species_summary_n_over_nt_final.pdf \
    03_PEPC_species_summary_n_over_nt_final.svg \
    03_PEPC_species_summary_n_over_nt_final.png \
    FORMAL_QC.tsv \
    > SHA256SUMS.txt

  sha256sum -c SHA256SUMS.txt
)
