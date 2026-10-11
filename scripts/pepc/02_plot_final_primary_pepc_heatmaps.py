#!/usr/bin/env python3

from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle


INPUT_DIR = Path(
    "08_focused_families/validation/step25kz15d/work/input"
)

OUTDIR = Path(
    "08_focused_families/validation/step25kz15d/release"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


COLORS = {
    "leaf": "#4daf4a",
    "stem": "#377eb8",
    "cotyledon": "#984ea3",
    "drought": "#e41a1c",
    "well_watered": "#377eb8",
    "08H": "#fee8c8",
    "16H": "#fdbb84",
    "24H": "#e34a33",
    "DAY": "#fdd835",
    "NIGHT": "#3949ab",
}


def read_tsv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_matrix(path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)

    genes = [row[0] for row in rows]
    samples = header[1:]

    matrix = []

    for row in rows:
        values = []

        for value in row[1:]:
            try:
                values.append(float(value))
            except ValueError:
                values.append(0.0)

        matrix.append(values)

    return genes, samples, np.asarray(matrix, dtype=float)


def sort_po(samples, annotation):
    tissue_order = {
        "leaf": 0,
        "stem": 1,
    }

    time_order = {
        "08H": 0,
        "16H": 1,
        "24H": 2,
    }

    condition_order = {
        "well_watered": 0,
        "drought": 1,
    }

    return sorted(
        samples,
        key=lambda sample: (
            tissue_order.get(annotation[sample]["tissue"], 99),
            time_order.get(annotation[sample]["time"], 99),
            condition_order.get(
                annotation[sample]["condition"],
                99,
            ),
            sample,
        ),
    )


def sort_ss(samples, annotation):
    tissue_order = {
        "cotyledon": 0,
        "leaf": 1,
    }

    period_order = {
        "DAY": 0,
        "NIGHT": 1,
    }

    return sorted(
        samples,
        key=lambda sample: (
            tissue_order.get(annotation[sample]["tissue"], 99),
            period_order.get(
                annotation[sample]["period"],
                99,
            ),
            sample,
        ),
    )


def draw_legend(axis, field, values, x, y):
    axis.text(
        x,
        y,
        field.capitalize(),
        ha="left",
        va="center",
        fontsize=8,
        fontweight="bold",
    )

    x += 1.7

    for value in values:
        color = COLORS[value]

        axis.add_patch(
            Rectangle(
                (x, y - 0.16),
                0.40,
                0.32,
                facecolor=color,
                edgecolor="black",
                linewidth=0.8,
            )
        )

        axis.text(
            x + 0.55,
            y,
            value.replace("_", " "),
            ha="left",
            va="center",
            fontsize=7,
        )

        x += max(2.0, 0.12 * len(value) + 0.9)


def draw_annotation_bars(axis, annotation, samples, fields):
    n_samples = len(samples)
    n_fields = len(fields)

    colors = []

    for field in fields:
        row = []

        for sample in samples:
            value = annotation[sample][field]
            row.append(
                matplotlib.colors.to_rgb(
                    COLORS[value]
                )
            )

        colors.append(row)

    colors = np.asarray(colors)

    # Use imshow for reliable per-cell colors.
    axis.imshow(
        colors,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
        extent=(0, n_samples, 0, n_fields),
    )

    # Solid black borders around every annotation cell.
    for x in range(n_samples + 1):
        axis.plot(
            [x, x],
            [0, n_fields],
            color="black",
            linewidth=0.75,
            solid_capstyle="butt",
        )

    for y in range(n_fields + 1):
        axis.plot(
            [0, n_samples],
            [y, y],
            color="black",
            linewidth=0.75,
            solid_capstyle="butt",
        )

    axis.set_xlim(0, n_samples)
    axis.set_ylim(0, n_fields + 2.2)
    axis.axis("off")

    # Labels correspond to the actual annotation rows.
    for index, field in enumerate(fields):
        y = n_fields - index - 0.5

        axis.text(
            -0.02,
            y,
            field.capitalize(),
            transform=axis.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=8,
            clip_on=False,
        )

def plot_species(species):
    prefix = f"{species}_PEPC"

    genes, samples, matrix = read_matrix(
        INPUT_DIR / f"{prefix}_expression_zscore.tsv"
    )

    # Keep only the formally classified primary PEPC genes.
    classification_path = Path(
        "08_focused_families/validation/step25kz16.6/release/"
        "01_pepc_candidate_classification.tsv"
    )

    primary_classes = {
        "PRIMARY_FULL_LENGTH_PEPC",
        "PRIMARY_RELAXED_PEPC",
    }

    primary_ids = set()

    with classification_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if (
                row.get("species") == species
                and row.get("classification")
                in primary_classes
            ):
                gene_id = (
                    row.get("tree_id")
                    or row.get("tip_id")
                    or row.get("sequence_id")
                    or ""
                )

                gene_id = gene_id.replace("OUR__", "")
                primary_ids.add(gene_id)
                primary_ids.add(gene_id.removesuffix(".t1"))
                primary_ids.add(gene_id.removesuffix(".t2"))

    excluded_ids = {
        "SS__g24599",
        "SS__g24599.t1",
        "g24599",
        "g24599.t1",
    }

    keep_indices = []

    for index, gene in enumerate(genes):
        normalized = gene.replace("OUR__", "")

        variants = {
            normalized,
            normalized.removesuffix(".t1"),
            normalized.removesuffix(".t2"),
        }

        if species == "SS" and variants & excluded_ids:
            continue

        if variants & primary_ids:
            keep_indices.append(index)

    genes = [genes[index] for index in keep_indices]
    matrix = matrix[keep_indices, :]

    print(
        f"{species}: retained {len(genes)} primary PEPC genes"
    )
    print("Retained genes:", ", ".join(genes))

    annotation_rows = read_tsv(
        INPUT_DIR / f"{prefix}_sample_annotation.tsv"
    )

    annotation = {
        row["sample_id"]: row
        for row in annotation_rows
    }

    if species == "PO":
        fields = [
            "tissue",
            "condition",
            "time",
            "period",
        ]

        legend_values = {
            "tissue": ["leaf", "stem"],
            "condition": ["drought", "well_watered"],
            "time": ["08H", "16H", "24H"],
            "period": ["DAY", "NIGHT"],
        }

        for row in annotation.values():
            row["period"] = row["day_night"]

        samples = sort_po(samples, annotation)

    else:
        # SS has only tissue and diel contrasts.
        fields = [
            "tissue",
            "period",
        ]

        legend_values = {
            "tissue": ["leaf", "cotyledon"],
            "period": ["DAY", "NIGHT"],
        }

        for row in annotation.values():
            row["period"] = row["day_night"]

        samples = sort_ss(samples, annotation)

    sample_indices = [
        list(read_matrix(
            INPUT_DIR / f"{prefix}_expression_zscore.tsv"
        )[1]).index(sample)
        for sample in samples
    ]

    matrix = matrix[:, sample_indices]
    matrix = np.clip(matrix, -2.5, 2.5)

    n_genes = len(genes)
    n_samples = len(samples)
    n_fields = len(fields)

    figure_width = 16.0

    figure_height = max(
        9.0,
        0.30 * n_genes + 4.5,
    )

    fig = plt.figure(
        figsize=(figure_width, figure_height),
        facecolor="white",
    )

    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[36, 1.0],
        height_ratios=[
            max(1.6, 0.55 * n_fields),
            max(5.0, 0.34 * n_genes),
        ],
        hspace=0.04,
        wspace=0.10,
    )

    annotation_axis = fig.add_subplot(grid[0, 0])
    heatmap_axis = fig.add_subplot(grid[1, 0])
    colorbar_axis = fig.add_subplot(grid[1, 1])

    draw_annotation_bars(
        annotation_axis,
        annotation,
        samples,
        fields,
    )

    if species == "PO":
        draw_legend(
            annotation_axis,
            "tissue",
            legend_values["tissue"],
            0,
            5.0,
        )

        draw_legend(
            annotation_axis,
            "condition",
            legend_values["condition"],
            10,
            5.0,
        )

        draw_legend(
            annotation_axis,
            "time",
            legend_values["time"],
            23,
            5.0,
        )

        draw_legend(
            annotation_axis,
            "period",
            legend_values["period"],
            0,
            4.45,
        )

    else:
        draw_legend(
            annotation_axis,
            "tissue",
            legend_values["tissue"],
            0,
            3.1,
        )

        draw_legend(
            annotation_axis,
            "period",
            legend_values["period"],
            0,
            2.55,
        )

    x_edges = np.arange(n_samples + 1)
    y_edges = np.arange(n_genes + 1)

    heatmap_axis.pcolormesh(
        x_edges,
        y_edges,
        matrix,
        cmap="RdBu_r",
        vmin=-2.5,
        vmax=2.5,
        edgecolors="black",
        linewidth=0.45,
        shading="flat",
    )

    heatmap_axis.set_xlim(0, n_samples)
    heatmap_axis.set_ylim(0, n_genes)

    heatmap_axis.set_yticks(
        np.arange(n_genes) + 0.5
    )

    heatmap_axis.set_yticklabels(
        [
            gene.replace("PO__", "").replace("SS__", "")
            for gene in genes
        ],
        fontsize=8,
    )

    sample_labels = []

    for sample in samples:
        row = annotation[sample]

        if species == "PO":
            label = (
                f"{row['tissue']}-"
                f"{row['condition'].replace('_', '-')}-"
                f"{row['time']}-"
                f"{sample.split('_')[-1]}"
            )
        else:
            label = (
                f"{row['tissue']}-"
                f"{row['period']}-"
                f"{sample.split('_')[-1]}"
            )

        sample_labels.append(label)

    heatmap_axis.set_xticks(
        np.arange(n_samples) + 0.5
    )

    heatmap_axis.set_xticklabels(
        sample_labels,
        rotation=90,
        ha="center",
        va="top",
        fontsize=13.0,
    )

    heatmap_axis.tick_params(
        axis="both",
        length=0,
        pad=7,
    )

    heatmap_axis.set_xlabel(
        "Samples",
        fontsize=9,
        labelpad=62,
    )

    heatmap_axis.set_ylabel(
        "PEPC genes",
        fontsize=9,
    )

    for spine in heatmap_axis.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)

    colorbar = fig.colorbar(
        ScalarMappable(
            norm=Normalize(vmin=-2.5, vmax=2.5),
            cmap="RdBu_r",
        ),
        cax=colorbar_axis,
    )

    colorbar.set_label(
        "Expression z-score",
        rotation=270,
        labelpad=13,
        fontsize=9,
    )

    output_prefix = OUTDIR / f"{prefix}_Figure3B_like_final_primary_PEPC"

    for extension in ["pdf", "png", "svg"]:
        output = output_prefix.with_suffix(
            f".{extension}"
        )

        fig.savefig(
            output,
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

        print(f"Created: {output}")

    plt.close(fig)


plot_species("PO")
plot_species("SS")

print("Completed Step25kz15 compact v4 heatmaps")
