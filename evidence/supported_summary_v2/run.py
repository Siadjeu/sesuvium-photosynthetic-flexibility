#!/usr/bin/env python3
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(os.environ["PROJECT_ROOT"]).resolve()
source = root / (
    "08_focused_families/downstream_supported/evidence_correction/"
    "release_21182796_20261023T225841120561Z"
)
stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
out = root / (
    "08_focused_families/downstream_supported/summary_figures_v2"
) / ("release_" + os.environ.get("SLURM_JOB_ID", "local") + "_" + stamp)

def require(condition, message):
    if not condition:
        raise SystemExit("STOP: " + message)

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write(name, records, fields):
    with (out / name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)

require(source.is_dir(), "The reviewed correction release is missing.")
for line in (source / "SHA256SUMS.txt").read_text().splitlines():
    expected, name = line.split("  ", 1)
    path = (source / name).resolve()
    require(path.is_relative_to(source.resolve()), "Unsafe checksum path.")
    require(sha(path) == expected, "Checksum mismatch: " + name)

summary = read(source / "02_family_species_summary_corrected.tsv")
species = ["AT", "BV", "MC", "PO", "SLP_Hap1", "SP_Hex", "SP_Tet", "SS"]
families = list(dict.fromkeys(r["wang_family"] for r in summary))
lookup = {(r["wang_family"], r["species"]): r for r in summary}
require(len(lookup) == len(summary) == len(families) * len(species),
        "Incomplete or duplicate family/species summary.")

metrics = [
    ("WGD_segmental", "WGD/segmental"),
    ("positive_selection", "Positive selection"),
    ("CAFE_primary_pipeline_expansion", "CAFE5 primary"),
    ("DE_any_contrast", "DE: any contrast"),
    ("CAFE_robust_pipeline_expansion", "CAFE5 both models"),
]

out.mkdir(parents=True, exist_ok=False)
for folder in ["tables", "figures", "methods", "inputs"]:
    (out / folder).mkdir()

for name in [
    "01_supported_locus_evidence_corrected.tsv",
    "02_family_species_summary_corrected.tsv",
    "09_family_species_gene_lists_corrected.tsv",
    "10_remaining_assessment_limits.tsv",
]:
    shutil.copyfile(source / name, out / "tables" / name)
shutil.copyfile(source / "METHODS_AND_LIMITATIONS.txt",
                out / "methods/SOURCE_CORRECTION_METHODS.txt")
shutil.copyfile(source / "QC.json", out / "inputs/CORRECTION_QC.json")
shutil.copyfile(source / "SHA256SUMS.txt",
                out / "inputs/CORRECTION_SOURCE_SHA256SUMS.txt")

plot_rows = []
cells = {}
for family in families:
    for sp in species:
        r = lookup[(family, sp)]
        for metric, label in metrics:
            a = r[metric + "_assessed"]
            p = r[metric + "_positive"]
            n = r["supported_loci"]
            if n == "NA":
                status, value, annotation = "DEFINITION_BLOCKED", np.nan, "NA"
            elif int(n) == 0:
                status, value, annotation = "NO_SUPPORTED_LOCI", np.nan, "—"
            elif a == "NA" or int(a) == 0 or p == "NA":
                status, value, annotation = "NO_ASSESSABLE_RESULT", np.nan, "NA"
            else:
                assessed, positive = int(a), int(p)
                require(0 <= positive <= assessed <= int(n),
                        "Invalid denominator: " + family + "/" + sp + "/" + metric)
                status = "ASSESSABLE"
                value = 100 * positive / assessed
                annotation = f"{positive}/{assessed}"
            cells[(family, sp, metric)] = (value, annotation)
            plot_rows.append(dict(
                wang_family=family, species=sp, metric=metric,
                supported_loci=n, assessed=a, positive=p,
                percent_of_assessed="NA" if np.isnan(value) else f"{value:.6f}",
                display_status=status, cell_label=annotation
            ))

write("tables/01_heatmap_values.tsv", plot_rows, [
    "wang_family", "species", "metric", "supported_loci",
    "assessed", "positive", "percent_of_assessed",
    "display_status", "cell_label",
])

plt.rcParams.update({
    "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none",
})
cmap = plt.get_cmap("YlGnBu").copy()
cmap.set_bad("#e5e5e5")
figure_index = []

def draw(values, annotations, xlabels, ylabels, title, stem):
    values = np.asarray(values, dtype=float)
    nr, nc = values.shape
    # Fixed physical cell dimensions within every heatmap.
    cell = 0.65
    left, right, bottom, top = 1.5, 1.3, 1.65, 0.7
    width, height = left + nc * cell + right, bottom + nr * cell + top
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([left/width, bottom/height,
                       nc*cell/width, nr*cell/height])
    image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap,
                      vmin=0, vmax=100, aspect="equal",
                      interpolation="nearest")
    ax.set_xticks(range(nc))
    ax.set_xticklabels(xlabels, rotation=45, ha="right")
    ax.set_yticks(range(nr))
    ax.set_yticklabels(ylabels)
    ax.set_title(title, pad=14)
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(nr):
        for j in range(nc):
            v = values[i, j]
            ax.text(j, i, annotations[i][j], ha="center", va="center",
                    fontsize=8, color="white" if np.isfinite(v) and v > 60
                    else "black")
    cbax = fig.add_axes([(left+nc*cell+0.25)/width, bottom/height,
                        0.16/width, nr*cell/height])
    fig.colorbar(image, cax=cbax, label="% of assessed loci")
    fig.text(
        0.03, 0.025,
        "Cells: positive/assessed. NA: unavailable assessment or blocked definition.\n"
        "—: no supported loci. Categories overlap; do not sum them.",
        fontsize=8,
    )
    for extension in ("pdf", "png", "svg"):
        filename = "figures/" + stem + "." + extension
        fig.savefig(out / filename, dpi=300)
        figure_index.append(dict(figure=filename, title=title))
    plt.close(fig)

for metric, label in metrics:
    values = [[cells[(f, s, metric)][0] for s in species] for f in families]
    labels = [[cells[(f, s, metric)][1] for s in species] for f in families]
    draw(values, labels, species, families, label,
         "all_families_" + metric)

for family in families:
    if not any(lookup[(family, s)]["supported_loci"] not in ("NA", "0")
               for s in species):
        continue
    values = [[cells[(family, s, m)][0] for m, _ in metrics] for s in species]
    labels = [[cells[(family, s, m)][1] for m, _ in metrics] for s in species]
    draw(values, labels, [label for _, label in metrics], species,
         family + ": supported-locus evidence",
         family + "_species_evidence")

write("FIGURE_INDEX.tsv", figure_index, ["figure", "title"])

methods = """SUPPORTED-FAMILY DOWNSTREAM EVIDENCE: METHODS

Scope and unit of analysis
The analysis was restricted to the curated set of 657 supported species/locus
pairs. Gene-family assignment does not establish a demonstrated C4 or CAM
function. Other evidence tiers were excluded from these summaries. MEP and PPA
remained definition-blocked. Species codes were AT, BV, MC, PO, SLP_Hap1,
SP_Hex, SP_Tet and SS; AT denotes Amaranthus tricolor.

Existing duplication, gene-family evolution, positive-selection and differential-
expression results were integrated using the established gene identifiers and
orthogroup mappings. No statistical analyses were rerun in this reporting step.

WGD/segmental duplication
Supported loci with an interpretable duplication classification formed the
assessment denominator. Loci classified as WGD_SEGMENTAL formed the numerator.
This classification combines whole-genome and segmental duplication evidence;
it does not distinguish the two mechanisms.

CAFE5 primary expansion
The primary model was M4_separate_sesuvium_terminals, error-corrected,
replicate 6. The sensitivity model was M2_tetraploid_hexaploid, without error
correction, replicate 3. Family probabilities and significance flags were read
from their existing Base_family_results.txt outputs. Recorded significance
flags were retained; probabilities were not retested or adjusted within the
supported-gene subset.

A positive primary expansion outcome required both a significant orthogroup-
level result and a significant expansion on the relevant terminal species
branch. A non-significant family result was classified as a negative outcome
for this combined criterion, even when no branch probability was reported.
For a significant family, a non-significant branch or a significant branch
without expansion was also negative. A significant family lacking an
interpretable branch result remained unresolved.

The denominator comprised supported loci with a decidable outcome under this
combined criterion. It therefore differs from the number of loci with explicit
branch probabilities. Family-level and branch-level evidence were retained
separately in the locus table. Missing model results were not converted to
negative outcomes.

Robust expansion across models
A robust positive required the expansion criterion to be satisfied in both
models. A negative outcome in either model ruled out this conjunction. If
neither model established a negative outcome and one remained unresolved,
the combined result remained unresolved. A decidable robust negative does not
mean that both models supplied a branch-level test.

CAFE5 results concern orthogroups and their branches. Reported gene counts are
the supported loci belonging to those orthogroups; they are not independent
statistical tests or counts of independent expansion events.

Positive selection
Existing valid aBSREL branch results were retained with their original global
Benjamini-Hochberg q-values. A locus was assessed when it had at least one valid
linked test with an interpretable global q-value, and positive when at least
one such result met global q < 0.05. The global correction was not recalculated
using only the supported-family subset.

BV and SLP_Hap1 were background-only under the original test design. AT was a
secondary foreground. Lack of a linked valid test was not interpreted as
absence of positive selection. For other unassessed loci, precise candidate-
selection or test-eligibility reasons remain unresolved unless documented.

Differential expression
DE evidence was restricted to the 14 previously authorized contrasts: 12 in
PO and two in SS. A gene/contrast result was assessable when both the within-
contrast adjusted p-value and the apeglm-shrunken log2 fold change were finite.
A positive result required adjusted p < 0.05 and absolute shrunken log2 fold
change >= 1.

The displayed 'DE: any contrast' category counted a locus once if it was
positive in at least one authorized contrast. Its denominator included loci
assessable in at least one contrast, not necessarily all contrasts. DE evidence
for the other species was unavailable. This union across contrasts must not be
interpreted as a uniform drought-response or CAM-response classification.

Heatmaps and denominators
For each family, species and evidence category, the percentage was calculated
as 100 times the number of positive supported loci divided by the number of
assessable supported loci. Cells display positive/assessed; colour spans
0-100 percent on the same scale for every figure.

Grey NA cells indicate unavailable assessments or blocked family definitions.
An em dash indicates no supported loci. Zero percent means that at least one
locus was assessed but none met the relevant positive criterion. Separate
columns record supported totals, assessed totals and display status.

Categories overlap. A locus may have duplication, expansion, selection and DE
evidence simultaneously; category counts must not be summed as unique genes.
These summaries are descriptive and do not constitute enrichment tests.

Reproducibility
Inputs were taken from correction release
release_21182796_20261023T225841120561Z and verified against its SHA256 manifest.
The release includes corrected locus evidence, family/species summaries, gene
lists, unresolved assessment limits, plotting values and this plotting script.
The source release and published repository were not modified.
"""
(out / "methods/METHODS.txt").write_text(methods)

(out / "methods/RESULTS_AND_LIMITATIONS.txt").write_text(
    "The supported set contains 657 loci. The primary CAFE5 expansion criterion "
    "was assessable for 629 loci, with nine positive loci, all previously "
    "identified PO PEPC members of one orthogroup. The remaining 28 loci "
    "comprised 13 without an assigned orthogroup and 15 without a recovered "
    "result in the selected models. No locus satisfied the expansion criterion "
    "in both models; the combined outcome was decidable for 629 loci.\n\n"
    "Positive selection remained six positive loci among 296 assessed. "
    "The unassessed set comprised 124 background-only BV/SLP_Hap1 loci, "
    "63 AT loci without a matched valid secondary-foreground test, and "
    "174 other loci with unresolved eligibility reasons. No additional "
    "positive-selection calls were inferred from missing results.\n\n"
    "Original WGD/segmental and DE evidence was unchanged. The expanded CAFE5 "
    "denominator reflects recovery of family-level negative evidence, not new "
    "branch tests. Unequal assessment coverage limits comparisons among species "
    "and families. Older and new CAFE5 gene-list labels may identify the same "
    "loci and must not be added together.\n"
)

shutil.copyfile(Path(__file__), out / "run.py")
qc = dict(
    status="COMPLETE",
    family_species_rows=len(summary),
    plotted_value_rows=len(plot_rows),
    figures=len(figure_index)//3,
    output_formats=["pdf", "png", "svg"],
    matplotlib_version=matplotlib.__version__,
    numpy_version=np.__version__,
    biological_analyses_rerun=False,
    source_files_modified=False,
    github_modified=False,
)
(out / "QC.json").write_text(json.dumps(qc, indent=2) + "\n")
files = sorted(p for p in out.rglob("*") if p.is_file())
(out / "SHA256SUMS.txt").write_text("".join(
    sha(p) + "  " + str(p.relative_to(out)) + "\n" for p in files
))
bundle = out.parent / (out.name + "_bundle.tar.gz")
with tarfile.open(bundle, "w:gz") as tar:
    tar.add(out, arcname=out.name)
print(json.dumps(qc, indent=2))
print("\nUpload this bundle:\n" + str(bundle))
