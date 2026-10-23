# Reproducibility and available inputs

Publication preparation copies reviewed evidence, verifies hashes and compiles
Python sources for syntax. It does not rerun biological tests or install tools.
The paths below are relative to the publication root. Activate your own software
environment. Use fresh output directories outside the publication snapshot.

## Verify the publication

```bash
python3 scripts/verify_publication.py
```

This verifies the curated snapshot, not biological validity. It replaces using
historical freeze-verification entry points that expect omitted archive files.

## Rebuild the Wang-layout evidence tables

Python standard library:

```bash
python3 evidence/family_validation/scripts/build_family_tables.py \
  --release evidence/family_validation --out ../family_tables_rebuilt
```

The original ledger and required input checksums are preserved. Definition-blocked
families and evidence tiers must not be silently collapsed.

## Replot the accepted replicate-mean expression heatmaps

Requires NumPy and Matplotlib:

```bash
python3 scripts/replot_expression.py --out ../expression_replotted
```

This portable wrapper uses the published input matrices and the reviewed v4
plotting calculations/layout. It does not copy earlier trees, fit DE models,
or change underlying gene counts. For each gene it averages three normalized
count replicates per group, then computes `(group_mean - mean_across_groups) /
SD_across_groups`, with population SD (`ddof=0`). Zero-variance rows are masked.
The published original HPC plotter remains in the reproduction_package for
provenance; it expects the original presentation release and PROJECT_ROOT.
Use the portable wrapper above from this publication instead.

## Reproduce SS metabolite figures and descriptive concordance

Requires NumPy, pandas, SciPy and Matplotlib:

```bash
python3 scripts/replot_ss_metabolites.py \
  --component metabolite_comparisons --out ../ss_metabolites_reproduced
python3 scripts/replot_ss_metabolites.py \
  --component transcript_metabolite_comparison --out ../ss_concordance_reproduced
```

The wrapper creates a temporary `data` directory from the frozen `inputs`,
checks the recorded input hashes where supplied, and invokes the archived
analysis script. The first command also reruns the documented small statistical
comparisons. These commands are optional; they are not executed during publication
preparation. No cross-assay biological sample pairing is introduced.

## Integration and PEPC reconstruction scripts

`evidence/supported_downstream/integrate_supported.py` records the exact
integration algorithm and accepts `--source-release` and `--output`. Full
re-execution requires its original source collection, including the complete
representative maps, OrthoFinder, MCScanX, CAFE5, aBSREL and DE result tables.
Those upstream genome-scale inputs are not all duplicated in this curated
publication. `SOURCE_MANIFEST.tsv` and the source selection record their origin.
The frozen integrated tables can be used directly for supported-subset analyses.
No claim is made that this repository alone reproduces every upstream HPC run.

`evidence/pepc/scripts/reconstruction/` retains the corrected reference/tree
workflow. Full reconstruction expects its documented original PROJECT_ROOT
releases plus MAFFT and IQ-TREE. Published corrected proteins, alignments,
reference records, trees, per-tree evidence and locus mappings are retained for
inspection and independent reconstruction. Use main-tree iTOL files for the
manuscript and sensitivity-tree files for supplementary presentation.

## Historical scripts and manifests

Original component scripts are retained as provenance even when their execution
requires external HPC inputs or historical directory layouts. Their presence
does not make all of them portable entry points. Publication SHA256SUMS files
are regenerated for the curated selections; historical source-manifest hashes
in QC describe their original releases. No source release was modified.

Join downstream tables by **species and gene ID**. The PEPC tree has two pairs
of identical SS proteins represented once per pair; expand through the locus
mapping when counting or joining. Preserve supported/provisional/unresolved
subtype status, assessment denominators, NA and untested categories.
