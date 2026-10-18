# Supported-family downstream evidence v1

This separate release integrates existing WGD/segmental duplication, CAFE5,
aBSREL positive-selection and PO/SS differential-expression results for 657
supported family loci. Read METHODS_RESULTS.txt before using the tables.

- 02_family_species_summary.tsv: family/species counts and denominators.
- 09_family_species_gene_lists.tsv: genes grouped by evidence category.
- matrix_*.tsv: Wang-style species-by-family tables.
- 04_DE_locus_contrast_evidence.tsv and 05_DE_family_contrast_summary.tsv:
  individual authorized contrasts; pooled DE totals mix biological comparisons.
- 01_supported_locus_evidence.tsv: identifiers and assessment coverage.

NA is not zero. Nine PO PEPC loci belong to a primary-model expanded
orthogroup; sensitivity-model branch tests were unavailable for this subset.
The 13 orthogroup-coverage exceptions are not failed gene-ID mappings.
Supported family membership does not establish a C4/CAM-specific role.

SHA256SUMS.txt verifies the original integration payload.
FREEZE_SHA256SUMS.txt additionally covers this README and freeze script.
SOURCE_BUNDLE_SHA256SUMS.txt applies to the earlier source collection.
Source data remain in the existing HPC collection; source hashes and relative
paths are retained in SOURCE_MANIFEST.tsv. No new biological analysis was run.

## Reproduce the integration

Use Python >=3.9 and the collected source release:

```bash
python3 integrate_supported.py --source-release "$SOURCE_RELEASE" --output "$NEW_OUTPUT"
```

For SLURM, set PROJECT_ROOT, SOURCE_RELEASE and INTEGRATION_SCRIPT_DIR to the
directory containing the integration script, then submit 03_integrate_supported.slurm.
The separate freeze script uses PROJECT_ROOT, optionally GITHUB_REPO_DIR and
SUPPORTED_RELEASE. Run without arguments to prepare; add --push to publish.
The expected verified integration manifest and base commit are pinned.

The parent is the unchanged photosynthesis-family-evidence-v1.0.0 freeze.
This branch/tag does not update main. Version tags are immutable by convention;
future corrections require a new version, not replacement of this tag.
