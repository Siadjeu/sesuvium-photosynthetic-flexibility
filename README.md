# Sesuvium photosynthetic flexibility: curated evidence

## Corrected downstream summaries

Use [the corrected count/percentage heatmaps](evidence/supported_summary_v2/figures/)
and [figure index](evidence/supported_summary_v2/FIGURE_INDEX.tsv) for the current
supported-family downstream presentation.

- [Corrected locus evidence](evidence/supported_summary_v2/tables/01_supported_locus_evidence_corrected.tsv)
- [Corrected family/species summaries](evidence/supported_summary_v2/tables/02_family_species_summary_corrected.tsv)
- [Gene lists](evidence/supported_summary_v2/tables/09_family_species_gene_lists_corrected.tsv)
- [Exact plotted values](evidence/supported_summary_v2/tables/01_heatmap_values.tsv)
- [Methods](evidence/supported_summary_v2/methods/METHODS.txt)
- [Results and limitations](evidence/supported_summary_v2/methods/RESULTS_AND_LIMITATIONS.txt)

CAFE5 primary expansion is positive for nine of 629 assessable supported loci;
28 remain unresolved. Expansion supported by both models is zero of 629
assessable loci. Positive selection remains six of 296 assessed loci.
CAFE5 gene counts describe orthogroup members, not independent expansion events.

These corrected CAFE5 pipeline summaries supersede the earlier branch-only
denominators for this presentation. Earlier branch-level evidence, DE contrast
tables and original releases remain available under `evidence/supported_downstream/`.
Use the new `CAFE_*_pipeline_expansion_*` columns for corrected CAFE5 summaries.
Old and new expansion gene-list category names can identify the same loci;
do not add them together. Unresolved positive-selection eligibility remains
documented rather than treated as a negative result.


Comparative gene-family evidence, evolutionary analyses, expression summaries,
and exploratory metabolite concordance for Sesuvium and related species.
This publication snapshot replaces the earlier PEPC-only working presentation.

## Start here

| Component | Authoritative files |
|---|---|
| Family definitions, evidence tiers and Wang-layout tables | [Family evidence](evidence/family_validation/README.md) and [generated tables](evidence/family_validation/tables/generated/) |
| Supported-family duplication, CAFE5, positive selection and DE | [Methods/results](evidence/supported_downstream/METHODS_RESULTS.txt), [657-locus ledger](evidence/supported_downstream/01_supported_locus_evidence.tsv), [family/species summaries](evidence/supported_downstream/02_family_species_summary.tsv) |
| Per-contrast PO/SS differential expression | [Gene/contrast evidence](evidence/supported_downstream/04_DE_locus_contrast_evidence.tsv) and [family/contrast summaries](evidence/supported_downstream/05_DE_family_contrast_summary.tsv) |
| Positive-selection tests and CAFE5 membership | [Selection tests](evidence/supported_downstream/06_positive_selection_tests.tsv), [CAFE5 evidence](evidence/supported_downstream/07_CAFE_model_evidence.tsv) and [gene lists](evidence/supported_downstream/09_family_species_gene_lists.tsv) |
| Accepted expression heatmaps | [Replicate-mean z-score figures/tables](evidence/expression/expression_mean_zscore/), [figure index](evidence/expression/FIGURE_INDEX.tsv) |
| Corrected broad PEPC phylogeny | [Main tree](evidence/pepc/itol/main/01_PEPC.treefile), [iTOL instructions](evidence/pepc/ITOL_START_HERE.txt), [51-locus subtype table](evidence/pepc/tables/10_locus_subtypes.tsv) |
| SS transcript/metabolite comparison | [Methods](evidence/ss_metabolite/docs/METHODS.txt), [results](evidence/ss_metabolite/docs/RESULTS.txt), [limitations](evidence/ss_metabolite/docs/LIMITATIONS.txt), [figures/evidence](evidence/ss_metabolite/evidence/) |

Read [reproducibility instructions](docs/REPRODUCIBILITY.md) and
[results and interpretation limits](docs/RESULTS_AND_LIMITATIONS.md).

## Scope and interpretation

The supported subset contains **657 distinct species/locus pairs**, not an
exhaustive inventory of functionally demonstrated C4/CAM genes. The wider family
ledger retains other evidence tiers. **MEP and PPA remain definition-blocked**;
unavailable values must remain NA, not zero. Species codes are AT, BV, MC, PO,
SLP_Hap1, SP_Hex, SP_Tet and SS. AT denotes Amaranthus tricolor, not Arabidopsis.

The integrated evidence records 465 WGD/segmental loci, nine PO PEPC loci in one
primary-model significantly expanded orthogroup, and six positive-selection loci
among 296 assessed. DE is restricted to PO and SS and its registered contrasts.
These evidence categories overlap and must not be summed as distinct genes.

The corrected PEPC tree has 81 tips: 49 protein representatives for 51 loci and
32 references. Subtype assignments comprise nine supported, 27 provisional and
15 unresolved loci. Maize Q84KR7 positions are **519, 780 and 890**. Historical
nucleotide references were corrected using annotated exon-aware CDS translations.
Earlier PEPC trees containing translated introns are superseded. The validated
51-locus count did not change. The main and sensitivity trees have identical
unrooted topology, not identical support values or proof of function.

Expression displays use the mean of three replicates per sample group followed
by row z-scores across group means. They use square 0.30-inch cells. Log2FC
heatmaps are not part of the accepted expression presentation; log2FC values
remain in DE evidence tables and the separate metabolite-concordance comparison.
Other-family presentation trees are provided as context, not new subtype or
functional validation; their older PEPC panel is omitted in favour of the
corrected PEPC release above.

SS transcriptomics and metabolomics were measured in **different plants**.
Their comparison is descriptive, not a paired gene–metabolite correlation.
Laboratory metabolite preprocessing metadata remain incompletely documented.

## Provenance and history

[Source selection](provenance/SOURCE_SELECTION.tsv) records the original relative
paths and hashes. [Supersession manifest](provenance/SUPERSESSION_MANIFEST.tsv)
records replacements for the old manuscript-facing files. Original HPC releases
remain unchanged. The following historical tags are preserved:

- `photosynthesis-family-evidence-v1.0.0`
- `supported-family-downstream-evidence-v1.0.0`

Historical tagged releases may contain superseded PEPC conclusions; use this
snapshot for current PEPC interpretation. No history rewriting is performed.
Publication checksums describe this curated layout and are regenerated because
redundant archives/checkpoints and superseded files are omitted.

Verify the publication from its root with `python3 scripts/verify_publication.py`.

## AI assistance

It has been partly written with the help of AI, and the scripts have been improved using AI.
