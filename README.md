# Sesuvium Photosynthetic Flexibility

Reproducible comparative genomics and transcriptomics analyses of photosynthetic flexibility in Sesuvium and related Caryophyllales.

This repository contains the validated PEPC analysis, including candidate classification, corrected translated phylogenies, Wang et al. (2023) Supplementary Table S7 sequences, maize-coordinate alignment windows, iTOL datasets, evidence integration, and PEPC expression heatmaps.

## Final PEPC decisions

- Arabidopsis references were restricted to canonical PPC1–PPC4.
- Wang S7 nucleotide records were translated before alignment.
- Maize PEPC `Q84KR7_MAIZE` was used as the coordinate reference.
- Alignment windows contain two amino acids before and after maize positions 519, 780, and 890.
- `SS__g8156` was translated from its DNA-like sequence.
- `SS__g24599` was excluded from the final PEPC tree after independent homology assessment.
- Identical protein copies were deduplicated for the tree but retained where appropriate in expression analyses.

The current expression figures use row-wise expression z-scores. TMM-normalized log2 fold-change heatmaps are a separate planned analysis.

## Frozen photosynthesis-family evidence release

The [18 October 2026 evidence release](releases/photosynthesis_family_evidence_v1_2026-10-18/README.md) integrates 21 of the 23 Wang Table S6 family categories across eight datasets. It contains 657 supported family–locus records, with model/panel caveats and unresolved records preserved separately; MEP and PPA remain NA. These are evidence counts, not exhaustive C4/CAM-functional gene counts.

- [Methods](releases/photosynthesis_family_evidence_v1_2026-10-18/docs/methods.md)
- [Results](releases/photosynthesis_family_evidence_v1_2026-10-18/docs/results.md)
- [Count and gene-list workbook](releases/photosynthesis_family_evidence_v1_2026-10-18/tables/Wang_family_counts_and_gene_lists.xlsx)
- [Reproducible reporting script](releases/photosynthesis_family_evidence_v1_2026-10-18/scripts/build_family_tables.py)

The earlier PEPC analysis remains unchanged. Use the frozen integrated ledger for the family evidence counts reported in this release.
