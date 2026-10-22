# Sesuvium Photosynthetic Flexibility

Reproducible comparative genomics and transcriptomics analyses of photosynthetic flexibility in Sesuvium and related Caryophyllales.

This repository contains the validated PEPC analysis, including candidate classification, corrected translated phylogenies, Wang et al. (2023) Supplementary Table S7 sequences, maize-coordinate alignment windows, iTOL datasets, evidence integration, and PEPC expression heatmaps. It has been partly written with the help of AI, and the scripts have been improved using AI.

## Final PEPC decisions

- Arabidopsis references were restricted to canonical PPC1–PPC4.
- Wang S7 nucleotide records were translated before alignment.
- Maize PEPC `Q84KR7_MAIZE` was used as the coordinate reference.
- Alignment windows contain two amino acids before and after maize positions 519, 780, and 890.
- `SS__g8156` was translated from its DNA-like sequence.
- `SS__g24599` was excluded from the final PEPC tree after independent homology assessment.
- Identical protein copies were deduplicated for the tree but retained where appropriate in expression analyses.

The current expression figures use row-wise expression z-scores. TMM-normalized log2 fold-change heatmaps are a separate planned analysis.
