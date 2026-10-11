# Methods: PEPC comparative analysis

## Species and inputs

The analysis used eight project assemblies represented by AT, BV, MC, PO, SLP_Hap1, SP_Hex, SP_Tet, and SS. SLP_Hap1 corresponds to *Selenicereus polyrhizus*. Hap2 was used only for validation and was not treated as an additional species.

## PEPC candidate recovery

Primary proteomes were screened using PEPC-related annotation keywords and similarity searches against curated PEPC references. Candidate classification combined sequence length, similarity, coverage, and the validated PEPC registry. PEPC kinase was treated as a separate family.

## Reference curation

The Arabidopsis reference set was audited and restricted to canonical plant-type PEPC proteins PPC1–PPC4. PPCK, PEPK, CPPM, and ZHD proteins were excluded from the PEPC reference set.

## Wang S7 and maize sequences

PEPC sequences from Wang et al. (2023) Supplementary Table S7 were retained. Nucleotide records were translated in the forward reading frames, with the best-supported frame selected using similarity evidence. The maize PEPC sequence `Q84KR7_MAIZE` was retained as the coordinate anchor.

## Alignment and phylogeny

Protein sequences were aligned with MAFFT. Identical ungapped amino-acid sequences were collapsed where required for the expanded tree. Maximum-likelihood phylogenies were reconstructed with IQ-TREE using bootstrap and SH-aLRT support.

## Coordinate windows

For maize positions 519, 780, and 890, five-residue windows were extracted using exact residue-to-alignment-column mapping. The displayed windows correspond to maize positions 517–521, 778–782, and 888–892. Gaps are reported as alignment gaps, ambiguous residues as `X`, and unavailable regions as `NO_COV`.

## Evidence integration

PEPC tips were annotated with direct WGD/segmental-retention evidence, robust CAFE5 evidence, positive-selection evidence, and differential-expression evidence. Missing evidence was retained as unavailable rather than interpreted as negative evidence.

## Expression heatmaps

PO and SS expression matrices were filtered to the final primary PEPC sets. Heatmap values are row-wise expression z-scores calculated across samples:

`z = (x - mean(x)) / standard_deviation(x)`

Positive values indicate expression above the gene-specific sample mean, whereas negative values indicate expression below that mean. The current heatmaps are z-score figures; TMM-normalized log2 fold-change figures remain a separate planned analysis.
