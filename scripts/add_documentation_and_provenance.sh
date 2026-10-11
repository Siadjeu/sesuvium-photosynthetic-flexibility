#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p methods results manifests

cat > README.md <<'EOF'
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
EOF

cat > methods/PEPC_analysis_methods.md <<'EOF'
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
EOF

cat > results/PEPC_analysis_results.md <<'EOF'
# Results: PEPC comparative analysis

The final corrected PEPC release passed checksum validation across the classification, translation, tree, evidence, and heatmap checkpoints.

The final translated PEPC tree contained 87 terminal sequences after sequence deduplication and exclusion of the independently unsupported `SS__g24599` candidate. The tree included the eight project assemblies, corrected Arabidopsis PPC references, Wang et al. (2023) Supplementary Table S7 sequences, and the maize PEPC coordinate reference.

The primary expression sets contained 13 PO PEPC genes and 9 SS PEPC genes. The SS set included the restored copies `SS__g2566.t1` and `SS__g4086.t1`, which were absent from an earlier incomplete expression matrix. The duplicate copies `SS__g29856.t1` and `SS__g7519.t1` were retained separately for expression interpretation where their evidence was inherited from their representative copies.

The detailed evidence table reports WGD/segmental retention, CAFE5, positive selection, and differential-expression fields by gene. `NA` indicates that the corresponding evidence was unavailable or not assessed; it does not imply gene absence.

Maize-centered coordinate windows were generated for positions 519, 780, and 890 using two residues on either side of each reference position. These data are supplied in tree order for direct use with iTOL or publication figures.
EOF

cat > manifests/final_release_manifest.tsv <<'EOF'
fileclassificationstatus
tables/01_pepc_candidate_classification.tsv candidate classificationcomplete
tables/01_all_species_PEPC_detailed_evidence.tsv detailed evidencecomplete
tables/03_tree_tip_metadata.tsv final tree metadatacomplete
figures/PEPC_corrected_translated.treefile final PEPC treecomplete
figures/03_iTOL_maize_coordinate_text.txt iTOL coordinate datasetcomplete
figures/PO_PEPC_Figure3B_like_final_primary_PEPC.pdf PO expression heatmapcomplete
figures/SS_PEPC_Figure3B_like_final_primary_PEPC.pdf SS expression heatmapcomplete
EOF

cat > manifests/software_versions.tsv <<'EOF'
softwarepurposeversion
MAFFTprotein alignment7.526
IQ-TREEphylogeny3.1.3
DIAMONDsimilarity searches2.2.5
HMMERhomology audit3.4
Pythonscript execution3.12
Biopythonsequence/tree handling1.88
EOF

chmod +x scripts/add_documentation_and_provenance.sh

git add \
  README.md \
  methods \
  results \
  manifests \
  scripts/add_documentation_and_provenance.sh

git diff --cached --check -- \
  '*.md' '*.sh'

if git diff --cached --quiet; then
  echo "No new changes to commit."
else
  git commit -m "Document PEPC methods results and provenance"
  git push origin main
fi

echo "Documentation and provenance step completed."
