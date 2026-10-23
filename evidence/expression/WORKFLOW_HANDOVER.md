# Mean-expression presentation handover
This release supersedes log2FC heatmaps for presentation. It exports the unchanged
21 gene-labelled tree packages and correspondence from the preceding presentation
release. The log2FC directory and its plotting scripts are omitted. The prior
release remains intact for provenance; this is not an in-place deletion.

For each gene, arithmetic means of DESeq2-normalized counts are calculated across
exactly three replicates of each sample group. PO groups are tissue x condition x
time x period (12 groups); SS groups are tissue x period (3 groups). No statistical
contrasts are averaged and no log transformation is applied. Row z-scores are
calculated across those group means: (mean_group - mean_across_groups) / population
SD_across_groups (ddof=0). Negative values mean below that gene's mean, not DE.
Zero-variance genes are grey; no missing values are imputed. Colour limits are
+/-2.5, while the exported numerical z-scores are uncapped. Group membership,
mean counts and z-scores are supplied. All 127 PO and 99 SS supported genes remain.

The old adaptive plotting layout did not enforce physically square cells.
This revision uses PO-PEPC as the common style reference and explicitly sets all
cells to 0.30 x 0.30 inches across all 36 panels, including PO-PEPC. Figure size
therefore varies with gene and group count. Do not resize exported panels by
different factors when assembling a composite figure. SVG/PDF preserve dimensions.
The established tissue/condition/time/period colour palette is retained.

Family validation remains 657 supported loci across eight species and 18 families.
AT means Amaranthus tricolor; Arabidopsis reference tips are separate. Five other
Wang framework families are not represented as supported sets: MEP/PPA definition
blocked, PEPCK/MPC/PPT zero supported (not proven absence). Tree memberships,
branches and assignments are unchanged. ALMT reference Q3E9Z9 is gap-rich; retain
that limitation. Neither normalized expression nor tree placement establishes
C4/CAM-specific function. The original family-validation GitHub freeze is preserved.
No publishing is performed. Next: inspect these figures, export labelled iTOL
trees, then prepare a separate publication freeze. Carry this bundle and the
prior handover into a new chat; automatic chat-capacity alerts are unavailable.

## Presentation revision v4
Only displayed formatting changed: solid black borders on all heatmap cells,
annotation cells and legend boxes; aligned legend rows headed Tissue, Condition,
Time and Period (applicable categories only for SS). Day and night appear in
lowercase in legends and sample labels; well-watered uses a hyphen. Original
machine-readable sample/group IDs remain unchanged for provenance. Group means,
z-scores, gene order, colour scale and 0.30-inch cell sizes are unchanged.
