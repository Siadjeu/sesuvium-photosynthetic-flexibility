Wang-compatible reconciliation, DOI 10.1093/plphys/kiad451
This is an evidence reconciliation release, NOT a final biological count release.
No analyses rerun; project inputs are read only. All new files are in this release.
01 contains all 15 species x 23 categories transcribed from the supplied Table S6.
03 preserves historical values without endorsing them.
04 has source-level evidence; 05 is the locus/count review queue.
06 counts unique exact candidate identifiers, not proven distinct loci.
Beta counts include supported and partial beta screen candidates; alpha/no-support rows
are retained in 05 but omitted from 06. Old mixed-reference reviewer acceptances are
not promoted. Tree exclusions are not automatic family exclusions.
07 is the required per-species/family completeness review.
08 intentionally contains NA until membership, locus identity and coverage are approved.
09 compares PO and BV only; AT is Amaranthus tricolor, not Arabidopsis.
SLP_Hap1 is not Wang's H. undatus. Do not substitute species or haplotypes.
10 compares PEPC labels after explicit evm.model./evm.TU. normalization; verify sequences.
12 inventories other possible sources without silently selecting newest or merging releases.
source_checkpoints contains small complete selected input files for review.
Never treat empty output as biological zero. Preserve separate gene loci even if
proteins are identical; collapse isoforms only using a documented locus mapping.
AT two-domain beta-CA protein and alternative isoform count as one locus once mapped.

Next: review 11 and 07, inspect source schemas and additional releases in 12.
In a COPY of 05, supply locus_id, mapping evidence and count_decision:
ACCEPT_FAMILY / EXCLUDE / UNRESOLVED, with reviewer_reason.
For pooled transporter records, assign exactly one supported Wang category or UNRESOLVED.
The script does not consume edited templates on rerun: retain edited copies separately.
A subsequent finalization must validate decisions, deduplicate (species,category,locus),
and require approved completeness for each reported cell. Only an explicitly verified
empty family can be reported as zero. C4/CAM functional roles remain separate annotations.
