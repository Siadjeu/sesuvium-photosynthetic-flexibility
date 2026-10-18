Primary-proteome integrity audit. Existing files are never modified.
01 reports every record, including ordinal indices to preserve duplicate IDs.
Nucleotide alphabet/composition flags are screening criteria, not universal proof
of DNA: amino-acid sequences can use the same letters. No records are deleted.
04 tests frame +1, standard genetic code, allowing N -> X and one terminal stop.
Exact translation matches support nucleotide/protein duplication; they do not
prove a unique locus mapping when multiple proteins share a sequence.
Same-stem matches are hints only. No automatic suffix stripping for final counts.
No match does not establish that the sequence is a protein or unusable CDS.
06 traces exact metadata-cell ID matches. Missing matches are explicitly in 07.
Producer scripts and input inventory are included where present. Source sequence
files are not automatically opened or copied based on metadata paths.
Absolute cluster paths in context/metadata are redacted; originals stay unchanged.
No corrected proteome is released until source provenance and locus mappings
are reviewed. Do not resume family-count finalization from composition alone.
All audit tables are included; only illustrative FASTA records are bundled.
