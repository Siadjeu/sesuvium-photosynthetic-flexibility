Provenance-based correction of step16.7 combined input.
Only records attributed to the confirmed SS candidate.cds.fa source are removed.
Every retained original identifier is looked up in its metadata-declared protein file.
The original sequence must reproduce the old record under its historical letter-only
normalization, or the script stops. Source checksums record the inputs used.
source_verified_records.fasta preserves source sequence characters.
corrected_protein_candidates.faa uppercases and removes at most one terminal stop.
Internal stops, invalid symbols and empty proteins are quarantined, never silently
stripped. The complete ledger records each action and sequence checksum.
Composition flags, including AT N-rich proteins, are retained for explicit review.
No identical protein sequences or transcript suffixes are collapsed. No CDS translation
is added. 287197 is a source record count, not an approved nonredundant locus count.
Unmatched CDS remains separate. Existing pipelines and inputs are not overwritten.
Do not automatically launch searches: review QC, quarantine, AT flags and source
selection before approval. Source provenance is not validation of biological annotation.
Full FASTAs stay in full_outputs on cluster. Their hashes are in 08; upload bundle
contains review tables and unmatched CDS but omits the large full_outputs directory.
