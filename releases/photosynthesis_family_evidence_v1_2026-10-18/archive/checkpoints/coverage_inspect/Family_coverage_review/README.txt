Family coverage and definition evidence audit.
No searches, alignments or trees are run. Existing files are read only.
01 is file-discovery coverage, NOT evidence that all genes were searched.
02 records original byte hashes. Source snapshots redact absolute cluster paths;
snapshot checksums are in SHA256SUMS.txt and therefore differ when redaction occurs.
04 summarizes existing filter/status labels without treating PASS as final membership.
05 retains distinct identifiers even when protein sequences are identical.
06 gives filter context; complete selected scripts are also included.
07 lists additional paths for follow-up; file-name hints are not family assignments.
09 is a reference curation template; no family identity is inferred from names alone.

Review order:
1. PPDK primary vs sensitivity vs legacy sets: identify exactly why candidates were lost.
2. V_ATPASE reference subunits and length filters: evaluate subunit-specific requirements.
3. PPDK_RP references and rejected candidates: do not interpret empty output as zero.
4. TDT_DTC_DIC references and candidates: split only with family-specific evidence.
5. Remaining Wang categories: distinguish missing search from uncollected evidence.
Do not force counts to match Wang. Final locus mapping and functional C4/CAM annotation
remain separate tasks. This audit uses the focused-family directory; absence of files
here does not exclude useful annotations elsewhere in the project.
Collection limits: 8 MB per source, 200 MB total. Skipped sources are explicit.
