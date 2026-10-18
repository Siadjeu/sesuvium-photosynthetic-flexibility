Coordinates are original genomic 1-based inclusive, except BED uses 0-based half-open.
Only target-overlapping plus-strand hits are audited; flagged queries remain controls.
CIGAR events describe alignment-implied disruptions, not validated variants or pseudogenes.
Terminal stop codons outside the PAF aligned interval are excluded. Split codons have explicit genomic_positions; BED spans their full interval.
Check 04_parser_checks.tsv before interpreting exact event positions. Any discrepancy requires parser review.
Splice dinucleotides and contexts are descriptive; canonical motifs do not validate an exon.
Upstream FASTA is existing 5 kb context only; no new search or gene rescue is claimed.
No BAM/CRAM evidence is inspected and no gene count or annotation is changed.
