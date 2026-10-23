SS transcript-metabolite evidence consolidation v1.0.0

This package provides manuscript draft methods, results, limitations and
figure legends, plus a standard-library-only HPC consolidation script.
Manually activate your preferred environment. No installation or activation
is included. Extract under PROJECT_ROOT and submit freeze.slurm.

The script verifies and copies these exact completed source releases:
08_focused_families/downstream_supported/metabolomics/release_20681419_20261021T134025Z
08_focused_families/downstream_supported/metabolomics/concordance_20732199_20261021T155959Z
No newest-release guessing, statistical rerun or source overwrite occurs.
A new timestamped release, checksum manifest, archive and archive checksum
are written under metabolomics/frozen_releases. Provenance paths are relative.

Read docs/METHODS.txt, RESULTS.txt, LIMITATIONS.txt and FIGURE_LEGENDS.txt.
The final concordance v3 figures are under evidence/transcript_metabolite_comparison/figures.
Supporting combined figures are under evidence/metabolite_comparisons/figures.
Use NADP-ME and NAD-ME as reaction-associated examples; retain other families
as supplementary/context displays. Each source retains inputs, script, QC,
tables and its original checksum manifest. No figures are changed here.

The primary-DE rule was verified in the prior integration script:
adjusted p < 0.05 AND absolute apeglm-shrunken log2FC >= 1.
There are 42 primary-DE genes: 21 up and 21 down in leaves; 83 assessable
of 99 supported. All genes remain in the tables, 83 in the final plots.

This freeze is archival, not an operating-system immutable lock. Any future
revision must receive a new version and manifest. No GitHub publication is
performed. Laboratory acquisition, normalization and detection-limit details
remain open and must be added before final manuscript submission.
