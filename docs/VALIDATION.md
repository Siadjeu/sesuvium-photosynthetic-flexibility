# Publication-preparation validation

The preparation workflow was exercised against copies of the reviewed source
files and a local Git test repository. Test-only commit IDs were used for the
fixture; the distributed script remains pinned to the audited HPC main commit.

- 1,734 selected source files were copied only after exact SHA256 verification.
- The staged Git blobs matched the assembled publication bytes.
- The original test checkout remained clean; main and historical refs were unchanged.
- The prepared tree passed the configured private-HPC-path scan, including
  decodable Flate-compressed PDF streams. This is not an audit of historical
  commits or a complete credential/OCR scan.
- Python scripts passed syntax parsing.
- The family table builder completed with 3,193 ledger records, 657 supported
  records, 184 family/species cells and MEP/PPA retained as blocked definitions.
- The portable expression wrapper generated 36 panels; all 72 group-mean and
  z-score tables matched the accepted numerical values within 1e-12 tolerance.
- Both SS reproduction commands completed. Nine metabolite tables and five
  descriptive concordance tables matched the frozen values within numerical
  tolerances (relative 1e-9, absolute 1e-10).

Figure binary files can differ between software versions; these checks establish
numeric reproduction and successful figure generation, not binary-identical plots.
The HPC preparation step itself performs no plotting or biological reanalysis.
