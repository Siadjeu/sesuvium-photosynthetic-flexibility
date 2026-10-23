Final consolidation of corrected PEPC job 20837873.
Python >=3.9 standard library only; no conda activation or package installation.
Place this package in PROJECT_ROOT and submit run.slurm with PROJECT_ROOT exported.
The source release is pinned in run.py; no newest-directory guessing.
Copies verified evidence to a new final_pepc_v1_JOB_TIMESTAMP directory beside
the source release. Generates methods, results, limitations, figure legend,
handover, species/status gene lists, summary tables and a compact iTOL ZIP.
Does not rerun MAFFT/IQ-TREE or alter subtype decisions, source files, family
ledgers, earlier freezes or GitHub. Existing unresolved categories remain.
Final reports refer specifically to the reviewed job 20837873; do not repoint
this package to another source without reviewing and updating the report text.
Output paths and the final bundle are printed in the job output log.
The main final SHA256SUMS covers every file, including the iTOL ZIP.
A separate checksum file covers the outer tar.gz bundle.
