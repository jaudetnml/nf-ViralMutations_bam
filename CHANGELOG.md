# CHANGELOG

## Release 1.0.0 - 2025-01-28

### Major changes
- Switched to containers (using Singularity at the moment) instead of Conda
  - Had to split some processes (mostly separate alignments from samtools usage)
- Re-structured the Preflight and Alignment workflows to:
  - Allow de-hosting of MinION data
  - Remove target reference pre-indexing with BWA for MinION data

### Other changes
- Moved resource assignment from processes to configs/resources.config
- Removed two processes that were leftovers from using LoFreq3

## Release 0.1.0 - 2024-11-07

Initial commit
Add CHANGELOG