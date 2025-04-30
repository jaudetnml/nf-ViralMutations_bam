# CHANGELOG

## Release 1.X.X - 2025-XX-XX

### Changes
- Switched read trimming to fastp/fastplong. These softwares also provide the pre- and post-trimming QC reports.
  - Since all arguments are named (rather than positional), a single process is sufficient to add custom features to the trimming
  - Since fastp/fastplong already provide QC reports (better ones than FastQC), the FastQC process was removed as well.

## Release 1.0.1 - 2025-04-01

### Fixes
- Fixed the issue related to Fastqc running out of memory
- Fixed an issue where the minimum number of lines output from FreeBayes was too low and allowed empty files (header only) to be pushed to the next processes.
- Fixed an issue collecting a pre-indexed host genome didn't work due to too many output files in the output tuple (one was duplicated)
  - Fixed an issue where the subworkflow was doing an `ifEmpty` check on the GetIndex output and the single value passed was triggering bad behaviour

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