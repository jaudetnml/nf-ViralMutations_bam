# nf-ViralMutations

## Pipeline intent

This pipeline is intended to be used to align Illumina or MinION reads to a reference and call the mutations and consensus.
The mutations called include non-consensus mutations, which is useful when looking at the dynamics of the viral population over time.

**This version of the pipeline is meant to start from 2 bam files that have been processed (deduplicated, filtered, primerclipped) and combine them.**
The first processing step is to merge the BAM files and get the depth, everything downstream follows the original [nf-ViralMutations](phac-nml/nf-ViralMutations).
  
The consensus will include IUPAC alternates if multiple alternate bases create a situation where neither an alternate nor the reference form consensus.

## Installation

The two installation requirements are [Nextflow](https://www.nextflow.io/docs/latest/install.html) and a container platform (the pipeline has profiles for [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) and [Apptainer](https://apptainer.org/)).

A profile is also set up to use `slurm`, remember to set the `SLURM_Queue` parameter.

## Running the pipeline

Assuming all the parameters are stored in a JSON or YAML file, navigate to the folder where you want the `nextflow` and `work` directories to be located and run the following command:

```shell
nextflow run PHAC-NML/nf-ViralMutations -r v1.1.0 -params-file /loc/of/Params_Exp1.yml -profile singularity,slurm -with-report -with-dag
```

## Parameters

The full set of parameters are described in the comments in the `nextflow.config` file.
Some less simple ones are explained below.
It is easier to assemble the pipeline parameters in a YAML or JSON file to be passed with the command line argument `-params-file`.

### snpEff configuration

The pipeline uses snpEff to infer the effect of mutations.
The snpEff database is "rebuilt" for every run.
In order for the pipeline to run you must pass a name for the pathogen and a folder containing the sequence and annotations.
The entry is added to a barebones snpEff.config file saved in the pipeline directory which passed along to the snpEff process.
The entry is built in the folder specified, make sure it is write-accessible.
The folder should have the same name as the name passed to snpEff and contain two files: 

- `sequence.fasta` contains the fasta sequence (sequences for segmented genomes)
- `genes.gbk` contains the annotations in GenBank format (concatenated entries for segmented genomes)

The easiest way to assemble this folder is to download the entries (selecting all of them for segmented genomes) from GenBank (or other repository) directly.
If the annotations were modified and exported in a Windows software (e.g. DNASTAR's SeqBuilder), make sure the files have LF (and not CRLF) line endings.

### Input reads

The `input` parameter takes a csv file that must have columns `sample`, `bam_1`, `bam_2`.
Absolute paths are generally preferred.

```CSV
sample,fastq_1,fastq_2
Samp1,/home/user/Data/Experiment3/Fastq/Samp1_S1_R1_001.fastq,/home/user/Data/Experiment3/Fastq/Samp1_S1_R2_001.fastq
Virus3,/home/user/Data/Experiment3/Fastq/Samp2_S5_R1_001.fastq,/home/user/Data/Experiment3/Fastq/Samp2_S5_R2_001.fastq
```

### Other parameters

- `outdir` defaults to `${launchDir}/Results`, but can be changed.
- `Seq_Tech` should be either `Illumina` or `MinION`.
- `Target_Reference` is the path to the viral genome you are aligning to. Since viral genomes are generally small, it is always re-indexed.
- `SLURM_Queue` to specify the SLURM queue or partition to be used. Only used with the `slurm` profile.
- `GenePos` (optional) an Excel file used to annotate the graph of SNPs with the following columns:
  - `CHR` The reference name.
  - `CDS_Name` The name to show on the graph. Assumes most annotation will be CDSs but doesn't check, so you can annotate any feature you want.
  - `Start` The 1-based position for the start of the CDS.
  - `Stop` The 1-based position for the end of the CDS.

## Outputs

Most intermediate files are saved into `outdir` to allow QC of the different steps and to diagnose unexpected results.

### Key Outputs

- `{SampleName}_variants.pdf` Graphical presentation of all the SNPs that pass the depth and minimum frequency. If there are too many, some of the labels will be missing. The grey line represents the depth at each position (left axis) and the red bars position the SNPs along the genome; the red dots show the frequency of the SNP (right axis) along with error bars (based on binomial distribution using frequency and depth).
- `{SampleName}_variants_annot_filtered.tsv`Tab-separated summary of the SNPs/indels along with snpEff annotations. Can be opened in Excel or imported in R for further processing.
- `{SampleName}_variants_annot.vcf` Cleaned-up FreeBayes output (not filtered) with snpEff annotations.
- `{SampleName}_variants.vcf` Raw FreeBayes output.
- `{SampleName}_variants_annot.html` snpEff summary of the effect of mutations.
- `Alignments/{SampleName}_Aligned_pe.bam` Initial alignment to the `Target_Reference`. All derivatives created during clean-up are in the `Alignments` sub-folder, the order they are created is: _dd (deduplicated; Illumina-only), _pc (primer-clipped), _noSplit (filtered), _ds (downsampled; if max-depth was specified).
- `{SampleName}_consensus.fasta` The consensus sequence.

## Example DAGs

Given the different combinations, the DAGs below are not exhaustive but they reflect 3 common scenarios.

### Illumina PCR sequencing

![DAG showing the processes involved in analysing Illumina tiling amplicon data](images/DAG_Illumina_PCR.png)

### Illumina Shotgun sequencing with pre-indexed host sequence

![DAG showing the processes involved in analysing Illumina shotgun data, including de-hosting with a pre-indexed host sequence](images/DAG_Illumina_Shotgun.png)

### MinION PCR sequencing

![DAG showing the processes involved in analysing MinION tiling amplicon data](images/DAG_MinION_PCR.png)

## Support
Should any issues arise when running this pipeline, please create an issue for the moment.

## Legal
Copyright Government of Canada 2024

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
