# RNA_seq_Snakemake
A snakemake pipeline for the analysis of RNA-seq data with the use of hisat2 and DESeq2

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

# Aim
Map, count, normalize and get differential expressions of paired-end Illumina RNA-seq data.

# Description
Pipeline that analyses raw RNA-seq data to a file containing normalized counts, differential expression and functions of transcripts. The raw fastq files will be trimmed for adaptors and quality with trimmomatic. next, necessery genome and transcriptome will be downloaded and the trimmed reads will be mapped using hisat2. with stringtie and a refference annotation a new annotation will be created. This new annotation will be used to obtain the raw counts and do a local blast to a transcriptome fasta containing predicted functions. The counts are normalized and differential expressions are calculated using DESeq2. This data is combined with the predicted functions to get the final results table.


# Content
- Snakefile containing the targeted output and the rules to generate them from the input files.
- config/ , folder containing the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules.
- data/, folder containing samples.txt (sample descriptions) and subsetted paired-end fastq files used to test locally the pipeline. Generated using [Seqtk](https://github.com/lh3/seqtk):
`seqtk sample -s100 {inputfile(can be gzipped)} 250000 > {output(always gunzipped)}`
- envs/, folder containing the environments needed for the Snakefile to run. Need to make one specifically for MACS2 as MACS2 uses python 2.7 following the information found [here](https://groups.google.com/forum/#!searchin/snakemake/macs%7Csort:relevance/snakemake/60txGSq81zE/NzCUTdJ_AQAJ).


# Usage

## Conda environment
First, you need to install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html).  
1. Create a virtual environment named "chipseq" from the `environment.yaml` file with the following command:
`conda env create --name RNA-Seq --file ~/envs/DCM.yaml`
2. Then, activate this virtual environment with: `source activate RNA-Seq`    
Now, all the basic softwares and packages versions in use are the one listed in the `DCM.yaml` file.
The other environments (hisat2, subRead etc) will automatically be created and activated when requested by a rule.

## Configuration file
The `~/configs/config.yaml` file specifies the sample data file, the genomic and transcriptomic reference fasta files to use, the parameters for the rules to use, etc. This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order.
It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Dry run
Use the command `snakemake --use-conda -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `Snakemake --use-conda` and provide the number of cores with `--cores 10` for ten cores for instance.  
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).

# Main outputs
the alignment files (*.bam), fastqc report files, raw countsfile (counts.txt), results file containing normalized differential expressions (results.tsv)
