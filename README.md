[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2393540.svg)](https://doi.org/10.5281/zenodo.2393540)

# RNA_seq_Snakemake
A snakemake pipeline for the analysis of RNA-seq data
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)  
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

# Description

## Aim
To align, count and normalize counts using single or paired-end Illumina RNA-seq data.

## Description
This pipeline analyses the raw RNA-seq data and produces two files containing the raw and normalized counts. 

1. The raw fastq files will be trimmed for adaptors and quality checked with `fastp`.  
2. The genome sequence FASTA file will be used for the mapping step of the trimmed reads using `STAR`. 
3. A GTF annotation file will be used to obtain the raw counts using `subread featureCounts`. 
4. The raw counts will be scaled using `DESeq2` to generate the scaled ("normalized") counts. 

## Prerequisites: what you should know before using this pipeline
- Some command of the Unix Shell to connect to a remote server where you will execute the pipeline (e.g. SURF Lisa Cluster). You can find a good tutorial from the Software Carpentry Foundation [here](https://swcarpentry.github.io/shell-novice/) and another one from Berlin Bioinformatics [here](http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/unix_for_beginners/README.html).
- Some command of the Unix Shell to transfer datasets to and from a remote server (to transfer sequencing files and retrieve the results/). The Berlin Bioinformatics Unix begginer guide available [here] should be sufficient for that (check the `wget` and `scp` commands).
- An understanding of the steps of a canonical RNA-Seq analysis (trimming, alignment, etc.). You can find some info [here](https://bitesizebio.com/13542/what-everyone-should-know-about-rna-seq/).

## Content of this GitHub repository
- `Snakefile`: a master file that contains the desired outputs and the rules to generate them from the input files.
- `config/samples.tsv`:  a file containing information about the names, the paths and the conditions of the samples used as input for the pipeline. **This file has to be adapted to your sample names before running the pipeline**.
- `config/config.yaml`: the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules.
- `config/refs/`: a folder containing
  - a genomic reference in fasta format. The `S_lycopersicum_chromosomes.4.00.chrom1.fa` is placed for testing purposes.
  - a GTF annotation file. The `ITAG4.0_gene_models.sub.gtf` for testing purposes.
- `fastq/`: a folder containing subsetted paired-end fastq files used to test locally the pipeline. Generated using [Seqtk](https://github.com/lh3/seqtk):
`seqtk sample -s100 {inputfile(can be gzipped)} 250000 > {output(always gunzipped)}`
This folder should contain the `fastq` of the paired-end RNA-seq data, you want to run.
- `envs/`: a folder containing the environments needed for the conda package manager. If run with the `--use-conda` command, Snakemake will install the necessary softwares and packages using the conda environment files.
- `Dockerfile`: a Docker file used to build the docker image that, once run using `docker run rnaseq:dockerfile Snakemake --cores N --use-conda` will trigger installation of the necessary softwares and run the Snakemake pipeline.


## Pipeline depencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [fastp](https://github.com/OpenGene/fastp)
* [STAR](https://github.com/alexdobin/STAR)   
* [subread](http://subread.sourceforge.net/)  
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).  


# Usage (local machine)

## Download or clone the Github repository
You will need a local copy of the `snakemake_rnaseq_to_counts` on your machine. You can either:
- use git in the shell: `git clone git@github.com:KoesGroup/Snakemake_hisat-DESeq.git`
- click on "Clone or download" and select `download`

## Installing the required softwares and packages 

### Option 1: conda
:round_pushpin: Option 1: using the conda package manager :one:  
First, you need to create an environment where core softwares such as `Snakemake` will be installed. Second, Snakemake itself will use conda to install the required softwares in each rule.
1. Install the [Miniconda3 distribution (Python 3.7 version)](https://docs.conda.io/en/latest/miniconda.html) for your OS (Windows, Linux or Mac OS X).  
2. Inside a Shell window (command line interface), create a virtual environment named `rnaseq` using the `environment.yaml` file with the following command: `conda env create --name rnaseq --file envs/environment.yaml`
3. Then, before you run the Snakemake pipeline, activate this virtual environment with `source activate rnaseq`.

The Snakefile will then take care of installing and loading the packages and softwares required by each step of the pipeline.

### Option 2: Docker 
:round_pushpin: Option 2: using a Docker container :two:  
1. Install Docker desktop for your operating system.
2. Open a Shell window and type: `docker pull mgalland/rnaseq` to retrieve a Docker container that includes the pipeline required softwares (Snakemake and conda and many others).
3. Run the pipeline on your system with:
` docker run --mount source=<data folder on host machine>,target=<folder inside container that will corresponds to data folder> mgalland/rnaseq Snakemake --use-conda --cores N`.

For instance, in a Shell window, go inside the `snakemake_rnaseq_to_counts/` directory and type: `docker run -it -v $PWD:/home/ mgalland/snakemake /bin/bash`.  
This will link your working directory (`snakemake_rnaseq_to_counts/`) to a directory called `/rnaseq` inside the container. Then, the ` -it` option will have you to enter the container where you can run Snakemake commands and retrieve your data folder.    
Finally, to export your results outside of the container, you can use the Docker `cp` command. See it there: https://docs.docker.com/engine/reference/commandline/cp/.

## Configuration 
:pencil2: :clipboard:  
You'll need to change a few things to accomodate this pipeline to your needs.
Make sure you have changed the parameters in the `config/config.yaml` file that specifies where to find the sample data file, the genomic and transcriptomic reference fasta files to use and the parameters for certains rules etc.  
This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order. It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Dry run
From the folder containing the `Snakefile`, use the command `snakemake --use-conda -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `Snakemake --use-conda` and provide the number of cores with `--cores 10` for ten cores for instance.  
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).

# Main outputs
- the fastp report files __\*.html__
- the unscaled RNA-Seq read counts: __counts.txt__

# Usage (HPC cluster)
singularity + docker image


# Directed Acyclic Graph of jobs
![dag](./dag.png)
