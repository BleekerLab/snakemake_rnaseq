#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
configfile: "config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

# fetch URL to transcriptome multi fasta from configfile
genome_url = config["refs"]["genome"]
transcriptome_fasta_url = config["refs"]["transcriptome_fasta"]
transcriptome_gtf_url= config["refs"]["transcriptome_gtf"]

########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions 
SAMPLES = units.index.get_level_values('sample').unique().tolist()
CONDITIONS = list(pd.read_table(config["units"])["condition"])
fwd        = dict(zip(list(pd.read_table(config["units"])["sample"]), list(pd.read_table(config["units"])["fq1"])))
rev        = dict(zip(list(pd.read_table(config["units"])["sample"]), list(pd.read_table(config["units"])["fq2"])))

def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_forward_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1"]].dropna()

def get_reverse_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq2"]].dropna()


#################
# Desired outputs
#################
rule all:
    input:
        FASTQC = expand(RESULT_DIR + "fastqc/{sample}.{step}.html", sample = SAMPLES,step=["original","trimmed"])
    message:
        "Job done! Removing temporary directory"
    shell:
        "rm -r {WORKING_DIR}"

#######
# Rules
#######

# download genome with the use of the URL
rule get_genome_fasta:
    output:
        "genome/genome.fasta"
    message:
        "downloading the required genomic fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_url}"

# download transcriptome fasta's with the use of the URL
rule get_transcriptome_fasta:
    output:
        "genome/ref_transcriptome.fasta"
    message:
        "downloading the required transcriptome fasta file"
    shell:
        "wget -O {output} {transcriptome_fasta_url}"

# download transcriptome gtf/gff's with the use of the URL
rule get_transcriptome_gtf:
    output:
        "genome/ref_transcriptome.gff"
    message:
        "downloading required transcriptome gtf file"
    shell:
        "wget -O {output} {transcriptome_gtf_url}"

# create transcriptome index, for blasting
rule get_ref_transcriptome_index:
    input:
        "genome/ref_transcriptome.fasta"
    output:
        ["genome/ref_transcriptome.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

# create quality reports of original reads
rule fastqc_before_trimming:
    input:
        fwd = get_forward_fastq,
        rev = get_reverse_fastq
    output:
        RESULT_DIR + "fastqc/{sample}.original.html"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastp.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "Quality check of trimmed {wildcards.sample} sample with fastp"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input.fwd} -I {input.rev} -h {output} 2>{log}"

# trim and quality filter of the reads
rule trimmomatic:
    input:
        get_fastq
    output:
        fw_reads        = "trimmed/{sample}_fw.fq",
        rev_reads       = "trimmed/{sample}_rev.fq",
        forwardUnpaired = "trimmed/{sample}_forward_unpaired.fastq",
        reverseUnpaired = "trimmed/{sample}_reverse_unpaired.fastq"
    message:
        "trimming reads"
    conda:
        "envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        seedMisMatches         =  str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =  str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold   =  str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual        =  str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual       =  str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize             =  str(config['trimmomatic']['windowSize']),
        avgMinQual             =  str(config['trimmomatic']['avgMinQual']),
        minReadLen             =  str(config['trimmomatic']['minReadLength']),
        phred                  =  str(config["trimmomatic"]["phred"]),
        adapters               = str(config["trimmomatic"]["adapters"])
    threads: 1
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input} "
        "{output.fw_reads} "
        "{output.forwardUnpaired} "
        "{output.rev_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

# create quality reports of trimmed reads
rule fastqc_after_trimming:
    input:
        fwd = "trimmed/{sample}_fw.fq",
        rev = "trimmed/{sample}_rev.fq"
    output:
        RESULT_DIR + "fastqc/{sample}.trimmed.html"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastp.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "Quality check of trimmed {wildcards.sample} sample with fastp"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input.fwd} -I {input.rev} -h {output} 2>{log}"

