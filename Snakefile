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
configfile: "config/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
#CONDITIONS = list(pd.read_table(config["units"])["condition"])
samplefile = config["units"]


###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return [WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"]


#################
# Desired outputs
#################
rule all:
    input:
        RESULT_DIR + "counts.txt"
    message:
        "Job done! Removing temporary directory"
    shell:
        "rm -r {WORKING_DIR}"

#######
# Rules
#######


##########
# trimming
##########

rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")

#########################
# RNA-Seq read alignement
#########################

rule index:
    input:
        config["refs"]["genome"]
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule hisat_mapping:
    input:
        get_trimmed,
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        sum   = RESULT_DIR + "logs/{sample}_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -U {input} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -1 {input[0]} -2 {input[1]} | samtools view -Sb -F 4 -o {output.bams}")



#########################################
# Get table containing the raw counts
#########################################

rule create_counts_table:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = config["refs"]["gff"]
    output:
        RESULT_DIR + "counts.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"

