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
samplefile = config["units"]


def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1","fq2"]].dropna()




###########################
# Input functions for rules
###########################
def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_forward_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1"]].dropna()

def get_reverse_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq2"]].dropna()

def get_bams():
    bams = [b for bam in glob("mapped/*.bam")]
    return bams

#################
# Desired outputs
#################
rule all:
    input:
        FASTQC = expand(RESULT_DIR + "fastqc/{sample}.{step}.html", sample = SAMPLES,step=["original","trimmed"]),
        GTF    = WORKING_DIR + "genome/stringtie_transcriptome.gtf",
        COUNTS = WORKING_DIR + "results/counts.txt"
    message:
        "Job done! Removing temporary directory"

#######
# Rules
#######


#####################
# Download references
#####################
rule get_genome_fasta:
    output:
        WORKING_DIR + "genome/genome.fasta"
    message:
        "downloading the required genomic fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_url}"

rule get_transcriptome_fasta:
    output:
        WORKING_DIR + "genome/ref_transcriptome.fasta"
    message:
        "downloading the required transcriptome fasta file"
    shell:
        "wget -O {output} {transcriptome_fasta_url}"

rule get_transcriptome_gtf:
    output:
        WORKING_DIR + "genome/ref_transcriptome.gff"
    message:
        "downloading required transcriptome gtf file"
    shell:
        "wget -O {output} {transcriptome_gtf_url}"

# create transcriptome index, for blasting
rule get_ref_transcriptome_index:
    input:
        WORKING_DIR + "genome/ref_transcriptome.fasta"
    output:
        [WORKING_DIR + "genome/ref_transcriptome.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

##################################
# Fastq QC before / after trimming
##################################
# create quality reports of original reads
rule fastqc_before_trimming:
    input:
        fwd = get_forward_fastq,
        rev = get_reverse_fastq
    output:
        html = RESULT_DIR + "fastqc/{sample}.original.html",
        json = RESULT_DIR + "fastqc/{sample}.original.json"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastp.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "Quality check of trimmed {wildcards.sample} sample with fastp"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input.fwd} -I {input.rev} -h {output.html} -j {output.json} 2>{log}"

# trim and quality filter of the reads
rule trimmomatic:
    input:
        get_fastq
    output:
        fw_reads        = WORKING_DIR + "trimmed/{sample}_fw.fq",
        rev_reads       = WORKING_DIR + "trimmed/{sample}_rev.fq",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq"
    message:
        "trimming {wildcards.sample} reads"
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
        fwd = WORKING_DIR + "trimmed/{sample}_fw.fq",
        rev = WORKING_DIR + "trimmed/{sample}_rev.fq"
    output:
        html = RESULT_DIR + "fastqc/{sample}.trimmed.html",
        json = RESULT_DIR + "fastqc/{sample}.trimmed.json",
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastp.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "Quality check of trimmed {wildcards.sample} sample with fastp"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input.fwd} -I {input.rev} -h {output.html} -j {output.json} 2>{log}"


#########################
# RNA-Seq read alignement
#########################
rule index:
    input:
        WORKING_DIR + "genome/genome.fasta"
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        "genome/genome"
    conda:
        "envs/hisat2.yaml"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule hisat_mapping:
    input:
        fwd   = WORKING_DIR + "trimmed/{sample}_fw.fq",
        rev   = WORKING_DIR + "trimmed/{sample}_rev.fq",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq",
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bamp  = temp(WORKING_DIR + "mapped/{sample}.pairs.bam"),
        bamf  = temp(WORKING_DIR + "mapped/{sample}.fwd.bam"),
        bamr  = temp(WORKING_DIR + "mapped/{sample}.rev.bam"),
        bams  = WORKING_DIR + "mapped/{sample}.bam"
    params:
        indexName = WORKING_DIR + "genome/genome"
    message:
        "mapping reads to genome to bam files."
    conda:
        "envs/hisat2_mapping.yaml"
    threads: 10
    shell:
        """
        hisat2 -p {threads} --dta -x {params.indexName} -1 {input.fwd} -2 {input.rev} | samtools view -Sb -F 4 -o {output.bamp}
        hisat2 -p {threads} --dta -x {params.indexName} -U {input.forwardUnpaired} | samtools view -Sb -F 4 -o {output.bamf}
        hisat2 -p {threads} --dta -x {params.indexName} -U {input.reverseUnpaired} | samtools view -Sb -F 4 -o {output.bamr}
        samtools merge {output.bams} {output.bamp} {output.bamf} {output.bamr}
        """

###########################################
# Create a de novo transcriptome annotation
###########################################
rule merge_bams:
    input:
        expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES)
    output:
        merged = temp(WORKING_DIR + "merged.bam"),
        bam_sorted = WORKING_DIR + "merged_sorted.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools merge {output.merged} {input}
        samtools sort {output.merged} -o {output.bam_sorted}
        """

rule create_stringtie_transcriptome:
    input:
        bam = WORKING_DIR + "merged_sorted.bam",
        Rtc = WORKING_DIR + "genome/ref_transcriptome.gff"
    output:
        WORKING_DIR + "genome/stringtie_transcriptome.gtf"
    #params:
        # some parameters
    conda:
        "envs/Stringtie.yaml"
    message:
        "creating transcriptome to stringtie_transcriptome.gtf."
    shell:
        "stringtie -G {input.Rtc} -o {output} {input.bam}"

rule gtf_to_fasta:
    input:
        Ntx  = WORKING_DIR + "genome/stringtie_transcriptome.gtf",
        gen  = WORKING_DIR + "genome/genome.fasta"
    output:
        WORKING_DIR + "genome/stringtie_transcriptome.fasta"
    conda:
        "envs/Tophat.yaml"
    shell:
        "gtf_to_fasta {input.Ntx} {input.gen} {output}"

rule blast_for_funtions:
    input:
        newTct = WORKING_DIR + "genome/stringtie_transcriptome.fasta",
        refTct = WORKING_DIR + "genome/ref_transcriptome.fasta",
        indexFiles = [WORKING_DIR + "genome/ref_transcriptome.fasta." + i for i in ("psq", "phr", "pin")]
    output:
        WORKING_DIR + "results/stringtie_transcriptome_blast.txt"
    params:
        evalue = "10",
        threads= "5"
    conda:
        "envs/blast.yaml"
    shell:
        "blastx -query {input.newTct} -db {input.refTct} -outfmt \"6 qseqid qlen slen evalue salltitles\" -out {output} -max_target_seqs 1 -num_threads {params.threads}"


rule create_counts_table:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/stringtie_transcriptome.gtf"
    output:
        WORKING_DIR + "results/counts.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t transcript -g gene_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"

rule DESeq2_analysis:
    input:
        counts      = WORKING_DIR + "results/counts.txt",
        samplefile  = samplefile
    output:
        WORKING_DIR + "results/result.csv"
    message:
        "normalizing read counts en creating differential expression table"
    params:
        maxfraction = float(config["DESeq2"]["maxfraction"])
    conda:
        "envs/deseq.yaml"
    shell:
        "Rscript scripts/DESeq2.R -c {input.counts} -s {input.samplefile} -o {output} -m {params.maxfraction}"
        
rule results_function:
    input:
        fa    = WORKING_DIR + "genome/stringtie_transcriptome.fasta",
        blast = WORKING_DIR + "results/stringtie_transcriptome_blast.txt",
        deseq = WORKING_DIR + "results/result.csv"
    output:
        final = "results/final.txt"
    shell:
        "python scripts/DE_with_Function.py {input.fa} {input.blast} {input.deseq} {output.final}"
