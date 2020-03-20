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
        "Job done! Removing temporary directory and copying config files."
    shell:
        "rm -r {WORKING_DIR} "
        "cp config/config.yaml {RESULT_DIR} "
        "cp config/samples.tsv {RESULT_DIR} "

#######
# Rules
#######


###########################
# Genome reference indexing
###########################

rule star_index:
    input:
        fasta = config["refs"]["genome"],
        gff =   config["refs"]["gff"]
    output:
         genome_index = [WORKING_DIR + "genome/" + f for f in ["chrLength.txt","chrNameLength.txt","chrName.txt","chrStart.txt","Genome","genomeParameters.txt","SA","SAindex"]]
    message:
        "generating STAR genome index"
    params:
        genome_dir = WORKING_DIR + "genome/",
        sjdb_overhang = config["star_index"]["sjdbOverhang"]
    threads: 10
    shell:
        "mkdir -p {params.genome_dir}; " # if directory not created STAR will ask for it
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.genome_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gff} "
        "--sjdbOverhang {params.sjdb_overhang}"

#######################
# RNA-seq read trimming
#######################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html",
        json = RESULT_DIR + "fastp/{sample}.json"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")



#########################
# RNA-Seq read alignement
#########################

rule map_to_genome_using_STAR:
    input:
        genome_index = rules.star_index.output,
        forward = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        reverse = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
    output:
        temp(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam"),
        RESULT_DIR +      "star/{sample}_Log.final.out",
    message:
        "mapping {wildcards.sample} reads to genome"
    params:
        prefix = RESULT_DIR + "star/{sample}_",
        maxmismatches = config["star"]["mismatches"],
        unmapped = config["star"]["unmapped"]   ,
        multimappers = config["star"]["multimappers"],
        matchNminoverLread = config["star"]["matchminoverlengthread"],
        outSamType = config["star"]["samtype"],
        outSAMattributes = config["star"]["samattributes"],
        intronmax = config["star"]["intronmax"],
        matesgap =  config["star"]["matesgap"],
        genome_index = WORKING_DIR + "genome/"
    threads: 10
    shell:
        "STAR --genomeDir {params.genome_index} "
        "--readFilesIn {input.forward} {input.reverse} "
        "--readFilesCommand zcat "
        "--outFilterMultimapNmax {params.multimappers} "
        "--outFilterMismatchNmax {params.maxmismatches} "
        "--alignMatesGapMax {params.matesgap} "
        "--alignIntronMax {params.intronmax} "
        "--outFilterMatchNminOverLread {params.matchNminoverLread} "
        "--alignEndsType EndToEnd "
        "--runThreadN {threads} "
        "--outReadsUnmapped {params.unmapped} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype {params.outSamType} "
        "--outSAMattributes {params.outSAMattributes} "


#########################################
# Get table containing the raw counts
#########################################

rule create_counts_table:
    input:
        bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        gff  = config["refs"]["gff"]
    output:
        RESULT_DIR + "counts.txt"
    threads: 10
    shell:
        "featureCounts -T {threads} -O -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"

