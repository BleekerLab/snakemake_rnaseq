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
FASTQC = expand(RESULT_DIR + "fastp/{sample}.html", sample = SAMPLES)
BAM_FILES = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES)


rule all:
    input:
        FASTQC,
        BAM_FILES, 
        RESULT_DIR + "raw_counts.tsv",
        RESULT_DIR + "scaled_counts.tsv"
    message:
        "Pipeline run complete!"
    shell:
        "cp config/config.yaml {RESULT_DIR};"
        "cp config/samples.tsv {RESULT_DIR}"

#######
# Rules
#######


###########################
# Genome reference indexing
###########################

rule star_index:
    input:
        fasta = config["refs"]["genome"],
        gtf =   config["refs"]["gtf"]
    output:
         genome_index = [WORKING_DIR + "genome/" + f for f in ["chrLength.txt","chrNameLength.txt","chrName.txt","chrStart.txt","Genome","genomeParameters.txt","SA","SAindex"]]
    message:
        "generating STAR genome index"
    params:
        genome_dir = WORKING_DIR + "genome/",
        sjdb_overhang = config["star_index"]["sjdbOverhang"],
        limit_genome_generate_ram = config["star_index"]["limitGenomeGenerateRAM"]
    threads: 10
    shell:
        "mkdir -p {params.genome_dir}; " # if directory not created STAR will ask for it
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.genome_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--limitGenomeGenerateRAM {params.limit_genome_generate_ram}"

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
        RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam",
        RESULT_DIR +      "star/{sample}_Log.final.out"
    message:
        "mapping {wildcards.sample} reads to genome"
    params:
        sample_name = "{sample}",
        prefix =             RESULT_DIR + "star/{sample}_",
        maxmismatches =      config["star"]["mismatches"],
        unmapped =           config["star"]["unmapped"]   ,
        multimappers =       config["star"]["multimappers"],
        matchNminoverLread = config["star"]["matchminoverlengthread"],
        outSamType =         config["star"]["samtype"],
        outSAMattributes =   config["star"]["samattributes"],
        intronmax =          config["star"]["intronmax"],
        matesgap =           config["star"]["matesgap"],
        genome_index =       WORKING_DIR + "genome/"
    threads: 10
    run:
        if sample_is_single_end(params.sample_name):
            shell("STAR --genomeDir {params.genome_index} --readFilesIn {input.forward} --readFilesCommand zcat --outFilterMultimapNmax {params.multimappers} \
            --outFilterMismatchNmax {params.maxmismatches} --alignMatesGapMax {params.matesgap} --alignIntronMax {params.intronmax}  \
            --outFilterMatchNminOverLread \ {params.matchNminoverLread} --alignEndsType EndToEnd --runThreadN {threads}  --outReadsUnmapped {params.unmapped} \
            --outFileNamePrefix {params.prefix} --outSAMtype {params.outSamType}  --outSAMattributes {params.outSAMattributes}")  
        else:
            shell("STAR --genomeDir {params.genome_index} --readFilesIn {input.forward} {input.reverse} --readFilesCommand zcat --outFilterMultimapNmax {params.multimappers} --outFilterMismatchNmax {params.maxmismatches} --alignMatesGapMax {params.matesgap} --alignIntronMax {params.intronmax} \
--outFilterMatchNminOverLread {params.matchNminoverLread} --alignEndsType EndToEnd  --runThreadN {threads}  --outReadsUnmapped {params.unmapped} --outFileNamePrefix {params.prefix}  --outSAMtype {params.outSamType} --outSAMattributes {params.outSAMattributes}")         


##################################
# Produce table of raw gene counts
##################################

rule create_counts_table:
    input:
        bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        gff  = config["refs"]["gtf"]
    output:
        RESULT_DIR + "raw_counts.tsv"
    message: "Producing the table of raw counts"
    threads: 10
    shell:
        "featureCounts -T {threads} -O -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"


rule parse_raw_counts:
    input:
        RESULT_DIR + "raw_counts.tsv"
    output:
        WORKING_DIR + "raw_counts.parsed.tsv"
    message: 
        "Parsing the raw counts file for scaling (removal of comment and header renaming)"
    params:
        tmp_file =             WORKING_DIR + "tmp.txt",
        star_result_dir_name = RESULT_DIR + "star/"
    shell:
        "tail -n +2 {input} | "
        "sed 's/_Aligned.sortedByCoord.out.bam//g' | "               # remove the file extension name in header
        "sed 's#{params.star_result_dir_name}##g'  > {output} "      # remove the file path in header


#########################################
# Produce table of normalised gene counts
#########################################

rule normalise_raw_counts:
    input:
        raw = WORKING_DIR + "raw_counts.parsed.tsv"
    output:
        norm = RESULT_DIR + "scaled_counts.tsv"
    message: 
        "Normalising raw counts the DESeq2 way"
    shell:
        "Rscript --vanilla scripts/deseq2_normalization.R {input.raw} {output.norm}" 
