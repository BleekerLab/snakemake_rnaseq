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
    """This function checks if the sample has paired end or single end reads and returns 1 or 2 names of the fastq files"""
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trim_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end
      2. Returns the correct input and output trimmed file names. 
    """
    if sample_is_single_end(wildcards.sample):
        inFile = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return "--in1 " + inFile[0] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz" 
    else:
        inFile = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz --out2 "  + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"

def get_star_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end.
      2. Returns the correct input file names for STAR mapping step.
    """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"     
    else:
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
BAM_FILES = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES)
MAPPING_REPORT = RESULT_DIR + "mapping_summary.csv"


rule all:
    input:
        MULTIQC,
        BAM_FILES, 
        MAPPING_REPORT,
        RESULT_DIR + "raw_counts.parsed.tsv",
        RESULT_DIR + "scaled_counts.tsv"
    message:
        "RNA-seq pipeline run complete!"
    shell:
        "cp config/config.yaml {RESULT_DIR};"
        "cp config/samples.tsv {RESULT_DIR};"
        "rm -r {WORKING_DIR}"

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
        limit_genome_generate_ram = config["star_index"]["limitGenomeGenerateRAM"],
        genome_sa = config["star_index"]["genomeSAindexNbases"],
        genome_chr_bin_n_bits = config["star_index"]["genomeChrBinNbits"]
    threads: 10
    resources: mem_mb=100000
    shell:
        "mkdir -p {params.genome_dir}; " # if directory not created STAR will ask for it
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.genome_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--limitGenomeGenerateRAM {params.limit_genome_generate_ram} "
        "--genomeSAindexNbases {params.genome_sa} "
        "--genomeChrBinNbits {params.genome_chr_bin_n_bits}"

#######################
# RNA-seq read trimming
#######################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = temp(WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz"),
        fq2  = temp(WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"),
        html = WORKING_DIR + "fastp/{sample}_fastp.html",
        json = WORKING_DIR + "fastp/{sample}_fastp.json"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        in_and_out_files =  get_trim_names,
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    resources: cpus=10
    shell:
        "touch {output.fq2};\
        fastp --thread {threads}  --html {output.html} --json {output.json} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        {params.in_and_out_files} \
        2>{log}"

rule multiqc:
    input:
        expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "

#########################
# RNA-Seq read alignement
#########################

rule map_to_genome_using_STAR:
    input:
        genome_index = rules.star_index.output,
        forward_read = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        reverse_read = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
    output:
        RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam",
        RESULT_DIR + "star/{sample}_Log.final.out"
    message:
        "mapping {wildcards.sample} reads to genome"
    params:
        sample_name           =  "{sample}",
        star_input_file_names =  get_star_names,
        prefix                =  RESULT_DIR + "star/{sample}_",
        maxmismatches         =  config["star"]["mismatches"],
        unmapped              =  config["star"]["unmapped"]   ,
        multimappers          =  config["star"]["multimappers"],
        matchNminoverLread    =  config["star"]["matchminoverlengthread"],
        outSamType            =  config["star"]["samtype"],
        outSAMattributes      =  config["star"]["samattributes"],
        intronmax             =  config["star"]["intronmax"],
        matesgap              =  config["star"]["matesgap"],
        genome_index          =  WORKING_DIR + "genome/"
    threads: 10
    resources: cpus=10
    shell:
        "STAR --genomeDir {params.genome_index} --readFilesIn {params.star_input_file_names} --readFilesCommand zcat --outFilterMultimapNmax {params.multimappers} \
        --outFilterMismatchNmax {params.maxmismatches} --alignMatesGapMax {params.matesgap} --alignIntronMax {params.intronmax}  \
        --outFilterMatchNminOverLread  {params.matchNminoverLread} --alignEndsType EndToEnd --runThreadN {threads}  --outReadsUnmapped {params.unmapped} \
        --outFileNamePrefix {params.prefix} --outSAMtype {params.outSamType}  --outSAMattributes {params.outSAMattributes}"


rule generate_mapping_summary:
    input:
        expand(RESULT_DIR + "star/{sample}_Log.final.out", sample = SAMPLES)
    output:
        RESULT_DIR + "mapping_summary.csv"
    message:
        "Concatenating STAR mapping report and generating .csv mapping summary."
    params:
        directory_with_mapping_reports = RESULT_DIR + "star/" 
    shell:
        "python scripts/generate_mapping_summary.py {params.directory_with_mapping_reports} {output}"

##################################
# Produce table of raw gene counts
##################################

if config["count_multimappers"]:
    rule create_counts_table:
        input:
            bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
            gff  = config["refs"]["gtf"]
        output:
            WORKING_DIR + "raw_counts.tsv"
        message: "Producing the table of raw counts (counting read multimappers)"
        threads: 10
        shell:
            "featureCounts -T {threads} -M -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"
else:
    rule create_counts_table:
        input:
            bams = expand(RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam", sample = SAMPLES),
            gff  = config["refs"]["gtf"]
        output:
            WORKING_DIR + "raw_counts.tsv"
        message: "Producing the table of raw counts (not counting multimappers)"
        threads: 10
        shell:
            "featureCounts -T {threads} -t exon -g transcript_id -F 'gtf' -a {input.gff} -o {output} {input.bams}"


rule parse_raw_counts:
    input:
        WORKING_DIR + "raw_counts.tsv"
    output:
        RESULT_DIR + "raw_counts.parsed.tsv"
    message: 
        "Parsing the raw counts file for scaling (removal of comment and header renaming)"
    params:
        star_result_dir_name = RESULT_DIR + "star/"
    shell:
        "python scripts/parse_raw_counts.py {input} {params.star_result_dir_name} {output}"

#########################################
# Produce table of normalised gene counts
#########################################

rule normalise_raw_counts:
    input:
        raw = RESULT_DIR + "raw_counts.parsed.tsv"
    output:
        norm = RESULT_DIR + "scaled_counts.tsv"
    message: 
        "Normalising raw counts the DESeq2 way"
    shell:
        "Rscript --vanilla scripts/deseq2_normalization.R {input.raw} {output.norm}" 
