# define configfile
configfile: "configs/configs.yaml"
# fetch URL to transcriptome multi fasta from configfile
genome_URL = config["refs"]["genome"]
transcriptome_fasta_URL = config["refs"]["transcriptomeFas"]
transcriptome_gtf_URL = config["refs"]["transcriptomeGtf"]
# create lists containing samplenames and conditions from the file: data/sampls.txt
import pandas as pd
SAMPLES = list(pd.read_table(config["units"])["sample"])
conditions = list(pd.read_table("data/samples.txt")["condition"])

rule all:
    input:
        countsTable = "results/counts.txt",
        resultsFile = "results/result.csv",
        final       = "results/final.txt",
        fwd         = expand("results/fastqc/{sample}_fw_fastqc.zip", sample = SAMPLES),
        rev         = expand("results/fastqc/{sample}_rev_fastqc.zip", sample = SAMPLES),
        fas         = "genome/stringtie_transcriptome.fasta",
        funcsList   = "results/stringtie_transcriptome_blast.txt"
    message:
        "Job done!"
# download genome with the use of the URL
rule get_genome_fasta:
    output:
        "genome/genome.fasta"
    message:"downloading required genomic fasta file"
    shell: "wget -O {output} {genome_URL}"

# download transcriptome fasta's with the use of the URL
rule get_transcriptome_fasta:
    output:
        "genome/ref_transcriptome.fasta"
    message:"downloading required transcriptome fasta file"
    shell: "wget -O {output} {transcriptome_fasta_URL}"

# download transcriptome gtf/gff's with the use of the URL
rule get_transcriptome_gtf:
    output:
        "genome/ref_transcriptome.gff"
    message:"downloading required transcriptome gtf file"
    shell: "wget -O {output} {transcriptome_gtf_URL}"

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

# trim and quality filter of the reads
rule trimmomatic:
    input:
        fq1 = "data/{SAMPLES}_R1.sub.fq",
        fq2 = "data/{SAMPLES}_R2.sub.fq",
        adapters = config["adapters"]
    output:
        fw_reads = "trimmed/{SAMPLES}_fw.fq",
        rev_reads = "trimmed/{SAMPLES}_rev.fq",
        forwardUnpaired = "trimmed/{SAMPLES}_forward_unpaired.fastq",
        reverseUnpaired = "trimmed/{SAMPLES}_reverse_unpaired.fastq"
#    message: "trimming reads"
#        "logs/trimmomatic/{SAMPLES}.log"
    params:
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 1
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.fq1} "
        "{input.fq2} "
        "{output.fw_reads} "
        "{output.forwardUnpaired} "
        "{output.rev_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen}" #" 2>{log}"

# create quality reports of trimmed reads
rule fastqc:
    input:
        fwd = "trimmed/{SAMPLES}_fw.fq",
        rev = "trimmed/{SAMPLES}_rev.fq",
    output:
        fwd="results/fastqc/{SAMPLES}_fw_fastqc.zip",
        rev="results/fastqc/{SAMPLES}_rev_fastqc.zip"
    log:
        "results/logs/fastqc/{SAMPLES}.fastqc.log"
    params:
        "results/fastqc/"
    message:
        "Quality check of trimmed samples with FASTQC" 		#removed, it was not working
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev}"


rule index:
    input:
        "genome/genome.fasta"
    output:
        ["genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:"indexing genome"
    params:
        "genome/genome"
    conda:
        "envs/Hisat.yaml"
    threads: 10
    shell:"hisat2-build -p {threads} {input} {params}"

rule hisat_mapping:
    input:
        fwd   = "trimmed/{SAMPLES}_fw.fq",
        rev   = "trimmed/{SAMPLES}_rev.fq",
        forwardUnpaired = "trimmed/{SAMPLES}_forward_unpaired.fastq",
        reverseUnpaired = "trimmed/{SAMPLES}_reverse_unpaired.fastq",
        indexFiles = ["genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bamp  = temp("mapped/{SAMPLES}pairs.bam"),
        bamf  = temp("mapped/{SAMPLES}fwd.bam"),
        bamr  = temp("mapped/{SAMPLES}rev.bam"),
        bams  = "mapped/{SAMPLES}.bam"
    params:
        indexName = "genome/genome"
    message:
        "mapping reads to genome to bam files."
    conda:
        "envs/hisat2.yaml"
    threads: 10
    shell:
        """
        hisat2 -p {threads} --dta -x {params.indexName} -1 {input.fwd} -2 {input.rev} | samtools view -Sb -F 4 -o {output.bamp}
        hisat2 -p {threads} --dta -x {params.indexName} -U {input.forwardUnpaired} | samtools view -Sb -F 4 -o {output.bamf}
        hisat2 -p {threads} --dta -x {params.indexName} -U {input.reverseUnpaired} | samtools view -Sb -F 4 -o {output.bamr}
        samtools merge {output.bams} {output.bamp} {output.bamf} {output.bamr}
        """

rule merge_bams:
    input:
        expand("mapped/{sample}.bam", sample=SAMPLES),
    output:
        merged = temp("mapped/merged.bam"),
        sorted = "mapped/merged_sorted.bam"
    conda:
        "envs/Samtools.yaml"
    shell:
        """
        samtools merge {output.merged} {input}
        samtools sort {output.merged} -o {output.sorted}
        """

rule create_stringtie_transcriptome:
    input:
        bam = "mapped/merged_sorted.bam",
        Rtc = "genome/ref_transcriptome.gff"
    output:
        "genome/stringtie_transcriptome.gtf"
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
        Ntx  = "genome/stringtie_transcriptome.gtf",
        gen  = "genome/genome.fasta"
    output:
        "genome/stringtie_transcriptome.fasta"
    conda:
        "envs/Tophat.yaml"
    shell:
        "gtf_to_fasta {input.Ntx} {input.gen} {output}"

rule blast_for_funtions:
    input:
        newTct = "genome/stringtie_transcriptome.fasta",
        refTct = "genome/ref_transcriptome.fasta",
        indexFiles = ["genome/ref_transcriptome.fasta." + i for i in ("psq", "phr", "pin")]
    output:
        "results/stringtie_transcriptome_blast.txt"
    params:
        evalue = "10"

    conda:
        "envs/blast.yaml"
    shell:
        "blastx -query {input.newTct} -db {input.refTct} -outfmt \"6 qseqid qlen slen evalue salltitles\" -out {output} -max_target_seqs 1"


rule create_counts_table:
    input:
        bams = expand("mapped/{sample}.bam", sample=SAMPLES),
        Ntx  = "genome/stringtie_transcriptome.gtf"
    output:
        "results/counts.txt"
    #params:
        # some parameters
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t transcript -g gene_id -F 'gff' -a {input.Ntx} -o {output} {input.bams}"

rule DESeq2_analysis:
    input:
        counts = "results/counts.txt",
        #functions = "results/stringtie_transcriptome_blast.txt"
    output:
        "results/result.csv"
    message:
        "normalizing read counts en creating differential expression table"
    params:
        con = expand("{cons}", cons=conditions)
    conda:
        "envs/deseq.yaml"
    shell:
        "Rscript scripts/DESeq2.R {input.counts} {params.con}"

rule results_function:
    input:
        fa    = "genome/ref_transcriptome.fasta",
        blast = "results/stringtie_transcriptome_blast.txt",
        deseq = "results/result.csv"
    output:
        final = "results/final.txt"
    shell:
        "python DE_with_Function.py {input.fa} {input.blast} {input.deseq} {output.final}"
