#!/usr/bin/env Rscript
library(DESeq2)
library("optparse")

# arguments to provide
option_list = list(
  make_option(c("-c", "--counts"), type="character", default="results/counts.txt", help="counts tabulated file from Feature Counts", metavar="character"),
  make_option(c("-s", "--samplefile"), type="character", default=NULL, help="sample files used to get conditions for DESEq2 model fit", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="results/deseq/diff.tsv", help="where to place differential expression files", metavar="character")
) 

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

######## Import and wrangle counts produced by FeatureCounts #######
countdata <- read.table(opt$counts, header=TRUE, row.names=1,stringsAsFactors = F)
countdata <- countdata[ ,6:ncol(countdata)] # Remove first five columns (chr, start, end, strand, length)

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

######### Extract experimental conditions from the file specifying the samples
samplefile = read.delim(file = opt$samplefile,header = T,stringsAsFactors = F)
col.names = colnames(samplefile)                                                # extract all column names
columns.to.discard = c("sample","fq1","fq2")                                    # column we don't need to extract all conditions
conditions = col.names[! col.names %in% columns.to.discard]                     # only keeping the column of interest

# one condition
if (length(conditions) == 1){
  condition <- factor(samplefile[,"conditions"])
# two conditions
} else if (length(conditions == 2)){
  # two conditions --> make a third column that is the product of the two
  samplefile$condition = with(samplefile,paste(condition[1],condition[2],sep = "-"))
} else if (length(conditions > 2)){
   print("too many conditions to compare, skipping")
}

(condition <- factor(args[2:length(args)]))
#(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value and mrge to normalized counts
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="results/result.csv")