# retrieving conditions from command line arguments
args <- commandArgs(trailingOnly=TRUE)
# Import data from featureCounts

countdata <- read.table("results/counts.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("\\mapped.Pehyb_Quatt", "", colnames(countdata))
# Convert to matrix
countdata <- as.matrix(countdata)

# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(args[2:length(args)]))
#(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

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