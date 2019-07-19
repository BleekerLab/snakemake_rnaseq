#!/usr/bin/env Rscript
library(DESeq2)
library(optparse)

# arguments to provide
option_list = list(
  make_option(c("-c", "--counts"), type="character", default="results/counts.txt", help="counts tabulated file from Feature Counts", metavar="character"),
  make_option(c("-s", "--samplefile"), type="character", default="data/samples2.tsv", help="sample files used to get conditions for DESEq2 model fit", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="results/deseqNew.csv", help="where to place differential expression files", metavar="character"),
  make_option(c("-f", "--helperfile"), type="character", default="results/helperfile.csv", help="helper file needed for the creation of the clustering and the heatmaps", metavar="character"),
  make_option(c("-m", "--maxfraction"), type="double", default=1.0, help="maximum fraction of the total number of genes to be allowed to be differential between two conditions to be included (number between 0 and 1)", metavar="double")
) 

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

######## Import and wrangle counts produced by FeatureCounts #######
countdata <- read.table(opt$counts, header=TRUE, skip=1, row.names=1,stringsAsFactors = F, check.names = F)
countdata <- countdata[ ,6:ncol(countdata)] # Remove first five columns (chr, start, end, strand, length)

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

######### Extract experimental conditions from the file specifying the samples
samplefile = read.delim(file = opt$samplefile,header = T,stringsAsFactors = F)
col.names = colnames(samplefile)                                                # extract all column names
columns.to.discard = c("sample","fq1","fq2")                                    # column we don't need to extract all conditions
colsForConditions = col.names[! col.names %in% columns.to.discard]                     # only keeping the column of interest

# one condition
if (length(colsForConditions) == 1){
  condition <- factor(samplefile[,colsForConditions])
  # two conditions
} else if (length(colsForConditions == 2)){
  # two conditions --> make a third column that is the product of the two
  samplefile$conditions = paste(samplefile[,colsForConditions[1]],samplefile[,colsForConditions[2]],sep = ".")
  condition <- factor(x = samplefile[,"conditions"],levels = unique(samplefile[,"conditions"]))
} else if (length(conditions > 2)){
  print("too many conditions to compare, skipping")
}

# create helper file for downstream use in the pipeline.
helperFile <- as.data.frame(condition)
rownames(helperFile) <- samplefile$sample
colnames(helperFile) <- NULL
print(helperFile)
write.csv(helperFile, file = opt$helperfile, quote = F, col.names = F, row.names = T)

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

# create dataframe containing normalized counts, to which the differential expression values will be added
resdata <- as.data.frame(counts(dds, normalized=TRUE))

## function for Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(padj), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(padj), labs=colnames(res), cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

# open file that will be containing a vulcano plot for each of the combinations of conditions
#pdf(file="vulcanoplots.pdf")

# iterate through the list of conditions to create differential expression (DE) values for all possible pairs
x <- 1
for(i in levels(condition)){
  x <- x + 1
  if (x <= length(levels(condition))){
    for(j in x:length(levels(condition))){
      res <- results(dds, contrast=c("condition", i,  levels(condition)[j]))                      # i = first in pair, levels(condition)[j] is the second in pair.
      d <- paste(i, levels(condition)[j], sep="&")                                                # paste the two conditions in one character, to be used as the pair name
      resP <- as.data.frame(table(res$padj<0.05))                                                 # get number of DE values with P-value < 0.05
      if(resP[1,2]< opt$maxfraction*nrow(resdata)){                                               # only continue with the pair if it is less then the maximum fraction(set by user in commandline) differentially expressed 
        #volcanoplot(res, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-10, 10), main=d)
        print(c(d,"Number of differentials is within accepted limit"))
        colnames(res) = paste(d,c(colnames(res)),sep = "-")                                       # paste the pair name to the column name
        resdata <- merge(as.data.frame(resdata), as.data.frame(res), by="row.names", sort=FALSE)  # merge the DE values to the matrix
        rownames(resdata) <- resdata$Row.names                                                    # redifine rownames as they have disapeared?? in the merging??
        resdata$Row.names <- NULL                                                                 # delete column containing the rownames, to avoid the same colname in the next iteration
      }
      else{
        print(c(d,"More differentials then allowed"))
      }
    }
  }
}
#dev.off()
# create first column containing the genenames.
resdata$genes <- rownames(resdata)
resdata <- merge(as.data.frame(resdata["genes"]), as.data.frame(resdata[,2:ncol(resdata)-1]), by="row.names", sort=FALSE)
resdata$Row.names <- NULL

# write the data to a file.
write.table(resdata, file=opt$outdir,sep = "\t",quote=F,row.names=F)
