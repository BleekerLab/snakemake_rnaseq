#!/usr/bin/env Rscript

# Rscript --vanilla deseq2_normalization.R raw_counts.txt scaled_counts.txt
args = commandArgs(trailingOnly=TRUE)

# import custom function
source("deseq2_normalization_function.R")

df <- read.csv(
	x = args[1], 
	header = TRUE, 
	row.names = 1)

normalised_data <- mor_normalization(df)

write.table(
	x = normalised_data, 
	file = args[2], 
	sep = "\t", 
	quote = FALSE, 
	row.names = TRUE)

