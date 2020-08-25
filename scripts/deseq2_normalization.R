#!/usr/bin/env Rscript

# Rscript --vanilla deseq2_normalization.R raw_counts.txt scaled_counts.txt
args = commandArgs(trailingOnly=TRUE)

# import custom function and libraries
source("scripts/deseq2_normalization_function.R")

if ("tidyverse" %in% installed.packages()){
	library(tidyverse)	
} else {
	install.packages("tidyverse")
	library("tidyverse")
}


# Read parsed counts and remove unnecessary columns
df <- read.delim(
	file = args[1], 
	header = TRUE)

df <- df %>% 
  dplyr::select(- Chr, - Start, - End, - Strand, - Length) %>% 
  column_to_rownames("Geneid")

normalised_data <- mor_normalization(df)

normalised_data = rownames_to_column(normalised_data, 
                                     var = "gene")

write.table(
	x = normalised_data, 
	file = args[2], 
	sep = "\t", 
	quote = FALSE, 
	row.names = FALSE)

