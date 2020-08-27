#!/usr/bin/env Rscript

# Rscript --vanilla deseq2_normalization.R raw_counts.txt scaled_counts.txt
args = commandArgs(trailingOnly=TRUE)

# import custom function and libraries
source("scripts/deseq2_normalization_function.R")

if ("dplyr" %in% installed.packages()){
	library("dplyr")	
} else {
	install.packages("dplyr", repos = "http://cran.us.r-project.org")
	library("dplyr")
}

if ("tibble" %in% installed.packages()){
	library("tibble")	
} else {
	install.packages("tibble", repos = "http://cran.us.r-project.org")
	library("tibble")
}

# Read parsed counts and remove unnecessary columns
df <- read.delim(
	file = args[1], 
	header = TRUE)

df <- df %>% 
  tibble::column_to_rownames("Geneid")

normalised_data <- mor_normalization(df)

normalised_data = tibble::rownames_to_column(normalised_data, var = "gene")

write.table(
	x = normalised_data, 
	file = args[2], 
	sep = "\t", 
	quote = FALSE, 
	row.names = FALSE)

