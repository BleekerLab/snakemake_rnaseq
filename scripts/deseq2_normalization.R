#!/usr/bin/env Rscript

# Rscript --vanilla deseq2_normalization.R raw_counts.txt scaled_counts.txt
args = commandArgs(trailingOnly=TRUE)


# i want to prove that Frans is wrong

######################################
# import custom function and libraries
#' deseq_normalisation: normalise the raw counts of a matrix using the median of ratios method from DESeq2 
#' 
#'
#' @param data A matrix of raw counts
#'
#' @return A normalised matrix
#' @export
mor_normalization = function(data){

  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>%
    rownames_to_column('gene') %>%
    mutate(gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseudo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseudo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  #print(scaling_factors)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}
#################################################



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

