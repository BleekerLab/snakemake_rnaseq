#!/usr/bin/env python
# coding: utf-8

# In[19]:


import sys
import pandas as pd


raw_counts_filepath = sys.argv[1]
star_results_directory_name = sys.argv[2]
out_parsed_counts_filepath = sys.argv[3]

############################################
# Parse raw counts to have clean column names
#############################################
raw_counts = pd.read_csv(raw_counts_filepath, skiprows=1, sep="\t")
raw_counts_parsed = raw_counts.drop(["Chr","Start","End","Strand","Length"], axis=1)

old_col_names = raw_counts_parsed.columns.values.tolist()

new_col_names = [col.replace(star_results_directory_name, '') for col in old_col_names]


new_col_names = [col.replace('_Aligned.sortedByCoord.out.bam', '') for col in new_col_names]
new_col_names

raw_counts_parsed.columns = new_col_names

######################################
# Write parsed raw counts to .csv file
######################################
raw_counts_parsed.to_csv(out_parsed_counts_filepath, sep="\t", index=False)







