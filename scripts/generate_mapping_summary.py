#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import os
from functools import reduce
import sys
from yaml import safe_load

directory_with_mapping_reports = sys.argv[1]
config_file_path = sys.argv[2] # to extract the name of the result directory from the config file
mapping_summary = sys.argv[3]

############################################################
# Reads each file. Add sample name in the column with values
############################################################

list_of_logs = [directory_with_mapping_reports + f for f in os.listdir(directory_with_mapping_reports) if f.endswith("Log.final.out")]
list_of_logs.sort()


sample_names = [log.replace("_Log.final.out","") for log in list_of_logs]
list_of_dfs = [pd.read_csv(log, sep = "\t", names=["attribute", str(sample)]) 
  for log,sample in zip(list_of_logs, sample_names)
]

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['attribute'], how='outer'), list_of_dfs)


########################################################################
# parse column names to get rid of the directory before the sample names
#########################################################################

if os.path.exists(config_file_path):
    with open(config_file_path, 'r') as f:
        # load configfile info
        config_info = safe_load(f)
else:
    print("Config file does not exist. Please check that you have a config/config.yaml '{}'".format(config_file_path), file=sys.stderr)

RESULT_DIR = config_info["result_dir"]

df_merged_parsed_colnames = df_merged.rename(columns = lambda x:x.replace(RESULT_DIR + "star/", ""))

####################
# Write to .csv file
####################
df_merged_parsed_colnames.to_csv(mapping_summary, sep=",", index=False)







