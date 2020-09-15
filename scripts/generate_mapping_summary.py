#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import os
from functools import reduce
import sys

directory_with_mapping_reports = sys.argv[1]
mapping_summary = sys.argv[2]

############################################################
# Reads each file. Add sample name in the column with values
############################################################

list_of_logs = [f for f in os.listdir(directory_with_mapping_reports) if f.endswith("Log.final.out")]
list_of_logs.sort()


sample_names = [log.replace("_Log.final.out","") for log in list_of_logs]
list_of_dfs = [pd.read_csv(log, sep = "\t", names=["attribute", str(sample)]) 
  for log,sample in zip(list_of_logs, sample_names)
]


df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['attribute'],
                                            how='outer'), list_of_dfs)
df_merged.to_csv(mapping_summary, sep="\t")







