#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
from functools import reduce
from optparse import OptionParser


####
#CLI
####
parser = OptionParser()
parser.add_option("-d", "--directory", default=None, dest="directory",
                  help="directory with the STAR mapping reports (*Log.final.out).")
parser.add_option("-o", "--outfile", default="None", dest="outfile", 
	              help="Tabulated separated file where mapping summaries for all samples will be written.") 
(options, args) = parser.parse_args()



############################################################
# Reads each file. Add sample name in the column with values
############################################################

list_of_logs = [f for f in os.listdir(options.directory) if f.endswith("Log.final.out")]
list_of_logs.sort()




sample_names = [log.replace("_Log.final.out","") for log in list_of_logs]
list_of_dfs = [pd.read_csv(
	log, 
	sep = "\t", 
	names=["attribute", str(sample)]) 
  for log,sample in zip(list_of_logs,sample_names)
]
list_of_dfs[0].head()


# In[17]:


df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['attribute'],
                                            how='outer'), list_of_dfs)

df_merged.to_csv(options.outfile,sep="\t")







