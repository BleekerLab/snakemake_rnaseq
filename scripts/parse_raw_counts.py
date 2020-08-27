#!/usr/bin/env python
# coding: utf-8

# In[19]:


import sys
import pandas as pd


# In[26]:


raw_counts_filepath = sys.argv[1]
star_results_directory_name = sys.argv[2]
out_parsed_counts_filepath = sys.argv[3]

raw_counts = pd.read_csv(raw_counts_filepath, skiprows=1, sep="\t")
raw_counts_parsed = raw_counts.drop(["Chr","Start","End","Strand","Length"], axis=1)


# In[27]:


old_col_names = raw_counts_parsed.columns.values.tolist()


# In[28]:


new_col_names = [col.replace(star_results_directory_name, '') for col in old_col_names]


# In[29]:


new_col_names = [col.replace('_Aligned.sortedByCoord.out.bam', '') for col in new_col_names]
new_col_names


# In[32]:


raw_counts_parsed.columns = new_col_names
raw_counts_parsed.to_csv(out_parsed_counts_filepath, sep="\t", index=False)


# In[ ]:




