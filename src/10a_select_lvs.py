#!/usr/bin/env python
# coding: utf-8

# # Select LVs
# Having 200 LVs is way more manageable than tens of thousands of genes, but that's still a lot of hypotheses to test. Inspired by the discussion section of the MultiPLIER paper, we'll filter out the LVs with low variance in the target dataset under the hypothesis that they are more likely to be reflective of non-brain biology.

# In[1]:


import pandas as pd


# In[2]:


lv_df = pd.read_csv('../output/NAc_PFC_VTA_LVs.txt', sep='\t')
# Drop information about the individual samples and shuffle to avoid subconcious bias
lv_df = lv_df.drop('sample', axis='columns')
lv_df = lv_df.sample(frac=1).reset_index(drop=True)


# In[3]:


lv_df


# ### Distribution of LV variances

# In[4]:


lv_df.var().hist(bins=50)


# There seems to be a cut point between the body of the distribution and its tail at a variance of .01. We'll keep the most variable LVs (those above .01) for the analysis.

# In[5]:


lvs_to_keep = lv_df.loc[:,lv_df.var() > .01].columns


# ## Save results

# In[6]:


lv_df = pd.read_csv('../output/NAc_PFC_VTA_LVs.txt', sep='\t')
lv_df = lv_df.set_index('sample')
lv_df = lv_df.loc[:,lvs_to_keep]
lv_df.to_csv('../output/high_variance_lvs.tsv', sep='\t')

