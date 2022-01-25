#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import anndata as ad
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


notes = pd.read_table('../Model_Runs/model_run_notes.txt')


# In[4]:


notes.head()


# In[5]:


cn = notes.iloc[13:37,0].to_list()
cn = [x.split(" ")[1] for x in cn]


# In[6]:


cn


# In[ ]:





# In[7]:


file_names = [x for x in os.listdir('../Model_Runs') if not "txt" in x]


# In[8]:


file_names


# In[9]:


a = pd.read_csv(os.path.join('../Model_Runs',file_names[0]), header=None, names=cn)


# In[10]:


a.head()


# In[11]:


a.insert(0, "time", a.index*0.1)


# In[12]:


a.head()


# In[ ]:





# In[13]:


plt.plot(a.iloc[0:10000,]['time'],a.iloc[:10000,]['CycE'])
plt.title('CycE')
plt.xlabel('Time in hours')
plt.ylabel('Num of molec')
plt.show()


# ## cellXgene
# For analysis using cellXgene (a project by the Chan-Zuckerberg Initiative), must get data in `anndata` format. See https://anndata.readthedocs.io/en/latest/

# In[ ]:





# In[ ]:





# In[ ]:




