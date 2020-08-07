#!/usr/bin/env python3
# coding: utf-8

# In[1]:

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import logging
import pandas as pd
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d


# In[2]:

#cluster15 = pd.read_csv("special_cluster15.csv")


# In[6]:


#cluster15.head()
#cluster15["cellname"] = cluster15["loom_cellid"]


# In[31]:


#set(cluster15["sample"])


# In[5]:


#cluster15["cellname"]


# In[19]:


cellid = pd.read_csv(snakemake.input[2])


# In[20]:


cellid.head()


# In[32]:


cellid.head()
cellid.index = cellid["final_cellnames"]


# In[47]:


# this is where I do massive filtering
respondervec = [ "NKD180900301", "NKD180900302", "NKD180900305", "NKD180900307","NKD180900310", "NKD180900311"]
ab1 = ["NKD180900302", "NKD180900303", "NKD180900304","NKD180900308", "NKD180900311", "NKD180900307"]
cellid["strain"] = np.where(cellid["sample"].isin(ab1), "ab1", "renca")
cellid["response"] = np.where(cellid["sample"].isin(respondervec), "responder", "nonresponder")



# In[48]:


#cellid["phenotype"] = cellid["sample"].map(typeofsample)


# In[49]:


cellid.head()


# In[37]:


renca_df = cellid[cellid["strain"] == "renca"]
ab1_df = cellid[cellid["strain"] == "ab1"]
renca_df.shape
#ab1_df.shape


# In[10]:


#cellid["CellID"] = cellid["sample"].map(str) + "_" + cellid["cellnames"].map(str)


# In[38]:


#cluster15mono = cellid['final_cellnames'][cellid["cellnames"].isin(cluster15["cellname"])]
monos = pd.read_csv(snakemake.input[3])
cluster15mono = monos["loom_cellid"]
#cluster15mono = cellid[cellid["cell_seurat_cluster"] == 14].index


# In[39]:


renca = vcy.load_velocyto_hdf5(snakemake.input[1])
a = pd.DataFrame(renca.delta_embedding)
#a["Irf"] = renca.velocity[94]
a.index = renca.ca["CellID"]
ly6 = a.reindex(cluster15mono)
ly6 = ly6.dropna()
ly6.columns = ["deltax", "deltay"]
ly6["vel"] = np.square(ly6["deltax"]) + np.square(ly6["deltay"])
ly6 = ly6.dropna()
ly6renca = pd.concat([ly6, cellid], axis = 1, join = "inner")
#ly6renca.to_csv("cluster15_renca.csv")
#ly6renca
del(renca)


# In[41]:


ab1 = vcy.load_velocyto_hdf5(snakemake.input[0])
b = pd.DataFrame(ab1.delta_embedding)
#a["Irf"] = renca.velocity[94]
b.index = ab1.ca["CellID"]
ly6 = b.reindex(cluster15mono)
ly6 = ly6.dropna()
ly6.columns = ["deltax", "deltay"]
ly6["vel"] = np.square(ly6["deltax"]) + np.square(ly6["deltay"])
ly6 = ly6.dropna()
ly6ab1 = pd.concat([ly6, cellid], axis = 1, join = "inner")
#ly6ab1.to_csv("cluster15_ab1.csv")
ly6ab1
del(ab1)


# In[45]:





# In[54]:


import seaborn as sns

ly6 = pd.concat([ly6renca, ly6ab1], axis = 0, join = "inner")

g = sns.FacetGrid(ly6, col="strain", hue = "response", sharex=True, sharey = False, col_wrap=2)
g.map(sns.kdeplot, "vel")
g.add_legend();
g.savefig(snakemake.output[0], dpi = 400)


# In[70]:


import seaborn as sns

ly6 = pd.concat([ly6renca, ly6ab1], axis = 0, join = "inner")

g = sns.FacetGrid(ly6, col="response", hue = "strain", sharex=True, sharey = False, col_wrap=2)
g.map(sns.kdeplot, "vel")
g.add_legend();
g.savefig(snakemake.output[1], dpi = 400)


