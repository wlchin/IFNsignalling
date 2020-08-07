#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import velocyto as vcy
import pandas as pd
import seaborn as sns
import numpy as np


# In[ ]:


cellid = pd.read_csv(snakemake.input[2])
cellid.index = cellid["final_cellnames"]


# In[ ]:


respondervec = [ "NKD180900301", "NKD180900302", "NKD180900305", "NKD180900307","NKD180900310", "NKD180900311"]
ab1 = ["NKD180900302", "NKD180900303", "NKD180900304","NKD180900308", "NKD180900311", "NKD180900307"]
cellid["strain"] = np.where(cellid["sample"].isin(ab1), "ab1", "renca")
cellid["response"] = np.where(cellid["sample"].isin(respondervec), "responder", "nonresponder")


# In[ ]:


targetmono = pd.read_csv(snakemake.input[3])
monos = targetmono["loom_cellid"]
ind = cellid["cellnames"].isin(monos)
cluster15mono = cellid.index[ind]


# In[ ]:


# extract projections from renca velocyto object and calculate momentum
renca = vcy.load_velocyto_hdf5(snakemake.input[1])
delta_embedding_renca = pd.DataFrame(renca.delta_embedding)
delta_embedding_renca.index = renca.ca["CellID"]
renca_ly6_cells_deltas = delta_embedding_renca.reindex(cluster15mono).dropna()
renca_ly6_cells_deltas.columns = ["deltax", "deltay"]
renca_ly6_cells_deltas["vel"] = np.square(renca_ly6_cells_deltas["deltax"]) + np.square(renca_ly6_cells_deltas["deltay"])
ly6crenca = pd.concat([renca_ly6_cells_deltas, cellid], axis = 1, join = "inner") # merge with cell metadata
del(renca)


# In[ ]:


# extract projections from ab1 velocyto object and calculate momentum
ab1 = vcy.load_velocyto_hdf5(snakemake.input[0])
delta_embedding_ab1 = pd.DataFrame(ab1.delta_embedding)
delta_embedding_ab1.index = ab1.ca["CellID"]
ab1_ly6_cells_deltas = delta_embedding_ab1.reindex(cluster15mono).dropna()
ab1_ly6_cells_deltas.columns = ["deltax", "deltay"]
ab1_ly6_cells_deltas["vel"] = np.square(ab1_ly6_cells_deltas["deltax"]) + np.square(ab1_ly6_cells_deltas["deltay"])
ly6cab1 = pd.concat([ab1_ly6_cells_deltas, cellid], axis = 1, join = "inner") # merge with cell metadata 
del(ab1)


# In[ ]:


ly6 = pd.concat([ly6crenca, ly6cab1], axis = 0, join = "inner")


# In[ ]:


g = sns.FacetGrid(ly6, col="strain", hue = "response", sharex=True, sharey = False, col_wrap=2)
g.map(sns.kdeplot, "vel")
g.add_legend();
g.savefig(snakemake.output[0], dpi = 400)


# In[ ]:


g = sns.FacetGrid(ly6, col="response", hue = "strain", sharex=True, sharey = False, col_wrap=2)
g.map(sns.kdeplot, "vel")
g.add_legend();
g.savefig(snakemake.output[1], dpi = 400)

