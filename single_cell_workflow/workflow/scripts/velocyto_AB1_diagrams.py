#!/usr/bin/env python3
# coding: utf-8

# In[3]:


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


# In[4]:


vlm = vcy.load_velocyto_hdf5(snakemake.input[0])


# In[ ]:





# In[ ]:


## color condes consistent with main text

from collections import defaultdict
monocytelabels = defaultdict(lambda: 'Grey')

#ice_cream = {}
monocytelabels['4'] = '#9382ae'
monocytelabels['3'] = '#fdceb8'
monocytelabels['12'] = '#df5974'

testmap = pd.Series(vlm.cluster_labels).map(str).map(monocytelabels)
#testmap = vers.map(ice_cream)
#testmap  = seuratlabels

#print(testmap.head)

from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors


def gaussian_kernel(X, mu = 0, sigma=1):
    return np.exp(-(X - mu)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)


## plot mini figure

plt.figure(None,(10,10))

steps = 45, 45
grs = []
for dim_i in range(vlm.embedding.shape[1]):
    m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
    m = m - 0.025 * np.abs(M - m)
    M = M + 0.025 * np.abs(M - m)
    gr = np.linspace(m, M, steps[dim_i])
    grs.append(gr)

meshes_tuple = np.meshgrid(*grs)
gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T
gridpoints_coordinates = gridpoints_coordinates + norm.rvs(loc=0, scale=0.15, size=gridpoints_coordinates.shape)

nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(gridpoints_coordinates, 20)
ix_choice = ixs[:,0].flat[:]
ix_choice = np.unique(ix_choice)

nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(vlm.embedding[ix_choice], 20)
density_extimate = gaussian_kernel(dist, mu=0, sigma=0.5).sum(1)
bool_density = density_extimate > np.percentile(density_extimate, 10)
ix_choice = ix_choice[bool_density]

plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c=testmap, alpha=0.2, s=120, edgecolor="")
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c=testmap[ix_choice], alpha=1, s=120, edgecolor="k")

quiver_kwargs=dict(scale=5, headaxislength=10, headlength=10, headwidth=4,linewidths=0.4, edgecolors="k", color="k", alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           **quiver_kwargs)

plt.xlim(-7.5,10) 
plt.ylim(-5, 15)
plt.savefig(snakemake.output[0], dpi = 400)


# In[ ]:


## plot large quiver plot

plt.figure(None,(14,14))
quiver_scale = 10

plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c="0.8", alpha=0.2, s=10, edgecolor="")

ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/8), replace=False)
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=10,linewidths=0.25, width=0.00045,edgecolors="k", color=vlm.colorandum[ix_choice], alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           scale=quiver_scale, **quiver_kwargs)
plt.axis("off")
plt.savefig(snakemake.output[1], dpi = 400)


# In[ ]:


#plot interferon velocity across digrams 
plt.figure(None,(14,14))

vlm.plot_velocity_as_color("Irf1")

plt.savefig(snakemake.output[2], dpi = 400)

