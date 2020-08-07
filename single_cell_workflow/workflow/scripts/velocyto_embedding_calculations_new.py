#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import velocyto as vcy
import pandas as pd


# In[ ]:


vlm = vcy.VelocytoLoom(snakemake.input[0])


# In[ ]:


seurat_object_cells = pd.read_csv(snakemake.input[1])
seurat_object_cells.index = seurat_object_cells["final_cellnames"]

loom_object_cells = pd.DataFrame(vlm.ca)
loom_object_cells.index = loom_object_cells['CellID']

cells_present_in_both = pd.concat([seurat_object_cells,loom_object_cells], axis=1, join='inner')


# In[ ]:


#create a boolean mask for use with vlm.filter_cells
loom_cells = loom_object_cells.index
seurat_cells = seurat_object_cells.index
ind = np.isin(loom_cells, seurat_cells)
vlm.filter_cells(ind)


# In[ ]:


#reindex the common_cells dataframe 
cells_present_in_both = cells_present_in_both.reindex(vlm.ca["CellID"]).dropna()


# In[ ]:


# assign clusters and UMAP coordinates using reindexed dataframe
vlm.set_clusters(cells_present_in_both["cell_seurat_cluster"])
vlm.ts = np.column_stack([cells_present_in_both["UMAP_1"], cells_present_in_both["UMAP_2"]])


# In[ ]:


# filter genes and normalize
vlm.score_cv_vs_mean(3000, plot=False, max_expr_avg=35)
vlm.score_detection_levels(min_expr_counts=5, min_cells_express=30)
vlm.score_cluster_expression()
vlm.filter_genes(by_detection_levels=True, by_cluster_expression = True, by_cv_vs_mean=True)
vlm.normalize()


# In[ ]:


vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=1000, b_maxl=1500, n_jobs=12)


# In[ ]:


vlm.fit_gammas()
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift()
vlm.extrapolate_cell_at_t(delta_t=1.)
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform = "sqrt", n_neighbors=1000, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.2)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=500)


# In[ ]:


vlm.to_hdf5(snakemake.output[0])

