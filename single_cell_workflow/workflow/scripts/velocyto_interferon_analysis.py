#!/usr/bin/env python
# coding: utf-8


import velocyto as vcy
import pandas as pd
import seaborn as sns
import numpy as np


vlm = vcy.load_velocyto_hdf5(snakemake.input[0])

# create gene list of fast isgs and their TFs
genes = pd.read_csv(snakemake.input[1], header = None)
#fastgenes_list = genes[0].to_list()
fastgenes_list = genes[0].tolist() # depending on pandas version
list2 = ["Irf1", "Irf7", "Irf9", "Stat1", "Stat2"]
fastgenes = fastgenes_list + list2

# subset velocyto object information -  cluster 12 cells only and their ISG velocities
ab1_pheno = pd.read_csv(snakemake.input[2])
ab1_pheno.index = ab1_pheno["final_cellnames"]
cluster12cellsname = ab1_pheno[ab1_pheno["cell_seurat_cluster"] == 12]["final_cellnames"]

geneindex = np.isin(vlm.ra["Gene"],fastgenes)
genenames = vlm.ra["Gene"][geneindex]
cellindex = np.isin(vlm.ca["CellID"],cluster12cellsname)
cellindex_names = vlm.ca["CellID"][cellindex]

# create pandas df for seaborn clustermap plot
genemat = vlm.velocity[geneindex][:,cellindex]
genemat = np.absolute(genemat)
cluster12cells = pd.DataFrame(genemat)
cluster12cells.index = genenames
cluster12cells.columns = cellindex_names


# create response df for seaborn clustermap plot
df = ab1_pheno.reindex(cellindex_names).dropna()
lut = dict(zip(df["response"].unique(), "rbg"))
response_map = df["response"]
rcolors = response_map.map(lut)

a = sns.clustermap(cluster12cells, standard_scale = 0, col_colors = rcolors.tolist(), xticklabels = False, yticklabels = True) # convert
a.savefig(snakemake.output[0], dpi = 400)

# create df for facetgrid plot of individual gene kdes
widedf = np.absolute(cluster12cells.transpose())
combination = pd.concat([ab1_pheno, widedf], axis = 1, join = "inner")
small = combination.loc[:,["Irf1", "Isg20", "Irf7", "Gbp2", "response"]]
testplot = pd.melt(small, id_vars=['response'], value_vars=['Irf1', "Isg20", "Irf7", "Gbp2"])


g = sns.FacetGrid(testplot, col="variable", hue = "response", sharex=False, sharey = False, col_wrap=2)
g.map(sns.kdeplot, "value")
g.add_legend();
g.savefig(snakemake.output[1], dpi = 400)

