#!/usr/bin/env python3
# coding: utf-8

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
import seaborn as sns



vlm = vcy.load_velocyto_hdf5(snakemake.input[0])
ab1_pheno = pd.read_csv(snakemake.input[1])

embedding = pd.DataFrame({'CellID':vlm.ca["CellID"], 'xdirect':vlm.delta_embedding[:,0], 'ydirect':vlm.delta_embedding[:,1]})
embedding["mag"] = np.square(embedding["xdirect"]) + np.square(embedding["ydirect"])
embedding.index = embedding["CellID"]
embedding

ab1_pheno.index = ab1_pheno["CellID"]
comb = pd.concat([ab1_pheno, embedding], axis = 1, join = "inner")
abbr = comb[(comb["cell_seurat_cluster"] == 12) | (comb["cell_seurat_cluster"] == 3) | (comb["cell_seurat_cluster"] == 5)]

g = sns.FacetGrid(abbr, col="cell_seurat_cluster", hue = "response")
g.map(sns.kdeplot, "mag")
g.add_legend();
g.savefig(snakemake.output[0], dpi = 400)

