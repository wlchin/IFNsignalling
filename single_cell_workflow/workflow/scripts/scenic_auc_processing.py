if __name__ == "__main__":

    import os
    import glob
    import pickle
    import pandas as pd
    import numpy as np

    #from dask.diagnostics import ProgressBar

    #from arboreto.utils import load_tf_names
    #from arboreto.algo import grnboost2

    from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies, load_motifs
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    from pyscenic.binarization import binarize


    
    with open(snakemake.input[0], "rb") as f:
        regulons = pickle.load(f)

    ex_matrix = pd.read_csv(snakemake.input[1], sep='\t', header=0, index_col=0).T
            
    print("mtx print")
    
    auc_mtx = aucell(ex_matrix, regulons)

    thresholds, mat = binarize(auc_mtx)

    print("binarize done")

    print("binarise save")
    thresholds.to_csv(snakemake.output[0])
