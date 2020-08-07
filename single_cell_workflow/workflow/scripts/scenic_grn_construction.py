if __name__ == "__main__":

    import os
    import glob
    import pickle
    import pandas as pd
    import numpy as np

    #from dask.diagnostics import ProgressBar

    from arboreto.utils import load_tf_names
    from arboreto.algo import grnboost2

    #from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
    #from pyscenic.utils import modules_from_adjacencies, load_motifs
    #from pyscenic.prune import prune2df, df2regulons
    #from pyscenic.aucell import aucell
    #from pyscenic.binarization import binarize

    

    ex_matrix = pd.read_csv(snakemake.input[0], sep='\t', header=0, index_col=0).T
    tf_names = load_tf_names("resources/network_analysis/mm_mgi_tfs.txt")

    print(ex_matrix.shape)

    print("finish loading, now do grn")
    
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

    print("grn_complete, now dumping")
    
    with open(snakemake.output[0], "wb") as f:
        pickle.dump(adjacencies, f)
        
