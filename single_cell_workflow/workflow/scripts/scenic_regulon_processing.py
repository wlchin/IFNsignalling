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
    from pyscenic.prune import prune2df
    from pyscenic.transform import df2regulons
    from pyscenic.aucell import aucell
    from pyscenic.binarization import binarize

    with open(snakemake.input[0], "rb") as f:
        adjacencies = pickle.load(f)

    ex_matrix = pd.read_csv(snakemake.input[1], sep='\t', header=0, index_col=0).T
        
    ## load the ex_matrix and the regulons
    
    DATABASES_GLOB = os.path.join("resources/network_analysis", "mm10_*.mc9nr.feather")    
    MOTIF_ANNOTATIONS_FNAME = os.path.join("resources/network_analysis", "motifs-v9-nr.mgi-m0.001-o0.0.tbl")

    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    print(dbs)
    
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, num_workers = 4)

    print("prune2df done, now saving")
    with open(snakemake.output[0], "wb") as f:
        pickle.dump(df, f)
    
    print("df2regulons carrying out")
    regulons = df2regulons(df)

    print("prunedone, now saving")
    with open(snakemake.output[1], "wb") as f:
        pickle.dump(regulons, f)
