
[![DOI](https://zenodo.org/badge/282825263.svg)](https://zenodo.org/badge/latestdoi/282825263)

# Dynamic interferon beta signalling underlies treatment response to immune checkpoint therapy in cancer

### Introduction
This repository contains the snakemake workflows for the RNAseq analysis for the manuscript. Each folder contains a snakemake workflow for the three major parts of the analysis, namely bulk RNAseq data, single RNAseq data and analysis of the human dataset from Griffiths <em>et al.</em> (2020).

### Attributions
This workflow uses [custom code](https://gist.github.com/Vessy/6562505) from Vesna Memisevic for constructing the HiveR object used for hiveplot visualisations.

The repository also contains java libraries from the *[CellRouter](https://github.com/edroaldo/cellrouter)* package (Lummertz da Rocha et al. (2018) since this software package is not available via the Conda package manager. In addition, some older dependencies used by this package are not compatible with newer R enviroments, so snakemake will run this analysis in a singularity container.

### References

1. Griffiths JI, Wallet P, Pflieger LT, et al. Circulating immune cell phenotype dynamics reflect the strength of tumor-immune cell interactions in patients during immunotherapy. Proc Natl Acad Sci U S A. 2020;117(27):16072-16082. doi:10.1073/pnas.1918937117

2. Lummertz da Rocha E, Rowe RG, Lundin V, et al. Reconstruction of complex single-cell trajectories using CellRouter. Nat Commun. 2018;9(1):892. Published 2018 Mar 1. doi:10.1038/s41467-018-03214-y
