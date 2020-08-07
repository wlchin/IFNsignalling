
source('workflow/scripts/CellRouter_Class.R')
library('dplyr')
library(Seurat)

load(snakemake@input[[1]])

diffvec <- cellrouter@clusters[[1]][[3]]

saveRDS(diffvec, file = snakemake@output[[7]])
names <- unique(names(cellrouter@pathsinfo$distr))
plotClusterHeatmap(cellrouter, names, 10, 10, 2, snakemake@output[[6]])

ind <- length(table(diffvec))

for(i in 1:ind){
    filenm <- snakemake@output[[i]]
    genes <- names(x[x ==i])
    write(genes, filenm) 
}


