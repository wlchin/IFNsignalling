
source('workflow/scripts/CellRouter_Class.R')
library('dplyr')

load(snakemake@input[[1]])

# get genelists for each component
diffvec <- cellrouter@clusters[[1]][[3]]
ind <- length(table(diffvec))

for(i in 1:ind){
    filenm <- snakemake@output[[i]]
    genes <- names(diffvec[diffvec ==i])
    write(genes, filenm) 
}

# store RDS and plot diffusion components
saveRDS(diffvec, file = snakemake@output[[6]])
names <- unique(names(cellrouter@pathsinfo$distr))
plotClusterHeatmap(cellrouter, names, 10, 10, 2, snakemake@output[[7]])



