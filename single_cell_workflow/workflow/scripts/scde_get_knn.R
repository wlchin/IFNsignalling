
library(scde)

cd <- readRDS(snakemake@input[[1]])
groups <- readRDS(snakemake@input[[2]])

knn <- knn.error.models(cd, groups = groups, k = ncol(cd)/4, n.cores = 5, save.model.plots = FALSE, verbose = 0)

saveRDS(knn, file = snakemake@output[[1]])
