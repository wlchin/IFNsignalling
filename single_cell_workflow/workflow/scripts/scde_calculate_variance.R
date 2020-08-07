
library(scde)

cd <- readRDS(snakemake@input[[1]])
knn <- readRDS(snakemake@input[[2]])

varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 6, plot = TRUE)

saveRDS(varinfo, file = snakemake@output[[1]])
