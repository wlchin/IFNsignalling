
library(scde)

df <- readRDS(snakemake@input[[1]])

test <- clean.counts(y[,rownames(df)], min.lib.size = 1800,  min.detected = 500, min.reads = 10)

df_clean <- df[rownames(df) %in% colnames(test),]

saveRDS(df_clean, snakemake@output[[1]])

knn <- knn.error.models(test, k = ncol(cd)/2, n.cores = 8, max.model.plots = 10)

saveRDS(knn, snakemake@output[[2]])

varinfo <- pagoda.varnorm(knn, counts = test, trim = 3/ncol(test), max.adj.var = 5, n.cores = 4, plot = TRUE)

saveRDS(varinfo, snakemake@output[[3]])
