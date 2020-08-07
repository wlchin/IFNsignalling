library(GENIE3)

set.seed(123)

mat <- readRDS(snakemake@input[[1]])
mat <- mat[rowSums(mat) > 0,]

weightMat <- GENIE3(mat, nCores = 8, verbose = TRUE)

saveRDS(weightMat, file = snakemake@output[[1]])
