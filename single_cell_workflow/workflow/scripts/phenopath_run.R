
library(phenopath)

pheno <- readRDS(snakemake@input[[2]])
mat <- readRDS(snakemake@input[[1]])

fit <- phenopath(t(mat), pheno, elbo_tol = 1e-6, thin = 40)

saveRDS(fit, file = snakemake@output[[1]])
