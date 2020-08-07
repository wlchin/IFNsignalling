library(reshape2)
library(data.table)

# this script creates a TF to ISG matrix from the direct TF to DE network
# load gene lists
# alpha and gamma are hallmark gene sets

alpha <- readRDS("resources/biomart_alpha_mm10.rds")
gamma <- readRDS("resources/biomart_gamma_mm10.rds")
comb <- unique(c(alpha$MGI.symbol,gamma$MGI.symbol))

# load the direct network of TF to DE genes

edgelist <- readRDS(snakemake@input[[1]])

# filter to keep only TFs as source nodes and ISGs as sink nodes

edgelist_filtered <- edgelist[edgelist$targetGene %in% comb,]

# reshape to matrix - TFs on y axis and ISGs on x axis
# GENIE3 weights as values

TFtoISGweights <- reshape2::dcast(edgelist_filtered, regulatoryGene ~ targetGene, value.var = "weight")

# NAs (no direct interactions) assign 0
TFtoISGweights[is.na(TFtoISGweights)] <- 0

# save matrix of weights for plots
ISGmat <- as.matrix(TFtoISGweights[,-c(1)]) #remove first column of TFs
rownames(ISGmat) <- TFtoISGweights$regulatoryGene   
colnames(ISGmat) <- colnames(TFtoISGweights)[-1]

saveRDS(ISGmat, file = snakemake@output[[1]])


