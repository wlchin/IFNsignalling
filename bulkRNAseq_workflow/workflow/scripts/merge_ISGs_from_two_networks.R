library(reshape2)
library(data.table)
library(GENIE3)
library(dplyr)

## this script creates a TF to ISG matrix from the direct TF to DE network
## load gene lists
## alpha and gamma are hallmark gene sets

alpha <- readRDS("resources/biomart_alpha_mm10.rds")
gamma <- readRDS("resources/biomart_gamma_mm10.rds")

comb <- unique(c(alpha$MGI.symbol,gamma$MGI.symbol))

ab1 <- readRDS(snakemake@input[[1]])
renca <- readRDS(snakemake@input[[2]])

ab1_full_GENIE3 <- readRDS(snakemake@input[[3]])
renca_full_GENIE3 <- readRDS(snakemake@input[[4]])

## get common between two columns

comgenes <- c(ab1$targetGene, renca$targetGene)
com_genes_ifn <- intersect(comgenes, comb)

TFs_in_direct_network <- unique(c(ab1$regulatoryGene, renca$regulatoryGene))

ab1_abbr_mat <- ab1_full_GENIE3[rownames(ab1_full_GENIE3) %in% TFs_in_direct_network, com_genes_ifn]
ab1_full_link_isg <- getLinkList(ab1_abbr_mat)

renca_abbr_mat <- renca_full_GENIE3[rownames(renca_full_GENIE3) %in% TFs_in_direct_network, com_genes_ifn]
renca_full_link_isg <- getLinkList(renca_abbr_mat)

### get the direct connections union set

ab1_ISGtarget <- ab1[ab1$targetGene %in% comb,]
renca_ISGtarget <- renca[renca$targetGene %in% comb,]
intergratedlist <- rbind(ab1_ISGtarget[,1:2], renca_ISGtarget[,1:2])
dum <- intergratedlist[!duplicated(intergratedlist),] 

rem_ab1 <- left_join(dum, ab1_full_link_isg)
rem_renca <- left_join(dum, renca_full_link_isg)

TFtoISGweights_ab1 <- reshape2::dcast(rem_ab1, regulatoryGene ~ targetGene, value.var = "weight")
TFtoISGweights_renca <- reshape2::dcast(rem_renca, regulatoryGene ~ targetGene, value.var = "weight")

TFtoISGweights_ab1[is.na(TFtoISGweights_ab1)] <- 0
ISGmat_ab1 <- as.matrix(TFtoISGweights_ab1[,-c(1)]) #remove first column of TFs
rownames(ISGmat_ab1) <- TFtoISGweights_ab1$regulatoryGene   
colnames(ISGmat_ab1) <- colnames(TFtoISGweights_ab1)[-1]

TFtoISGweights_renca[is.na(TFtoISGweights_renca)] <- 0
ISGmat_renca <- as.matrix(TFtoISGweights_renca[,-c(1)]) #remove first column of TFs
rownames(ISGmat_renca) <- TFtoISGweights_renca$regulatoryGene   
colnames(ISGmat_renca) <- colnames(TFtoISGweights_renca)[-1]

saveRDS(ISGmat_renca, file = snakemake@output[[1]])
saveRDS(ISGmat_ab1, file = snakemake@output[[2]])
