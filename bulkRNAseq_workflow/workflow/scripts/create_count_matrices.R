#!/usr/bin/env Rscript

library("tximport")
library("dplyr")

## read phenodata table and gene to transcript map 
tot <- readRDS(snakemake@input[[1]])
phenodata <- readRDS(snakemake@input[[2]])

## append file path to hd5 file
phenodata$file_location <- file.path(phenodata$path,"abundance.h5")


## create two phenodata dfs - one for responders and one for non responders
responders <- dplyr::filter(phenodata, Group == "RS")
nonresponders <- dplyr::filter(phenodata, Group == "NR")

## use tximport to create R and NR count matrices as input to GENIE3

kallisto_countsR <- tximport(responders$file_location, type = "kallisto", tx2gene = tot, ignoreAfterBar = F)
kallisto_countsR_mat <- kallisto_countsR$counts
colnames(kallisto_countsR_mat) <- responders$sample
saveRDS(kallisto_countsR_mat, file = snakemake@output[[1]])

kallisto_countsNR <- tximport(nonresponders$file_location, type = "kallisto", tx2gene = tot, ignoreAfterBar = F)
kallisto_countsNR_mat <- kallisto_countsNR$counts
colnames(kallisto_countsNR_mat) <- nonresponders$sample
saveRDS(kallisto_countsNR_mat, file = snakemake@output[[2]])



