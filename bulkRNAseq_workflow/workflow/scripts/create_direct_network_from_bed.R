library(data.table)
library(dplyr)

## this script appends GENIE3 weights to the direct interactions between TFs and DE genes
## JASPAR TFBS bedfile has TFs in format e.g. FOS::JUN and SMAD2::SMAD3::SMAD4
## the first part of the script parses these formats, splitting into individual genes

processedbed <- read.table(snakemake@input[[1]], stringsAsFactors = F)

message("processing non dimeric and non trimeric TFs")

edges_single <- processedbed[!grepl("::", processedbed$V10),]
singles <- data.frame(regulatoryGene = edges_single[,10], targetGene = edges_single[,4], stringsAsFactors = F)

message("processing dimeric/trimeric TFs")

edges_multimers <- processedbed[grepl("::", processedbed$V10),]

multimerdf <- list()

for(i in 1:nrow(edges_multimers)){
    edges_multimers_regulatory <- as.character(edges_multimers[i,10])
    edges_multimers_target <- as.character(edges_multimers[i,4])
    edgepair <- data.frame(regulatoryGene = unlist(strsplit(edges_multimers_regulatory, "::")), targetGene = edges_multimers_target, stringsAsFactors = F)
    multimerdf[[i]] <- edgepair
}

multimer <- do.call(rbind, multimerdf)
directnetwork <- rbind(singles, multimer)

directnetwork <- directnetwork[!duplicated(directnetwork),] # remove duplicated entries
saveRDS(directnetwork, snakemake@output[[1]])
