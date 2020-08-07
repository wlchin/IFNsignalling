
library(tidyr)
library(dplyr)

# this script creates aggregate matrices for TF activation following SCENIC analysis

AUC_results <- read.csv(snakemake@input[[1]], stringsAsFactors = F)

irf_TFs <- AUC_results[,grepl("IRF", colnames(AUC_results))]
stat_TFs <- AUC_results[,grepl("STAT", colnames(AUC_results))]

TFs <- cbind(stat_TFs, irf_TFs)
TFs$cellnames <- AUC_results$Cell

metadata <- readRDS(snakemake@input[[2]])
metadata$cellnames <- paste0("X",rownames(metadata)) #clean cellnames

# per cell TF activation as rds object

comb <- dplyr::inner_join(irf_all, metadata)
saveRDS(comb, file = snakemake@input[[1]])

# aggregated TF activation

TF_agg <- comb %>% group_by(timepoint, response) %>%
    summarise_at(colnames(irf_all),
                 mean, na.rm = T)
write.table(TF_agg, quote = F, file = snakemake@output[[2]])
