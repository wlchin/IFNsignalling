library("sleuth")

#read phenodata table and append sample name and paths for kallisto

raw_phenodata <- read.csv(snakemake@input[[1]], stringsAsFactors = F) 

for(i in 1:nrow(raw_phenodata)){
    regx <- raw_phenodata$AGRF_ID[i] #pattern to find correct folder
    x <- list.files("results/kallisto_aligned", pattern = regx, full.names = T)
    y <- list.files("results/kallisto_aligned", pattern = regx, full.names = F)
    raw_phenodata$path[i] <- x #append path to df
    raw_phenodata$sample[i] <- y #append sample name to df
}

## select only relevant columns

metadata <- dplyr::select(raw_phenodata, sample, Timepoint, Group, path)

##load pre-processed gencode gene to transcript map

tot <- readRDS(snakemake@input[[2]])

## run sleuth with gene aggregation mode TRUE; LRT test for DE

so <- sleuth_prep(metadata, target_mapping = tot,
                  aggregation_column = 'gene_name', extra_bootstrap_summary = TRUE, num_cores = 1)

so <- sleuth_fit(so, ~Timepoint + Group + Timepoint:Group, 'full')
so <- sleuth_fit(so, ~Timepoint + Group, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

## create list of DE genes based on qval =<0.05 and save sleuth gene table output as rds

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene_filtered <- dplyr::filter(sleuth_table_gene, qval <= 0.05)

write(sleuth_table_gene_filtered$target_id, file = snakemake@output[[1]])
saveRDS(sleuth_table_gene, file = snakemake@output[[2]])

## save sleuth object and phenodata file for later analysis
saveRDS(so, file = snakemake@output[[3]])
saveRDS(metadata, file = snakemake@output[[4]])
