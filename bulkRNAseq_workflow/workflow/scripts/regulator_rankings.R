## this script loads the GENIE3 linklist
## it then measures regulator influence by summing edges based grouped_by regulatoryGenes

#library(data.table)

mat <- readRDS(snakemake@input[[1]])

regsums <- sort(rowSums(mat), decreasing = T)

## save regulators

saveRDS(regsums, file = snakemake@output[[1]])







