
library(scde)
library(Seurat)
library(dplyr)

# this script reads scde output into a dataframe for downstream visualisation

go.env <- readRDS(snakemake@input[[1]])
go.env.custom <- readRDS(snakemake@input[[2]]) 
varinfo <- readRDS(snakemake@input[[3]])

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 6)
pwpca_custom <- pagoda.pathway.wPCA(varinfo, go.env.custom, n.components = 1, n.cores = 6)

# incorporate the additional custom enrichment into the same list object

pwpca[["fast_isgs"]] <- pwpca_custom[[1]]

listofvals <- list()

for(i in 1:7){
vec_vals <- as.numeric(pwpca[[i]]$xv)
nameofterm <- names(pwpca)[i]
listofvals[[nameofterm]] <- vec_vals
}

df_new <- data.frame(do.call(cbind, listofvals))
rownames(df_new) <- colnames(pwpca[[i]]$xv)
df_new$CellID <- colnames(pwpca[[i]]$xv)

saveRDS(df_new, file = snakemake@output[[1]])

metadata <- readRDS(snakemake@input[[4]])
metadata$CellID <- rownames(metadata)

combined_df <- dplyr::inner_join(df_new, metadata, by = "CellID")
write.table(combined_df, sep = "," quote = F, file = snakemake@output[[2]])
