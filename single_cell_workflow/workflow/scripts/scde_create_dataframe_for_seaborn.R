
library(scde)
library(Seurat)
library(dplyr)

go.env <- readRDS(snakemake@input[[1]])
go.env.custom <- readRDS(snakemake@input[[2]]) 
varinfo <- readRDS(snakemake@input[[3]])

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 6)
pwpca_custom <- pagoda.pathway.wPCA(varinfo, go.env.custom, n.components = 1, n.cores = 6)

# incorporate the additional custom enrichment into the same list object
pwpca[["fast_isgs"]] <- pwpca_custom[[1]]

## Step 2 create df

listofvals <- list()

for(i in 1:7){
vec_vals <- as.numeric(pwpca[[i]]$xv)
nameofterm <- names(pwpca)[i]
listofvals[[nameofterm]] <- vec_vals
}

df_new <- data.frame(do.call(cbind, listofvals))
rownames(df_new) <- colnames(pwpca[[i]]$xv)
df_new$CellID <- colnames(pwpca[[i]]$xv)

## Step 3 create UMAP coordinates and response data from seurat object and merge
seuratobject <- readRDS(snakemake@input[[4]])
metadata <- cbind(seuratobject@meta.data, seuratobject@reductions$umap@cell.embeddings)
metadata$CellID <- rownames(metadata)
df_comb <- inner_join(df_new, metadata, by = "CellID")

write.csv(df_comb, quote = F, file = snakemake@output[[1]])

