
source('workflow/scripts/CellRouter_Class.R')
library('dplyr')
library(Seurat)

AB1 <- readRDS(snakemake@input[[1]])

(Idents(AB1) <- 'seurat_clusters')
ab1_monos <- subset(AB1, idents = c("4", "5", "12", "3"))# subset by monocyte clusters

exp <- as.data.frame(as.matrix(ab1_monos@assays$SCT@data)) # extract normalised count matrix for cellrouter object
cellrouter <- CellRouter(exp)
cellrouter <- scaleData(cellrouter)
cellrouter <- computePCA(cellrouter, num.pcs = 50, seed=42)

samples <- data.frame(clustering = paste0("cluster_", ab1_monos$seurat_clusters),
                      row.names = colnames(ab1_monos))

cellrouter <- addInfo(cellrouter,
                      metadata = samples,
                      colname = 'clustering',
                      metadata.column='clustering')

coords <- ab1_monos$umap@cell.embeddings
cellrouter <- customSpace(cellrouter, coords)

save(cellrouter, file = snakemake@output[[1]], version = 2)
