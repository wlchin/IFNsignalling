
source('scripts/CellRouter_Class.R')
library('dplyr')
library(Seurat)
library(RColorBrewer)

AB1 <- readRDS(snakemake@input[[1]])

(Idents(AB1) <- 'seurat_clusters')
ab1_monos <- subset(AB1, idents = c("4", "12", "3"))

ab1_mono_307 <- ab1_monos
exp <- as.data.frame(as.matrix(ab1_mono_307@assays$SCT@data))

cellrouter <- CellRouter(exp)
cellrouter <- scaleData(cellrouter)
cellrouter <- computePCA(cellrouter, num.pcs = 50, seed=42)
lab <- as.character(ab1_mono_307$seurat_clusters)

num <- length(lab)

for(i in 1:num){
if(lab[i] == "12"){
    lab[i] <- "Cluster 1 Monocytes"}
else if(lab[i] == "3"){
        lab[i] <- "Cluster 2 Monocytes"}
else if(lab[i] == "4"){
        lab[i] <- "Cluster 3 Monocytes"}
else{lab[i] <- "."}
}

samples <- data.frame(clustering = paste0("cluster_", ab1_mono_307$seurat_clusters, "_", ab1_mono_307$response),
                      row.names = colnames(ab1_mono_307))
labtot <- paste0(lab, " (", ab1_mono_307$response, ")") 

samples <- data.frame(clustering = labtot,
                      row.names = colnames(ab1_mono_307))

testcolours <- palette(brewer.pal(n = 6, name = "Dark2"))
names(testcolours) <- unique(labtot)
colorvec <- testcolours[labtot]

nicecolours <- data.frame(nicecolours = colorvec,
                      row.names = colnames(ab1_mono_307))

cellrouter <- addInfo(cellrouter,
                      metadata = samples,
                      colname = 'clustering',
                      metadata.column='clustering')

cellrouter <- addInfo(cellrouter,
                      metadata = nicecolours,
                      colname = 'nicecolours',
                      metadata.column='nicecolours')

coords <- ab1_mono_307$umap@cell.embeddings
cellrouter <- customSpace(cellrouter, coords)

markers <- findSignatures(cellrouter, column = "clustering", pos.only = TRUE, fc.threshold = 0.5)


top10 <- markers %>% group_by(population) %>% top_n(5, fc)

filename <- snakemake@output[[1]]

plotSignaturesHeatmap(cellrouter, markers = top10, genes.show = top10$gene, column.ann = 'clustering', column.color = 'nicecolours', num.cells = 200, threshold = 3, width = 20, height = 12, filename=filename)
