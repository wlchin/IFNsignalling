library(Seurat)

future::plan(strategy = 'multicore', workers = 1)
options(future.globals.maxSize = 10 * 1024 ^ 3)

filevec <- unlist(snakemake@input)

list_seurat <- list()

for(i in filevec){
    x <- readRDS(i)
    list_seurat[[i]] <- x
}
 

sample.list <- list_seurat
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features, verbose = FALSE)

sample.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", anchor.features = sample.features, verbose = FALSE)

sample.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT", verbose = FALSE)

sample.integrated <- RunPCA(sample.integrated,  verbose = FALSE)
sample.integrated <- RunUMAP(sample.integrated, reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, dims = 1:30, verbose = FALSE)
sample.integrated <- FindClusters(sample.integrated, verbose = FALSE)

saveRDS(sample.integrated, file = snakemake@output[[1]])
