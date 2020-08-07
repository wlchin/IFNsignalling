
library(Seurat)

x <- readRDS(snakemake@input[[1]])

umapAB1 <- data.frame(x$umap@cell.embeddings)
colnames(umapAB1) <- c("UMAP_1", "UMAP_2")

tot <- data.frame(cellnames = colnames(x), sample = x$sample)
tot$clean_cellnames <- gsub('.{2}$', '', tot$cellnames)
tot$final_cellnames <- paste0(tot$sample,":",tot$clean_cellnames,"x")

df <- cbind(tot, umapAB1)

if("final_ident" %in% colnames(x@meta.data)){
    df$celltype <- x$final_ident
}

if("analysis_ident" %in% colnames(x@meta.data)){
    df$celltype <- x$analysis_ident
}

if("seurat_clusters" %in% colnames(x@meta.data)){
    df$cell_seurat_cluster <- x$seurat_clusters
}

if("response" %in% colnames(x@meta.data)){
    df$response <- x$response
}

write.csv(df, file = snakemake@output[[1]])

