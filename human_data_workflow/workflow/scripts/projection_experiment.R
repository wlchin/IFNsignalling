
library(Seurat)
library(dplyr)

ab1 <- readRDS(snakemake@input[[1]])
pbmc <- readRDS(snakemake@input[[2]])

Idents(ab1) <- "seurat_clusters"

initvec <- rep("other", length(levels(ab1)))
initvec[levels(ab1) %in% "12"] <- "Mono1"
initvec[levels(ab1) %in% "3"] <- "Mono2"
initvec[levels(ab1) %in% "4"] <- "Mono3"

names(initvec) <- levels(ab1)

ab1 <- RenameIdents(ab1, initvec)
ab1[["tagged"]] <- Idents(ab1)

anchors <- FindTransferAnchors(ab1, pbmc, normalization.method = "SCT")

gc()

#save predictions and seurat object with added metadata

predictions <- TransferData(anchorset = anchors, refdata = ab1$tagged)
saveRDS(predictions, file = snakemake@output[[1]])
pbmc <- AddMetaData(pbmc, metadata = predictions)
saveRDS(pbmc@meta.data, file = snakemake@output[[3]])
