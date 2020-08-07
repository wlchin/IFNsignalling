

library(Seurat)
library(ggplot2)

x <- readRDS(snakemake@input[[1]])
park <- readRDS(snakemake@input[[2]])

DefaultAssay(park) <- "RNA"
park <- NormalizeData(park)
FeaturePlot(park, features = c("Fxyd2", "Lrp2", "Slc27a2"), ncol = 3)
ggsave("results/markers_park.jpeg", width = 21, height = 7)

DefaultAssay(park) <- "SCT"
tumour_cluster <- c("12", "5", "6", "4", "1", "2", "24", "3", "10", "9", "7")
initvec <- rep("NT", length(levels(park)))
initvec[levels(park) %in% tumour_cluster] <- "T"
names(initvec) <- levels(park)
park <- RenameIdents(park, initvec)
park[["tagged"]] <- Idents(park)

pancreas.anchors <- FindTransferAnchors(park, x, normalization.method = "SCT")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = park$tagged)

x <- AddMetaData(x, metadata = predictions)
saveRDS(x@meta.data, file = snakemake@output[[1]])

head(x@meta.data)

FeaturePlot(x, features = "prediction.score.T")
ggsave("results/prediction.jpeg")



