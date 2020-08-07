
library(Seurat)
library(cowplot)

x <- readRDS(snakemake@input[[1]]) 
y <- readRDS(snakemake@input[[2]])

forrep <- length(y@assays$RNA["HA-tag",] > 0)
ind <- rep("non_tumour", forrep)
ind[as.logical(y@assays$RNA["HA-tag",] > 0)] <- "tumour"

y[["tagged"]] <- ind
(Idents(y) <- "tagged")

pancreas.anchors <- FindTransferAnchors(y, x)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = y$tagged)
x <- AddMetaData(x, metadata = predictions)

saveRDS(x@meta.data, file = snakemake@output[[1]])

p1 <- FeaturePlot(x, features = "prediction.score.tumour")
p2 <- DimPlot(x)
plot_grid(p1,p2)
ggplot2::ggsave(snakemake@output[[2]], height = 7, width = 14)


