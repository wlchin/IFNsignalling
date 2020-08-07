
library(Seurat)
library(SingleR)
library(cowplot)

x <- readRDS(snakemake@input[[1]]) 
ref_mouse <- readRDS("resources/databases/mouse.rds")

DefaultAssay(x) <- "SCT" 
cluster_ids <- x$seurat_clusters  
raw_counts <- x@assays$SCT@data
cluster_ids <- as.character(cluster_ids)

singler_mouse = SingleR(method = "cluster",
                        test = raw_counts,
                        ref = ref_mouse$data,
                        labels = ref_mouse$main_types,
                        clusters = cluster_ids)

x[["cluster_input"]] <- as.character(cluster_ids)

Idents(x) <- "cluster_input"
clus <- singler_mouse@rownames
labs <- singler_mouse$labels
names(labs) <- clus
x <- RenameIdents(x, labs)

singlerlabels <- as.character(Idents(x))
saveRDS(singlerlabels, file = snakemake@output[[1]])

p1 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Idents(x) <- "seurat_clusters"
p2 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot_grid(p1, p2)
ggplot2::ggsave(snakemake@output[[2]], height = 7, width = 14)






