
library(Seurat)
library(SingleR)
library(cowplot)


x <- readRDS(snakemake@input[[1]]) 

load(snakemake@input[[3]])

DefaultAssay(x) <- "SCT" 
cluster_ids <- x$seurat_clusters
raw_counts_init <- readRDS(snakemake@input[[2]])
raw_counts <- raw_counts_init[,colnames(x)]


singler_mouse = SingleR(method = "cluster",
                        sc_data = raw_counts,
                        ref_data = hpca$data,
                        types = hpca$types,
                        clusters = cluster_ids,
  genes = "de", quantile.use = 0.8, p.threshold = 0.05,
  fine.tune = TRUE, fine.tune.thres = 0.05, sd.thres = 1,
  do.pvals = T, numCores = 2)

saveRDS(singler_mouse, file = snakemake@output[[3]])

(Idents(x) <- "seurat_clusters")
new.cluster.ids <- singler_mouse$labels1
names(new.cluster.ids) <- levels(x)

message("mapping")

x <- RenameIdents(x, new.cluster.ids)
singlerlabels <- as.character(Idents(x))
x[["human"]] <- singlerlabels

saveRDS(singlerlabels, file = snakemake@output[[1]])

message("plotting")

p1 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Idents(x) <- "seurat_clusters"
p2 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot_grid(p1, p2)
ggplot2::ggsave(snakemake@output[[2]], height = 7, width = 14)

message("saving singler")

saveRDS(x, file = snakemake@output[[4]])
