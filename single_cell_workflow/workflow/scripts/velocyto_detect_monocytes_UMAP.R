
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

#comb <- readRDS("../../results/test_combined_seurat_objects/renca_and_ab1_combined.rds")
#ab1 <- readRDS("../../results/test_combined_seurat_objects/ab1_combined.rds")

comb <- readRDS(snakemake@input[[2]])
ab1 <- readRDS(snakemake@input[[1]])

combmeta <- comb@meta.data
ab1meta <- ab1@meta.data
ab1meta$CellID <- paste0(ab1meta$sample, "_", rownames(ab1meta))
combmeta$CellID <- paste0(combmeta$sample, "_", rownames(combmeta))
combmeta$loom_cellid <- rownames(combmeta)

combmeta$big_clustering <- combmeta$seurat_clusters
filtered_ab1_match_clustering <- inner_join(ab1meta, combmeta, by = "CellID")
monos <- filtered_ab1_match_clustering[filtered_ab1_match_clustering$seurat_clusters.x == 12,]

tab <- table(monos$big_clustering)
highmono <- names(which(tab == max(tab)))
ly6mono <- combmeta[combmeta$big_clustering == highmono,]

write.csv(ly6mono, file = snakemake@output[[1]])



DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb)

p1 <- FeaturePlot(comb, features = c("Ly6c1"), order = F)
p2 <- DimPlot(comb, label = T) + NoLegend()
p3 <- DimPlot(comb, cells.highlight = rownames(ly6mono))

plot_grid(p1, p2, p3, ncol = 3)
ggplot2::ggsave(snakemake@output[[2]], height = 7, width = 21)
