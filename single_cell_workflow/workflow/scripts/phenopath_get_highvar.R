library(scater)
library(phenopath)
library(Seurat)
library(scde)
library(ggplot2)

x <- readRDS(snakemake@input[[1]])
(Idents(x) <- "seurat_clusters")
DefaultAssay(x) <- "SCT"

x_small <- subset(x, idents = c("12", "3"))
x_small <- FindVariableFeatures(x_small, selection.method = "vst", nfeatures = 3000)
highvar <- VariableFeatures(x_small)

nom <- x_small[["SCT"]]@data
nom <- as.matrix(nom)
nom_c <- nom[highvar,]

pheno <- x_small$response

message("highvargenes in 12 to 3:")
length(highvar)

saveRDS(nom_c, file = snakemake@output[[1]])
saveRDS(pheno, file = snakemake@output[[2]])

