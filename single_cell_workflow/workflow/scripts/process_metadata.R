
library(Seurat)

seurat_object <- readRDS(snakemake@input[[1]])

x <- readRDS(snakemake@input[[3]])
y <- readRDS(snakemake@input[[2]])

seurat_object[["singler"]] <- y 

tumour_info <- table(x$predicted.id, x$seurat_cluster)
NT <- tumour_info[1,]
T <- tumour_info[2,]
ratios <- T/(NT + T)

tumourident <- rep("NT", length(ratios))
tumourident[ratios > 0.1] <- "T"

names(tumourident) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, tumourident)

seurat_object[["tumourident"]] <- Idents(seurat_object) 

tumour_ident <- seurat_object$tumourident
final_ident <- seurat_object$singler

for(i in 1:length(y)){
    if(tumour_ident[i] == "T"){
        final_ident[i] <- "Tumour"
    }}

seurat_object[["analysis_ident"]] <- final_ident

vecres <- seurat_object$sample
ind <- length(vecres)

for(i in 1:ind){
    test <- vecres[i]
    if(grepl("301|305|309|302|307|311", test)){
        vecres[i] <- "responder"
        }else{ vecres[i] <- "nonresponder"}
}

seurat_object[["response"]] <- vecres
saveRDS(seurat_object, file = snakemake@output[[1]])

Idents(seurat_object) <- "analysis_ident"
DimPlot(seurat_object, label = T) + NoLegend()
ggplot2::ggsave(snakemake@output[[2]])

