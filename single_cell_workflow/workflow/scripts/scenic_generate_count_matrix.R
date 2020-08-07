

library(Seurat)

AB1 <- readRDS(snakemake@input[[1]])

#AB1 <- readRDS("../data/tumour_and_single_r_labelled_phenotypes.seurat.rds")

(Idents(AB1) <- "sample")

## this will go but I need to ensure the pipeline runs
#AB1 <- subset(AB1, idents = c("NKD180900302","NKD180900308")) # test with smaller

(Idents(AB1) <- "analysis_ident")
#AB1 <- subset(AB1, idents = "Monocytes")
DefaultAssay(AB1) <- "SCT"

#responder_seurat <- subset(AB1, ident = "responder")
#vect <- AB1$res_clus 

dat <- data.frame(CellID = colnames(AB1),
                  cluster_res = paste0(AB1$seurat_clusters, "_", AB1$response))

write.table(dat, sep = ",", row.names = F, quote = F, file = snakemake@output[[2]])

all_mat <- as.matrix(GetAssayData(object = AB1, slot = "counts"))

## subset delete later
#someTFs <- c("Irf1", "Irf7")
#ind <- c(rownames(all_mat)[1:1000], someTFs)
#all_mat <- all_mat[ind,] # selects only the first 1000 genes for speed - delete later

str(all_mat)

nom <- data.frame(cellid = rownames(all_mat), all_mat)
write.table(nom, row.names = F, sep = "\t", quote = F, file = snakemake@output[[1]])


