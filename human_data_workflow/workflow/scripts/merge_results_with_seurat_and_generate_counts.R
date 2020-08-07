
library(scde)
library(Seurat)

## load results and seurat object for metadata merge

x <- readRDS(snakemake@input[[1]])
df_counts <- readRDS(snakemake@input[[2]])
predictions <- readRDS(snakemake@input[[3]])
human <- readRDS(snakemake@input[[4]])

metadata <- read.csv(snakemake@input[[5]], stringsAsFactors = F)

samvec <- x$orig.ident
indvec <- metadata$timepoint
names(indvec) <- metadata$sample_id

indvecres <- metadata$response
names(indvecres) <- metadata$sample_id

timepointvec <- indvec[samvec]
responsevec <- indvecres[samvec]

x[["predictions"]] <- predictions$predicted.id
x[["timepoint"]] <- timepointvec
x[["response"]] <- responsevec
x[["human_labels"]] <- human

## visusalise predictions

monocytes <- names(human)[grepl("Monocyte", human)]
monos_only <- subset(x, cells = monocytes)

## save metadata and counts for later

saveRDS(x@meta.data, file = snakemake@output[[1]])
write.table(x@meta.data, row.names = F, sep = "\t", quote = F, file = snakemake@output[[2]])

all_mat <- df_counts[,monocytes]
nom <- data.frame(cellid = rownames(all_mat), all_mat)
saveRDS(nom, file = snakemake@output[[3]])
write.table(nom, row.names = F, sep = "\t", quote = F, file = snakemake@output[[4]])
