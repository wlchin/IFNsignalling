library(scde)
library(Seurat)

x <- readRDS(snakemake@input[[1]])

(Idents(x) <- "sample")

cd <- clean.counts(as.matrix(x@assays$SCT@counts))
cd <- apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x})

groups <- x$analysis_ident[colnames(cd)]

saveRDS(cd, snakemake@output[[1]])
saveRDS(groups, snakemake@output[[2]])
