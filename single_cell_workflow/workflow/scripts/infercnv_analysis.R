
library(infercnv)
library(Seurat)

x <- readRDS(snakemake@input[[1]])

genes_ <- read.table("resources/mm10.data")
genes <- data.frame(genes_[,2:4], row.names = genes_[,1])
annots <- data.frame(label = paste0(Idents(x),"_", x$seurat_clusters), row.names = colnames(x))

annots <- data.frame(label = x$analysis_ident, row.names = colnames(x))

data = x@assays$RNA@counts[,colnames(x)]

testvec <- as.character(unique(annots$label))[-1]

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data, 
                                               gene_order_file=genes,
                                               annotations_file=annots,
                                               ref_group_names=testvec)

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir=snakemake@output[[1]], 
                              cluster_by_groups=TRUE, 
                              denoise=TRUE,
                              HMM=FALSE,
                              num_threads=12,
                              no_plot=FALSE)

saveRDS(infercnv_obj, file = snakemake@output[[2]])
