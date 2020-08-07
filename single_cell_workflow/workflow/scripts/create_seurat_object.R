
library(Seurat)

data_dir <- snakemake@params$fpath

list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx

expression_matrix <- Read10X(data.dir = data_dir)

s_obj <- CreateSeuratObject(counts = expression_matrix, 
                            project = "scell", 
                            min.cells = 3, 
                            min.features = 200)

#str(s_obj)

s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^mt-")

s_obj <- SCTransform(s_obj, vars.to.regress = "percent.mt")

s_obj[["sample"]] <- snakemake@params$folder 

saveRDS(s_obj, file = snakemake@output[[1]])
