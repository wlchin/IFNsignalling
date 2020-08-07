

library(Seurat)

expr <- read.table(snakemake@input[[1]])

park <- CreateSeuratObject(counts = expr,
                   project = "Park",
                   min.cells = 3,
                   min.features = 200)

park[["percent.mt"]] <- PercentageFeatureSet(park, pattern = "^mt-")

park <- SCTransform(park, vars.to.regress = "percent.mt", verbose = FALSE)

park <- RunPCA(park, verbose = FALSE)
park <- RunUMAP(park, dims = 1:30, verbose = FALSE)

park <- FindNeighbors(park, dims = 1:30, verbose = FALSE)
park <- FindClusters(park, verbose = FALSE)

saveRDS(park, file = snakemake@output[[1]])

