library(dplyr)
library(pheatmap)

phenodata <- readRDS(snakemake@input[[1]])
Rcounts <- readRDS(snakemake@input[[2]])
NRcounts <- readRDS(snakemake@input[[3]])

## combine both responder and non responder matrices
## for loop to extract groups + timepoints
## average counts placed into list
## merge entries from list into a (average) count matrix

average_counts <- list()
combined_count_matrix <- cbind(Rcounts, NRcounts)

for(i in  c("RS", "NR")){
    for(j in c("d0", "d2", "d4", "d6")){
        samps <- dplyr::filter(phenodata, Timepoint == j, Group == i)
        phenolist <- samps$sample
        mat <- rowMeans(combined_count_matrix[,phenolist])
        pheno_label <- paste0(i,"_",j)
        average_counts[[pheno_label]] <- mat
    }
}

matrix_meancounts <- do.call(cbind, average_counts)

## scaling function

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

scaled <-scale_rows(matrix_meancounts)
scaled_with_NAs_removed <- scaled[complete.cases(scaled),]

saveRDS(scaled_with_NAs_removed, file = snakemake@output[[1]]) 
individual_counts <- list()
combined_count_matrix <- cbind(Rcounts, NRcounts)

for(i in  c("RS", "NR")){
    for(j in c("d0", "d2", "d4", "d6")){
        samps <- dplyr::filter(phenodata, Timepoint == j, Group == i)
        phenolist <- samps$sample
        mat <- combined_count_matrix[,phenolist]
        pheno_label <- paste0(i,"_",j)
        individual_counts[[pheno_label]] <- mat
    }
}

matrix_individualcounts <- do.call(cbind, individual_counts)
scaled <-scale_rows(matrix_individualcounts)
scaled_with_NAs_removed <- scaled[complete.cases(scaled),]
saveRDS(scaled_with_NAs_removed, file = snakemake@output[[2]]) #write out the rds for downstre

