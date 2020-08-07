
library(dplyr)
library(Seurat)
options(future.globals.maxSize = 8000 * 1024^2)

x <- readRDS(snakemake@input[[1]])
cellabels <- unique(x$analysis_ident)

DE_metadata <- paste0(x$analysis_ident, "_", x$response)
x$DEmeta <- DE_metadata

(Idents(x) <- 'DEmeta')

DE_output <- list()

for(i in cellabels){
    responder_type = paste0(i,"_","responder")
    nonresponder_type = paste0(i,"_","nonresponder")
    DElist <- FindMarkers(x, ident.1 = responder_type, ident.2 = nonresponder_type, verbose = FALSE)
    DE_output[[i]] <- DElist
}

saveRDS(DE_output, file = snakemake@output[[1]])

x <- DE_output
celltypes <- names(x)

listsofDE <- list()

for(i in celltypes){
df <- x[[i]]
df$genes <- rownames(df) 
df2 <- df %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_logFC) > 1)
listsofDE[[i]] <- df2$genes
}

comb <- data.frame(lapply(listsofDE, "length<-", max(lengths(listsofDE))), stringsAsFactors = F)
write.csv(comb, row.names = F, quote = F, na = "", file = snakemake@output[[2]])


