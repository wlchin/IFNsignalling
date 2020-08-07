
library(dplyr)
library(pheatmap)

message("loading binarised data")

binarised_auc <- readRDS(snakemake@input[[1]])

cellmetadata <- read.csv(snakemake@input[[2]],
                         colClasses=c("cluster_res"="character"))

cleanrownames <- gsub("\\.", "-", rownames(binarised_auc))
binarised_df_per_cell <- data.frame(binarised_auc, CellID = cleanrownames, stringsAsFactors = F) 

combined_df <- dplyr::inner_join(cellmetadata, binarised_df_per_cell)

message("calculate the aggregate")

TF_average_matrix <- combined_df %>%
    group_by(cluster_res) %>%
    summarise_at(colnames(binarised_auc), mean, na.rm = T)

nonres <- TF_average_matrix[grepl("_nonresponder$", TF_average_matrix$cluster_res),]
res <- TF_average_matrix[grepl("_responder$", TF_average_matrix$cluster_res),]

delta_matrix <- as.matrix(res[,-1] - nonres[,-1])
rownames(delta_matrix) <- gsub("_responder$", "", res$cluster_res)

ifn_reg <- c("Irf7", "Stat2", "Irf1", "Stat1", "Irf9")

ind <- rownames(delta_matrix) %in% c("12","3","4")

ind_gene <- colnames(delta_matrix) %in% ifn_reg

mattoplot_small <- delta_matrix[ind,ind_gene]

mattoplot_big <- delta_matrix[ind,]

rownames(mattoplot_small) <- c("Monocytes cluster 1",
                               "Monocytes cluster 2",
                               "Monocytes cluster 3")

rownames(mattoplot_big) <- c("Monocytes cluster 1",
                             "Monocytes cluster 2",
                             "Monocytes cluster 3")


pdf(snakemake@output[[1]], height = 7, width = 7)
pheatmap(mattoplot_small,
         cellwidth = 20,
         cellheight = 20,
         fontsize = 10,
         treeheight_row = 10,
         treeheight_column = 10,
         cluster_rows = F)
dev.off()


pdf(snakemake@output[[2]], height = 7, width = 30)
pheatmap(mattoplot_big,
         cellwidth = 5,
         cellheight = 5,
         fontsize = 4,
         treeheight_row = 10,
         treeheight_column = 10,
         cluster_rows = F)
dev.off()


