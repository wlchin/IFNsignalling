library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

ab1_ISG <- readRDS(snakemake@input[[1]])
ab1_counts <- readRDS(snakemake@input[[2]])

renca_ISG <- readRDS(snakemake@input[[3]])
renca_counts <- readRDS(snakemake@input[[4]])

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
       
)

expression_theme <- function(mat, labs, wide, name){
  d <- Heatmap(mat, cluster_columns = F, width = unit(wide*1.3, "cm"), name = name,
               height = unit(10, "cm"), show_row_names = F, column_names_gp = gpar(fontsize = 8), show_column_names = T,
               show_row_dend = F, col = col_fun2, heatmap_legend_param = list(title = "Gene Expression"), row_dend_side = "right", km = 2,
               show_heatmap_legend = TRUE, column_title = labs, column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  return(d)
}

weightmat_theme <- function(mat, labs, wide, name){
  b <- Heatmap(mat, 
               na_col = "blue", width = unit(wide*1.5, "cm"), 
               height = unit(10, "cm"), cluster_columns = F, column_title = labs, 
               column_names_gp = gpar(fontsize = 6), name = name,
               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
               show_column_dend = FALSE, show_column_names = T, show_row_names = F,
               show_row_dend = FALSE, col = col_fun, heatmap_legend_param = list(title = "GENIE3 Scores"))
  return(b)
}

col_fun = colorRamp2(c(0, 0.002, 0.004, 0.006, 0.008), rev(brewer.pal(5,"YlGnBu")))
col_fun2 = colorRamp2(c(-2, -1, 0, 1, 2), rev(brewer.pal(5,"RdBu")))

allness <- unique(c(colnames(ab1_ISG), colnames(renca_ISG)))

## filtering the matrices

common_column_index <- intersect(colnames(ab1_ISG), colnames(renca_ISG))

ab1_expr <- ab1_counts[common_column_index,1:4]
renca_expr <- renca_counts[common_column_index,1:4]

colnames(ab1_expr) <- c("Day 0", "Day 2", "Day 4", "Day 6")
colnames(renca_expr) <- c("Day 0", "Day 2", "Day 4", "Day 6")


ab1_edges_all <- t(ab1_ISG[,common_column_index])
renca_edges_all <- t(renca_ISG[,common_column_index])

top10ab1 <- names(sort(colSums(ab1_edges_all), decreasing = T)[1:10])
top10renca <- names(sort(colSums(renca_edges_all), decreasing = T)[1:10])

## extract top ten for simpler visualisation - full heatmap in supplementary

ht1 <- weightmat_theme(ab1_edges_all[,top10ab1], "AB1", 1.5, "ht1")
ht2 <- expression_theme(ab1_expr, "nom", 2.0, "ht2")

ht3 <- weightmat_theme(renca_edges_all[,top10renca], "Renca", 1.5, "ht3")
ht4 <- expression_theme(renca_expr, "nom", 2.0, "ht4")

htlist <- ht1+ht2+ht3+ht4

str(row_order(ht2))

length(rownames(ab1_expr)[row_order(ht2)[[1]]])
length(rownames(ab1_expr)[row_order(ht2)[[2]]])

## output pdf

pdf(snakemake@output[[1]])
draw(htlist, merge_legend = T,
     ht_gap = unit(c(0.2,0.8,0.2,0.8), "cm"), main_heatmap = 2, 
     row_sub_title_side = "left", annotation_legend_side = "bottom")
dev.off()
