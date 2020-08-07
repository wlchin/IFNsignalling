library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
       
)

col_fun = colorRamp2(c(0, 0.002, 0.004, 0.006, 0.008), rev(brewer.pal(5,"YlGnBu")))
col_fun2 = colorRamp2(c(-2, -1, 0, 1, 2), rev(brewer.pal(5,"RdBu")))

# Annotations
renca_ranks <- readRDS(snakemake@input[[1]])
ab1_ranks <- readRDS(snakemake@input[[2]])
ab1_counts <- readRDS(snakemake@input[[3]])
renca_counts <- readRDS(snakemake@input[[4]])

#process genes according to rank for rows

genelist_ab1 <- names(head(ab1_ranks, 1000))
genelist_renca <- names(head(renca_ranks, 1000))

#genelist_renca <- head(dplyr::arrange(renca_ranks, desc(V1)), 851)$regulatoryGene
genelist <- intersect(genelist_ab1, genelist_renca)

#filter from whole matrix
abr_mat_ab1 <- ab1_counts[rownames(ab1_counts) %in% genelist[1:99],]
abr_mat_renca <- renca_counts[rownames(renca_counts) %in% genelist[1:99],]

timevec = c(rep("Day0", 12),rep("Day2",8),rep("Day4",8),rep("Day6",8),rep("Day0", 12),rep("Day2",8),rep("Day4",8),rep("Day6",8)) 
phenovec = c(rep("responder",36),rep("non-responder",36))

hatest1 = HeatmapAnnotation(
  Timepoint = timevec,
  col = list(Timepoint = c("Day0" = "#DCDCDC","Day2" = "#A9A9A9","Day4" = "#778899","Day6" = "#000000")), 
  name = "hatest1",
  simple_anno_size = unit(0.2, "cm"), show_annotation_name = F)

hatest2 = HeatmapAnnotation(
  Phenotype = phenovec,
  col = list(
             Phenotype = c("responder" = "#1a9641", "non-responder" = "#FFA500")
  ), name = "hatest2",
  simple_anno_size = unit(0.2, "cm"), show_annotation_name = F
)

expression_theme <- function(mat, labs, wide, tall, name, show_legend = T){
  d <- Heatmap(mat, cluster_columns = F, width = unit(wide, "cm"), 
               height = unit(tall, "cm"), show_row_names = F, column_names_gp = gpar(fontsize = 9), show_column_names = F,
               show_row_dend = F, col = col_fun2, heatmap_legend_param = list(title = "Gene Expression"), row_dend_side = "right",
               show_heatmap_legend = show_legend, column_title = labs, column_title_gp = gpar(fontsize = 12, fontface = "bold"),
               top_annotation = hatest1, bottom_annotation = hatest2, name = name)
  return(d)
}

ht1 <- expression_theme(abr_mat_ab1, "AB1", 4.0, 5, "ht1", T)
ht2 <- expression_theme(abr_mat_renca, "Renca", 4.0, 5, "ht2", T)
htlist <- ht1+ht2

pdf(snakemake@output[[1]])
draw(htlist, merge_legend = T)
dev.off()
