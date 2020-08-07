library(HiveR)
library(igraph)
library(grid)
library(dplyr)
library(RColorBrewer)

source(snakemake@input[[2]]) # load custom hiveR constructor

## Step 1: Annotate high value edges (by quantile) with colour

ab1_resp_filt_1 <- readRDS(snakemake@input[[1]])

threshold <- quantile(ab1_resp_filt_1$weight, 0.95)

for(i in 1:nrow(ab1_resp_filt_1)){
  if(ab1_resp_filt_1$weight[i] > threshold){
      ab1_resp_filt_1$color[i] <- "red"
    } else{
      ab1_resp_filt_1$color[i] <- "grey"
    }
}

ifn_tfs <- c("Stat1","Stat2","Irf9","Irf7","Irf1")
non_ifn_tfs <- c("Klf5","Sp3","Spi1","Rara","Sp1",
                 "Rxra","Tfap2a","Klf16","Sp2","Plag1",
                 "Glis2", "Tcf12", "Pparg") 

# Step 2 create dataframe using modified hiveR constructor

graph_direct_edges <- graph_from_data_frame(ab1_resp_filt_1[,1:2])
nodes <- names(V(graph_direct_edges)) 

axis <- rep(2, length(nodes))
extrafast <- ab1_resp_filt_1$targetGene[ab1_resp_filt_1$regulatoryGene %in% ifn_tfs] 
ind4 <- which(nodes %in% extrafast)
axis[ind4] <- 4
ind1 <- which(nodes %in% non_ifn_tfs)
axis[ind1] <- 3
ind3 <- which(nodes %in% ifn_tfs)
axis[ind3] <- 1

degree <- igraph::degree(graph_direct_edges, V(graph_direct_edges), mode = "all")

nodestuff <- data.frame(nodes = nodes,
                        deg = degree,
                        axis = as.integer(axis),
                        colour = "black", size = 0)

# Step 3 - create the node labels metadatafile
book1 <- data.frame(node.lab = ifn_tfs, 
           node.text = ifn_tfs,
           angle = rep(0,5),
           radius = rep(0,5),
           offset = rep(-0.23,5),
           hjust = rep(0,5),
           vjust = rep(0,5))

write.table(book1, row.names = F, quote = F, file = "metadata_hive.csv", sep = ",")

# Step 4 draw hive plot

hive_alpha <- mod.edge2HPD(edge_df = data.frame(ab1_resp_filt_1[,1:2]), 
             edge.weight = ab1_resp_filt_1$weight * 200, 
             edge.color = ab1_resp_filt_1$color, 
             node.color = nodestuff[,c("nodes", "colour")], 
             node.size = nodestuff[,c("nodes", "size")], 
             node.radius = nodestuff[,c("nodes", "deg")], 
             node.axis = nodestuff[,c("nodes", "axis")])

pdf(snakemake@output[[1]])
plotHive(hive_alpha, method = "ranknorm", np = TRUE, ch = 0.1, 
         anNodes = "metadata_hive.csv", bkgnd = "white",
         anNode.gpar =  gpar(col = "black", fontsize = 14),
         axLabs = c("", "", "", ""),
         axLab.gpar = gpar(col = "orange", fontsize = 14))
dev.off()
