library(igraph)

network <- readRDS(snakemake@input[[1]])
ab1_scaled_counts <- readRDS(snakemake@input[[2]])

ab1_g <- graph_from_data_frame(network[,1:2])
top_ab1_nodes <- head(sort(degree(ab1_g, mode = "all"), decreasing = T), 10)

vertices <- as.integer(unlist(V(ab1_g)))
names(vertices) <- names(V(ab1_g))
top10 <- vertices[names(top_ab1_nodes)]
basegraph <- induced.subgraph(ab1_g, top10)

listvert <- list()

for(i in top10){
  nom <- unlist(ego(ab1_g, order = 1, nodes = i, mode = "out"))
  namevert <- names(V(ab1_g)[i])
  listvert[[namevert]] <- names(nom)
}

saveRDS(listvert, file = snakemake@output[[1]])

allnodes <- do.call(c, listvert)
net.bg <- induced.subgraph(ab1_g, allnodes)

comp <-components(net.bg)
ds <- as.integer(unlist(V(net.bg)))
names(ds) <- names(V(net.bg))
conn <- ds[which(membership(comp) == 1)]

net.bg2 <- net.bg

V(net.bg2)$size <- 7
V(net.bg2)$frame.color <- "black"
V(net.bg2)$color <- "orange"
V(net.bg2)$label <- "" 
E(net.bg2)$arrow.mode <- 0
E(net.bg2)$width <- 1
l <- layout_with_kk(net.bg2)

c_scale <- colorRamp(c('blue','white', 'red'))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

rn <- names(V(net.bg2))
scaledmat <- range01(ab1_scaled_counts)

plot_graph <- function(ind){
  V(net.bg2)$exp <- scaledmat[rn,ind]
  V(net.bg2)$color = apply(c_scale(V(net.bg2)$exp), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
  plot(net.bg2, layout = l)
  box()
}

pdf(snakemake@output[[2]])
par(mfrow=c(2,4), mar=c(0.5,0.5,0.5,0.5))
plot_graph(1)
plot_graph(2)
plot_graph(3)
plot_graph(4)
plot_graph(5)
plot_graph(6)
plot_graph(7)
plot_graph(8)
dev.off()
