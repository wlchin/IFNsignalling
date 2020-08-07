
source('workflow/scripts/CellRouter_Class.R')
library(igraph)

libdir <- paste0(getwd(), "/resources/CellRouter")
load(file = snakemake@input[[1]])

cellrouter <- findClusters(cellrouter, method="graph.clustering", num.pcs=15, k=20)
cellrouter <- buildKNN(cellrouter, k = 20, column.ann = 'clustering', num.pcs = 20, sim.type = 'jaccard')

sources <- c('cluster_12')
targets <- c('cluster_3')
methods <- c("euclidean", "maximum", "manhattan","canberra","binary")
filename <- "cell_edge_weighted_network.txt"

write.table(cellrouter@graph$edges, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE) #input network

cellrouter <- findPaths(cellrouter, column='clustering', libdir, paste0(getwd(),"/"), method="graph")
ranks <- c('path_cost', 'path_flow', 'rank', 'length')

cellrouter <- processTrajectories(cellrouter, rownames(cellrouter@ndata), path.rank=ranks[3], num.cells = 3, neighs = 3,column.ann = 'population', column.color = 'colors')

names <- unique(names(cellrouter@pathsinfo$distr))
cellrouter <- correlationPseudotime(cellrouter, type='spearman')
cellrouter <- topGenes(cellrouter, 0.8, 0.1)
cellrouter <- smoothDynamics(cellrouter, names)

cellrouter <- clusterGenesPseudotime(cellrouter, 5)

save(cellrouter, file = snakemake@output[[1]])
