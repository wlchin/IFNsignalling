library(scde)
library(org.Mm.eg.db)

## need to make sure rownames are here

cd <- readRDS(snakemake@input[[1]])

ids <- unlist(lapply(mget(rownames(cd), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids);
names(rids) <- ids 
terms <- as.list(org.Mm.egGO2ALLEGS)

testterm <- c("GO:0035456", "GO:0032608", "GO:0032609", "GO:0034341", "GO:0032606", "GO:0032481")

gos.interest <- testterm # selection of go terms from the list in org.MM
go.env <- lapply(mget(gos.interest, org.Mm.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

saveRDS(go.env, file = snakemake@output[[1]])


x <- readRDS(snakemake@input[[2]])
gocustom <- list(x)
names(gocustom) <- "fast_isg"
y <- list2env(gocustom)
saveRDS(y, file = snakemake@output[[2]])
