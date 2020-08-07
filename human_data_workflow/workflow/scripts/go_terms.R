library(scde)
library(org.Hs.eg.db)

cd <- readRDS(snakemake@input[[1]])

ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids)
names(rids) <- ids 
terms <- as.list(org.Hs.egGO2ALLEGS)

testterm <- c("GO:0035456", "GO:0032608", "GO:0032609", "GO:0034341", "GO:0032606", "GO:0032481")

gos.interest <- testterm # selection of go terms from the list in org.MM
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

saveRDS(go.env, file = snakemake@output[[1]])

fast_isg_list <- readRDS(snakemake@input[[2]])
gocustom <- list(toupper(fast_isg_list))
names(gocustom) <- "fast_isg"
fast_isg_env <- list2env(gocustom)
saveRDS(fast_isg_env, file = snakemake@output[[2]])
