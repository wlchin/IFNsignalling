
x <- readRDS(snakemake@input[[1]])

ind <- length(table(x))

for(i in 1:ind){
    filenm <- snakemake@output[[i]]
    genes <- names(x[x ==i])
    write(genes, filenm) 
}
