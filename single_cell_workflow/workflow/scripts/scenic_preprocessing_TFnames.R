
library(tidyr)

cleanTFnames <- function(df){
    init_vec <- colnames(df)
    substr(init_vec, 1, nchar(init_vec)-3) 
    colnames(df) <- substr(init_vec, 1, nchar(init_vec)-3) 
    return(df)
}

clean_rownames <- function(df){
    rownames(df) <- df[,1]
    df <- df[,-1]
    return(df)
}

x <- read.csv(snakemake@input[[1]])
test <- x  %>% cleanTFnames() %>% clean_rownames()
saveRDS(test, file = snakemake@output[[1]])

