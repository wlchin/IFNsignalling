library(GENIE3)
library(stringr)

## this transfers the weights from the genie3 matrix to the direct TF-to-DE matrix.

x_full <- readRDS(snakemake@input[[1]])
x_direct <- readRDS(snakemake@input[[2]])

nom <- str_to_title(x_direct$regulatoryGene)
nom_clean <- gsub("\\(Var\\.2\\)|\\(Var\\.3\\)", "", nom)
x_direct$regulatoryGene <- nom_clean

## remove the ones from the main df
## Spread it out 

ind <- grepl("Ewsr1-Fli1", x_direct$regulatoryGene)
stuff <- x_direct[!ind,]

df1 <- data.frame(regulatoryGene = "Ewsr1", targetGene = x_direct$targetGene[ind])
df2 <- data.frame(regulatoryGene = "Fli1", targetGene = x_direct$targetGene[ind])

df_tot <- rbind(df1, df2)
final <- rbind(stuff, df_tot)
final_no_dup <- final[!duplicated(final),]

regs <- unique(final_no_dup$regulatoryGene)
downstream <- unique(final_no_dup$targetGene)

gene_weights_mat <- x_full[rownames(x_full) %in% regs, colnames(x_full) %in% downstream]

abbr <- getLinkList(gene_weights_mat)

oinks <- dplyr::inner_join(abbr, final_no_dup)

saveRDS(oinks, file = snakemake@output[[1]])

