
#library(data.table)
library(ggplot2)

cutoff <- 0.0003125

## Filter AB1 direct network using cutoff

x <- readRDS(snakemake@input[[1]])
filt_x <- x[x$weight > cutoff,]

saveRDS(filt_x, file = snakemake@output[[2]])

## Filter AB1 direct network using cutoff

y <- readRDS(snakemake@input[[2]])
filt_y <- y[y$weight > cutoff,]

saveRDS(filt_y, file = snakemake@output[[3]])

## Filter AB1 direct network using cutoff

combined_edgelist <- rbind(x, y)

ggplot(combined_edgelist) +
  geom_density(aes(weight), adjust = 3) + ylim(0,1000) + theme_classic() + xlab("GENIE3 scores") +
    geom_vline(xintercept = cutoff, color = "red") +
     theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))
ggsave(snakemake@output[[1]])
