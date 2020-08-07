library(dplyr)

nom <- read.table(snakemake@input[[1]], stringsAsFactors = F)

importantgenes <- nom %>% filter(interaction_effect_size < 0 & significant_interaction == TRUE) %>% arrange(desc(interaction_effect_size))

vecneg <- importantgenes$feature

importantgenes <- nom %>% filter(interaction_effect_size > 0 & significant_interaction == TRUE) %>% arrange(desc(interaction_effect_size))

vecpos <- importantgenes$feature

write(vecneg, file = snakemake@output[[1]])
write(vecpos, file = snakemake@output[[2]])
