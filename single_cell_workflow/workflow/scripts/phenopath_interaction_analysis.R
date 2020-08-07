
library(phenopath)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scater)


fit <- readRDS(snakemake@input[[1]])
ints <- interactions(fit)

interaction_terms <- ints %>%
    dplyr::filter(significant_interaction == TRUE) %>%
    dplyr::arrange(desc(interaction_effect_size)) %>%
    select(c("feature", "interaction_effect_size"))

interaction_terms <- ints %>%
    dplyr::filter(significant_interaction == TRUE) 

write.table(interaction_terms, file = snakemake@output[[1]], row.names = T, quote = F)

# plot the top10 genes
data_abbr1 <- dplyr::arrange(interaction_terms, interaction_effect_size)[1:20,]
data_abbr2 <- dplyr::arrange(interaction_terms, desc(interaction_effect_size))[1:20,]
data_abbr <- rbind(data_abbr1, data_abbr2)

chi_cutoff <- sort(ints$chi)[10]

ggplot(ints, aes(x = interaction_effect_size, y = 1 / chi, 
                 color = significant_interaction)) +
  geom_point() +
  geom_text_repel(data = data_abbr, 
                  aes(label = feature)) +
  scale_colour_brewer(palette = "Set1")

ggsave(snakemake@output[[2]])


ggplot(ints, aes(x = pathway_loading, y = interaction_effect_size, 
                 color = significant_interaction)) +
  geom_point() +
  geom_text_repel(data = data_abbr, 
                  aes(label = feature), size = 5) +
  scale_colour_brewer(palette = "Set1")  +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

ggsave(snakemake@output[[3]])

trajectory_output <- trajectory(fit)

saveRDS(trajectory_output, file = snakemake@output[[4]])

