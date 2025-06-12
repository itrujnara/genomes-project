library(dplyr)
library(tidyr)
library(vroom)
library(data.table)
library(ggplot2)
library(plotly)
library(umap)

kmers <- vroom("../getkmers_dummy/remote/out/kmers/joined/all_5mers.tsv")
meta_sample <- vroom("~/Downloads/p10k_sample_metadata.csv")
meta_genome <- vroom("../getkmers_dummy/p10k_genome_metadata.csv")

metadata <- meta_genome %>% left_join(meta_sample, by = c("SampleID" = "Sample ID"))

kmers_sel <- kmers %>% select(-c("kmer"))

kmers_T <- kmers_sel %>% apply(2, function(row) {
  row / sum(row) * 10^6
}) %>% data.frame() %>% transpose()

# kmers_T <- kmers_sel %>% transpose()

colnames(kmers_T) <- kmers$kmer
rownames(kmers_T) <- colnames(kmers_sel)

kmers_T$id <- rownames(kmers_T)
kmers_T$assembly <- sapply(strsplit(kmers_T$id, "_"), tail, n = 1) 

kmers_clean <- kmers_T %>% select(-c("id", "assembly"))

pca <- prcomp(na.omit(kmers_clean), center = TRUE, scale = TRUE)

summary(pca)

pctable <- bind_cols(na.omit(kmers_T)$assembly, pca$x)
colnames(pctable)[1] <- "assembly"
pctable <- pctable %>% left_join(metadata, by = c("assembly" = "AssemblyID"))

pca_plot <- plot_ly()

pca_plot <- pca_plot %>%
  add_trace(
    type = "scatter",
    mode = "markers",
    x = pctable$PC1,
    y = pctable$PC2,
    color = pctable$AnnotLevel,
    text = pctable$Species.x,
    hovertemplate = paste('<b>%{text}</b>',
                          '<br>X: %{x}<br>',
                          'Y: %{y}')
  )

# pca_plot <- pca_plot %>% layout(xaxis = list(type="log"), yaxis = list(type="log"))

pca_plot

umap_basic <- umap(na.omit(kmers_clean))

umap_plot <- plot_ly()

umap_plot <- umap_plot %>% 
  add_trace(
    type = "scatter",
    mode = "markers",
    x = umap_basic$layout[,1],
    y = umap_basic$layout[,2],
    text = pctable$Species.x,
    hovertemplate = paste('<b>%{text}</b>',
                          '<br>X: %{x}<br>',
                          'Y: %{y}')
  )

umap_plot
