library(ggplot2)
library(dplyr)
library(data.table)

kmers <- read.delim("out_10/kmers/joined/all_5mers.tsv")

kmers$pfal_norm <- kmers$pfal / sum(kmers$pfal) * 10^6
kmers$tgon_norm <- kmers$tgon / sum(kmers$tgon) * 10^6
kmers$log_pfal <- log10(kmers$pfal)
kmers$log_tgon <- log10(kmers$tgon)

ggplot(kmers) +
  aes(x = pfal, y = tgon) +
  geom_point() +
  scale_x_log10() + 
  scale_y_log10()

ggplot(kmers) +
  aes(x = pfal_norm, y = tgon_norm) +
  geom_point() +
  scale_x_log10() + 
  scale_y_log10()

cor(log10(kmers$pfal), log10(kmers$tgon), method = "spearman")

cor.test(kmers$pfal, kmers$tgon, method = "spearman", exact = FALSE)

ggplot(kmers) +
  aes(x = pfal) +
  geom_histogram(aes(y = after_stat(density)), bins = 200) +
  stat_function(fun = dnorm, args = list(mean = mean(kmers$pfal), sd = sd(kmers$pfal)))

ggplot(kmers) +
  aes(x = log_pfal) +
  geom_histogram(aes(y = after_stat(density)), bins = 200) +
  stat_function(fun = dnorm, args = list(mean = mean(kmers$log_pfal), sd = sd(kmers$log_pfal)))

shapiro.test(kmers$log_pfal)

kmers_sel <- kmers %>% select(-c("kmer"))

kmers_T <- transpose(data.table(kmers_sel))

rownames(kmers_T) <- colnames(kmers_sel)
colnames(kmers_T) <- kmers$kmer

pca <- prcomp(na.omit(kmers_T), scale = T)
summary(pca)

ggplot(pca$x) +
  aes(x = PC1, y = PC2) +
  geom_point()
