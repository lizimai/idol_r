test <- read.csv("data/processed/significant_genes_extranidalActivity_20231109_093105.csv")

library(ggplot2)

ggplot(test, aes(x = logFC)) +
  geom_histogram()

test2 <- read.csv("data/processed/significant_genes_PC1_20231109_093208.csv")
ggplot(test2, aes(x = logFC)) +
  geom_histogram()
