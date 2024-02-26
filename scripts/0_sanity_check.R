rm(list = ls())

deg <- read.csv("data/processed/DEG_genes__20231109_155247.csv")
deg_without_p440 <- read.csv("data/processed/sanity_check/DEG_genes__20231116_132326.csv")

common <- semi_join(deg_without_p440, deg, by = ("gene" = "gene"))
144
