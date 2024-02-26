rm(list=ls())
library(dplyr)
library(tidyr)

DEG_all <- read.csv("data/processed/edgeR/DEG_genes__20231023_164959.csv", row.names = "X")
immune_gene <- read.csv("data/raw/Obir_immunity_share/immune_gene_ids.csv", row.names = "X")

DEG_with_immune_annotation <- left_join(DEG_all, immune_gene, by = join_by("gene" == "id_gcf")) %>% 
  drop_na(id_gca) 

write.csv(DEG_with_immune_annotation, file = "data/processed/edgeR/DEG_immune_genes.csv")