rm(list = ls())
library(dplyr)

obiroi_annotation_heiko <- read.table("data/raw/Ooceraea-biroi_rna-from-genomic_omicsbox_Generic-Annotation_with_gene_id.txt", header = T)
obiroi_annotation_heiko$gene_id <- as.character(obiroi_annotation_heiko$gene_id)
id_name <- unique(read.table("data/raw/gene2exon.txt"))

obiroi_annotation_heiko_gene_names <- left_join(id_name, obiroi_annotation_heiko, by = c("V1" = "gene_id")) %>% 
  dplyr::select(V2) %>% 
  distinct()

foo <- id_name %>% 
  dplyr::select(V2) %>% 
  distinct()

# all the gene ID: V1 (16395 without LOC),  gene names V2 (15511) from the genome can be found in obiroi_annotation_heiko, therefore we should be able to use it
obiroi_GO <- left_join(id_name, obiroi_annotation_heiko, by = c("V1" = "gene_id")) %>% 
  dplyr::select(V2, Annotation.GO.ID) %>% 
  unique() %>% 
  drop_na()

# Convert the list to a format suitable for saving
split_go_ids <- strsplit(obiroi_GO$Annotation.GO.ID, split = ";")
gene_to_go <- setNames(split_go_ids, obiroi_GO$V2)

mapping_lines <- sapply(names(gene_to_go), function(x) {
  paste(x, paste(gene_to_go[[x]], collapse = "\t"), sep = "\t")
}, USE.NAMES = FALSE)

# Write to a file
writeLines(mapping_lines, "data/raw/gene_to_GO_mappings_obiroi.txt")

