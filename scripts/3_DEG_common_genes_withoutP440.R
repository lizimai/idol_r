# load libraries
rm(list=ls())

library(UpSetR)
library(dplyr)
library(tidyr)
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Read the files as you provided
files <- list.files(path = "data/processed/sanity_check", pattern = "^significant_genes_", full.names = TRUE)

# Initialize vectors to store results and a list to store the actual data
behaviour_metrics <- character(length(files))
num_sig_genes <- numeric(length(files))
file_data <- vector("list", length(files))

# Loop through each file and count rows and also store the actual data
for (i in seq_along(files)) {
  data <- read.csv(files[i], header = TRUE)
  
  # Extract behavior metric from the filename
  filename_parts <- unlist(strsplit(basename(files[i]), "_"))
  behaviour_metrics[i] <- filename_parts[3]
  
  num_sig_genes[i] <- nrow(data)
  
  # Store the data for later use in UpSet plot
  file_data[[i]] <- data[[1]]  # Assuming the unnamed first column contains the gene names
}

# Create a binary matrix indicating the presence/absence of genes in each file
all_genes <- unique(unlist(file_data))
matrix_data <- matrix(0, nrow=length(all_genes), ncol=length(file_data))
rownames(matrix_data) <- all_genes
colnames(matrix_data) <- behaviour_metrics

for(i in seq_along(file_data)) {
  matrix_data[file_data[[i]], i] <- 1
}

# Generate the UpSet plot
desired_order <- rev(c("extranidalActivity", "pathLength", "spatialCoverage", "spatialEntropies", "meanVelocityOutside", "activity", "velBurstiness", "tortuosity", "nodeDegree", "contactBurstiness", "PC1", "PC2")) # replace with your desired order
upset_plot <- upset(as.data.frame(matrix_data), sets = desired_order, order.by = "degree", keep.order=TRUE, matrix.color="#D55E00", sets.bar.color="#009E73", text.scale = 1.3) 
pdf(file=paste0("output/plots/sanity_check/upset_plot_",current_time,".pdf"), width = 8, height = 6) # or other device
upset_plot
dev.off()


# Summary ----
# Convert the matrix into a long format dataframe
long_data <- as.data.frame(matrix_data) %>%
  rownames_to_column("gene") %>%
  gather(key = "behaviour_metric", value = "presence", -gene)

# Filter to keep only rows where gene is present
overlapping_genes <- long_data %>%
  filter(presence == 1) %>%
  group_by(gene) %>%
  summarize(behaviours = list(behaviour_metric), 
            n_behaviours = n())

# Specify desired metrics for overlap
obiroi_annotation_heiko <- read.table("data/raw/Ooceraea-biroi_rna-from-genomic_omicsbox_Generic-Annotation_with_gene_id.txt", header = T)
obiroi_annotation_heiko$gene_id <- as.character(obiroi_annotation_heiko$gene_id)
id_name <- unique(read.table("data/raw/gene2exon.txt"))

all_metrics_with_DEG <- c("extranidalActivity", "pathLength", "spatialCoverage", "spatialEntropies", "activity", "nodeDegree", "PC1")
empty_metrics <- c()
functional_annotation <- function(overlapping_genes, obiroi_annotation_heiko, id_name, specific_metrics){
  # Filter overlapping_genes to ensure all specific metrics are present in behaviours
  specific_overlap_genes <- overlapping_genes %>%
    filter(sapply(behaviours, function(x) all(specific_metrics %in% x)))
  
  # Convert the behaviours list column to a comma-separated string column
  specific_overlap_genes$behaviours <- sapply(specific_overlap_genes$behaviours, paste, collapse=", ")
  
  # Find gene annotation
  gene_list <- left_join(specific_overlap_genes, id_name, by = join_by(gene == V2)) %>% 
    rename(gene_id = V1)
  
  gene_with_annotation <- left_join(gene_list, obiroi_annotation_heiko, by = join_by(gene_id)) %>% 
    group_by(gene_id) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    select(gene, behaviours, Sequence.Description, InterPro.GO.ID, InterPro.GO.Term, InterPro.GO.Category)
  
  current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(gene_with_annotation, paste0("data/processed/sanity_check/DEG_genes_", paste0(specific_metrics, collapse = '_'), "_", current_time,".csv"))
}

functional_annotation(overlapping_genes, obiroi_annotation_heiko, id_name, all_metrics_with_DEG)
functional_annotation(overlapping_genes, obiroi_annotation_heiko, id_name, empty_metrics)
