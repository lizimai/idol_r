# Load Libraries ----
rm(list=ls())
library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel) 
library(UpSetR)

current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Data Importation and Preprocessing ----
# Import counting table
counts <- read.delim("data/raw/salmon.merged.gene_counts.tsv", check.names = F, stringsAsFactors = F)
colnames(counts)<-gsub("X","",colnames(counts))
# Remove gene with duplicate ids
# Create a table of counts for each gene_name
gene_counts <- table(counts$gene_name)
# Get gene names that appear only once
unique_genes <- names(gene_counts[gene_counts == 1])
# Subset the data frame to keep only rows with those unique gene names
counts <- counts[counts$gene_name %in% unique_genes, ]
# Remove redundant gene ID column
counts <- counts[, -1]
# Remove the suffix of the column names
rownames(counts) <- counts[,1]
# Remove P4-40
counts <- counts[, -c(70:78)]

# Import ant Table
behaviour <- read.csv("../idol_matlab/data/processed/updatedAntTable_20231023_162828.csv")
# Remove any underscores in the colony names
behaviour$colony <- gsub("-", "", behaviour$colony)
# Concatenate the colony and antColour columns to create the new antID
behaviour$antID <- paste(behaviour$colony, behaviour$antColour, sep = "_")
behaviour <- behaviour[, -c(1,2,4)]
row.names(behaviour) <- behaviour$antID
behaviour <- behaviour[order(behaviour$antID), ]
behaviour <- behaviour[, c("antID", setdiff(names(behaviour), "antID"))]
# center and scale parameters
scaled_vector <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# Apply the function to your data
behaviour[, -c(1,2)] <- lapply(behaviour[, -c(1,2)], scaled_vector)
# match antIDs to the ones in counts
group <- behaviour[behaviour$antID %in% colnames(counts),]


# Filtering: Remove Low Count Reads ----
# Convert raw counts to CPM
cpm_values <- cpm(DGEList(counts = counts))
# Define a threshold for lowly expressed genes; for example, genes with CPM < 1 in more than half the samples will be removed
lowly_expressed_genes <- rowSums(cpm_values < 1) > (ncol(counts)/2)
# Filter out the lowly expressed genes from your original raw counts
filtered_counts <- counts[!lowly_expressed_genes, ]

# Creating DGEList Object for Analysis ----
dge <- DGEList(counts = filtered_counts[,2:69], genes = filtered_counts[,1], group=group$PC1)
# Recompute lib sizes
dge$samples$lib.size <- colSums(dge$counts)
# Normalisation
dge <- calcNormFactors(dge)

# Calculate distance matrices ----
# Normalize the data and compute Counts Per Million (CPM)
cpm_values <- cpm(dge, log = TRUE, prior.count = 5)

# Export Euclidean Distance Matrix for Mantel test (Step 2 code)
expr_dist_matrix <- as.matrix(dist(t(cpm_values), method = "euclidean"))
behav_dist_matrix <- as.matrix(dist(group[, -c(1,2,13,14)], method = "euclidean"))
save(expr_dist_matrix, behav_dist_matrix, file = paste0("data/processed/sanity_check/dist_matrix_expr_behav_", current_time,".rds"))


# Define Differential Expression Analysis Function ----
diff_exp_analysis <- function(counts, group, design_col) {
  
  # Set up the design matrix using the column provided
  design <- model.matrix(as.formula(paste0("~ ", design_col)), data = group)
  
  # Estimate the correlation for the random effect (colony)
  correlation <- duplicateCorrelation(counts, design, block=group$colony)
  
  # Use the correlation estimate in the linear modeling
  v <- voom(counts, design, plot=TRUE, block=group$colony, correlation=correlation$consensus)
  
  # Fit the model
  fit <- lmFit(v, design, block=group$colony, correlation=correlation$consensus)
  fit <- eBayes(fit)
  
  # Extract top genes
  top_genes <- topTable(fit, coef=2, n=Inf, sort.by="p")
  sig_genes <- top_genes[top_genes$adj.P.Val < 0.05, ]
  
  # Visualize the results
  num_sig_genes <- nrow(subset(top_genes, adj.P.Val < 0.05))
  plot <- ggplot(top_genes, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(alpha=0.4) +  # Plotting all genes
    geom_point(data = subset(top_genes, adj.P.Val < 0.05), aes(x=logFC, y=-log10(adj.P.Val)), color="red", alpha=0.6) +  # Highlighting significant genes in red
    theme_minimal() +
    labs(title=paste0("Volcano Plot", ", Design metric: ", design_col),
         subtitle = paste("Number of genes significantly differentially expressed =", num_sig_genes),
         x="Log2 Fold Change", y="-Log10 adjusted P-value") +
    geom_hline(yintercept=-log10(0.05), color="red", linetype="dashed") + # FDR threshold line
    annotate("text", x = Inf, y = -log10(0.05) + 0.5, label = "FDR = 0.05", vjust = 3, hjust = 1)  # Text label for the FDR threshold
  
  # Save the plot and significant genes
  current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  pdf(paste0("output/plots/sanity_check/volcano_plot_", design_col, "_", current_time, ".pdf"))
  print(plot)
  dev.off()
  write.csv(sig_genes, paste0("data/processed/sanity_check/significant_genes_", design_col, "_", current_time, ".csv"), row.names=TRUE)
  write.csv(top_genes, paste0("data/processed/sanity_check/top_genes_", design_col, "_", current_time, ".csv"), row.names=TRUE)
  
  # Return the plot and significant genes for further exploration if needed
  list(plot = plot, sig_genes = sig_genes)
}

# Run Differential Expression Analysis ----
# Extract the column names after "colony" in the group dataframe
cols_to_analyze <- names(group)[which(names(group) == "colony"):length(names(group))][-1]
# Apply the function for each column
for (col in cols_to_analyze) {
  diff_exp_analysis(counts = dge$counts, group = group, design_col = col)
}

# Summary for all behaviour metrics ----
# Get a list of all CSV files that start with 'significant_genes_'
files <- list.files(path = "data/processed/sanity_check/", pattern = "^significant_genes_", full.names = TRUE)

# Initialize vectors to store results
behaviour_metrics <- character(length(files))
num_sig_genes <- numeric(length(files))

# Loop through each file and count rows
for (i in seq_along(files)) {
  data <- read.csv(files[i], header = TRUE)
  
  # Extract behavior metric from the filename
  filename_parts <- unlist(strsplit(basename(files[i]), "_"))
  behaviour_metrics[i] <- filename_parts[3]
  
  num_sig_genes[i] <- nrow(data)
}

# Combine the results into a dataframe
df_summary <- data.frame(Behaviour_Metrics = behaviour_metrics, Num_Significant_Genes = num_sig_genes)
write.table(df_summary, sep="\t", row.names=FALSE, quote=FALSE)




