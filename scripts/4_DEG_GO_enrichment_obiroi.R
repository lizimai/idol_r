rm(list = ls())
library(tidyverse)
library(topGO)
library(ALL)
library(pheatmap)

#### Load GO annotations
geneID2GO_Obir  <- readMappings(file = "data/raw/gene_to_GO_mappings_obiroi.txt") ## original annot

make_named_numeric_vector <- function(list_file_name) {
  full_list <-
    as.list(read.table(list_file_name, sep = ",", head = T))
  full_list_GL <- as.numeric(as.character(full_list$adj.P.Val))
  names(full_list_GL) <- as.character(full_list$X)
  return(full_list_GL)
}

# Function to find the latest file based on the variable name
find_latest_file <-
  function(variable_name, directory = "data/processed") {
    pattern <-
      sprintf("^top_genes_%s_\\d{8}_\\d{6}\\.csv$", variable_name)
    files <-
      list.files(directory, pattern = pattern, full.names = TRUE)
    
    if (length(files) == 0) {
      stop("No files found for the given variable.")
    }
    
    files_df <- data.frame(filename = files,
                           mtime = file.info(files)$mtime)
    
    latest_file <- files_df %>%
      arrange(desc(mtime)) %>%
      dplyr::slice(1) %>%
      pull(filename)
    
    return(latest_file)
  }

#### Read in the significant genes for PC1
#PC1_GL  <- make_named_numeric_vector("data/processed/edgeR/top_genes_PC1_20231016_210926.csv")
#PC1_GL_df  <- read.csv("data/processed/edgeR/top_genes_PC1_20231016_210926.csv")

## set the sig_for_go to the threshold for significant DE genes (This is used for the fishers test, the GSEA doesn't use it)
make_named_numeric_vector <- function(list_file_name) {
  full_list <- as.list(read.table(list_file_name, sep = ",", head = T))
  full_list_GL <- as.numeric(as.character(full_list$adj.P.Val))
  names(full_list_GL) <- as.character(full_list$X)
  return(full_list_GL)
}

find_latest_file <- function(variable_name, directory = "data/processed") {
  pattern <- sprintf("^top_genes_%s_\\d{8}_\\d{6}\\.csv$", variable_name)
  files <- list.files(directory, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) stop("No files found for the given variable.")
  
  files_df <- data.frame(filename = files, mtime = file.info(files)$mtime)
  latest_file <- files_df %>% arrange(desc(mtime)) %>% dplyr::slice(1) %>% pull(filename)
  
  return(latest_file)
}

run_enrichment <- function(genelist, ref, sig_for_GO) {
  topDiffGenes <- function(allScore) { allScore < sig_for_GO }
  
  # BP analysis
  GODATA_BP <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_BP_list = GODATA_BP@graph@nodes
  N_GO_term_use_BP <- length(GODATA_BP@graph@nodes)
  result_GSEA_BP <- runTest(GODATA_BP, algorithm = "elim", statistic = "ks")
  resultFisher_BP <- runTest(GODATA_BP, algorithm = "weight01", statistic = "fisher")
  allRes1_BP <- GenTable(GODATA_BP, GSEA = result_GSEA_BP, Fisher_w01 = resultFisher_BP,
                         ranksOf = "GSEA", topNodes = N_GO_term_use_BP, numChar = 200)
  sig_GSEA_BP_GO <- subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID
  sig_fisher_BP_GO <- subset(allRes1_BP, allRes1_BP$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_BP <- sum(GODATA_BP@feasible == TRUE)
  
  # MF analysis
  GODATA_MF <- new("topGOdata", ontology = "MF", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_MF_list = GODATA_MF@graph@nodes
  N_GO_term_use_MF <- length(GODATA_MF@graph@nodes)
  result_GSEA_MF <- runTest(GODATA_MF, algorithm = "elim", statistic = "ks")
  resultFisher_MF <- runTest(GODATA_MF, algorithm = "weight01", statistic = "fisher")
  allRes1_MF <- GenTable(GODATA_MF, GSEA = result_GSEA_MF, Fisher_w01 = resultFisher_MF,
                         ranksOf = "GSEA", topNodes = N_GO_term_use_MF, numChar = 200)
  sig_GSEA_MF_GO <- subset(allRes1_MF, allRes1_MF$GSEA < sig_for_GO)$GO.ID
  sig_fisher_MF_GO <- subset(allRes1_MF, allRes1_MF$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_MF <- sum(GODATA_MF@feasible == TRUE)
  
  # CC analysis
  GODATA_CC <- new("topGOdata", ontology = "CC", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_CC_list = GODATA_CC@graph@nodes
  N_GO_term_use_CC <- length(GODATA_CC@graph@nodes)
  result_GSEA_CC <- runTest(GODATA_CC, algorithm = "elim", statistic = "ks")
  resultFisher_CC <- runTest(GODATA_CC, algorithm = "weight01", statistic = "fisher")
  allRes1_CC <- GenTable(GODATA_CC, GSEA = result_GSEA_CC, Fisher_w01 = resultFisher_CC,
                         ranksOf = "GSEA", topNodes = N_GO_term_use_CC, numChar = 200)
  sig_GSEA_CC_GO <- subset(allRes1_CC, allRes1_CC$GSEA < sig_for_GO)$GO.ID
  sig_fisher_CC_GO <- subset(allRes1_CC, allRes1_CC$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_CC <- sum(GODATA_CC@feasible == TRUE)
  
  
  allRes1_BPMFCC <- rbind(allRes1_BP, allRes1_MF, allRes1_CC)
  allRes1_BPMFCC$GO_type <- c(rep("BP", length(allRes1_BP[, 1])), 
                              rep("MF", length(allRes1_MF[, 1])), 
                              rep("CC", length(allRes1_CC[, 1])))
  
  out_list <- list("N_GO_term_use_BP" = N_GO_term_use_BP, "GO_term_use_BP_list" = GO_term_use_BP_list, 
                   "allRes1_BP" = allRes1_BP, "sig_GSEA_BP_GO" = sig_GSEA_BP_GO, 
                   "sig_fisher_BP_GO" = sig_fisher_BP_GO, "GODATA_BP" = GODATA_BP, 
                   "N_GO_term_use_MF" = N_GO_term_use_MF, "GO_term_use_MF_list" = GO_term_use_MF_list, 
                   "allRes1_MF" = allRes1_MF, "sig_GSEA_MF_GO" = sig_GSEA_MF_GO, 
                   "sig_fisher_MF_GO" = sig_fisher_MF_GO, "GODATA_MF" = GODATA_MF, 
                   "N_GO_term_use_CC" = N_GO_term_use_CC, "GO_term_use_CC_list" = GO_term_use_CC_list, 
                   "allRes1_CC" = allRes1_CC, "sig_GSEA_CC_GO" = sig_GSEA_CC_GO, 
                   "sig_fisher_CC_GO" = sig_fisher_CC_GO, "GODATA_CC" = GODATA_CC, 
                   "allRes1_BPMFCC" = allRes1_BPMFCC, "N_genes_annot_BP" = N_genes_annot_BP, 
                   "N_genes_annot_MF" = N_genes_annot_MF, "N_genes_annot_CC" = N_genes_annot_CC)
  return(out_list)
}

get_trans_from_GOs <- function(GOLIST, GODATA) {
  out_me <- list()
  for (i in 1:length(GOLIST)) {
    mygenes <- genesInTerm(GODATA, GOLIST[i])
    out_me <- append(out_me, mygenes)
  }
  all_genes <- unique(unlist(out_me))
  return(all_genes)
}

add_GO_to_DE <- function(DE_df, GO_genes_list, GO_col_name) {
  with_go <- ifelse(DE_df[, 1] %in% GO_genes_list, 1, 0)
  DE_df$GO_in <- with_go
  colnames(DE_df) <- c(colnames(DE_df), GO_col_name)
  return(DE_df)
}


# Function to automate the enrichment analysis
automate_enrichment_analysis <- function(variable_name, geneID2GO, immune_terms) {
  gene_list_file <- find_latest_file(variable_name)
  output_folder <- dirname(gene_list_file)
  genelist <- make_named_numeric_vector(gene_list_file)
  
  enrichment_results <- run_enrichment(genelist, geneID2GO, 0.05)
  current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_path <- paste0(output_folder, "/enrich_obiroi_", variable_name, "_", current_time, ".csv")
  
  write.csv(enrichment_results$allRes1_BPMFCC[, c(1, 2, 7, 9)], results_path, row.names = FALSE)
  
  all_immune_genes <- get_trans_from_GOs(immune_terms, enrichment_results$GODATA_BP)
  
  genelist_df <- read.table(gene_list_file, sep = ",", head = TRUE)
  genelist_df <- genelist_df[genelist_df$adj.P.Val < 0.05, ]
  genelist_immune <- genelist_df[genelist_df$X %in% all_immune_genes, ]
  
  if (nrow(genelist_immune) >= 2) {
    genelist_immune_FC <- as.data.frame(genelist_immune$logFC)
    colnames(genelist_immune_FC) <- "logFC"
    rownames(genelist_immune_FC) <- genelist_immune$X
    breaksList <- seq(-4.25, 4.25, by = 0.5)
    plot_dims <- c(max(5, nrow(genelist_immune_FC) / 10), 6)
    
    heatmap_filename <- paste0("output/plots/heatmap_", variable_name, ".png")
    png(heatmap_filename, width = plot_dims[2], height = plot_dims[1], units = "in", pointsize = 12, bg = "white", res = 300)
    
    pheatmap(genelist_immune_FC, clustering_distance_rows = "euclidean", cluster_cols = FALSE, clustering_method = "ward.D2", 
             color = colorRampPalette(c("#08103A", "#08306B", "#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9DCBE1", "#9ECAE1", 
                                        "#FFFFFF", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(length(breaksList)), 
             breaks = breaksList, show_rownames = TRUE, border_color = NA)
    
    dev.off()
    
    genelist_with_immune <- add_GO_to_DE(genelist_df, all_immune_genes, "immune_GO")
    immune_path <- paste0("data/processed/top_genes_with_immune_terms_", variable_name, ".csv")
    write.csv(genelist_with_immune, immune_path, row.names = FALSE)
  }
}


run_enrichment_and_export <- function(variable_name, geneID2GO) {
  gene_list_file <- find_latest_file(variable_name)
  output_folder <- dirname(gene_list_file)
  genelist <- make_named_numeric_vector(gene_list_file)
  
  enrichment_results <- run_enrichment(genelist, geneID2GO, 0.05)
  current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_path <- paste0(output_folder, "/enrich_obiroi_", variable_name, "_", current_time, ".csv")
  
  write.csv(enrichment_results$allRes1_BPMFCC[, c(1, 2, 7, 9)], results_path, row.names = FALSE)
  
  return(list(enrichment_results = enrichment_results, output_folder = output_folder, gene_list_file = gene_list_file))
}

process_immune_genes_and_heatmap <- function(variable_name, immune_terms, enrichment_results, gene_list_file, output_folder) {
  all_immune_genes <- get_trans_from_GOs(immune_terms, enrichment_results$GODATA_BP)
  
  genelist_df <- read.table(gene_list_file, sep = ",", head = TRUE)
  genelist_df <- genelist_df[genelist_df$adj.P.Val < 0.05, ]
  genelist_immune <- genelist_df[genelist_df$X %in% all_immune_genes, ]
  
  if (nrow(genelist_immune) >= 2) {
    create_heatmap(genelist_immune, variable_name, output_folder)
  }
  
  genelist_with_immune <- add_GO_to_DE(genelist_df, all_immune_genes, "immune_GO")
  immune_path <- paste0("data/processed/top_genes_with_immune_terms_", variable_name, ".csv")
  write.csv(genelist_with_immune, immune_path, row.names = FALSE)
}

# Example usage
results_ex <- run_enrichment_and_export("extranidalActivity", geneID2GO_Obir)
results_pl <- run_enrichment_and_export("pathLength", geneID2GO_Obir)
results_sc <- run_enrichment_and_export("spatialCoverage", geneID2GO_Obir)
results_se <- run_enrichment_and_export("spatialEntropies", geneID2GO_Obir)
results_vo <- run_enrichment_and_export("meanVelocityOutside", geneID2GO_Obir)
results_ac <- run_enrichment_and_export("activity", geneID2GO_Obir)
results_vb <- run_enrichment_and_export("velBurstiness", geneID2GO_Obir)
results_tt <- run_enrichment_and_export("tortuosity", geneID2GO_Obir)
results_nd <- run_enrichment_and_export("nodeDegree", geneID2GO_Obir)
results_cb <- run_enrichment_and_export("contactBurstiness", geneID2GO_Obir)
results_pc1 <- run_enrichment_and_export("PC1", geneID2GO_Obir)
results_pc2 <- run_enrichment_and_export("PC2", geneID2GO_Obir)

process_immune_genes_and_heatmap("my_variable", immune_terms, results$enrichment_results, results$gene_list_file, results$output_folder)


# Run the function for a specific variable name
immune_terms = c("GO:0002376")
automate_enrichment_analysis("extranidalActivity", geneID2GO_nrDroso, immune_terms)

# Run GO enrichment only on DEG
all_deg <- read.csv("data/processed/DEG_genes__20231109_155420.csv")
all_deg_gene <- all_deg$gene
