rm(list = ls())
library(tidyverse)
library(topGO)

all_deg <- read.csv("data/processed/DEG_genes__20231109_155420.csv")
geneID2GO  <- readMappings(file = "data/raw/gene_to_GO_mappings_obiroi.txt")
geneID2GO_old<- readMappings(file = "data/raw/Obir.OGS5.2.longest.pepforTopGO.txt")

all_deg_gene <- all_deg$gene

deg_list <- tibble(gene = names(geneID2GO)) %>% 
  mutate(value = ifelse(gene %in% all_deg_gene, 1, 0))

deg_list_old <- tibble(gene = names(geneID2GO_old)) %>% 
  mutate(value = ifelse(gene %in% all_deg_gene, 1, 0))

geneList <- deg_list$value
names(geneList) <- deg_list$gene

geneList_old <- deg_list_old$value
names(geneList_old) <- deg_list_old$gene

topDiffGenes <- function(allScore) { allScore == 1 }

run_enrichment <- function(genelist, ref, sig_for_GO) {
  topDiffGenes <- function(allScore) { allScore > sig_for_GO }
  
  # BP analysis
  GODATA_BP <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_BP_list = GODATA_BP@graph@nodes
  N_GO_term_use_BP <- length(GODATA_BP@graph@nodes)
  resultFisher_BP <- runTest(GODATA_BP, algorithm = "weight01", statistic = "fisher")
  allRes1_BP <- GenTable(GODATA_BP, Fisher_w01 = resultFisher_BP,
                         ranksOf = "Fisher_w01", topNodes = N_GO_term_use_BP, numChar = 200)
  sig_GSEA_BP_GO <- subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID
  sig_fisher_BP_GO <- subset(allRes1_BP, allRes1_BP$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_BP <- sum(GODATA_BP@feasible == TRUE)
  
  # MF analysis
  GODATA_MF <- new("topGOdata", ontology = "MF", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_MF_list = GODATA_MF@graph@nodes
  N_GO_term_use_MF <- length(GODATA_MF@graph@nodes)
  resultFisher_MF <- runTest(GODATA_MF, algorithm = "weight01", statistic = "fisher")
  allRes1_MF <- GenTable(GODATA_MF, Fisher_w01 = resultFisher_MF,
                         ranksOf = "Fisher_w01", topNodes = N_GO_term_use_MF, numChar = 200)
  sig_GSEA_MF_GO <- subset(allRes1_MF, allRes1_MF$GSEA < sig_for_GO)$GO.ID
  sig_fisher_MF_GO <- subset(allRes1_MF, allRes1_MF$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_MF <- sum(GODATA_MF@feasible == TRUE)
  
  # CC analysis
  GODATA_CC <- new("topGOdata", ontology = "CC", allGenes = genelist, geneSel = topDiffGenes, 
                   annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  GO_term_use_CC_list = GODATA_CC@graph@nodes
  N_GO_term_use_CC <- length(GODATA_CC@graph@nodes)
  resultFisher_CC <- runTest(GODATA_CC, algorithm = "weight01", statistic = "fisher")
  allRes1_CC <- GenTable(GODATA_CC, Fisher_w01 = resultFisher_CC,
                         ranksOf = "Fisher_w01", topNodes = N_GO_term_use_CC, numChar = 200)
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

############## ############## ############## ############## ############## ############## ############## ##############
############## extracting expression data for GO terms
### Function to get transcript IDs annotated with particular GO terms
get_trans_from_GOs = function(GOLIST, GODATA) {
  #mygenes <- genesInTerm(GODATA, GOLIST)
  out_me <- list()
  
  for (i in 1:length(GOLIST))
  {
    #print(GOLIST[i])
    mygenes <- genesInTerm(GODATA, GOLIST[i])
    #print(mygenes)
    mygenes_GOs <- append(GOLIST[i], mygenes)
    #print(mygenes_GOs)
    out_me = append(out_me, mygenes)
  }
  
  #print(out_me)
  
  all_genes <- c()
  for (n in seq(1, length(out_me))) {
    print(n)
    all_genes <- c(all_genes, out_me[[n]])
  }
  
  all_genes_u <- unique(all_genes)
  return(all_genes_u)
}

###################################
### add immune genes column to DE data (or whatever your list is)
add_GO_to_DE <- function(DE_df, GO_genes_list, GO_col_name) {
  with_go <- c()
  for (i in 1:length(DE_df[, 1])) {
    gene_name <- DE_df[i, 1]
    gene_in <- ifelse(gene_name %in% GO_genes_list, 1, 0)
    with_go <- c(with_go, gene_in)
  }
  
  out_header <- c(colnames(DE_df), GO_col_name)
  
  DE_df$GO_in <- with_go
  colnames(DE_df) <- out_header
  #print(head(DE_df))
  
  return(DE_df)
  
}

run_enrichment_and_export <- function(genelist, geneID2GO, output_folder) {
  enrichment_results <- run_enrichment(genelist, geneID2GO, 0.05)
  current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_path <- paste0(output_folder, "/enrich_obiroi_all_DEG_obiroi", "_", current_time, ".csv")
  
  write.csv(enrichment_results$allRes1_BPMFCC[, c(1, 2, 6, 7)], results_path, row.names = FALSE)
  
  return(list(enrichment_results = enrichment_results, output_folder = output_folder))
}

enrichment_results <- run_enrichment_and_export(geneList, geneID2GO, "data/processed")





