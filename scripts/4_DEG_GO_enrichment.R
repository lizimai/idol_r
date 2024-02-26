rm(list = ls())
library(tidyverse)
library(topGO)
library(ALL)
library(pheatmap)

# library("VennDiagram")
# library(gridExtra)
# library(grid)
# library(ggplot2)
# library("SuperExactTest")
# library(cowplot)
# print (sessionInfo())
#### Load GO annotations

geneID2GO_Obir    <- readMappings(file = "../iDOLDataAnalysis/data/raw/GOAnnotation/Obir.OGS5.2.longest.pepforTopGO.txt") ## original annot
#geneID2GO_nrArth  <- readMappings(file = "data/raw/gcf_003672135_1_obir_v5_4_rna_from_genomic_longest_iso_to_nr_arth_with_interpro.annot_fortopgo_gene.txt") ## annotated to nr-Arthropod
geneID2GO_nrDroso <- readMappings(file = "data/raw/gcf_003672135_1_obir_v5_4_rna_from_genomic_longest_iso_to_nr_droso_with_interpro.annot_fortopgo_gene.txt") ## annotated to nr-Drosophila

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
run_enrichment <- function(genelist, ref, sig_for_GO) {
  ### make rule for classing sig / non-sig - note this rule is not used for the GSEA
  
  topDiffGenes <- function(allScore) {
    return(allScore < sig_for_GO)
  }
  # topDiffGenes <- function(allScore) {return(allScore < 1)} ## as a check - setting to one gives the same pvalues for the GSEA
  
  #### make GOdata object
  #### setting node size as 10 so at least 10 genes must be annot per GO terms
  #### do enrichment test
  
  ## BP
  GODATA_BP = new(
    "topGOdata",
    ontology = "BP",
    allGenes = genelist,
    geneSel = topDiffGenes,
    annot = annFUN.gene2GO,
    gene2GO = ref,
    nodeSize = 10
  )
  GO_term_use_BP_list = GODATA_BP@graph@nodes
  N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
  result_GSEA_BP     <-
    runTest(GODATA_BP, algorithm = "elim", statistic = "ks")
  resultFisher <-
    runTest(GODATA_BP, algorithm = "weight01", statistic = "fisher")
  allRes1_BP <-
    GenTable(
      GODATA_BP,
      GSEA = result_GSEA_BP,
      Fisher_w01 = resultFisher,
      ranksOf = "GSEA",
      topNodes = length(GODATA_BP@graph@nodes),
      numChar = 200
    )
  sig_GSEA_BP_GO     = subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID
  sig_fisher_BP_GO     = subset(allRes1_BP, allRes1_BP$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_BP = (sum(GODATA_BP@feasible == TRUE))
  
  ## MF
  GODATA_MF = new(
    "topGOdata",
    ontology = "MF",
    allGenes = genelist,
    geneSel = topDiffGenes,
    annot = annFUN.gene2GO,
    gene2GO = ref,
    nodeSize = 10
  )
  GO_term_use_MF_list = GODATA_MF@graph@nodes
  N_GO_term_use_MF = length(GODATA_MF@graph@nodes)
  result_GSEA_MF     <-
    runTest(GODATA_MF, algorithm = "elim", statistic = "ks")
  resultFisher <-
    runTest(GODATA_MF, algorithm = "weight01", statistic = "fisher")
  allRes1_MF <-
    GenTable(
      GODATA_MF,
      GSEA = result_GSEA_MF,
      Fisher_w01 = resultFisher,
      ranksOf = "GSEA",
      topNodes = length(GODATA_MF@graph@nodes),
      numChar = 200
    )
  sig_GSEA_MF_GO     = subset(allRes1_MF, allRes1_MF$GSEA < sig_for_GO)$GO.ID
  sig_fisher_MF_GO     = subset(allRes1_MF, allRes1_MF$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_MF = (sum(GODATA_MF@feasible == TRUE))
  
  ## CC
  GODATA_CC = new(
    "topGOdata",
    ontology = "CC",
    allGenes = genelist,
    geneSel = topDiffGenes,
    annot = annFUN.gene2GO,
    gene2GO = ref,
    nodeSize = 10
  )
  GO_term_use_CC_list = GODATA_CC@graph@nodes
  N_GO_term_use_CC = length(GODATA_CC@graph@nodes)
  result_GSEA_CC     <-
    runTest(GODATA_CC, algorithm = "elim", statistic = "ks")
  resultFisher <-
    runTest(GODATA_CC, algorithm = "weight01", statistic = "fisher")
  allRes1_CC <-
    GenTable(
      GODATA_CC,
      GSEA = result_GSEA_CC,
      Fisher_w01 = resultFisher,
      ranksOf = "GSEA",
      topNodes = length(GODATA_CC@graph@nodes),
      numChar = 200
    )
  sig_GSEA_CC_GO     = subset(allRes1_CC, allRes1_CC$GSEA < sig_for_GO)$GO.ID
  sig_fisher_CC_GO     = subset(allRes1_CC, allRes1_CC$Fisher_w01 < sig_for_GO)$GO.ID
  N_genes_annot_CC = (sum(GODATA_CC@feasible == TRUE))
  
  ## cat together
  
  allRes1_BPMFCC <- rbind(allRes1_BP, allRes1_MF, allRes1_CC)
  allRes1_BPMFCC$GO_type <-
    c(rep("BP", length(allRes1_BP[, 1])), rep("MF", length(allRes1_MF[, 1])), rep("CC", length(allRes1_CC[, 1])))
  1
  ## return everything!
  out_list = list(
    "N_GO_term_use_BP" = N_GO_term_use_BP,
    "GO_term_use_BP_list" = GO_term_use_BP_list,
    "allRes1_BP" = allRes1_BP,
    "sig_GSEA_BP_GO" = sig_GSEA_BP_GO,
    "sig_fisher_BP_GO" = sig_fisher_BP_GO,
    "GODATA_BP" = GODATA_BP,
    "N_GO_term_use_MF" = N_GO_term_use_MF,
    "GO_term_use_MF_list" = GO_term_use_MF_list,
    "allRes1_MF" = allRes1_MF,
    "sig_GSEA_MF_GO" = sig_GSEA_MF_GO,
    "sig_fisher_MF_GO" = sig_fisher_MF_GO,
    "GODATA_MF" = GODATA_MF,
    "N_GO_term_use_CC" = N_GO_term_use_CC,
    "GO_term_use_CC_list" = GO_term_use_CC_list,
    "allRes1_CC" = allRes1_CC,
    "sig_GSEA_CC_GO" = sig_GSEA_CC_GO,
    "sig_fisher_CC_GO" = sig_fisher_CC_GO,
    "GODATA_CC" = GODATA_CC,
    "allRes1_BPMFCC" = 	allRes1_BPMFCC,
    "N_genes_annot_BP" = N_genes_annot_BP,
    "N_genes_annot_MF" = N_genes_annot_MF,
    "N_genes_annot_CC" = N_genes_annot_CC
  )
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

# Function to automate the enrichment analysis
automate_enrichment_analysis <-
  function(variable_name, geneID2GO, immune_terms) {
    # Find the latest file for the variable
    gene_list_file <- find_latest_file(variable_name)
    output_folder <- dirname(gene_list_file)
    
    # Read in the significant genes
    genelist <- make_named_numeric_vector(gene_list_file)
    
    # Run enrichment analysis
    enrichment_results <- run_enrichment(genelist, geneID2GO, 0.05)
    # Export enrichment results
    current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(
      enrichment_results$allRes1_BPMFCC[, c(1, 2, 7, 9)],
      paste0(
        output_folder,
        "/enrich_nrDroso_",
        variable_name,
        "_",
        current_time,
        ".csv"
      ),
      row.names = FALSE
    )
    
    all_immune_genes <-
      get_trans_from_GOs(immune_terms, enrichment_results$GODATA_BP)
    #######################################################################################
    ### draw heatmap
    genelist_df <- read.table(gene_list_file, sep = ",", head = T)
    genelist_df <- genelist_df[genelist_df$adj.P.Val < 0.05, ]
    
    # Check if the number of rows in genelist_df is greater than or equal to 2
    genelist_immune    <-
      genelist_df[genelist_df$X %in% all_immune_genes, ]
    genelist_immune_FC <- as.data.frame(genelist_immune$logFC)
    if (nrow(genelist_immune_FC) >= 2) {
      colnames(genelist_immune_FC) <- "logFC"
      rownames(genelist_immune_FC) <- genelist_immune$X
      
      breaksList = c(
        -4.25,
        -3.75,
        -3.25,
        -2.75,
        -2.25,
        -1.75,
        -1.25,
        -0.75,
        -0.25,
        0.25,
        0.75,
        1.25,
        1.75,
        2.25,
        2.75,
        3.25,
        3.75,
        4.25
      )
      
      # Calculate plot dimensions
      num_rows <- nrow(genelist_immune_FC)
      plot_height <-
        max(5, num_rows / 10)  # Set a minimum height and adjust based on number of rows
      plot_width <- 6  # Fixed width, but you can adjust this as well
      
      # Save the plot
      png(
        filename = paste0("output/plots/heatmap_", variable_name, ".png"),
        width = plot_width,
        height = plot_height,
        units = "in",
        pointsize = 12,
        bg = "white",
        res = 300
      )
      
      pheatmap(
        genelist_immune_FC,
        clustering_distance_rows = "euclidean",
        cluster_cols = FALSE,
        clustering_method = "ward.D2",
        color = colorRampPalette(
          c(
            "#08103A",
            "#08306B",
            "#08417C",
            "#08519C",
            "#2171B5",
            "#4292C6",
            "#6BAED6" ,
            "#9DCBE1",
            "#9ECAE1",
            "#FFFFFF",
            "#FFFFFF"  ,
            "#FCBBA1",
            "#FCBBA1" ,
            "#FC9272" ,
            "#FB6A4A" ,
            "#EF3B2C" ,
            "#CB181D" ,
            "#A50F15" ,
            "#67000D",
            "#67000D",
            "#67000D"
          )
        )(length(breaksList)),
        breaks = breaksList,
        show_rownames = T,
        border_color = NA
      )
      
      dev.off()
      
      genelist_with_immune <-
        add_GO_to_DE(genelist_df, all_immune_genes, "immune_GO")
      # Save the updated genelist
      write.csv(
        genelist_with_immune,
        paste0(
          "data/processed/top_genes_with_immune_terms_",
          variable_name,
          ".csv"
        ),
        row.names = F
      )
    }
  }


# Run the function for a specific variable name
immune_terms = c("GO:0002376")
automate_enrichment_analysis("extranidalActivity", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("pathLength", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("spatialCoverage", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("spatialEntropies", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("meanVelocityOutside", geneID2GO_nrDroso, immune_terms) # no DEG
automate_enrichment_analysis("activity", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("velBurstiness", geneID2GO_nrDroso, immune_terms) # no DEG
automate_enrichment_analysis("tortuosity", geneID2GO_nrDroso, immune_terms) # no DEG
automate_enrichment_analysis("nodeDegree", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("contactBurstiness", geneID2GO_nrDroso, immune_terms) # no DEG
automate_enrichment_analysis("PC1", geneID2GO_nrDroso, immune_terms)
automate_enrichment_analysis("PC2", geneID2GO_nrDroso, immune_terms) # no DEG