rm(list = ls())
# Load necessary library
library(dplyr)
library(stringr)
library(Biostrings)

# Read the BLAST results into a data frame
blast_results <- read.table("data/raw/gene_blast/blast_for2GCF.txt", header = FALSE, sep = "\t")
index <- read.table("data/raw/refGenome/gene2transcript.txt")

# Name the columns for easier referencing
colnames(blast_results) <- c("query_id", "subject_id", "percent_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

# Group by query_id and filter for the best hit (highest bit_score) for each query
best_hits <- blast_results %>%
  group_by(query_id) %>%
  filter(bit_score == max(bit_score, na.rm = TRUE)) %>% 
  distinct(query_id, .keep_all = TRUE)

# Extract gene IDs
best_hits_ids <- best_hits %>%
  select(query_id, subject_id) %>% 
  unique()

# Read the genome fasta files
ref_genome <- readDNAStringSet("data/raw/obiroi_genome/cds_from_genomic_GCF.fna")

# Extract unique IDs from best_hits_ids
reference_ids <- unique(best_hits_ids$subject_id)

# Define a function to extract protein and gene IDs from the DNAStringSet object's names.
# Define a function to extract protein and gene IDs from the DNAStringSet object's names.
extract_ids <- function(ids, genome, pattern) {
  sapply(ids, function(id) {
    header <- names(genome)[which(startsWith(names(genome), id))]
    if (length(header) > 0) {
      match <- str_match(header, pattern)
      if(!is.na(match[1,2])) {  # Check if a match was found
        return(match[1,2])     # Return the captured group
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
}

# Define patterns to extract protein and gene IDs
reference_pattern<- "gene=([^\\]]+)"

# Extract IDs
ref_protein_ids <- extract_ids(reference_ids, ref_genome, query_pattern)

# Create data frames for the extracted IDs
reference_id_df <- data.frame(subject_id = reference_ids, reference_gene_id = ref_protein_ids)

# Combine best_hits_ids with the extracted IDs using left_join
result <- best_hits_ids %>%
  ungroup() %>% 
  left_join(reference_id_df, by = "subject_id") %>% 
  unique() 
