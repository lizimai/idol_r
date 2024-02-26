rm(list=ls())

# Loading libraries ---- 
# Optionally, you can visualize the distance matrices
library(ggplot2)
library(reshape2)
library(vegan)

# Load distance matrices ----
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
load(file = "data/processed/dist_matrix_expr_behav_20231109_151441.rds")
# Check dimensions to ensure they're square matrices
dim(expr_dist_matrix)
dim(behav_dist_matrix)

# Visualisation ----
# For gene expression
expr_melted <- melt(as.matrix(expr_dist_matrix))
ggplot(expr_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "ant ID", y = "ant ID", title = "Eucledian distance between individual gene expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("output/plots/expr_dist_matrix_", current_time, ".pdf"), width = 12, height = 9)

# For behavioral data
behav_melted <- melt(as.matrix(behav_dist_matrix))
ggplot(behav_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "ant ID", y = "ant ID", title = "Eucledian distance between individual behaviour matrics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("output/plots/behav_dist_matrix_", current_time, ".pdf"), width = 12, height = 9)

mantel_result = mantel(behav_dist_matrix, expr_dist_matrix, method="pearson", permutations=999)
print(mantel_result)

# do the same test without P4-40
expr_dist_matrix_sub <- expr_dist_matrix[-c(68:77), -c(68:77)] 
behav_dist_matrix_sub <- behav_dist_matrix[-c(68:77), -c(68:77)] 
mantel_result_sub = mantel(behav_dist_matrix_sub, expr_dist_matrix_sub, method="pearson", permutations=999)
print(mantel_result_sub)

