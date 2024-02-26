library(edgeR)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)

# Import counting table
counts <- read.delim("data/raw/nfcore/salmon.merged.gene_counts.tsv", check.names = F, stringsAsFactors = F)
# Remove redundant gene ID column
counts <- counts[, -1]
# Remove the suffix of the column names
colnames(counts)<-gsub("X","",colnames(counts))

#Read sample description
exp_desc <- read.csv("data/processed/DEseq2/meta.csv", sep = ",")
exp_desc_arr <- exp_desc %>% 
  arrange(X)
colony <- relevel(as.factor(exp_desc_arr$colony), ref = "18")

# creating the DGEList data class
y <- DGEList(counts = counts[,2:78], genes = counts[,1], group=colony)

# filtering: remove low count reads
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
# Recompute lib sizes
y$samples$lib.size <- colSums(y$counts)
# Normalisation
# Normalize
y <- calcNormFactors(y)

y <- estimateDisp(y, robust = T)
pdf("exp1_plotBCV.pdf", width = 10, height = 10)
plotBCV(y)
dev.off()

# plot MDS
p <- plotMDS(y)
d <- data.frame(p$x, p$y)
row.names(d) <- colnames(counts[,-1])
d <- left_join(rownames_to_column(d), exp_desc_arr, by=c("rowname" = "X"))
row.names(d) <- colnames(counts[,-1])

pdf("exp1_plotMDS.pdf", width = 10, height = 10)
ggplot(d, aes(p.x, p.y, label=rownames(d), color=exp_desc_arr$colony))+
  geom_point(size=2)+
  coord_cartesian(clip="off")+
  scale_shape_manual(values=c(23,24))+
  #scale_color_manual(values=c("blue", "red"))+
  geom_text_repel(xlim=c(-Inf,Inf), ylim=c(-Inf, Inf), box.padding = 0.3, min.segment.length = 0, seed = 42, show.legend = F)+
  xlab("Leading logFC dim 1 (64%)")+
  ylab("Leading logFC dim 2 (18%)")
dev.off()

