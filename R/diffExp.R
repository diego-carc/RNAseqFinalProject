#Name: diffExp.R
#Author: Diego Carmona Campos
#Last update:  04/02/2024
#Description:  Generates the SummarizedExperiment object for SRP165875 record used in this project

# Load libraries
library(limma)
library(ggplot2)
library(sessioninfo)
library("pheatmap")

# Load data
rse_gene_SRP165875_norm <- readRDS("../data/rse_gene_SRP165875_norm.RDS")

# Check age bias
ggplot(as.data.frame(colData(rse_gene_SRP165875_norm)), aes(y = assigned_gene_prop, x=sra_attribute.age, group =sra_attribute.age, color=factor(sra_attribute.age))) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  scale_x_discrete(breaks = c(6,12,18)) +
  scale_color_discrete(name="Mouse Age") +
  ylab("Assigned Gene Prop") +
  xlab("Age Group (months)")


# Check sex bias
ggplot(as.data.frame(colData(rse_gene_SRP165875_norm)), aes(y = assigned_gene_prop, x = sra_attribute.Sex, group = sra_attribute.Sex, color=factor(sra_attribute.Sex))) +
  geom_boxplot() +
  scale_color_discrete(name="Mouse Sex") +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Sex Group")

# Crate model
mod <- model.matrix(~ sra_attribute.age + sra_attribute.Sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP165875_norm))

# Model variance
vGene <- voom(assays(rse_gene_SRP165875_norm)$normcounts, mod, plot = TRUE)

# Compute model coefficients
eb_results <- eBayes(lmFit(vGene))

#
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP165875_norm),
  sort.by = "none"
)

# MA plot
limma::plotMA(eb_results, coef = 2)

# Volcano plot
positions <- match(rownames(eb_results$lods), rownames(rowData(rse_gene_SRP165875_norm)))
volcanoplot(eb_results, coef = 2, highlight = 3, names=rowData(rse_gene_SRP165875_norm)$gene_name[positions])


# Heatmap
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 30, ]
positions <- match(rownames(exprs_heatmap),  rownames(rowData(rse_gene_SRP165875_norm)))
rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP165875_norm)$gene_name[positions]
df <- as.data.frame(colData(rse_gene_SRP165875_norm)[, c("sra_attribute.age", "sra_attribute.Sex")])
colnames(df) <- c("AgeGroup", "Sex")

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)

# Save session info
writeLines(capture.output(session_info()), "../doc/diffExp_session_info.txt")

