#Name: dataPreProcessing.R
#Author: Diego Carmona Campos
#Last update:  04/02/2024
#Description:  Transforms the data computing read counts, applying normalization and
# performing filter by expression.

# Load libraries
library(recount3)
library(edgeR)
library(sessioninfo)

# Read data
rse_gene_SRP165875 <- readRDS("../rawData/rse_gene_SRP165875.Rds")

# Get read counts from base-pair counts
assay(rse_gene_SRP165875, "counts") <- compute_read_counts(rse_gene_SRP165875)
rse_gene_SRP165875$assigned_gene_prop <- rse_gene_SRP165875$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP165875$recount_qc.gene_fc_count_all.total

# Expand sra attributes
rse_gene_SRP165875 <- expand_sra_attributes(rse_gene_SRP165875)
rse_gene_SRP165875$sra_attribute.age <- as.factor(rse_gene_SRP165875$sra_attribute.age)
rse_gene_SRP165875$sra_attribute.Sex <- as.factor(rse_gene_SRP165875$sra_attribute.Sex)
rse_gene_SRP165875$sra_attribute.source_name <- as.factor(rse_gene_SRP165875$sra_attribute.source_name)
rse_gene_SRP165875$sra_attribute.generation  <- as.factor(rse_gene_SRP165875$sra_attribute.generation)
rse_gene_SRP165875$sra_attribute.tissue  <- as.factor(rse_gene_SRP165875$sra_attribute.tissue)

# Apply normalization log2(CPM + 0.5)
dge <- DGEList(
  counts = assay(rse_gene_SRP165875, "counts"),
  genes = rowData(rse_gene_SRP165875)
)
dge <- calcNormFactors(dge)

assays(rse_gene_SRP165875)$logcounts <- cpm(dge, log=TRUE, prior.count=0.5)
assays(rse_gene_SRP165875)$normcounts <- dge$counts

# Filter low expression values
keep <- filterByExpr(dge,
                     keep.lib.sizes=FALSE,
                     group = rse_gene_SRP165875$sra_attribute.age,
)

rse_gene_SRP165875_filt <- rse_gene_SRP165875[keep, ]

## Recovery percentage = 32.84 %
round(nrow(rse_gene_SRP165875_filt) / nrow(rse_gene_SRP165875) * 100, 2)

# Save normalized data
saveRDS(rse_gene_SRP165875_filt, "../data/rse_gene_SRP165875_norm.RDS")

# Save session info
writeLines(capture.output(session_info()), "../doc/dataPreProcessing_session_info.txt")


### Plot data distribution
data <- data.frame(counts = as.vector(assays(rse_gene_SRP165875_filt)$counts))
ggplot(data, aes(x = counts)) +
  geom_histogram(aes(y = ..density..), colour = "darkgray", fill = "lightgray") +
  binwi
  geom_density(fill = "#69b3a2", alpha = 0.3) +
  labs(x = "raw counts", y = "Frecuency") +
  theme(plot.margin = unit(c(2.5, 5, 2.5, 5), "cm"))

filt_logcounts_data <- data.frame(logcounts = as.vector(assays(rse_gene_SRP165875_filt)$logcounts))
ggplot(filt_logcounts_data, aes(x = logcounts)) +
  geom_histogram(aes(y = ..density..), colour = "darkgray", fill = "lightgray") +
  geom_density(fill = "#69b3a2", alpha = 0.3) +
  labs(x = "log(CPM+0.5)", y = "Frecuency") +
  theme(plot.margin = unit(c(2.5, 5, 2.5, 5), "cm"))
