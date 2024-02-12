#Name: getData.R
#Author: Diego Carmona Campos
#Last update:  04/02/2024
#Description:  Generates the SummarizedExperiment object for SRP165875 record used in this project

# Load libraries
library(recount3)
library(sessioninfo)

# Select the study of interest
rse_gene_SRP165875 <- create_rse_manual(
  project = "SRP165875",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

# Save data
saveRDS(rse_gene_SRP165875, "../rawData/rse_gene_SRP165875.Rds")

# Save metadata for reproducibility
capture.output(metadata(rse_gene_SRP165875), file="../doc/rse_gene_SRP165875.metadata")
writeLines(capture.output(session_info()), "../doc/getData_session_info.txt")
