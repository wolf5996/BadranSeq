# code to prepare `gsea_object` dataset goes here

# load libraries
library(clusterProfiler)
library(tidyverse)
library(ReactomePA)

# import GSEA results
gsea_go_bp_results <- read_rds(file = "inst/extdata/gsea_go_bp_results.rds")
gsea_go_reactome_results <- read_rds(file = "inst/extdata/gsea_reactome_results.rds")

# save GSEA object as rda file for package use
usethis::use_data(gsea_go_bp_results, gsea_go_reactome_results, overwrite = TRUE)
