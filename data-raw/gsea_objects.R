# Code to prepare GSEA results datasets goes here

# Check if required packages are available
required_packages <- c("clusterProfiler", "ReactomePA")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("The following packages are required to build package data:\n",
       paste(missing_packages, collapse = ", "), "\n",
       "Install with: BiocManager::install(c('",
       paste(missing_packages, collapse = "', '"), "'))")
}

# Import GSEA results from external data files
gsea_go_bp_results <- readr::read_rds(file = "inst/extdata/gsea_go_bp_results.rds")
gsea_reactome_results <- readr::read_rds(file = "inst/extdata/gsea_reactome_results.rds")

# Save GSEA objects as .rda files for package use
usethis::use_data(gsea_go_bp_results, gsea_reactome_results, overwrite = TRUE)
