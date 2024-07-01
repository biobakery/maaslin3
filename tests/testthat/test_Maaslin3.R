library(testthat)
library(maaslin3)

expected_results_run1 <- read.table("expected_results_run1.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")

taxa_table <- read.table(system.file(package="maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, sep="\t")
metadata <- read.table(system.file(package="maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, sep="\t")
rownames(taxa_table) <- taxa_table$ID; taxa_table$ID <- NULL
rownames(metadata) <- metadata$ID; metadata$ID <- NULL

#Prepare parameter lists 
param_list <- list(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = 'output', 
                   normalization = 'TSS', 
                   transform = 'LOG', 
                   formula = '~ diagnosis + dysbiosisUC + dysbiosisCD + antibiotics + age + reads_filtered + (1 | subject)', 
                   save_models = FALSE, 
                   plot_heatmap = T, 
                   plot_associations = T, 
                   max_significance = 0.1, 
                   augment = TRUE, 
                   median_comparison_abundance = TRUE, 
                   median_comparison_prevalence = FALSE, 
                   cores=1)

#Run MaAsLin3
fit_out <- maaslin3::maaslin3(param_list)

maaslin_results = read.table(file.path("output","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(expected_results_run1$N.not.0[1:50],equals(maaslin_results$N.not.0[1:50]))
expect_that(round(as.numeric(expected_results_run1$pval.value[1:50]),10),equals(round(as.numeric(maaslin_results$pval.value[1:50]),10)))
expect_that(round(as.numeric(expected_results_run1$qval.value[1:50]),10),equals(round(as.numeric(maaslin_results$qval.value[1:50]),10)))
