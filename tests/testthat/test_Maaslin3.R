library(testthat)
library(Maaslin3)

expected_results_run1 <- read.table("expected_results_run1.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")

features <- read.table(system.file(package="Maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, row.names = 1, sep="\t")
metadata <- read.table(system.file(package="Maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, row.names = 1, sep="\t")
rownames(features) <- taxa_table$ID; taxa_table$ID <- NULL
rownames(metadata) <- metadata$ID; metadata$ID <- NULL

# run with dysbiosis as a nested categorical variable
param_list <- list(input_data = taxa_table, input_metadata = metadata, min_abundance = 0, min_prevalence = 0, output = 'tmp/', 
                   min_variance = 0, normalization = 'TSS', transform = 'LOG', analysis_method = 'LM', 
                   formula = '~ diagnosis + dysbiosisUC + dysbiosisCD + antibiotics + age + (1 | subject)', 
                   save_scatter = FALSE, save_models = FALSE, plot_heatmap = T, plot_scatter = F, 
                   max_significance = 0.1, augment = TRUE, iterative_mode = TRUE, cores=1)
fit_out <- Maaslin3(param_list)

write.table(rbind(fit_out$fit_data_non_zero$results, 
                  fit_out$fit_data_binary$results), 
            'significant_results.tsv', sep = '\t', row.names = F)
maaslin_results = read.table(file.path("output","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(expected_results_run1$N.not.0[1:50],equals(maaslin_results$N.not.0[1:50]))
expect_that(round(as.numeric(expected_results_run1$pval.value[1:50]),10),equals(round(as.numeric(maaslin_results$pval.value[1:50]),10)))
expect_that(round(as.numeric(expected_results_run1$qval.value[1:50]),10),equals(round(as.numeric(maaslin_results$qval.value[1:50]),10)))
