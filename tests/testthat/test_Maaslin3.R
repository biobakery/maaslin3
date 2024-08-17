library(testthat)
library(maaslin3)

expected_results_run1 <- read.table("expected_results_run1.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")

taxa_table <- read.table(system.file(package="maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, sep="\t")
metadata <- read.table(system.file(package="maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, sep="\t")

metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
metadata$antibiotics <- factor(metadata$antibiotics, levels = c('No', 'Yes'))

# Run MaAsLin 3
output_tmp <- tempfile()
fit_out <- maaslin3::maaslin3(input_data = taxa_table, 
                            input_metadata = metadata, 
                            output = output_tmp, 
                            normalization = 'TSS', 
                            transform = 'LOG', 
                            formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                            save_models = FALSE, 
                            plot_summary_plot = T, 
                            plot_associations = T, 
                            max_significance = 0.1, 
                            augment = TRUE, 
                            median_comparison_abundance = TRUE, 
                            median_comparison_prevalence = FALSE, 
                            cores=1)

maaslin_results = read.table(file.path(output_tmp, "significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(expected_results_run1$N.not.0[1:50],equals(maaslin_results$N.not.0[1:50]))
expect_that(round(as.numeric(expected_results_run1$pval_individual[1:50]),10),equals(round(as.numeric(maaslin_results$pval_individual[1:50]),10)))
expect_that(round(as.numeric(expected_results_run1$qval_individual[1:50]),10),equals(round(as.numeric(maaslin_results$qval_individual[1:50]),10)))

unlink(output_tmp, recursive = T)
