library(testthat)
library(maaslin3)

expected_results_run1 <- read.table("expected_results_run1.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")

taxa_table <- read.table(system.file(package="maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, sep="\t", row.names = 1)
metadata <- read.table(system.file(package="maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, sep="\t", row.names = 1)

metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
metadata$antibiotics <- factor(metadata$antibiotics, levels = c('No', 'Yes'))

# Run MaAsLin 3
output_tmp <- tempfile()
set.seed(1)
fit_out <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

maaslin_results = read.table(file.path(output_tmp, "significant_results.tsv"), 
                             header = TRUE, 
                             stringsAsFactors = FALSE)

expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(round(as.numeric(expected_results_run1$pval_individual[1:50]),10),
            equals(round(as.numeric(maaslin_results$pval_individual[1:50]),10)))
expect_that(round(as.numeric(expected_results_run1$qval_individual[1:50]),10),
            equals(round(as.numeric(maaslin_results$qval_individual[1:50]),10)))

se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(taxa_table = t(taxa_table)),
    colData = metadata
)

output_tmp <- tempfile()
fit_out <- maaslin3(input_data = se, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

se_multi <- SummarizedExperiment::SummarizedExperiment(
    assays = list(taxa_table_junk = matrix(0, nrow = ncol(taxa_table), ncol = nrow(taxa_table)),
                  another_taxa_table = t(taxa_table)),
    colData = metadata
)

metadata_df <- as(metadata, "DataFrame")
output_tmp <- tempfile()
fit_out <- maaslin3(input_data = se_multi, 
                    input_metadata = metadata_df, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR',
                    assay.type = 'another_taxa_table')

output_tmp <- tempfile()
fit_out <- maaslin3(input_data = se_multi, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR',
                    assay.type = 2)

output_tmp <- tempfile()
set.seed(1)
fit_out <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

# Weird names
taxa_table <- read.table(system.file(package="maaslin3","extdata","HMP2_taxonomy.tsv"), header = TRUE, sep="\t", row.names = 1)
metadata <- read.table(system.file(package="maaslin3","extdata","HMP2_metadata.tsv"), header = TRUE, sep="\t", row.names = 1)

metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
metadata$antibiotics <- factor(metadata$antibiotics, levels = c('No', 'Yes'))

colnames(metadata)[colnames(metadata) == 'dysbiosis_state'] <- c("dysbiosis state")
colnames(metadata)[colnames(metadata) == 'red_meat'] <- c("red meat")
metadata$`red meat` <- factor(metadata$`red meat`)
rownames(metadata) <- gsub('_', ' ', rownames(metadata))
metadata$`dysbiosis state` <- gsub('_', ' ', as.character(metadata$`dysbiosis state`))
metadata$`dysbiosis state` <- factor(metadata$`dysbiosis state`, levels = c('none', 'dysbiosis UC', 'dysbiosis CD'))

rownames(taxa_table) <- gsub('_', ' ', rownames(taxa_table))
colnames(taxa_table) <- gsub('_', ' ', colnames(taxa_table))

output_tmp <- tempfile()
set.seed(1)
fit_out <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ ordered(`red meat`)', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = FALSE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

output_tmp <- tempfile()
fit_out2 <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    ordered_effects = 'red meat',
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = FALSE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

expect_that(fit_out$fit_data_abundance$results$coef,
            equals(fit_out2$fit_data_abundance$results$coef))


output_tmp <- tempfile()
fit_out <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ group(`red meat`)', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = FALSE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

output_tmp <- tempfile()
fit_out2 <- maaslin3(input_data = taxa_table, 
                     input_metadata = metadata, 
                     output = output_tmp, 
                     normalization = 'TSS', 
                     transform = 'LOG', 
                     group_effects = 'red meat',
                     save_models = FALSE, 
                     plot_summary_plot = FALSE, 
                     plot_associations = FALSE, 
                     max_significance = 0.1, 
                     augment = TRUE, 
                     median_comparison_abundance = FALSE, 
                     median_comparison_prevalence = FALSE, 
                     cores=1, 
                     verbosity = 'ERROR')

expect_that(fit_out$fit_data_abundance$results$coef,
            equals(fit_out2$fit_data_abundance$results$coef))

output_tmp <- tempfile()
fit_out <- maaslin3(input_data = taxa_table, 
                    input_metadata = metadata, 
                    output = output_tmp, 
                    normalization = 'TSS', 
                    transform = 'LOG', 
                    formula = '~ diagnosis + `dysbiosis state` + antibiotics + age + reads', 
                    save_models = FALSE, 
                    plot_summary_plot = FALSE, 
                    plot_associations = FALSE, 
                    max_significance = 0.1, 
                    augment = TRUE, 
                    median_comparison_abundance = TRUE, 
                    median_comparison_prevalence = FALSE, 
                    cores=1, 
                    verbosity = 'ERROR')

output_tmp <- tempfile()
expect_error(maaslin3(input_data = taxa_table, 
    input_metadata = metadata, 
    output = output_tmp, 
    normalization = 'TSS', 
    transform = 'LOG', 
    formula = '~ diagnosis + `dysbiosis state` + antibiotics + age + reads', 
    fixed_effects = 'diagnosis',
    save_models = FALSE, 
    plot_summary_plot = FALSE, 
    plot_associations = FALSE, 
    max_significance = 0.1, 
    augment = TRUE, 
    median_comparison_abundance = TRUE, 
    median_comparison_prevalence = FALSE, 
    cores=1, 
    verbosity = 'ERROR'))
