library(testthat)
library(maaslin3)

output_tmp <- tempfile()
maaslin_log_arguments(input_data = 'something1', 
                      input_metadata = 'something2',
                      output = output_tmp, 
                      formula = 'something3', 
                      fixed_effects = 'something4', 
                      reference = 'something5', 
                      random_effects = 'something6', 
                      group_effects = 'something7',
                      ordered_effects = 'something8', 
                      strata_effects = 'something9',
                      feature_specific_covariate = 'something10',
                      feature_specific_covariate_name = 'something11',
                      feature_specific_covariate_record = 'something12',
                      min_abundance = 0,
                      min_prevalence = 1,
                      zero_threshold = 2,
                      min_variance = 3,
                      max_significance = 4,
                      normalization = 'TSS',
                      transform = 'LOG',
                      correction = 'BH',
                      standardize = TRUE,
                      unscaled_abundance = 'something16',
                      median_comparison_abundance = FALSE,
                      median_comparison_prevalence = FALSE,
                      median_comparison_abundance_threshold = 5,
                      median_comparison_prevalence_threshold = 6,
                      subtract_median = FALSE,
                      warn_prevalence = TRUE,
                      augment = TRUE,
                      evaluate_only = NULL,
                      plot_summary_plot = TRUE,
                      summary_plot_first_n = 7,
                      coef_plot_vars = 'something18',
                      heatmap_vars = 'something19',
                      plot_associations = TRUE,
                      max_pngs = 8,
                      cores = 9,
                      save_models = FALSE,
                      verbosity = 'FINEST')

lines_in <- readLines(file.path(output_tmp, 'maaslin3.log'))
lines_in <- sub('.*::', '', lines_in)

lines_to_compare <- c("Writing function arguments to log file",
                        "Function arguments",
                        "Input data file: something1",
                        "Input metadata file: something2",
                        "Output folder:",
                        "Formula: something3",
                        "Fixed effects: something4",
                        "Reference: something5",
                        "Random effects: something6",
                        "Group effects: something7",
                        "Ordered effects: something8",
                        "Strata effects: something9",
                        "Feature specific covariate: something10",
                        "Feature specific covariate name: something11",
                        "Feature specific covariate include: something12",
                        "Min Abundance: 0",
                        "Min Prevalence: 1",
                        "Zero Threshold: 2",
                        "Min variance: 3",
                        "Max significance: 4",
                        "Normalization: TSS",
                        "Transform: LOG",
                        "Correction method: BH",
                        "Standardize: TRUE",
                        "Unscaled abundance: something16",
                        "Abundance median comparison: FALSE",
                        "Prevalence median comparison: FALSE",
                        "Abundance median comparison threshold: 5",
                        "Prevalence median comparison threshold: 6",
                        "Subtract median: FALSE",
                        "Warn prevalence: TRUE",
                        "Augment: TRUE",
                        "Evaluate only:",
                        "Cores: 9",
                        "Verifying options selected are valid")

expect_starts_with <- function(strings, prefixes) {
    expect_equal(sapply(seq_along(strings), function(i) {
        startsWith(strings[i], prefixes[i])
    }), rep(TRUE, length(strings)))
}

expect_starts_with(lines_in, lines_to_compare)

logging::logReset()
unlink(output_tmp, recursive = T)
