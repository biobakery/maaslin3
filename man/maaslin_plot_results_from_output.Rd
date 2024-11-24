\name{maaslin_plot_results_from_output}
\alias{maaslin_plot_results_from_output}
\title{
Plot the results from a MaAsLin 3 run.
}
\description{
Two types of plots are generated. First, the summary plot contains
sorted per-feature coefficients plotted with their standard errors for
key variables and a heatmap summarizing the remaining variables. Second,
for significant features, association plots (scatterplots, boxplots, or
tables depending on the association) are generated to visualize and
verify the model fits. The data are shown with their transformed values
in the association plots since this is the scale on which the models are
fit. In comparison to \code{maaslin_plot_results} that needs the entire
\code{maaslin_fit} list, only the parameter list and an outputs
directory containing a completed run are needed for
\code{maaslin_plot_results_from_output}.
}
\usage{
maaslin_plot_results_from_output(output,
                                metadata,
                                normalization,
                                transform,
                                feature_specific_covariate = NULL,
                                feature_specific_covariate_name = NULL,
                                feature_specific_covariate_record = NULL,
                                median_comparison_abundance = TRUE,
                                median_comparison_prevalence = FALSE,
                                max_significance = 0.1,
                                plot_summary_plot = TRUE,
                                summary_plot_first_n = 25,
                                coef_plot_vars = NULL,
                                heatmap_vars = NULL,
                                plot_associations = TRUE,
                                max_pngs = 30,
                                balanced = FALSE)
}
\arguments{
    \item{output}{The output folder to write results.}
    \item{metadata}{A data frame of per-sample metadata. It should be
    formatted with variables as columns and samples as rows. The column and
    row names should be the variable names and sample names respectively.}
    \item{normalization}{The normalization to apply to the features before
    transformation and analysis. The option \code{TSS} (total sum scaling)
    is recommended, but \code{CLR} (centered log ratio) and \code{NONE} can
    also be used.}
    \item{transform}{The transformation to apply to the features after
    normalization and before analysis. The option \code{LOG} (base 2) is
    recommended, but \code{PLOG} (pseudo-log) and \code{NONE} can also be
    used.}
    \item{feature_specific_covariate}{A table of feature-specific covariates
    or a filepath to a tab-delimited file with feature-specific covariates.
    It should be formatted with features as columns and samples as rows (or
    the transpose). The row names and column names should be the same as
    those of the \code{input_data}: the column and row names should be the
    feature names and sample names respectively.}
    \item{feature_specific_covariate_name}{The name for the feature-specific
    covariates when fitting the models.}
    \item{feature_specific_covariate_record}{Whether to keep the
    feature-specific covariates in the outputs when calculating p-values,
    writing results, and displaying plots.}
    \item{median_comparison_abundance}{Test abundance coefficients against a
    null value corresponding to the median coefficient for a metadata
    variable across the features. This is recommended for relative abundance
    data but should not be used for absolute abundance data.}
    \item{median_comparison_prevalence}{Test prevalence coefficients against
    a null value corresponding to the median coefficient for a metadata
    variable across the features. This is only recommended if the analyst is
    interested in how feature prevalence associations compare to each other
    or if there is likely strong compositionality-induced sparsity.}
    \item{max_significance}{The FDR corrected q-value threshold for
    significance used in selecting which associations to write as
    significant and to plot.}
    \item{plot_summary_plot}{Generate a summary plot of significant
    associations.}
    \item{summary_plot_first_n}{Include the top \code{summary_plot_first_n}
    features with significant associations.}
    \item{coef_plot_vars}{Vector of variable names to be used in the
    coefficient plot section of the summary plot. Continuous variables
    should match the metadata column name, and categorical variables should
    be of the form \code{"[variable] [level]"}.}
    \item{heatmap_vars}{Vector of variable names to be used in the heatmap
    section of the summary plot. Continuous variables should match the
    metadata column name, and categorical variables should be of the form
    \code{"[variable] [level]"}.}
    \item{plot_associations}{Whether to generate plots for significant
    associations.}
    \item{max_pngs}{The top \code{max_pngs} significant associations will be
    plotted.}
    \item{balanced}{If set to TRUE the summary plot will 
    show the top N features of each variable included in 
    \code{coef_plot_vars} where N is equal to: 
    \code{ceiling(summary_plot_first_n/length(coef_plot_vars))}. Will error
    if \code{coef_plot_vars} = \code{NULL}}
}

\value{
Results will be written to the \code{figures} folder within the folder
\code{output}. The list of individual association plots is returned if
\code{plot_associations=TRUE}. In the heatmap of the summary plot, one
star corresponds to the user-set \code{max_significance} and two stars
corresponds to the user-set \code{max_signifiance} divided by 10.
}
\author{
    William Nickols<willnickols@g.harvard.edu>,\cr
    Jacob Nearing<nearing@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{
    # Read features table
    taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package =
    "maaslin3")
    taxa_table <- read.csv(taxa_table_name, sep = '\t')
    
    # Read metadata table
    metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package =
    "maaslin3")
    metadata <- read.csv(metadata_name, sep = '\t')
    
    metadata$diagnosis <-
        factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
    metadata$dysbiosis_state <-
        factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC',
        'dysbiosis_CD'))
    metadata$antibiotics <-
        factor(metadata$antibiotics, levels = c('No', 'Yes'))
    
    #Run MaAsLin3
    maaslin3::maaslin_log_arguments(
        input_data = taxa_table,
        input_metadata = metadata,
        output = 'output',
        formula = '~ diagnosis + dysbiosis_state + antibiotics +
        age + reads')
    
    read_data_list <- maaslin3::maaslin_read_data(
        taxa_table,
        metadata)
    read_data_list <- maaslin3::maaslin_reorder_data(
        read_data_list$data,
        read_data_list$metadata)
    
    data <- read_data_list$data
    metadata <- read_data_list$metadata
    
    formulas <- maaslin3::maaslin_check_formula(
        data,
        metadata,
        input_formula = '~ diagnosis + dysbiosis_state + antibiotics +
        age + reads')
    
    formula <- formulas$formula
    random_effects_formula <- formulas$random_effects_formula
    
    normalized_data = maaslin3::maaslin_normalize(data,
                                output = 'output')
                                
    filtered_data = maaslin3::maaslin_filter(normalized_data,
                                output = 'output')
                                
    transformed_data = maaslin3::maaslin_transform(filtered_data,
                                output = 'output')
                                
    standardized_metadata = maaslin3::maaslin_process_metadata(
        metadata,
        formula = formula)
    
    maaslin_results = maaslin3::maaslin_fit(
        filtered_data,
        transformed_data,
        standardized_metadata,
        formula,
        random_effects_formula,
        warn_prevalence = FALSE)
    
    maaslin3::maaslin_write_results(
        output = 'output',
        maaslin_results$fit_data_abundance,
        maaslin_results$fit_data_prevalence,
        random_effects_formula)
    
    maaslin3::maaslin_plot_results_from_output(
        output = 'output',
        metadata,
        normalization = "TSS",
        transform = "LOG")
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
