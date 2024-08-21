\name{maaslin_fit}
\alias{maaslin_fit}
\title{
Fit MaAsLin 3 models.
}
\description{
Fit the abundance data with abundance and prevalence models to discover
feature-metadata associations.
}
\usage{
maaslin_fit(filtered_data,
            transformed_data,
            metadata,
            formula,
            random_effects_formula,
            feature_specific_covariate = NULL,
            feature_specific_covariate_name = NULL,
            feature_specific_covariate_record = NULL,
            zero_threshold = 0,
            max_significance = 0.1,
            correction = 'BH',
            median_comparison_abundance = TRUE,
            median_comparison_prevalence = FALSE,
            median_comparison_abundance_threshold = 0.25,
            median_comparison_prevalence_threshold = 0.25,
            subtract_median = FALSE,
            augment = TRUE,
            evaluate_only = NULL,
            cores = 1,
            save_models = FALSE)
}
\arguments{
    \item{filtered_data}{A data frame of filtered feature abundances. It
    should be formatted with features as columns and samples as rows. The
    column and row names should be the feature names and sample names
    respectively.}
    \item{transformed_data}{A data frame of transformed feature abundances.
    It should be formatted with features as columns and samples as rows. The
    column and row names should be the feature names and sample names
    respectively.}
    \item{metadata}{A data frame of per-sample metadata. It should be
    formatted with variables as columns and samples as rows. The column and
    row names should be the variable names and sample names respectively.}
    \item{formula}{A formula in \code{lme4} format as from
    \code{\link{maaslin_check_formula}}.}
    \item{random_effects_formula}{A formula in \code{lme4} format as from
    \code{\link{maaslin_check_formula}}.}
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
    \item{zero_threshold}{Abundances less than or equal to
    \code{zero_threshold} will be treated as zeros. This is primarily to be
    used when the abundance table has likely low-abundance false positives.}
    \item{max_significance}{The FDR corrected q-value threshold for
    significance used in selecting which associations to write as
    significant and to plot.}
    \item{correction}{The correction to obtain FDR-corrected q-values from
    raw p-values. Any valid options for \code{p.adjust} can be used.}
    \item{median_comparison_abundance}{Test abundance coefficients against a
    null value corresponding to the median coefficient for a metadata
    variable across the features. This is recommended for relative abundance
    data but should not be used for absolute abundance data.}
    \item{median_comparison_prevalence}{Test prevalence coefficients against
    a null value corresponding to the median coefficient for a metadata
    variable across the features. This is only recommended if the analyst is
    interested in how feature prevalence associations compare to each other
    or if there is likely strong compositionality-induced sparsity.}
    \item{median_comparison_abundance_threshold}{Coefficients within
    \code{median_comparison_abundance_threshold} of the median association
    will automatically be counted as insignificant (p-value set to 1) since
    they likely represent compositionality-induced associations. This
    threshold will be divided by the metadata variable's standard deviation
    if the metadatum is continuous to ensure the threshold applies to the
    right scale.}
    \item{median_comparison_prevalence_threshold}{Same as
    \code{median_comparison_abundance_threshold} but applied to the
    prevalence associations.}
    \item{subtract_median}{Subtract the median from the coefficients.}
    \item{augment}{Add extra lowly-weighted 0s and 1s to avoid linear
    separability.}
    \item{evaluate_only}{Whether to evaluate just the abundnace
    ("abundance") or prevalence ("prevalence") models}
    \item{cores}{How many cores to use when fitting models. (Using multiple
    cores will likely be faster only for large datasets or complex models.}
    \item{save_models}{Whether to return the fit models and save them to an
    RData file.}
}
\value{
A list containing the following named items:
    \describe{
    \item{(1)}{\code{fit_data_abundance}: The results from the fit abundance
    models.}
    \item{(2)}{\code{fit_data_prevalence}: The results from the fit
    prevalence models.}
    }
    
The \code{fit_data_abundance} and \code{fit_data_prevalence} items have
the same structure. They are both lists with the following named items:
\describe{
    \item{(1)}{\code{results}: A results table with the modeled associations
    (see below).}
    \item{(2)}{\code{residuals}: A features (rows) by samples (columns)
    dataframe of residuals from the models.}
    \item{(3)}{\code{fitted}: A features (rows) by samples (columns)
    dataframe of fitted values from the models.}
    \item{(4)}{\code{ranef}: A features (rows) by random effect (columns)
    dataframe of random effects from the models. If multiple random effects
    are specified, this is a dataframe of dataframes.}
    \item{(5)}{\code{fits}: If \code{save_models=TRUE}, this is a list of
    the fit models.}
    }
    
The \code{results} tables contain the following columns for each
association (row):
\describe{
    \item{(1)}{\code{feature}: The feature involved in the association.}
    \item{(2)}{\code{metadata}: The metadata variable involved in the
    association.}
    \item{(3)}{\code{value}: The value of the metadata variable: the
    metadata variable itself if continuous or the level if categorical.}
    \item{(4)}{\code{name}: The name of the model component involved in the
    association: the metadata variable itself if continuous or a
    concatenated version of the metadata variable and level if categorical.}
    \item{(5)}{\code{coef}: The coefficient of the association: the slope
    coefficient in the abundance model and the change in log odds in the
    prevalence model.}
    \item{(6)}{\code{stderr}: The standard error of the coefficient.}
    \item{(7)}{\code{pval_individual}: The (uncorrected) p-value of the
    association.}
    \item{(8)}{\code{error}: Any error produced by the model during fitting.
    NA otherwise.}
    \item{(9)}{\code{qval_individual}: The FDR corrected q-value of the
    association. FDR correction is performed over all associations in the
    abundance and prevalence modeling separately.}
    \item{(10)}{\code{model}: \code{linear} for the abundance models and
    \code{logistic} for the prevalence models.}
    \item{(11)}{\code{N}: The number of data points for the association's
    feature.}
    \item{(12)}{\code{N.not.zero}: The number of non-zero data points for
    the association's feature.}
    \item{(13)}{\code{pval_joint}: The p-value of the overall association
    (combining abundance and prevalence) by taking the minimum of the
    abundance and logistic p-values and applying the Beta(1,2) CDF. These
    will be the same in the abundance and prevalence results for an
    association.}
    \item{(14)}{\code{qval_joint}: The FDR corrected q-value of the
    association. FDR correction is performed over all joint p-values.}
    }
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
        age + reads',
        plot_summary_plot = FALSE,
        plot_associations = FALSE)
    
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
        random_effects_formula)
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
