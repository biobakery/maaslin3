\name{maaslin3}
\alias{maaslin3}
\title{
MaAsLin 3: A multivariable statistical framework for finding abundance
and
prevalence associations between metadata and high-dimensional microbial
multi-omics data.
}
\description{
This wrapper for all MaAsLin 3 steps finds abundance and prevalence
    associations between microbiome meta-omics features and
    complex metadata in population-scale epidemiological studies.
    The software includes multiple analysis methods (including
    support for multiple covariates, repeated measures, and ordered
    predictors),
    filtering, normalization, and transform options to customize analysis
    for
    your specific study.
}
\usage{
maaslin3(input_data,
    input_metadata,
    output,
    formula = NULL,
    fixed_effects = NULL,
    reference = NULL,
    random_effects = NULL,
    group_effects = NULL,
    ordered_effects = NULL,
    strata_effects = NULL,
    feature_specific_covariate = NULL,
    feature_specific_covariate_name = NULL,
    feature_specific_covariate_record = NULL,
    min_abundance = 0,
    min_prevalence = 0.1,
    zero_threshold = 0,
    min_variance = 0,
    max_significance = 0.1,
    normalization = 'TSS',
    transform = 'LOG',
    correction = 'BH',
    standardize = TRUE,
    unscaled_abundance = NULL,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    median_comparison_abundance_threshold = 0.25,
    median_comparison_prevalence_threshold = 0.25,
    augment = TRUE,
    evaluate_only = NULL,
    plot_summary_plot = TRUE,
    summary_plot_first_n = 25,
    coef_plot_vars = NULL,
    heatmap_vars = NULL,
    plot_associations = TRUE,
    max_pngs = 30,
    cores = 1,
    save_models = FALSE)
}
\arguments{
    \item{input_data}{A data frame of feature abundances or read counts or a
    filepath to a tab-delimited file with abundances. It should be formatted
    with features as columns and samples as rows (or the transpose). The
    column
    and row names should be the feature names and sample names
    respectively.}
    \item{input_metadata}{A data frame of per-sample metadata or a filepath
    to
    a tab-delimited file with metadata. It should be formatted with
    variables
    as columns and samples as rows (or the transpose). The column and row
    names
    should be the variable names and sample names respectively.}
    \item{output}{The output folder to write results.}
    \item{formula}{A formula in \code{lme4} format. Random effects,
    interactions, and functions of the metadata can be included (note that
    these functions will be applied after standardization if
    \code{standardize=TRUE}). Group, ordered, and strata variables can be
    specified as: \code{group(grouping_variable)},
    \code{ordered(ordered_variable)} and \code{strata(strata_variable)}. The
    other variable options below will not be considered if a formula is
    set.}
    \item{fixed_effects}{A vector of variable names to be included as fixed
    effects.}
    \item{reference}{For a variable with more than two levels supplied with
    \code{fixed_effects}, the factor to use as a reference provided as a
    string
    of 'variable,reference' semi-colon delimited for multiple variables.}
    \item{random_effects}{A vector of variable names to be included as
    random
    intercepts.}
    \item{group_effects}{A factored categorical variable to be included for
    group testing. An ANOVA-style test will be performed to assess whether
    any
    of the variable's levels are significant, and no coefficients or
    individual
    p-values will be returned.}
    \item{ordered_effects}{A factored categorical variable to be included.
    Consecutive levels will be tested for significance against each other,
    and
    the resulting associations will correspond to effect sizes, standard
    errors, and significances of each level versus the previous.}
    \item{strata_effects}{A vector with one variable name to be included as
    the
    strata variable in case-control studies. Strata cannot be combined with
    random effects.}
    \item{feature_specific_covariate}{A table of feature-specific covariates
    or
    a filepath to a tab-delimited file with feature-specific covariates. It
    should be formatted with features as columns and samples as rows (or the
    transpose). The row names and column names should be the same as those
    of
    the \code{input_data}: the column and row names should be the feature
    names
    and sample names respectively.}
    \item{feature_specific_covariate_name}{The name for the feature-specific
    covariates when fitting the models.}
    \item{feature_specific_covariate_record}{Whether to keep the
    feature-specific covariates in the outputs when calculating p-values,
    writing results, and displaying plots.}
    \item{min_abundance}{Features with abundances more than
    \code{min_abundance} in \code{min_prevalence} of the samples will be
    included for analysis. The threshold is applied after normalization and
    before transformation.}
    \item{min_prevalence}{See \code{min_abundance}.}
    \item{zero_threshold}{Abundances less than or equal to
    \code{zero_threshold} will be treated as zeros. This is primarily to be
    used when the abundance table has likely low-abundance false positives.}
    \item{min_variance}{Features with abundance variances less than or equal
    to
    \code{min_variance} will be dropped. This is primarily used for dropping
    features that are entirely zero.}
    \item{max_significance}{The FDR corrected q-value threshold for
    significance used in selecting which associations to write as
    significant
    and to plot.}
    \item{normalization}{The normalization to apply to the features before
    transformation and analysis. The option \code{TSS} (total sum scaling)
    is
    recommended, but \code{CLR} (centered log ratio) and \code{NONE} can
    also
    be used.}
    \item{transform}{The transformation to apply to the features after
    normalization and before analysis. The option \code{LOG} (base 2) is
    recommended, but \code{PLOG} (pseudo-log) and \code{NONE} can also be
    used.}
    \item{correction}{The correction to obtain FDR-corrected q-values from
    raw
    p-values. Any valid options for \code{p.adjust} can be used.}
    \item{standardize}{Whether to apply z-scores to continuous metadata
    variables so they are on the same scale. This is recommended in order to
    compare coefficients across metadata variables, but note that functions
    of the metadata specified in the \code{formula} will apply after
    standardization.}
    \item{unscaled_abundance}{A data frame with a single column of absolute
    abundances or a filepath to such a tab-delimited file. The row names
    should match the names of the samples in \code{input_data} and
    \code{input_metadata}. When using spike-ins, the single column should
    have the same name as one of the features in \code{input_data}, and the
    \code{unscaled_abundance} should correspond to the absolute quantity of
    the spike-in. When using total abundance scaling, the single column
    should have the name 'total', and the \code{unscaled_abundance} should
    correspond to the total abundance of each sample.}
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
    \item{augment}{Add extra lowly-weighted 0s and 1s to avoid linear
    separability.}
    \item{evaluate_only}{Whether to evaluate just the abundnace
    ("abundance") or prevalence ("prevalence") models}
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
    \item{cores}{How many cores to use when fitting models. (Using multiple
    cores will likely be faster only for large datasets or complex models.}
    \item{save_models}{Whether to return the fit models and save them to an
    RData file.}
}
\value{
    A list containing the following items:
    \describe{
    \item{(1)}{\code{data}: A dataframe of feature abundances with the
    retained samples for fitting.}
    \item{(2)}{\code{normalized_data}: A dataframe of normalized feature
    abundances.}
    \item{(3)}{\code{filtered_data}: A dataframe of feature abundances on
    the original scale after normalization and filtering.}
    \item{(4)}{\code{transformed_data}: A dataframe of feature abundances
    after filtering, normalization, and transformation.}
    \item{(5)}{\code{metadata}: A dataframe of metadata with the retained
    samples for fitting.}
    \item{(6)}{\code{standardized_metadata}: A dataframe of metadata after
    scaling (if selected).}
    \item{(7)}{\code{formula}: Checked or constructed formula(s) specifying
    the model to be fit.}
    \item{(8)}{\code{fit_data_abundance}: The results from the fit abundance
    models (see \code{\link{maaslin_fit}}).}
    \item{(9)}{\code{fit_data_prevalence}: The results from the fit
    prevalence models (see \code{\link{maaslin_fit}}).}
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
    fit_out <- maaslin3::maaslin3(input_data = taxa_table,
                                input_metadata = metadata,
                                output = 'output',
                                formula = '~ diagnosis + dysbiosis_state +
                                antibiotics + age + reads',
                                plot_summary_plot = FALSE,
                                plot_associations = FALSE)
                                
    unlink('output', recursive=TRUE)
    logging::logReset()
}
