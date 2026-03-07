\name{maaslin_contrast_test}
\alias{maaslin_contrast_test}
\title{
Perform a contrast from a fit MaAsLin 3 model.
}
\description{
Perform a contrast test (\code{lmerTest::contest} for mixed effects linear; 
\code{multcomp::glht} for all others) using a named contrast matrix and 
right hand side. One contrast test is applied per row of the matrix.
}
\usage{
maaslin_contrast_test(maaslin3_fit,
                        contrast_mat,
                        rhs = NULL,
                        max_significance = 0.1,
                        correction = 'BH',
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        subtract_median = FALSE,
                        small_random_effects = FALSE,
                        bypass_small_group_warning = FALSE,
                        evaluate_only = NULL)
}
\arguments{
    \item{maaslin3_fit}{The output of \code{\link{maaslin_fit}} with 
    \code{save_models = TRUE}.}
    \item{contrast_mat}{A matrix with one row per contrast test to run. The
    columns will be matched to the coefficients of the model by name. 
    Contrast vector coefficients need not be specified if they would be zero. 
    If row names are provided, they will be used to label the test in the
    results.}
    \item{rhs}{The right hand size of the contrast test. The length should
    be the same as the number of rows in the \code{contrast_mat}. This will
    default to 0 or the median comparison if \code{median_comparison=TRUE}.}
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
    \item{subtract_median}{Subtract the median from the coefficients.}
    \item{small_random_effects}{Automatically replace random effects with 
    fixed effects in the logistic prevalence model to handle low numbers
    of observations per group.}
    \item{bypass_small_group_warning}{When \code{TRUE}, bypass the warning
    that is normally added to results when there are fewer than 4 average
    observations per random effect group.}
    \item{evaluate_only}{Whether to evaluate just the abundnace
    ("abundance") or prevalence ("prevalence") models}
}
\value{
A dataframe with the following columns:
\describe{
    \item{(1)}{\code{feature}: The feature involved in the association.}
    \item{(2)}{\code{test}: The contrast test name.}
    \item{(3)}{\code{coef}: The coefficient of the association: the slope
    coefficient in the abundance model and the change in log odds in the
    prevalence model.}
    \item{(4)}{\code{null_hypothesis}: The value of the null hypothesis
    against which the coefficients are tested (zero or the per-metadatum 
    median).}
    \item{(5)}{\code{stderr}: The standard error of the coefficient.}
    \item{(6)}{\code{pval_individual}: The (uncorrected) p-value of the
    association.}
    \item{(7)}{\code{qval_individual}: The FDR corrected q-value of the
    association. FDR correction is performed over all associations in the
    abundance and prevalence modeling without errors together.}
    \item{(8)}{\code{pval_joint}: The p-value of the overall association
    (combining abundance and prevalence) by taking the minimum of the
    abundance and logistic p-values and applying the Beta(1,2) CDF. These
    will be the same in the abundance and prevalence results for an
    association.}
    \item{(9)}{\code{qval_joint}: The FDR corrected q-value of the
    association. FDR correction is performed over all joint p-values
    without errors.}
    \item{(10)}{\code{error}: Any error produced by the model during fitting.
    NA otherwise.}
    \item{(11)}{\code{model}: \code{linear} for the abundance models and
    \code{logistic} for the prevalence models.}
    \item{(12)}{\code{N}: The number of data points for the association's
    feature.}
    \item{(13)}{\code{N_not_zero}: The number of non-zero data points for
    the association's feature.}
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
    taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1)
    
    # Read metadata table
    metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package =
    "maaslin3")
    metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)
    
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
                                
    contrast_mat <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 
        ncol = 4, nrow = 2, byrow = TRUE)
        
    colnames(contrast_mat) <- c("diagnosisUC",
                                "diagnosisCD",
                                "dysbiosis_statedysbiosis_UC",
                                "dysbiosis_statedysbiosis_CD")
                                
    rownames(contrast_mat) <- c("diagnosis_test", "dysbiosis_test")
                                
    maaslin_contrast_test(maaslin3_fit = fit_out,
                            contrast_mat = contrast_mat)
                                
    unlink('output', recursive=TRUE)
    logging::logReset()
}
