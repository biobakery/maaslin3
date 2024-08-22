\name{maaslin_check_arguments}
\alias{maaslin_check_arguments}
\title{
Check parameter arguments to ensure a successful MaAsLin 3 run.
}
\description{
Check the arguments provided are valid for further MaAsLin 3 use.
}
\usage{
maaslin_check_arguments(feature_specific_covariate = NULL,
                        feature_specific_covariate_name = NULL,
                        feature_specific_covariate_record = NULL,
                        zero_threshold = 0,
                        normalization = 'TSS',
                        transform = 'LOG',
                        correction = 'BH',
                        evaluate_only = NULL,
                        unscaled_abundance = NULL,
                        median_comparison_abundance = TRUE)
}
\arguments{
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
    \item{normalization}{The normalization to apply to the features before
    transformation and analysis. The option \code{TSS} (total sum scaling)
    is recommended, but \code{CLR} (centered log ratio) and \code{NONE} can
    also be used.}
    \item{transform}{The transformation to apply to the features after
    normalization and before analysis. The option \code{LOG} (base 2) is
    recommended, but \code{PLOG} (pseudo-log) and \code{NONE} can also be
    used.}
    \item{correction}{The correction to obtain FDR-corrected q-values from
    raw p-values. Any valid options for \code{p.adjust} can be used.}
    \item{warn_prevalence}{Warn when prevalence associations are likely
    induced by abundance associations. This requires re-fitting the linear
    models on the TSS log-transformed data.}
    \item{evaluate_only}{Whether to evaluate just the abundnace
    ("abundance") or prevalence ("prevalence") models}
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
}
\value{
    No value is returned, but incompatibile arguments will produce an error.
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
    
    # Prepare parameter lists
    maaslin3::maaslin_check_arguments(zero_threshold = 0,
                                    normalization = 'TSS',
                                    transform = 'LOG',
                                    correction = 'BH',
                                    median_comparison_abundance = TRUE)
                                    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
