\name{maaslin_reorder_data}
\alias{maaslin_reorder_data}
\title{
Reorder the abundance data and metadata.
}
\description{
Reorder the abundance data and metadata to ensure samples are rows and
remove any samples without abundances or metadata.
}
\usage{
maaslin_reorder_data(data,
                    metadata,
                    feature_specific_covariate = NULL,
                    unscaled_abundance = NULL)
}
\arguments{
    \item{data}{A data frame of feature abundances or read counts. It should
    be formatted with features as columns and samples as rows (or the
    transpose). The column and row names should be the feature names and
    sample names respectively.}
    \item{metadata}{A data frame of per-sample metadata. It should be
    formatted with variables as columns and samples as rows (or the
    transpose). The column and row names should be the variable names and
    sample names respectively.}
    \item{feature_specific_covariate}{A table of feature-specific
    covariates. It should be formatted with features as columns and samples
    as rows (or the transpose). The row names and column names should be the
    same as those of the \code{input_data}: the column and row names should
    be the feature names and sample names respectively.}
    \item{unscaled_abundance}{A data frame with a single column of absolute
    abundances. The row names should match the names of the samples in
    \code{input_data} and \code{input_metadata}. When using spike-ins, the
    single column should have the same name as one of the features in
    \code{input_data}, and the \code{unscaled_abundance} should correspond
    to the absolute quantity of the spike-in. When using total abundance
    scaling, the single column should have the name 'total', and the
    \code{unscaled_abundance} should correspond to the total abundance of
    each sample.}
}
\value{
    A list containing the following items:
    \describe{
    \item{(1)}{\code{data}: A data frame of feature abundances.}
    \item{(2)}{\code{metadata}: A data frame of metadata.}
    \item{(3)}{\code{feature_specific_covariate}: A data frame of feature
    specific covariates.}
    \item{(4)}{\code{unscaled_abundance}: A data frame of unscaled
    abundances.}
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
        taxa_table,
        metadata)
    
    unlink('output')
    logging::logReset()
}
