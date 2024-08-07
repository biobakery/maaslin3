\name{maaslin_filter}
\alias{maaslin_filter}
\title{
Filter abundance data before MaAsLin 3 model fitting.
}
\description{
Set abundances below \code{zero_threshold} to zero, remove features
without abundances more than \code{min_abundance} in
\code{min_prevalence} of the samples, and remove features with variances
less than or equal to \code{min_variance}.
}
\usage{
maaslin_filter(normalized_data,
            output,
            min_abundance = 0,
            min_prevalence = 0.1,
            zero_threshold = 0,
            min_variance = 0)
}
\arguments{
    \item{normalized_data}{A data frame of normalized feature abundances. It
    should be formatted with features as columns and samples as rows. The
    column and row names should be the feature names and sample names
    respectively.}
    \item{output}{The output folder to write results.}
    \item{min_abundance}{Features with abundances more than
    \code{min_abundance} in \code{min_prevalence} of the samples will be
    included for analysis. The threshold is applied after normalization and
    before transformation.}
    \item{min_prevalence}{See \code{min_abundance}.}
    \item{zero_threshold}{Abundances less than or equal to
    \code{zero_threshold} will be treated as zeros. This is primarily to be
    used when the abundance table has likely low-abundance false positives.}
    \item{min_variance}{Features with abundance variances less than or equal
    to \code{min_variance} will be dropped. This is primarily used for
    dropping features that are entirely zero.}
}
\value{
A dataframe of filtered features (features are columns; samples are
rows).
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
    
    normalized_data = maaslin3::maaslin_normalize(data,
                                    output = 'output')
                                    
    filtered_data = maaslin3::maaslin_filter(normalized_data,
                                    output = 'output')
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
