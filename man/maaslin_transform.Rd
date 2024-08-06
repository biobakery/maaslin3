\name{maaslin_transform}
\alias{maaslin_transform}
\title{
Transform abundance data for MaAsLin 3 modeling.
}
\description{
Transform the abundance data according to the \code{transform}
parameter.
}
\usage{
maaslin_transform(filtered_data,
                output,
                transform = 'LOG')
}
\arguments{
    \item{filtered_data}{A data frame of filtered feature abundances. It
    should be formatted with features as columns and samples as rows. The
    column and row names should be the feature names and sample names
    respectively.}
    \item{output}{The output folder to write results.}
    \item{transform}{The transformation to apply to the features after
    normalization and before analysis. The option \code{LOG} (base 2) is
    recommended, but \code{PLOG} (pseudo-log) and \code{NONE} can also be
    used.}
}
\value{
A dataframe of transformed features (features are columns; samples are
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
                                    
    transformed_data = maaslin3::maaslin_filter(filtered_data,
                                    output = 'output')
    
    unlink('output')
    logging::logReset()
}
