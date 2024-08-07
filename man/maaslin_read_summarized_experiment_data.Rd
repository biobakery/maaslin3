\name{maaslin_read_summarized_experiment_data}
\alias{maaslin_read_summarized_experiment_data}
\title{
Read abundance data and metadata from a SummarizedExperiment.
}
\description{
Read abundance data and metadata from a SummarizedExperiment.
}
\usage{
maaslin_read_summarized_experiment_data(summarized_experiment)
}
\arguments{
    \item{summarized_experiment}{A \code{SummarizedExperiment} object containing
    the features as the \code{assays} and metadata as \code{colData}.}
}
\value{
    A list containing the following items:
    \describe{
    \item{(1)}{\code{data}: A data frame of feature abundances.}
    \item{(2)}{\code{metadata}: A data frame of metadata.}
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

    # Create example SE
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = t(read_data_list$data), 
        colData = read_data_list$metadata)
    
    # This replaces maaslin_read_data above
    read_data_list <- maaslin_read_summarized_experiment_data(se)
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
