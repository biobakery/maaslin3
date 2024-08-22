\name{maaslin_write_results_lefse_format}
\alias{maaslin_write_results_lefse_format}
\title{
Write the results from a MaAsLin 3 run in LEfSe format.
}
\description{
Write the results from a MaAsLin 3 run to the output folder in LEfSe
format.
}
\usage{
maaslin_write_results_lefse_format(output,
                                fit_data_abundance,
                                fit_data_prevalence)
}
\arguments{
    \item{output}{The output folder to write results.}
    \item{fit_data_abundance}{The abundance outputs of
    \code{\link{maaslin_fit}}.}
    \item{fit_data_prevalence}{The prevalence outputs of
    \code{\link{maaslin_fit}}.}
}
\value{
Results will be written to the \code{lefse_style_results_abundance.res}
file in the folder \code{output}. No value is returned.
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
        random_effects_formula,
        warn_prevalence = FALSE)
    
    maaslin3::maaslin_write_results_lefse_format(
        output = 'output',
        maaslin_results$fit_data_abundance,
        maaslin_results$fit_data_prevalence)
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
