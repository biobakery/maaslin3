\name{maaslin_check_formula}
\alias{maaslin_check_formula}
\title{
Check a MaAsLin 3 formula to ensure a proper MaAsLin 3 run.
}
\description{
Ensure that the formula provided is valid. Only one of
\code{maaslin_compute_formula} or \code{maaslin_check_formula} should be
used.
}
\usage{
maaslin_check_formula(data,
                    metadata,
                    input_formula = NULL,
                    feature_specific_covariate_name = NULL)
}
\arguments{
    \item{data}{A data frame of feature abundances. It should be formatted
    with features as columns and samples as rows. The column and row names
    should be the feature names and sample names respectively.}
    \item{metadata}{A data frame of per-sample metadata. It should be
    formatted with variables as columns and samples as rows. The column and
    row names should be the variable names and sample names respectively.}
    \item{input_formula}{A formula in \code{lme4} format. Random effects,
    interactions, and functions of the metadata can be included (note that
    these functions will be applied after standardization if
    \code{standardize=TRUE}). Group, ordered, and strata variables can be
    specified as: \code{group(grouping_variable)},
    \code{ordered(ordered_variable)} and \code{strata(strata_variable)}. The
    other variable options below will not be considered if a formula is
    set.}
    \item{feature_specific_covariate_name}{The name for the feature-specific
    covariates when fitting the models.}
}
\value{
A list containing the following named items:
    \describe{
    \item{(1)}{\code{formula}: The constructed formula.}
    \item{(2)}{\code{random_effects_formula}: A formula for the random
    effects.}
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
    
    unlink('output', recursive=TRUE)
    logging::logReset()
}
