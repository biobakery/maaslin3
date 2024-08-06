\name{maaslin_compute_formula}
\alias{maaslin_compute_formula}
\title{
Compute a formula for a MaAsLin 3 run based on the specified effects.
}
\description{
Compute a formula using variables provided through \code{fixed_effects},
\code{random_effects}, \code{group_effects}, \code{ordered_effects}, and
\code{strata_effects}. Only one of \code{maaslin_compute_formula} or
\code{maaslin_check_formula} should be used.
}
\usage{
maaslin_compute_formula(data,
                        metadata,
                        fixed_effects = NULL,
                        random_effects = NULL,
                        group_effects = NULL,
                        ordered_effects = NULL,
                        strata_effects = NULL,
                        feature_specific_covariate_name = NULL)
}
\arguments{
    \item{data}{A data frame of feature abundances. It should be formatted
    with features as columns and samples as rows. The column and row names
    should be the feature names and sample names respectively.}
    \item{metadata}{A data frame of per-sample metadata. It should be
    formatted with variables as columns and samples as rows. The column and
    row names should be the variable names and sample names respectively.}
    \item{fixed_effects}{A vector of variable names to be included as fixed
    effects.}
    \item{random_effects}{A vector of variable names to be included as
    random intercepts.}
    \item{group_effects}{A factored categorical variable to be included for
    group testing. An ANOVA-style test will be performed to assess whether
    any of the variable's levels are significant, and no coefficients or
    individual p-values will be returned.}
    \item{ordered_effects}{A factored categorical variable to be included.
    Consecutive levels will be tested for significance against each other,
    and the resulting associations will correspond to effect sizes, standard
    errors, and significances of each level versus the previous.}
    \item{strata_effects}{A vector with one variable name to be included as
    the strata variable in case-control studies. Strata cannot be combined
    with random effects.}
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
        fixed_effects = c('diagnosis', 'dysbiosis_state', 'antibiotics',
                        'age', 'reads'),
        random_effects = c('participant_id'),
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
    
    formulas <- maaslin3::maaslin_compute_formula(
        data,
        metadata,
        fixed_effects = c('diagnosis', 'dysbiosis_state', 'antibiotics',
                        'age', 'reads'),
        random_effects = c('participant_id'))
    
    unlink('output')
    logging::logReset()
}
