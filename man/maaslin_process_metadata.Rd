\name{maaslin_process_metadata}
\alias{maaslin_process_metadata}
\title{
Process metadata before MaAsLin 3 model fitting.
}
\description{
Check that references are set properly if the metadata variables are categorical and provided through \code{fixed_effects}. Standardize the continuous metadata variables as a z-score (subtract the mean, divide by the standard deviation) if \code{standardize} is set.
}
\usage{
maaslin_process_metadata(metadata,
                         formula = NULL,
                         fixed_effects = NULL,
                         reference = NULL,
                         feature_specific_covariate_name = NULL,
                         standardize = TRUE)
}
\arguments{
    \item{metadata}{A data frame of per-sample metadata. It should be formatted with variables as columns and samples as rows. The column and row names should be the variable names and sample names respectively.}
    \item{formula}{A formula in \code{lme4} format. Random effects, interactions, and functions of the metadata can be included (note that these functions will be applied after standardization if \code{standardize=TRUE}). Group, ordered, and strata variables can be specified as: \code{group(grouping_variable)}, \code{ordered(ordered_variable)} and \code{strata(strata_variable)}. The other variable options below will not be considered if a formula is set.}
    \item{fixed_effects}{A vector of variable names to be included as fixed effects.}
    \item{reference}{For a variable with more than two levels supplied with \code{fixed_effects}, the factor to use as a reference provided as a string of 'variable,reference' semi-colon delimited for multiple variables.}
    \item{feature_specific_covariate_name}{The name for the feature-specific covariates when fitting the models.}
    \item{standardize}{Whether to apply z-scores to continuous metadata variables so they are on the same scale. This is recommended in order to compare coefficients across metadata variables, but note that functions of the metadata specified in the \code{formula} will apply after standardization.}
}
\value{
The processed metadata.
}
\author{
    William Nickols<willnickols@g.harvard.edu>,\cr
    Jacob Nearing<nearing@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{
  # Read features table 
  taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package = "maaslin3")
  taxa_table <- read.csv(taxa_table_name, sep = '\t')

  # Read metadata table
  metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")
  metadata <- read.csv(metadata_name, sep = '\t')

  metadata$diagnosis <- 
    factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
  metadata$dysbiosis_state <- 
    factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
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
                                      
  standardized_metadata = maaslin3::maaslin_process_metadata(
    metadata,
    formula = formula)

  unlink('output')
  logging::logReset()
}
