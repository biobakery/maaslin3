\name{maaslin_check_formula}
\alias{maaslin_check_formula}
\title{
Check a MaAsLin 3 formula to ensure a proper MaAsLin 3 run.
}
\description{
Ensure that the formula provided as part of the \code{param_list} is valid. Only one of \code{maaslin_compute_formula} or \code{maaslin_check_formula} should be used.
}
\usage{
maaslin_check_formula(params_and_data)
}
\arguments{
    \item{params_and_data}{
    The results of \code{maaslin_reorder_data} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: The parameter list with \code{unscaled_abundance} trimmed of samples with missing abundances or metadata and a \code{formula} specified.}
      \item{(2)}{\code{data}: A data frame of feature abundances trimmed of samples with missing abundances or metadata.}
      \item{(3)}{\code{metadata}: A data frame of metadata trimmed of samples with missing abundances or metadata.}
    }
}
}
\value{
A list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: The same input parameter list.}
      \item{(2)}{\code{data}: The same input data data frame.}
      \item{(3)}{\code{metadata}: The same input metadata data frame.}
      \item{(4)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
    }
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
  rownames(taxa_table) <- taxa_table$ID; taxa_table$ID <- NULL
  
  # Read metadata table
  metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")
  metadata <- read.csv(metadata_name, sep = '\t')
  rownames(metadata) <- metadata$ID; metadata$ID <- NULL
  
  metadata$diagnosis <- 
    factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
  metadata$dysbiosis_state <- 
    factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
  metadata$antibiotics <- 
    factor(metadata$antibiotics, levels = c('No', 'Yes'))
  
  # Prepare parameter lists 
  param_list <- list(input_data = taxa_table, 
                     input_metadata = metadata, 
                     output = 'output', 
                     formula = '~ diagnosis + dysbiosis_state + 
                        antibiotics + age + reads', 
                     plot_summary_plot = FALSE, 
                     plot_associations = FALSE)
  
  #Run MaAsLin3
  param_list <- 
    maaslin3::maaslin_log_arguments(param_list)
  params_and_data <- 
    maaslin3::maaslin_read_data(param_list)
  params_and_data <- 
    maaslin3::maaslin_reorder_data(params_and_data)
  params_and_data_and_formula <- 
    maaslin3::maaslin_check_formula(params_and_data)
}
