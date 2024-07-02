\name{maaslin_filter_and_standardize}
\alias{maaslin_filter_and_standardize}
\title{
Filter abundance data and standardize metadata before MaAsLin 3 model fitting.
}
\description{
Set abundances below \code{min_abundance} to zero, remove features without abundances of at least \code{min_abundance} in  \code{min_prevalence} of the samples, and remove features with variances less than or equal to \code{min_variance}. Check that references are set properly if the metadata variables are categorical and provided through \code{fixed_effects}. Standardize the continuous metadata variables as a z-score (subtract the mean, divide by the standard deviation) if \code{standardize} is set.
}
\usage{
maaslin_filter_and_standardize(params_and_data_and_formula)
}
\arguments{
    \item{params_and_data_and_formula}{
    The results of \code{maaslin_check_formula} or \code{maaslin_compute_formula} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: An input parameter list with \code{min_abundance}, \code{min_prevalence}, \code{min_variance}, and \code{standardize}.}
      \item{(2)}{\code{data}: A feature abundance dataframe.}
      \item{(3)}{\code{metadata}: A metadata dataframe.}
      \item{(4)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
    }
}
}
\value{
A list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: The same input parameter list.}
      \item{(2)}{\code{data}: The same input data dataframe.}
      \item{(3)}{\code{filtered_data}: The abundance data dataframe after filtering.}
      \item{(4)}{\code{unfiltered_metadata}: The same input metadata dataframe.}
      \item{(5)}{\code{metadata}: The metadata dataframe after filtering.}
      \item{(6)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
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
  
  # Prepare parameter lists 
  param_list <- list(input_data = taxa_table, 
                     input_metadata = metadata, 
                     output = 'output', 
                     formula = '~ diagnosis + dysbiosisUC + dysbiosisCD + antibiotics + 
                        age + reads_filtered + (1 | subject)', 
                     plot_heatmap = FALSE, 
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
  params_and_data_and_formula <- 
    maaslin3::maaslin_filter_and_standardize(params_and_data_and_formula)
}
