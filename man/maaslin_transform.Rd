\name{maaslin_transform}
\alias{maaslin_transform}
\title{
Transform abundance data for MaAsLin 3 modeling.
}
\description{
Transform the abundance data according to the \code{transform} parameter.
}
\usage{
maaslin_transform(params_and_data_and_formula)
}
\arguments{
    \item{params_and_data_and_formula}{
    The results of \code{maaslin_normalize} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: An input parameter list with \code{transform}.}
      \item{(2)}{\code{data}: The input data dataframe.}
      \item{(3)}{\code{filtered_data}: The abundance data dataframe after filtering.}
      \item{(4)}{\code{filtered_data_norm}: The abundance data after filtering and normalization.}
      \item{(5)}{\code{unfiltered_metadata}: The input metadata dataframe.}
      \item{(6)}{\code{metadata}: The metadata dataframe after filtering.}
      \item{(7)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
    }
}
}
\value{
A list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: An input parameter list with \code{normalization}.}
      \item{(2)}{\code{data}: The input data dataframe.}
      \item{(3)}{\code{filtered_data}: The abundance data dataframe after filtering.}
      \item{(4)}{\code{filtered_data_norm_transformed}: The abundance data after filtering, normalization, and transformation.}
      \item{(5)}{\code{unfiltered_metadata}: The input metadata dataframe.}
      \item{(6)}{\code{metadata}: The metadata dataframe after filtering.}
      \item{(7)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
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
  
  metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
  metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
  metadata$antibiotics <- factor(metadata$antibiotics, levels = c('No', 'Yes'))
  
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
  params_and_data_and_formula <- 
    maaslin3::maaslin_filter_and_standardize(params_and_data_and_formula)
  params_and_data_and_formula <- 
    maaslin3::maaslin_normalize(params_and_data_and_formula)
  params_and_data_and_formula <- 
    maaslin3::maaslin_transform(params_and_data_and_formula)
}
