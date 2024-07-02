\name{maaslin_reorder_data}
\alias{maaslin_reorder_data}
\title{
Reorder the abundance data and metadata.
}
\description{
Reorder the abundance data and metadata to ensure samples are rows and remove any samples without abundances or metadata.
}
\usage{
maaslin_reorder_data(params_and_data)
}
\arguments{
    \item{params_and_data}{
    The results of \code{maaslin_read_data} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: The parameter list with \code{unscaled_abundance} replaced by a data frame of absolute abundances if specified.}
      \item{(2)}{\code{data}: A data frame of feature abundances.}
      \item{(3)}{\code{metadata}: A data frame of metadata.}
    }
}
}
\value{
    A list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: The parameter list with \code{unscaled_abundance} trimmed of samples with missing abundances or metadata.}
      \item{(2)}{\code{data}: A data frame of feature abundances trimmed of samples with missing abundances or metadata.}
      \item{(3)}{\code{metadata}: A data frame of metadata trimmed of samples with missing abundances or metadata.}
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
                     plot_summary_plot = FALSE, 
                     plot_associations = FALSE)
  
  #Run MaAsLin3
  param_list <- 
    maaslin3::maaslin_log_arguments(param_list)
  params_and_data <- 
    maaslin3::maaslin_read_data(param_list)
  params_and_data <- 
    maaslin3::maaslin_reorder_data(params_and_data)
}
