\name{maaslin_normalize}
\alias{maaslin_normalize}
\title{
Normalize abundance data for MaAsLin 3 model fitting.
}
\description{
Normalize the abundance data according to the \code{normalization} parameter. If \code{unscaled_abundance} is specified, compute the absolute abundances.
}
\usage{
maaslin_normalize(data,
                  output,
                  zero_threshold = 0,
                  normalization = 'TSS',
                  unscaled_abundance = NULL)
}
\arguments{
    \item{data}{A data frame of feature abundances. It should be formatted with features as columns and samples as rows. The column and row names should be the feature names and sample names respectively.}
    \item{output}{The output folder to write results.}
    \item{zero_threshold}{Abundances less than or equal to \code{zero_threshold} will be treated as zeros. This is primarily to be used when the abundance table has likely low-abundance false positives.}
    \item{normalization}{The normalization to apply to the features before transformation and analysis. The option \code{TSS} (total sum scaling) is recommended, but \code{CLR} (centered log ratio) and \code{NONE} can also be used.}
    \item{unscaled_abundance}{A data frame with a single column of absolute abundances. The row names should match the names of the samples in \code{input_data} and \code{input_metadata}. When using spike-ins, the single column should have the same name as one of the features in \code{input_data}, and the \code{unscaled_abundance} should correspond to the absolute quantity of the spike-in. When using total abundance scaling, the single column should have the name 'total', and the \code{unscaled_abundance} should correspond to the total abundance of each sample.}
}
\value{
A dataframe of normalized features (features are columns; samples are rows).
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
      
  normalized_data = maaslin3::maaslin_normalize(data,
                                      output = 'output')

  unlink('output')
  logging::logReset()
}
