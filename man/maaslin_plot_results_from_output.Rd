\name{maaslin_plot_results_from_output}
\alias{maaslin_plot_results_from_output}
\title{
Plot the results from a MaAsLin 3 run.
}
\description{
Two types of plots are generated. First, the summary plot contains sorted per-feature coefficients plotted with their standard errors for key variables and a heatmap summarizing the remaining variables. Second, for significant features, association plots (scatterplots, boxplots, or tables depending on the association) are generated to visualize and verify the model fits. The data are shown with their transformed values in the association plots since this is the scale on which the models are fit. In comparison to \code{maaslin_plot_results} that needs the entire \code{maaslin_fit} list, only the parameter list and an outputs directory containing a completed run are needed for \code{maaslin_plot_results_from_output}.
}
\usage{
maaslin_plot_results_from_output(param_list)
}
\arguments{
    \item{param_list}{
    The \code{param_list} that was used to fit the MaAsLin 3 model originally (see \code{\link{maaslin_parse_param_list}}). The relevant plotting variables can be changed from their original settings: \code{plot_summary_plot}, \code{summary_plot_first_n}, \code{heatmap_vars}, \code{plot_associations}, \code{coef_plot_vars}, \code{max_pngs}. The function will automatically look for the \code{all_results.tsv} and \code{filtered_data_norm_transformed.tsv} files in the \code{output} directory.
  }
}
\value{
Results will be written to the \code{figures} folder within the folder \code{output}. The list of individual association plots is returned if \code{plot_associations=TRUE}. In the heatmap of the summary plot, one star corresponds to the user-set \code{max_significance} and two stars corresponds to the user-set \code{max_signifiance} divided by 10.
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
                     formula = '~ diagnosis + dysbiosis_state + antibiotics + 
                        age + reads', 
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
  params_data_formula_fit <- 
    maaslin3::maaslin_fit(params_and_data_and_formula)
  maaslin3::maaslin_plot_results_from_output(param_list)
}
