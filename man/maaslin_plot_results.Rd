\name{maaslin_plot_results}
\alias{maaslin_plot_results}
\title{
Plot the results from a MaAsLin 3 run.
}
\description{
Two types of plots are generated. First, the summary plot contains sorted per-feature coefficients plotted with their standard errors for key variables and a heatmap summarizing the remaining variables. Second, for significant features, association plots (scatterplots, boxplots, or tables depending on the association) are generated to visualize and verify the model fits. The data are shown with their transformed values in the association plots since this is the scale on which the models are fit.
}
\usage{
maaslin_plot_results(params_data_formula_fit)
}
\arguments{
    \item{params_data_formula_fit}{
    The results of \code{maaslin_fit} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: An input parameter list. The parameters \code{plot_summary_plot} and \code{plot_associations} determine whether to plot the summary plot and individual associations respectively. The parameter \code{summary_plot_first_n} determines how many features to include in the summary plot. The parameter \code{max_pngs} determines how many individual associations will be plotted. The parameters \code{heatmap_vars} and \code{coef_plot_vars} are vectors of names that specify which variables should be included in the heatmap and coefficient plot portion of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form \code{"[variable] [level]"}.}
      \item{(2)}{\code{data}: The input data dataframe.}
      \item{(3)}{\code{filtered_data}: The abundance data dataframe after filtering.}
      \item{(4)}{\code{filtered_data_norm}: The abundance data after filtering and normalization.}
      \item{(5)}{\code{unfiltered_metadata}: The input metadata dataframe.}
      \item{(6)}{\code{metadata}: The metadata dataframe after filtering.}
      \item{(7)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
      \item{(8)}{\code{fit_data_abundance}: The results from the fit abundance models.}
      \item{(9)}{\code{fit_data_prevalence}: The results from the fit prevalence models.}
    }
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
  maaslin3::maaslin_plot_results(params_data_formula_fit)
}
