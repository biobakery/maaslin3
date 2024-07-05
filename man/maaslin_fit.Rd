\name{maaslin_fit}
\alias{maaslin_fit}
\title{
Fit MaAsLin 3 models.
}
\description{
Fit the abundance data with abundance and prevalence models to discover feature-metadata associations.
}
\usage{
maaslin_fit(params_and_data_and_formula)
}
\arguments{
    \item{params_and_data_and_formula}{
    The results of \code{maaslin_normalize} or a list containing the following named items:
    \describe{
      \item{(1)}{\code{param_list}: An input parameter list.}
      \item{(2)}{\code{data}: The input data dataframe.}
      \item{(3)}{\code{filtered_data}: The abundance data dataframe after filtering.}
      \item{(4)}{\code{filtered_data_norm_transformed}: The abundance data after filtering, normalization, and transformation.}
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
      \item{(4)}{\code{filtered_data_norm}: The abundance data after filtering and normalization.}
      \item{(5)}{\code{unfiltered_metadata}: The input metadata dataframe.}
      \item{(6)}{\code{metadata}: The metadata dataframe after filtering.}
      \item{(7)}{\code{formula}: A list with the named items \code{formula} as a computed formula object and \code{random_effects_formula} as the random effects component of the formula.}
      \item{(8)}{\code{fit_data_abundance}: The results from the fit abundance models.}
      \item{(9)}{\code{fit_data_prevalence}: The results from the fit prevalence models.}
    }
    
The \code{fit_data_abundance} and \code{fit_data_prevalence} items have the same structure. They are both lists with the following named items:
  \describe{
      \item{(1)}{\code{results}: A results table with the modeled associations (see below).}
      \item{(2)}{\code{residuals}: A features (rows) by samples (columns) dataframe of residuals from the models.}
      \item{(3)}{\code{fitted}: A features (rows) by samples (columns) dataframe of fitted values from the models.}
      \item{(4)}{\code{ranef}: A features (rows) by random effect (columns) dataframe of random effects from the models. If multiple random effects are specified, this is a dataframe of dataframes.}
      \item{(5)}{\code{fits}: If \code{save_models=TRUE}, this is a list of the fit models.}
    }
    
The \code{results} tables contain the following columns for each association (row):
  \describe{
      \item{(1)}{\code{feature}: The feature involved in the association.}
      \item{(2)}{\code{metadata}: The metadata variable involved in the association.}
      \item{(3)}{\code{value}: The value of the metadata variable: the metadata variable itself if continuous or the level if categorical.}
      \item{(4)}{\code{name}: The name of the model component involved in the association: the metadata variable itself if continuous or a concatenated version of the metadata variable and level if categorical.}
      \item{(5)}{\code{coef}: The coefficient of the association: the slope coefficient in the abundance model and the change in log odds in the prevalence model.}
      \item{(6)}{\code{stderr}: The standard error of the coefficient.}
      \item{(7)}{\code{pval_individual}: The (uncorrected) p-value of the association.}
      \item{(8)}{\code{error}: Any error produced by the model during fitting. NA otherwise.}
      \item{(9)}{\code{qval_individual}: The FDR corrected q-value of the association. FDR correction is performed over all associations in the abundance and prevalence modeling separately.}
      \item{(10)}{\code{model}: \code{LM} for the abundance models and \code{logistic} for the prevalence models.}
      \item{(11)}{\code{N}: The number of data points for the association's feature.}
      \item{(12)}{\code{N.not.zero}: The number of non-zero data points for the association's feature.}
      \item{(13)}{\code{pval_joint}: The p-value of the overall association (combining abundance and prevalence) by taking the minimum of the abundance and logistic p-values and applying the Beta(1,2) CDF. These will be the same in the abundance and prevalence results for an association.}
      \item{(14)}{\code{qval_joint}: The FDR corrected q-value of the association. FDR correction is performed over all joint p-values.}
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
}
