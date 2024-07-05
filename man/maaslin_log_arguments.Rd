\name{maaslin_log_arguments}
\alias{maaslin_log_arguments}
\title{
Log MaAsLin 3 parameters.
}
\description{
Parse the \code{param_list} provided, check that it is valid for further MaAsLin 3 use, and open a logger to log the parameters.
}
\usage{
maaslin_log_arguments(param_list)
}
\arguments{
    \item{param_list}{
    The parameter list containing the following required named items:
    \describe{
      \item{(1)}{\code{input_data}: A data frame of feature abundances or read counts or a filepath to a tab-delimited file with abundances. It should be formatted with features as columns and samples as rows (or the transpose). The column and row names should be the feature names and sample names respectively.}
      \item{(2)}{\code{input_metadata}: A data frame of per-sample metadata or a filepath to a tab-delimited file with metadata. It should be formatted with variables as columns and samples as rows (or the transpose). The column and row names should be the variable names and sample names respectively.}
      \item{(3)}{\code{output}: The output folder to write results.}
    }
    A model should also be specified using one of the options below. If no formula is specified, at least one of  \code{fixed_effects}, \code{group_effects}, or \code{ordered_effects} should be specified.
    \itemize{
      \item \code{formula}: A formula in \code{lme4} format. Random effects, interactions, and functions of the metadata can be included (note that these functions will be applied after standardization if \code{standardize=TRUE}). Group and ordered variables can be specified as: \code{group(grouping_variable)} and \code{ordered(ordered_variable)}. The other variable options below will not be considered if a formula is set.
      \item \code{fixed_effects}: A vector of variable names to be included as fixed effects.
      \item \code{reference}: For a variable with more than two levels supplied with \code{fixed_effects}, the factor to use as a reference provided as a string of 'variable,reference' semi-colon delimited for multiple variables.
      \item \code{random_effects}: A vector of variable names to be included as random intercepts.
      \item \code{group_effects}: A factored categorical variable to be included for group testing. An ANOVA-style test will be performed to assess whether any of the variable's levels are significant, and no coefficients or individual p-values will be returned.
      \item \code{ordered_effects}: A factored categorical variable to be included. Consecutive levels will be tested for significance against each other, and the resulting associations will correspond to effect sizes, standard errors, and significances of each level versus the previous.
    }
    The parameter list can also contain the following analysis options:
    \itemize{
      \item \code{min_abundance} (default \code{0}): Features with abundances of at least \code{min_abundance} in \code{min_prevalence} of the samples will be included for analysis. The threshold is applied before normalization and transformation. The options \code{min_abundance} and \code{min_prevalence} should usually be 0 since MaAsLin 3 is designed to model sparsity.
      \item \code{min_prevalence} (default \code{0}): See above.
      \item \code{zero_threshold} (default \code{0}): Abundances less than or equal to \code{zero_threshold} will be treated as zeros. This is primarily to be used when the abundance table has likely low-abundance false positives.
      \item \code{min_variance} (default \code{0}): Features with abundance variances less than or equal to \code{min_variance} will be dropped. This is primarily used for dropping features that are entirely zero.
      \item \code{max_significance} (default \code{0.1}): The FDR corrected q-value threshold for significance used in selecting which associations to write as significant and to plot.
      \item \code{normalization} (default \code{TSS}): The normalization to apply to the features before transformation and analysis. The option \code{TSS} is recommended, but \code{CLR}, \code{CSS}, \code{NONE}, and \code{TMM} can also be used.
      \item \code{transform} (default \code{LOG}): The transformation to apply to the features after normalization and before analysis. The option \code{LOG} is recommended, but \code{LOGIT}, \code{AST}, and \code{NONE} can also be used.
      \item \code{correction} (default \code{BH}): The correction to obtain FDR-corrected q-values from raw p-values. Any valid options for \code{p.adjust} can be used.
      \item \code{standardize} (default \code{TRUE}): Whether to apply z-scores to continuous metadata variables so they are on the same scale. This is recommended in order to compare coefficients across metadata variables, but note that functions of the metadata specified in the \code{formula} will apply after standardization.
       }
    The parameter list may also specify how to account for compositionality:
    \itemize{
       \item \code{unscaled_abundance}: A data frame with a single column of absolute abundances or a filepath to such a tab-delimited file. The row names should match the names of the samples in \code{input_data} and \code{input_metadata}. When using spike-ins, the single column should have the same name as one of the features in \code{input_data}, and the \code{unscaled_abundance} should correspond to the absolute quantity of the spike-in. When using total abundance scaling, the single column should have the name 'total', and the \code{unscaled_abundance} should correspond to the total abundance of each sample.
       \item \code{median_comparison_abundance} (default \code{TRUE}): Test abundance coefficients against a null value corresponding to the median coefficient for a metadata variable across the features. This is recommended for relative abundance data but should not be used for absolute abundance data.
       \item \code{median_comparison_prevalence} (default \code{FALSE}): Test prevalence coefficients against a null value corresponding to the median coefficient for a metadata variable across the features. This is only recommended if the analyst is interested in how feature prevalence associations compare to each other or if there is likely strong compositionality-induced sparsity.
       \item \code{median_comparison_abundance_threshold} (default \code{0.25}): Coefficients within \code{median_comparison_abundance_threshold} of the median association will automatically be counted as insignificant (p-value set to 1) since they likely represent compositionality-induced associations. This threshold will be divided by the metadata variable's standard deviation if the metadatum is continuous to ensure the threshold applies to the right scale.
       \item \code{median_comparison_prevalence_threshold} (default \code{0.25}): Same as \code{median_comparison_abundance_threshold} but applied to the prevalence associations.
       \item \code{augment} (default \code{TRUE}): Add extra lowly-weighted 0s and 1s to avoid linear separability.
       }
    The parameter list can also include plotting parameters:
    \itemize{
       \item \code{plot_summary_plot} (default \code{TRUE}): Generate a summary plot of significant associations.
       \item \code{summary_plot_first_n} (default \code{25}): Include the top \code{summary_plot_first_n} features with significant associations.
       \item \code{coef_plot_vars}: Vector of variable names to be used in the coefficient plot section of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form \code{"[variable] [level]"}.
       \item \code{heatmap_vars}: Vector of variable names to be used in the heatmap section of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form \code{"[variable] [level]"}.
       \item \code{plot_associations} (default \code{TRUE}): Whether to generate plots for significant associations.
       \item \code{max_pngs} (default \code{30}): The top \code{max_pngs} associations will be plotted.
    }
    The parameter list can also include the following technical parameters:
    \itemize{
       \item \code{cores} (default \code{1}): How many cores to use when fitting models. (Using multiple cores will likely be faster only for large datasets or complex models.)
       \item \code{save_models} (default \code{FALSE}): Whether to return the fit models and save them to an RData file.
    }
}
}
\value{
    The same \code{param_list} input.
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
  param_list <- maaslin3::maaslin_log_arguments(param_list)
}
