# MaAsLin 3 #

[**MaAsLin 3**](http://huttenhower.sph.harvard.edu/MaAsLin3) is the next generation of **M**a**A**s**L**in (**M**icrobiome **M**ultivariable **A**ssociations with **L**inear Models). This comprehensive R package efficiently determines multivariable associations between clinical metadata and microbial meta-omics features. Relative to MaAsLin 2, MaAsLin 3 introduces the ability to quantify and test for **both abundance and prevalence associations** while **accounting for compositionality**. By incorporating generalized linear models, MaAsLin 3 accomodates most modern epidemiological study designs including cross-sectional and longitudinal studies.

If you use the MaAsLin 3 software, please cite our manuscript: 

> William A. Nickols, Jacob T. Nearing, Kelsey N. Thompson, Jiaxian Shen, Curtis Huttenhower 
MaAsLin 3: Refining and extending generalized multivariate linear models for meta-omic association discovery. (In progress).

### Support ###
Check out the [MaAsLin 3 tutorial](https://github.com/biobakery/biobakery/wiki/MaAsLin3) for an overview of analysis options and some example runs.

If you have questions, please direct them to the [MaAsLin 3 Forum](https://forum.biobakery.org/c/Downstream-analysis-and-statistics/maaslin). 

--------------------------------------------

## Contents ##
- [Introduction](#introduction)
  - [Support](#support)
- [Contents](#contents)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running MaAsLin 3](#running-maaslin-3)
  - [Input data](#input-data)
  - [Output files](#output-files)
  - [Run a demo](#run-a-demo)
    - [In R](#in-r)
    - [Command line](#command-line)
  - [Options](#options)
    - [Required parameters](#required-parameters)
    - [Model formula](#model-formula)
    - [Analysis options](#analysis-options)
    - [Compositionality corrections](#compositionality-corrections)
    - [Plotting parameters](#plotting-parameters)
    - [Technical parameters](#technical-parameters)
- [Troubleshooting](#troubleshooting)

## Requirements ##

MaAsLin3 is an R package that can be run on the command line or as an R function.

## Installation ##

#### Install using GitHub and devtools

The latest development version of MaAsLin 3 can be installed from GitHub using the `devtools` package.

```
# Install devtools if not present
if (!require('devtools', character.only = TRUE)) {
  install.packages('devtools')
}

# Install MaAsLin 3
library("devtools")
install_github("biobakery/maaslin3")
library(maaslin3)
```

## Running MaAsLin 3 ##

MaAsLin3 can be run from the command line or as an R function. Both methods require the same arguments, have the same options, and use the same default settings. To run MaAsLin 3, the user must supply a table of per-sample feature abundances (with zeros still included), a table of per-sample metadata, and a model specifying how the metadata should relate to the feature prevalence (how likely the feature is to be present or absent) and abundance (how much of the feature is there if it's there). MaAsLin 3 will return a table of associations including an effect size and p-value for each feature-metadatum association and a folder of visuals including a summary plot and diagnostic plots for significant associations.

### Input data ###

MaAsLin3 requires two input files.

1. Feature abundance data frame
    * Formatted with features as columns and samples as rows.
    * The transpose of this format is also okay.
    * Possible features include taxonomy or genes. These can be relative abundances or counts.
    * This can be a filepath to a tab-delimited file.
2. Metadata data frame
    * Formatted with variables as columns and samples as rows.
    * The transpose of this format is also okay.
    * Possible metadata include gender or age.
    * This can be a filepath to a tab-delimited file.

The data file can contain samples not included in the metadata file (along with the reverse case). For both cases, those samples not included in both files will be removed from the analysis. Also, the samples do not need to be in the same order in the two files.

To run MaAsLin 3, it is also necessary to specify a model. The model can come from a formula or vectors of terms. In either case, variable names should not have spaces or unusual characters.

* *Formula*: The `formula` parameter should be set to any formula that satisfies the `lme4` specifications: fixed effects, random effects, interaction terms, polynomial terms, and more can all be included. If categorical variables are included as fixed effects, each level will be tested against the first factor level. In addition, ordered predictors and group predictors can be included by including `group(variable_name)` and `ordered(variable_name)` in the formula.
* *Vectors*: Alternatively, a vector of variable names can be supplied to the parameters `fixed_effects`, `random_effects`, `group_effects`, and `ordered_effects`. Categorical variables should either be ordered as factors beforehand, or `reference` should be provided as a string of 'variable,reference' semi-colon delimited for multiple variables (e.g., `variable_1,reference_1;variable_2,reference_2`). NOTE: adding a space between the variable and level might result in the wrong reference level being used.

**Because MaAsLin 3 identifies prevalence (presence/absence) associations, sample read depth (number of reads) should be included as a covariate if available. Deeper sequencing will likely increase feature detection in a way that could spuriously correlate with metadata of interest when read depth is not included in the model.**


### Output files ###

MaAsLin 3 generates two types of output files explained below: data and visualizations. In addition, the object returned from `maaslin3` contains all the data and results (see `?maaslin_fit`).

1. Data output files
    * ``all_results.tsv``
        * This file contains all results ordered by increasing q-value.
        * `feature` and `metadata` are the feature and metadata names.
        * `value` and `name` are the value of the metadata and variable name from the model.
        * `coef` and `stderr` are the fit coefficient and standard error from the model. In abundance models, a one-unit change in the metadatum variable corresponds to a $2^{\textrm{coef}}$ fold change in the relative abundance of the feature. In prevalence models, a one-unit change in the metadatum variable corresponds to a $\textrm{coef}$ change in the log-odds of a feature being present.
        * `pval_individual` is the p-value of the individual association.
        * `error` lists any errors from the model fitting.
        * `qval_individual` is the false discovery rate (FDR) corrected q-value of the individual association. FDR correction is performed over all non-NA p-values in the abundance and prevalence modeling separately.
        * `model` specifies whether the association is abundance or prevalence.
          * `N` and `N.not.zero` are the total number of data points and the total number of non-zero data points for the feature.
        * `pval_joint` and `qval_joint` are the p-value and q-value of the joint prevalence and abundance association. The p-value comes from plugging in the minimum of the association's abundance and prevalence p-values into the Beta(1,2) CDF. It is interpreted as the probability that either the abundance or prevalence association would be as extreme as observed if there was neither an abundance nor prevalence association between the feature and metadatum.
    * ``significant_results.tsv``
        * This file is a subset of the results in the first file.
        * It only includes associations with q-values <= to the threshold and no errors.
    * ``features``
        * This folder includes the filtered, normalized, and transformed versions of the input feature table.
        * These steps are performed sequentially in the above order.
        * If an option is set such that a step does not change the data, the resulting table will still be output.
    * ``models_LM.rds`` and ``models_logistic.rds``
        * These files contain a list with every model fit object (`LM` for linear models, `logistic` for logistic models).
        * It will only be generated if `save_models` is set to TRUE.
    * ``residuals_LM.rds`` and ``residuals_logstic.rds``
        * These files contain a data frame with residuals for each feature.
    * ``fitted_LM.rds`` and ``fitted_logistic.rds``
        * These files contain a data frame with fitted values for each feature.
    * ``ranef_LM.rds`` and ``ranef_logistic.rds``
        * These files contain a data frame with extracted random effects for each feature (when random effects are specified).
    * ``maaslin3.log``
        * This file contains all log information for the run.
        * It includes all settings, warnings, errors, and steps run.
2. Visualization output files
    * ``summary_plot.pdf``
        * This file contain a combined coefficient plot and heatmap of the most significant associations. In the heatmap, one star indicates the individual q-value is below the parameter `max_significance`, and two stars indicate the individual q-value is below `max_significance` divided by 10.
    * ``association_plots/[metadatum]_[feature]_[association].pdf``
        * A plot is generated for each significant association up to `max_pngs`.
        * Scatter plots are used for continuous metadata abundance associations.
        * Box plots are used for categorical data abundance associations.
        * Box plots are used for continuous data prevalence associations.
        * Grids are used for categorical data prevalence associations.
        * Data points plotted are after filtering, normalization, and transformation so that the scale in the plot is the scale that was used in fitting.

At the top right of each association plot is the name of the significant association in the results file, the FDR corrected q-value for the individual association, the number of samples in the dataset, and the number of samples with non-zero abundances for the feature. In the plots with categorical metadata variables, the reference category is on the left, and the significant q-values and coefficients in the top right are in the order of the values specified above. Because the displayed coefficients correspond to the full fit model with (possibly) scaled metadata variables, the marginal association plotted might not match the coefficient displayed. However, the plots are intended to provide an interpretable visual while usually agreeing with the full model.

#### Diagnostics

There are a few common issues to check for in the results:

1. When warnings or errors are thrown during the fitting process, they are recorded in the `error` column of `all_results.tsv`. Often, these indicate model fitting failures or poor fits that should not be trusted, but sometimes the warnings can be benign, and the model fit might still be reasonable. Users should check associations of interest if they produce errors.
2. Despite the error checking, significant results could still result from poor model fits. These can usually be diagnosed with the visuals in the `association_plots` directory. 
    * Any significant abundance associations with a categorical variable should usually have **at least 10 observations in each category**.
    * Significant prevalence associations with categorical variables should also have **at least 10 samples in which the feature was present and at least 10 samples in which it was absent for each significant category**. 
    * Significant abundance associations with continuous metadata should be checked visually for influential outliers.

### Run a demo ###

Example input files can be found in the ``inst/extdata`` folder of the MaAsLin 3 source. The files provided were generated from the Human Microbiome Project 2 (HMP2) data which can be downloaded from https://ibdmdb.org/.

* ``HMP2_taxonomy.tsv``: a tab-delimited file with samples as rows and species as columns. It is a subset of the full HMP2 taxonomy that includes just some of the the species abundances.
* ``HMP2_metadata.tsv``: a tab-delimited file with samples as rows and metadata as columns. It is a subset of the full HMP2 metadata that includes just some of the fields.

#### In R ####

The following code identifies associations between patient metadata and microbial species in the HMP2 cohort.

```
# Read abundance table
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

# Prepare parameter lists 
param_list <- list(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = 'hmp2_output', 
                   formula = '~ diagnosis + dysbiosis_state + antibiotics + age + reads',
                   normalization = 'TSS', 
                   transform = 'LOG', 
                   augment = TRUE, 
                   standardize = TRUE,
                   max_significance = 0.1, 
                   median_comparison_abundance = TRUE, 
                   median_comparison_prevalence = FALSE, 
                   max_pngs = 100,
                   cores = 1)

fit_out <- maaslin3(param_list)
```

The `maaslin3` wrapper function runs all steps at once, but equivalent results would be produced from running one step at a time:

```
param_list <- maaslin_log_arguments(param_list)
params_and_data <- maaslin_read_data(param_list)
params_and_data <- maaslin_reorder_data(params_and_data)
params_and_data_and_formula <- maaslin_check_formula(params_and_data)
params_and_data_and_formula <- maaslin_filter_and_standardize(params_and_data_and_formula)
params_and_data_and_formula <- maaslin_normalize(params_and_data_and_formula)
params_and_data_and_formula <- maaslin_transform(params_and_data_and_formula)
fit_out <- maaslin_fit(params_and_data_and_formula)
maaslin_write_results(fit_out)
invisible(maaslin_plot_results(fit_out)) # invisible to avoid dumping plots into Rmd
```

#### Command line ####

MaAsLin 3 can also be run with a command line interface. For example, the first HMP2 analysis can be performed with:

```
./R/maaslin3.R inst/extdata/HMP2_taxonomy.tsv inst/extdata/HMP2_metadata.tsv command_line_output --formula='~ diagnosis + dysbiosis_state + antibiotics + age + reads' --reference='diagnosis,nonIBD;dysbiosis_state,none;antibiotics,No'
```

* Make sure to provide the full path to the MaAsLin3 executable (i.e. `./R/maaslin3.R`).
* In the demo command:
    * ``inst/extdata/HMP2_taxonomy.tsv`` is the path to your data (or features) file
    * ``inst/extdata/HMP2_metadata.tsv`` is the path to your metadata file
    * ``command_line_output`` is the path to the folder to write the output

### Options ###

From the command line, the following command will print the list of MaAsLin 3 options and default settings:

```
$ ./R/maaslin3.R --help
```

When running MaAsLin 3 in R, the manual page for each function (e.g., `?maaslin3`) will show the available options and default settings. For both, the options and settings are as follows:

#### Required parameters ####

* `input_data`: A data frame of feature abundances or read counts or a filepath to a tab-delimited file with abundances. It should be formatted with features as columns and samples as rows (or the transpose). The column and row names should be the feature names and sample names respectively.
* `input_metadata`: A data frame of per-sample metadata or a filepath to a tab-delimited file with metadata. It should be formatted with variables as columns and samples as rows (or the transpose). The column and row names should be the variable names and sample names respectively.
* `output`: The output folder to write results.

#### Model formula ####

* `formula`: A formula in `lme4` format. Random effects, interactions, and functions of the metadata can be included (note that these functions will be applied after standardization if `standardize = TRUE`). Group, ordered, and strata variables can be specified as: `group(grouping_variable)`, `ordered(ordered_variable)` and `strata(strata_variable)`.  Ordered and group predictors should stand alone in the formula (i.e., no group predictors in random effects). The other variable options below will not be considered if a formula is set.
* `fixed_effects`: A vector of variable names to be included as fixed effects.
  * Fixed effects models are fit with `lm` (linear) or `glm` (logistic).
* `reference`: For a variable with more than two levels supplied with `fixed_effects`, the factor to use as a reference provided as a string of 'variable,reference' semi-colon delimited for multiple variables.
* `random_effects`: A vector of variable names to be included as random intercepts. **Random intercept models may produce poor model fits when there are fewer than 5 observations per group.** In these scenarios, per-group fixed effects should be used and subsequently filtered out.
  * Random effects models are fit with `lmer` (linear) and `glmer` (logistic), and the significance tests come from `lmerTest` and `glmer` respectively.
* `group_effects`: A factored categorical variable to be included for group testing. An ANOVA-style test will be performed to assess whether any of the variable's levels are significant, and no coefficients or individual p-values will be returned.
  * Tests are performed with the `anova` function's `LRT` option (logistic fixed and mixed effects), the `anova` function's F test (linear fixed effects), or `lmerTest::contest` (linear mixed effects).
* `ordered_effects`: A factored categorical variable to be included. Consecutive levels will be tested for significance against each other with contrast tests, and the resulting associations will correspond to effect sizes, standard errors, and significances of each level versus the previous.
  * Contrast tests are performed with `multcomp::glht` (fixed effects and logistic mixed effects) and `lmerTest::contest` (linear mixed effects).
* `strata_effects`: A single grouping variable to be included in matched case-control studies. If a strata variable is included, no random effects can be included. When a strata variable is included, a conditional logistic regression will be run to account for the strata. The abundance model will be run with a random intercept in place of the strata. Strata can include more than two observations per group. Only variables that differ within the groups can be tested.

#### Feature specific covariates ####
Particularly for use in metatranscriptomics workflows, a table of feature-specific covariates can be included. A feature's covariates will be included like a fixed effect metadatum when fitting the model for that feature. The covariate's name does not need to be included in the formula.

* `feature_specific_covariate`: A table of feature-specific covariates or a filepath to a tab-delimited file with feature-specific covariates. It should be formatted with features as columns and samples as rows (or the transpose). The row names and column names should be the same as those of the `input_data`: the column and row names should be the feature names and sample names respectively.
* `feature_specific_covariate_name`: The name for the feature-specific covariates when fitting the models.
* `feature_specific_covariate_record`: Whether to keep the feature-specific covariates in the outputs when calculating p-values, writing results, and displaying plots.


#### Analysis options ####

* `min_abundance` (default `0`): Features with abundances of at least `min_abundance` in `min_prevalence` of the samples will be included for analysis. The threshold is applied before normalization and transformation. **The options `min_abundance` and `min_prevalence` should usually be 0 since MaAsLin 3 is designed to model sparsity.**
* `min_prevalence` (default `0`): See above.
* `zero_threshold` (default `0`): Abundances less than or equal to `zero_threshold` will be treated as zeros. This is primarily to be used when the abundance table has likely low-abundance false positives.
* `min_variance` (default `0`): Features with abundance variances less than or equal to `min_variance` will be dropped. This is primarily used for dropping features that are entirely zero.
* `max_significance` (default `0.1`): The FDR corrected q-value threshold for significance used in selecting which associations to write as significant and to plot.
* `normalization` (default `TSS`): The normalization to apply to the features before transformation and analysis. The option `TSS` is recommended, but `CLR`, `CSS`, `NONE`, and `TMM` can also be used.
* `transform` (default `LOG`): The transformation to apply to the features after normalization and before analysis. The option `LOG` is recommended, but `LOGIT`, `AST`, and `NONE` can also be used.
* `correction` (default `BH`): The correction to obtain FDR-corrected q-values from raw p-values. Any valid options for `p.adjust` can be used.
* `standardize` (default `TRUE`): Whether to apply z-scores to continuous metadata variables so they are on the same scale. This is recommended in order to compare coefficients across metadata variables, but note that functions of the metadata specified in the `formula` will apply after standardization.

#### Compositionality corrections ####

##### Absolute abundance

Most microbiome methodologies have historically focused on relative abundances (proportions out of 1). However, some experimental protocols can enable estimation of absolute abundances (cell count/concentration). MaAsLin 3 can be used with two types of absolute abundance estimation: spike-ins and total abundance scaling. In a spike-in procedure, a small, known quantity of a microbe that otherwise would not be present in the sample is added, and the sequencing procedure is carried out as usual. Then, the absolute abundance of a microbe already in the community is estimated as:
$$\textrm{Absolute abundance other microbe}=\frac{\textrm{Relative abundance other microbe}}{\textrm{Relative abundance spike-in microbe}}\cdot (\textrm{Absolute abundance spike-in microbe})$$
Alternatively, the total microbial abundance of a sample can be determined (e.g., with qPCR of a marker gene or by cell counting). Then, the absolute abundance of a microbe in the community is estimated as:
$$\textrm{Absolute abundance microbe}=(\textrm{Total absolute abundance})\cdot(\textrm{Relative abundance microbe})$$

##### Compositionality corrections continued

* `unscaled_abundance`: A data frame with a single column of absolute abundances or a filepath to such a tab-delimited file. The row names should match the names of the samples in `input_data` and `input_metadata`. When using spike-ins, the single column should have the same name as one of the features in `input_data`, and the values should correspond to the absolute quantity of the spike-in. For example, if the same spike-in quantity is used in each sample, the entire column can be set to 1. When using total abundance scaling, the single column should have the name 'total', and the values should correspond to the total abundance of each sample. In both cases, `median_comparison_abundance` should be set to `FALSE` since the spike-in or total abundance normalization accounts for compositionality.

  Alternatively, if the `input_data` abundances have already been scaled to be absolute abundances, the user should set `normalization = NONE` and `median_comparison_abundance = FALSE` and not include anything for `unscaled_abundance`. Then, the absolute abundances will be log transformed, and models will be fit on those values directly.
  
##### Median comparisons

When `median_comparison_abundance` or `median_comparison_prevalence` are on, the coefficients for a metadatum will be tested against the median coefficient for that metadatum (median across the features). Otherwise, the coefficients will be tested against 0. For abundance associations, this is designed to account for compositionality, the idea that if only one feature has a positive association with a metadatum on the absolute scale (cell count), the other features will have apparent negative associations with that metadatum on the relative scale (proportion of the community) because relative abundances must sum to 1. More generally, associations on the relative scale are not necessarily the same as the associations on the absolute scale in magnitude or sign, so **testing against zero on the relative scale is not equivalent to testing against zero on the absolute scale**. When testing associations on the relative scale, the coefficients should be tested against 0 (median comparison off). However, since these tests do not correspond to tests for associations on the absolute scale, we instead use a test against the median, which can enable some inference on the absolute scale. There are two interpretations of this test for absolute abundance associations:

1. In linear models, if two features' associations with a particular metadatum on the log *absolute* scale differ by some value $d$, the features' associations with that metadatum on the log *relative* scale (total-sum scaling) will also differ by $d$. That is, the absolute and relative coefficients for a particular feature-metadatum association are different, but **the ordering of the relative coefficients is the same as the ordering of the absolute coefficients for a metadatum, and the difference between two coefficients on the relative scale is the same as the difference between the corresponding coefficients on the absolute scale**. Therefore, the test against the relative coefficient median can always be interpreted as "a test of whether a particular association is significantly different from the typical (median) association on the absolute scale."
2. Under the assumption that at least half the features are not changing on the absolute scale, the median true absolute coefficient is 0, so this can be interpreted as a test of whether the feature has a non-zero association on the absolute scale.

By contrast, sparsity should be less affected by compositionality since a feature should still be present even if another increases or decreases in abundance. (Note that, because the read depth is finite, this might not always be true in practice.) Therefore, `median_comparison_prevalence` is off by default, but it can be turned on if the user is interested in testing whether a particular prevalence association is significantly different from the typical prevalence association.

In both cases, if the tested coefficient is within `median_comparison_[abundance/prevalence]_threshold` of the median, it will automatically receive a p-value of 1. This is based on the idea that the association might be statistically significantly different but not substantially different from the median and therefore is likely still a result of compositionality.

To conclude:

* `median_comparison_abundance` is `TRUE` by default and should be used to make inference on the absolute scale when using relative abundance data. When `median_comparison_abundance` is `TRUE`, only the p-values and q-values change; the coefficients returned are still the relative abundance coefficients.
* `median_comparison_abundance` should be `FALSE` when (1) testing for significant relative associations, (2) testing for absolute abundance associations under the assumption that the total absolute abundance is not changing, or (3) testing for significant absolute associations when supplying spike-in or total abundances with `unscaled_abundance`.
* `median_comparison_prevalence` is `FALSE` by default.

##### Compositionality corrections continued
  
* `median_comparison_abundance` (default `TRUE`): Test abundance coefficients against a null value corresponding to the median coefficient for a metadata variable across the features. Otherwise, test against 0. This is recommended for relative abundance data but should not be used for absolute abundance data.
* `median_comparison_prevalence` (default `FALSE`): Test prevalence coefficients against a null value corresponding to the median coefficient for a metadata variable across the features. Otherwise, test against 0. This is only recommended if the analyst is interested in how feature prevalence associations compare to each other or if there is likely strong compositionality-induced sparsity.
* `median_comparison_abundance_threshold` (default `0.25`): Coefficients within `median_comparison_abundance_threshold` of the median association will automatically be counted as insignificant (p-value set to 1) since they likely represent compositionality-induced associations. This threshold will be divided by the metadata variable's standard deviation if the metadatum is continuous to ensure the threshold applies to the right scale.
* `median_comparison_prevalence_threshold` (default `0.25`): Same as `median_comparison_abundance_threshold` but applied to the prevalence associations.
* `augment` (default `TRUE`): To avoid linear separability in the logistic regression, at each input data point, add an extra 0 and an extra 1 observation weighted as the number of predictors divided by two times the number of data points. This is almost always recommended to avoid linear separability while having a minor effect on fit coefficients otherwise.

#### Plotting parameters ####

* `plot_summary_plot` (default `TRUE`): Generate a summary plot of significant associations.
* `summary_plot_first_n` (default `25`): Include the top `summary_plot_first_n` features with significant associations.
* `coef_plot_vars`: Vector of variable names to be used in the coefficient plot section of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form `"[variable] [level]"`.

  By default, the (up to) two metadata variables with the most significant associations will be plotted in the coefficient plot, and the rest will be plotted in the heatmap. Because predicting the output variable names can be tricky, it is recommended to first run `maaslin3` without setting `coef_plot_vars` or `heatmap_vars`, look at the names of the variables in the summary plot, and then rerun with `maaslin_plot_results_from_output` after updating `coef_plot_vars` and `heatmap_vars` in `param_list` with the desired variables.
* `heatmap_vars`: Vector of variable names to be used in the heatmap section of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form `"[variable] [level]"`.
* `plot_associations` (default `TRUE`): Whether to generate plots for significant associations.
* `max_pngs` (default `30`): The top `max_pngs` significant associations will be plotted.

#### Technical parameters ####

* `cores` (default `1`): How many cores to use when fitting models. (Using multiple cores will likely be faster only for large datasets or complex models.)
* `save_models` (default `FALSE`): Whether to return the fit models and save them to an RData file.
		
## Troubleshooting ##

1. Question: When I run from the command line I see the error ``maaslin3.R: command not found``. How do I fix this? 
    * Answer: Provide the full path to the executable when running maaslin3.R
2. Question: When I run as a function I see the error ``Error in library(maaslin3): there is no package called 'maaslin3'``. How do I fix this? 
    * Answer: Install the R package and then try loading the library again.
