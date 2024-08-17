#!/usr/bin/env Rscript

###############################################################################
# MaAsLin 3

# Copyright (c) 2024 Harvard School of Public Health

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

###############################################################
# If running on the command line, load other Maaslin3 modules #
###############################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
    # source all R in Maaslin3 package, relative to this folder
    # same method as original maaslin
    script_options <- commandArgs(trailingOnly = FALSE)
    script_path <-
        sub("--file=", "", script_options[grep("--file=", script_options)])
    script_dir <- dirname(script_path)
    script_name <- basename(script_path)
    R_files <- c("fit.R", "utility_scripts.R", "viz.R")
    `%>%` <- dplyr::`%>%`
    
    # load in the required libraries, report an error if they are not installed
    # for (lib in c(
    #     'optparse',
    #     'logging',
    #     'data.table',
    #     'dplyr',
    #     'pbapply',
    #     'lmerTest',
    #     'parallel',
    #     'lme4',
    #     'multcomp',
    #     'ggplot2',
    #     'viridis',
    #     "grid",
    #     'RColorBrewer',
    #     'patchwork',
    #     'scales'
    # )) {
    #     suppressPackageStartupMessages(require(lib, character.only = TRUE))
    # }
    for (R_file in R_files)
    {
        if (!(R_file == script_name))
            source(file.path(script_dir, R_file))
    }
}

#### Set the default options ####

normalization_choices <- c("TSS", "CLR", "NONE")
transform_choices <- c("LOG", "PLOG", "NONE")

correction_choices <-
    c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY")

# set the default run options
args <- list()
args$input_data <- NULL
args$input_metadata <- NULL
args$output <- NULL
args$min_abundance <- 0.0
args$zero_threshold <- 0.0
args$min_prevalence <- 0.0
args$min_variance <- 0.0
args$max_significance <- 0.1
args$normalization <- normalization_choices[1]
args$transform <- transform_choices[1]
args$correction <- correction_choices[1]
args$formula <- NULL
args$random_effects <- NULL
args$group_effects <- NULL
args$ordered_effects <- NULL
args$strata_effects <- NULL
args$fixed_effects <- NULL
args$feature_specific_covariate <- NULL
args$feature_specific_covariate_name <- NULL
args$feature_specific_covariate_record <- NULL
args$standardize <- TRUE
args$median_comparison_abundance <- TRUE
args$median_comparison_prevalence <- FALSE
args$median_comparison_abundance_threshold <- 0.25
args$median_comparison_prevalence_threshold <- 0.25
args$subtract_median <- FALSE
args$augment <- TRUE
args$evaluate_only <- NULL
args$unscaled_abundance <- NULL
args$plot_summary_plot <- TRUE
args$summary_plot_first_n <- 25
args$coef_plot_vars <- NULL
args$heatmap_vars <- NULL
args$plot_associations <- TRUE
args$max_pngs <- 30
args$cores <- 1
args$save_models <- FALSE
args$reference <- NULL

#### end ####

#### Add command line arguments ####

options <-
    optparse::OptionParser(usage = paste(
        "%prog",
        "<data.tsv>",
        "<metadata.tsv>",
        "<output_folder>",
        "[options]"
    ))
options <-
    optparse::add_option(
        options,
        c("--formula"),
        type = "character",
        dest = "formula",
        default = args$formula,
        help = paste("The formula for the model",
                    "[ Default: all variables fixed ]")
    )
options <-
    optparse::add_option(
        options,
        c("--fixed_effects"),
        type = "character",
        dest = "fixed_effects",
        default = args$fixed_effects,
        help = paste(
            "The fixed effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: all ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--reference"),
        type = "character",
        dest = "reference",
        default = args$reference,
        help = paste(
            "The factor to use as a reference for",
            "a variable with more than two levels",
            "provided as a string of 'variable,reference'",
            "semi-colon delimited for multiple variables [ Default: NA ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--random_effects"),
        type = "character",
        dest = "random_effects",
        default = args$random_effects,
        help = paste(
            "The random effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--group_effects"),
        type = "character",
        dest = "group_effects",
        default = args$group_effects,
        help = paste(
            "The group effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--ordered_effects"),
        type = "character",
        dest = "ordered_effects",
        default = args$ordered_effects,
        help = paste(
            "The ordered effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--strata_effects"),
        type = "character",
        dest = "strata_effects",
        default = args$ordered_effects,
        help = paste(
            "The strata effects for the model.",
            "Only one is allowed.",
            "[ Default: none ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--feature_specific_covariate"),
        type = "character",
        dest = "feature_specific_covariate",
        default = args$feature_specific_covariate,
        help = paste(
            "The table to use for feature-specific",
            "covariates. Row and column names should",
            "match the data input."
        )
    )
options <-
    optparse::add_option(
        options,
        c("--feature_specific_covariate_name"),
        type = "character",
        dest = "feature_specific_covariate_name",
        default = args$feature_specific_covariate,
        help = paste("The name of the feature-specific covariate")
    )
options <-
    optparse::add_option(
        options,
        c("--feature_specific_covariate_record"),
        type = "character",
        dest = "feature_specific_covariate_record",
        default = args$feature_specific_covariate_record,
        help = paste("Whether to include the feature-specific 
                    covariate in the outputs")
    )
options <-
    optparse::add_option(
        options,
        c("--min_abundance"),
        type = "double",
        dest = "min_abundance",
        default = args$min_abundance,
        help = paste(
            "The minimum abundance for each feature 
            (before normalization and transformation)",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--min_prevalence"),
        type = "double",
        dest = "min_prevalence",
        default = args$min_prevalence,
        help = paste(
            "The minimum proportion of samples for which",
            "a feature is detected at minimum abundance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--zero-threshold"),
        type = "double",
        dest = "zero_threshold",
        default = args$zero_threshold,
        help = paste(
            "The minimum abundance to be considered non-zero",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--min_variance"),
        type = "double",
        dest = "min_variance",
        default = args$min_variance,
        help = paste(
            "Keep features with variances",
            "greater than value",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--max_significance"),
        type = "double",
        dest = "max_significance",
        default = args$max_significance,
        help = paste("The q-value threshold for significance",
                    "[ Default: %default ]")
    )
options <-
    optparse::add_option(
        options,
        c("--normalization"),
        type = "character",
        dest = "normalization",
        default = args$normalization,
        help = paste(
            "The normalization method to apply",
            "[ Default: %default ] [ Choices:",
            toString(normalization_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--transform"),
        type = "character",
        dest = "transform",
        default = args$transform,
        help = paste(
            "The transform to apply [ Default: %default ] [ Choices:",
            toString(transform_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--correction"),
        type = "character",
        dest = "correction",
        default = args$correction,
        help = paste(
            "The correction method for computing",
            "the q-value [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--standardize"),
        type = "logical",
        dest = "standardize",
        default = args$standardize,
        help = paste(
            "Apply z-score so continuous metadata are on",
            "the same scale [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--unscaled_abundance"),
        type = "character",
        dest = "unscaled_abundance",
        default = args$unscaled_abundance,
        help = paste(
            "The table to use as an unscaled",
            "abundance reference (the single column",
            "name must be the same as one of the",
            "features or 'total')"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--median_comparison_abundance"),
        type = "logical",
        dest = "median_comparison_abundance",
        default = args$median_comparison_abundance,
        help = paste(
            "Test abundance coefficients against the median",
            "association [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--median_comparison_prevalence"),
        type = "logical",
        dest = "median_comparison_prevalence",
        default = args$median_comparison_prevalence,
        help = paste(
            "Test prevalence coefficients against the median",
            "association [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--median_comparison_abundance_threshold"),
        type = "double",
        dest = "median_comparison_abundance_threshold",
        default = args$median_comparison_abundance_threshold,
        help = paste(
            "Radius within which the median adjustment",
            "gives a p-value of 1 [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--median_comparison_prevalence_threshold"),
        type = "double",
        dest = "median_comparison_prevalence_threshold",
        default = args$median_comparison_prevalence_threshold,
        help = paste(
            "Radius within which the median adjustment",
            "gives a p-value of 1 [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--subtract_median"),
        type = "logical",
        dest = "subtract_median",
        default = args$subtract_median,
        help = paste(
            "Subtract the median from coefficients when",
            "doing median comparisons [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--augment"),
        type = "logical",
        dest = "augment",
        default = args$augment,
        help = paste(
            "Add weighted extra 0s and 1s to avoid linear",
            "separability [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--evaluate_only"),
        type = "character",
        dest = "evaluate_only",
        default = args$evaluate_only,
        help = paste(
            "Whether to evaluate just the abundnace or 
            prevalence models [ Default: %default ] [ Choices:",
            toString(c("abundance", "prevalence"),),
            "]"
        )
        
    )

options <-
    optparse::add_option(
        options,
        c("--plot_summary_plot"),
        type = "logical",
        dest = "plot_summary_plot",
        default = args$plot_summary_plot,
        help = paste(
            "Generate a summary plot for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--summary_plot_first_n"),
        type = "double",
        dest = "summary_plot_first_n",
        default = args$summary_plot_first_n,
        help = paste(
            "In summary plot, plot top N features with significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--coef_plot_vars"),
        type = "character",
        dest = "coef_plot_vars",
        default = args$coef_plot_vars,
        help = paste(
            "The variables to use in the coefficient plot",
            "section of the summary plot provided as a",
            "comma-separated string. Continuous variables",
            "should match the metadata column name, and",
            "categorical variables should be of the form:",
            "[metadata] [level]. [ Default: NA ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--heatmap_vars"),
        type = "character",
        dest = "heatmap_vars",
        default = args$coef_plot_vars,
        help = paste(
            "The variables to use in the heatmap",
            "section of the summary plot provided as a",
            "comma-separated string. Continuous variables",
            "should match the metadata column name, and",
            "categorical variables should be of the form:",
            "[metadata] [level]. [ Default: NA ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--plot_associations"),
        type = "logical",
        dest = "plot_associations",
        default = args$plot_associations,
        help = paste(
            "Generate associations plots for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--max_pngs"),
        type = "double",
        dest = "max_pngs",
        default = args$max_pngs,
        help = paste(
            "The maximum number of association plots for signficant",
            "associations to save as png files [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--cores"),
        type = "double",
        dest = "cores",
        default = args$cores,
        help = paste(
            "The number of R processes to",
            "run in parallel [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("--save_models"),
        type = "logical",
        dest = "save_models",
        default = args$save_models,
        help = paste(
            "Return the full model outputs",
            "and save to an RData file [ Default: %default ]"
        )
    )

option_not_valid_error <- function(message, valid_options) {
    logging::logerror(paste(message, ": %s"), toString(valid_options))
    stop("Option not valid", call. = FALSE)
}

#### end ####

####################################
# Check valid options are selected #
####################################

maaslin_check_arguments <-
    function(feature_specific_covariate = NULL,
            feature_specific_covariate_name = NULL,
            feature_specific_covariate_record = NULL,
            zero_threshold = 0,
            normalization = 'TSS',
            transform = 'LOG',
            correction = 'BH',
            evaluate_only = NULL,
            unscaled_abundance = NULL,
            median_comparison_abundance = TRUE) {
        normalization <- toupper(normalization)
        transform <- toupper(transform)
        
        # Match variable ignoring case then set 
        # correctly as required for p.adjust
        correction <- correction_choices[match(toupper(correction),
                                            toupper(correction_choices))]
        
        if (!is.null(unscaled_abundance) &&
            median_comparison_abundance) {
            stop(
                "`median_comparison_abundance` usually should not be 
                TRUE (default) with unscaled abundances. To bypass this 
                check, run the maaslin steps individually, 
                skipping `maaslin_check_arguments`."
            )
        }
        
        if (!((
            is.null(feature_specific_covariate) +
            is.null(feature_specific_covariate_name) +
            is.null(feature_specific_covariate_record)
        ) %in% c(0, 3))) {
            stop(
                "`feature_specific_covariate`, 
                `feature_specific_covariate_name`, 
                and `feature_specific_covariate_record`
                should all be null or all be non-null"
            )
        }
        
        # Check valid normalization option selected
        logging::loginfo("Verifying options selected are valid")
        if (!normalization %in% normalization_choices) {
            option_not_valid_error(
                paste(
                    "Please select a normalization",
                    "from the list of available options"
                ),
                toString(normalization_choices)
            )
        }
        
        if (!is.null(unscaled_abundance) & normalization != 'TSS') {
            stop("Normalization must be TSS if using unscaled abundance")
        }
        
        # check valid transform option selected
        if (!transform %in% transform_choices) {
            option_not_valid_error(
                "Please select a transform from the list of available options",
                toString(transform_choices)
            )
        }
        
        # check valid correction method selected
        if (!correction %in% correction_choices) {
            option_not_valid_error(
                paste(
                    "Please select a correction method",
                    "from the list of available options"
                ),
                toString(correction_choices)
            )
        }
        
        if (normalization == 'CLR' & transform != 'NONE') {
            stop("normalization = CLR can only be used with transform NONE")
        }
        
        if (transform == 'PLOG' & zero_threshold >= 0) {
            stop(
                "transform set to PLOG, but zero_threshold is >= 0. 
                Set zero_threshold to -1 to count all 
                features as present and apply PLOG."
            )
        }
        
        if (!is.null(evaluate_only) &&
            !evaluate_only %in% c("abundance", "prevalence")) {
            stop("evaluate_only must be NULL, 'abundance', or 'prevalence'")
        }
        
        if (transform == 'PLOG' &&
            (is.null(evaluate_only) ||
            evaluate_only != "abundance")) {
            stop("PLOG should only be used for 
                abundance-only modeling (evaluate_only)")
        }
    }

#####################################
# Create log file and log arguments #
#####################################

maaslin_log_arguments <- function(input_data,
                                input_metadata,
                                output,
                                formula = NULL,
                                fixed_effects = NULL,
                                reference = NULL,
                                random_effects = NULL,
                                group_effects = NULL,
                                ordered_effects = NULL,
                                strata_effects = NULL,
                                feature_specific_covariate = NULL,
                                feature_specific_covariate_name = NULL,
                                feature_specific_covariate_record = NULL,
                                min_abundance = 0,
                                min_prevalence = 0.0,
                                zero_threshold = 0,
                                min_variance = 0,
                                max_significance = 0.1,
                                normalization = 'TSS',
                                transform = 'LOG',
                                correction = 'BH',
                                standardize = TRUE,
                                unscaled_abundance = NULL,
                                median_comparison_abundance = TRUE,
                                median_comparison_prevalence = FALSE,
                                median_comparison_abundance_threshold = 0.25,
                                median_comparison_prevalence_threshold = 0.25,
                                subtract_median = FALSE,
                                augment = TRUE,
                                evaluate_only = NULL,
                                plot_summary_plot = TRUE,
                                summary_plot_first_n = 25,
                                coef_plot_vars = NULL,
                                heatmap_vars = NULL,
                                plot_associations = TRUE,
                                max_pngs = 30,
                                cores = 1,
                                save_models = FALSE) {
    # Allow for lower case variables
    normalization <- toupper(normalization)
    transform <- toupper(transform)
    
    # Match variable ignoring case then set correctly as required for p.adjust
    correction <- correction_choices[match(toupper(correction),
                                        toupper(correction_choices))]
    
    # If formula is a formula object, convert it back to a string
    if (methods::is(formula, "formula")) {
        formula <- paste0(trimws(deparse(formula)), collapse = " ")
    }
    
    # create an output folder
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    
    # create log file (write info to stdout and debug level to log file)
    # set level to finest so all log levels are reviewed
    log_file <- file.path(output, "maaslin3.log")
    
    # remove log file if already exists (to avoid append)
    if (file.exists(log_file)) {
        logging::logwarn(paste(
            "Warning: Deleting existing log file:", log_file))
        unlink(log_file)
    }
    
    logging::logReset()
    logging::basicConfig(level = 'FINEST')
    logging::addHandler(logging::writeToFile,
                        file = log_file, level = "DEBUG")
    logging::setLevel(20, logging::getHandler('basic.stdout'))
    
    logging::loginfo("Writing function arguments to log file")
    logging::logdebug("Function arguments")
    if (is.character(input_data)) {
        logging::logdebug("Input data file: %s", input_data)
    }
    if (is.character(input_metadata)) {
        logging::logdebug("Input metadata file: %s", input_metadata)
    }
    logging::logdebug("Output folder: %s", output)
    logging::logdebug("Min Abundance: %f", min_abundance)
    logging::logdebug("Zero Threshold: %f", zero_threshold)
    logging::logdebug("Min Prevalence: %f", min_prevalence)
    logging::logdebug("Normalization: %s", normalization)
    logging::logdebug("Transform: %s", transform)
    logging::logdebug("Max significance: %f", max_significance)
    logging::logdebug("Random effects: %s", random_effects)
    logging::logdebug("Fixed effects: %s", fixed_effects)
    logging::logdebug("Group effects: %s", group_effects)
    logging::logdebug("Ordered effects: %s", ordered_effects)
    logging::logdebug("Strata effects: %s", strata_effects)
    logging::logdebug("Formula: %s", formula)
    logging::logdebug("Correction method: %s", correction)
    logging::logdebug("Standardize: %s", standardize)
    logging::logdebug("Augment: %s", augment)
    logging::logdebug("Evaluate only: %s", evaluate_only)
    logging::logdebug("Cores: %d", cores)
    logging::logdebug("Abundance median comparison: %s",
                    median_comparison_abundance)
    logging::logdebug("Prevalence median comparison: %s",
                    median_comparison_prevalence)
    logging::logdebug(
        "Abundance median comparison threshold: %s",
        median_comparison_abundance_threshold
    )
    logging::logdebug(
        "Prevalence median comparison threshold: %s",
        median_comparison_prevalence_threshold
    )
    logging::logdebug(
        "Subtract median: %s",
        subtract_median
    )
    if (is.character(unscaled_abundance)) {
        logging::logdebug("Unscaled abundance: %s", unscaled_abundance)
    }
    if (!is.null(feature_specific_covariate)) {
        if (is.character(feature_specific_covariate)) {
            logging::logdebug("Feature specific covariate: %s",
                            feature_specific_covariate)
        }
    }
    if (!is.null(feature_specific_covariate_name)) {
        if (is.character(feature_specific_covariate_name)) {
            logging::logdebug("Feature specific covariate name: %s",
                            feature_specific_covariate_name)
        }
    }
    if (!is.null(feature_specific_covariate_record)) {
        if (is.character(feature_specific_covariate_record)) {
            logging::logdebug(
                "Feature specific covariate include: %s",
                feature_specific_covariate_record
            )
        }
    }
    
    maaslin_check_arguments(
        feature_specific_covariate,
        feature_specific_covariate_name,
        feature_specific_covariate_record,
        zero_threshold,
        normalization,
        transform,
        correction,
        evaluate_only,
        unscaled_abundance,
        median_comparison_abundance
    )
}

#################################
# Read in the data and metadata #
#################################

maaslin_read_data <- function(input_data,
                            input_metadata,
                            feature_specific_covariate = NULL,
                            unscaled_abundance = NULL) {
    # if a character string then this is a file name, else it
    # is a data frame
    if (is.character(input_data) && file.exists(input_data)) {
        data <- read.table(input_data,
                        header = TRUE,
                        row.names = 1)
    } else if (is.data.frame(input_data)) {
        if (!tibble::has_rownames(input_data)) {
            stop("If supplying input_data as a data frame, 
                it must have appropriate rownames!")
        }
        data <-
            as.data.frame(input_data) # in case it's a tibble or something
    } else if (is.matrix(input_data)) {
        logging::logwarn("Input is a matrix, 
                        passing through as.data.frame() .")
        data <- as.data.frame(input_data)
    } else {
        stop("input_data is neither a file nor a data frame!")
    }
    
    if (is.character(input_metadata) &&
        file.exists(input_metadata)) {
        metadata <- read.table(input_metadata,
                            header = TRUE,
                            row.names = 1)
    } else if (is.data.frame(input_metadata)) {
        if (!tibble::has_rownames(input_metadata)) {
            stop(
                "If supplying input_metadata as a data frame, 
                it must have appropriate rownames!"
            )
        }
        metadata <-
            as.data.frame(input_metadata) # in case it's a tibble or something
    } else {
        stop("input_metadata is neither a file nor a data frame!")
    }
    
    if (is.character(unscaled_abundance) &&
        file.exists(unscaled_abundance)) {
        unscaled_abundance <-
            read.table(unscaled_abundance,
                    header = TRUE,
                    row.names = 1)
    } else if (is.data.frame(unscaled_abundance)) {
        if (!tibble::has_rownames(unscaled_abundance)) {
            stop(
                "If supplying unscaled_abundance as a data frame, 
                it must have appropriate rownames!"
            )
        }
        unscaled_abundance <-
            as.data.frame(unscaled_abundance)
    } else if (!is.null(unscaled_abundance)) {
        stop("unscaled_abundance is not a file or data frame!")
    }
    
    if (is.character(feature_specific_covariate) &&
        file.exists(feature_specific_covariate)) {
        feature_specific_covariate <-
            read.table(feature_specific_covariate,
                    header = TRUE,
                    row.names = 1)
    } else if (is.data.frame(feature_specific_covariate)) {
        if (!tibble::has_rownames(feature_specific_covariate)) {
            stop(
                "If supplying feature_specific_covariate as a data frame, 
                it must have appropriate rownames!"
            )
        }
        feature_specific_covariate <-
            as.data.frame(feature_specific_covariate) 
    } else if (!is.null(feature_specific_covariate)) {
        stop("feature_specific_covariate is not a file or data frame!")
    }
    
    return(
        list(
            "data" = data,
            "metadata" = metadata,
            "feature_specific_covariate" = feature_specific_covariate,
            "unscaled_abundance" = unscaled_abundance
        )
    )
}

###############################################################
# Determine orientation of data in input and reorder to match #
###############################################################

maaslin_reorder_data <- function(data,
                                metadata,
                                feature_specific_covariate = NULL,
                                unscaled_abundance = NULL) {
    logging::loginfo("Determining format of input files")
    samples_row_row <- intersect(rownames(data), rownames(metadata))
    if (length(samples_row_row) > 0) {
        # this is the expected formatting so do not modify data frames
        logging::loginfo(paste(
            "Input format is data samples",
            "as rows and metadata samples as rows"
        ))
    } else {
        samples_column_row <- intersect(colnames(data), rownames(metadata))
        
        if (length(samples_column_row) == 0) {
            # modify possibly included special chars in sample names in metadata
            rownames(metadata) <- make.names(rownames(metadata))
            
            samples_column_row <-
                intersect(colnames(data), rownames(metadata))
        }
        
        if (length(samples_column_row) > 0) {
            logging::loginfo(
                paste(
                    "Input format is data samples",
                    "as columns and metadata samples as rows"
                )
            )
            # transpose data frame so samples are rows
            data <- as.data.frame(t(data))
            logging::logdebug("Transformed data so samples are rows")
        } else {
            samples_column_column <-
                intersect(colnames(data), colnames(metadata))
            if (length(samples_column_column) > 0) {
                logging::loginfo(
                    paste(
                        "Input format is data samples",
                        "as columns and metadata samples as columns"
                    )
                )
                data <- as.data.frame(BiocGenerics::t(data))
                metadata <- type.convert(
                    as.data.frame(BiocGenerics::t(metadata)))
                logging::logdebug("Transformed data and metadata 
                                so samples are rows")
            } else {
                samples_row_column <-
                    intersect(rownames(data), colnames(metadata))
                
                if (length(samples_row_column) == 0) {
                    # modify possibly included special chars 
                    # in sample names in data
                    rownames(data) <- make.names(rownames(data))
                    
                    samples_row_column <-
                        intersect(rownames(data), colnames(metadata))
                }
                
                if (length(samples_row_column) > 0) {
                    logging::loginfo(
                        paste(
                            "Input format is data samples",
                            "as rows and metadata samples as columns"
                        )
                    )
                    metadata <-
                        type.convert(as.data.frame(BiocGenerics::t(metadata)))
                    logging::logdebug("Transformed metadata so 
                                    samples are rows")
                } else {
                    logging::logerror(
                        paste(
                            "Unable to find samples in data and",
                            "metadata files.",
                            "Rows/columns do not match."
                        )
                    )
                    logging::logdebug("Data rows: %s",
                                    paste(rownames(data), collapse = ","))
                    logging::logdebug("Data columns: %s",
                                    paste(colnames(data), collapse = ","))
                    logging::logdebug("Metadata rows: %s",
                                    paste(rownames(metadata), collapse = ","))
                    logging::logdebug("Metadata columns: %s",
                                    paste(colnames(data), collapse = ","))
                    stop()
                }
            }
        }
    }
    
    if (!is.null(feature_specific_covariate)) {
        samples_row_row <-
            intersect(rownames(data),
                    rownames(feature_specific_covariate))
        samples_col_col <-
            intersect(colnames(data),
                    colnames(feature_specific_covariate))
        if (length(samples_row_row) > 0 &
            length(samples_col_col) > 0) {
            # this is the expected formatting so do not modify data frames
            logging::loginfo(
                paste(
                    "Input format is data samples",
                    "as rows and feature_specific_covariate samples as rows"
                )
            )
        } else {
            samples_column_row <-
                intersect(colnames(data),
                        rownames(feature_specific_covariate))
            samples_row_column <-
                intersect(rownames(data),
                        colnames(feature_specific_covariate))
            
            if (length(samples_column_row) == 0 |
                length(samples_row_column) == 0) {
                # modify possibly included special 
                # chars in sample names in metadata
                rownames(feature_specific_covariate) <-
                    make.names(rownames(feature_specific_covariate))
                rownames(data) <- make.names(rownames(data))
                
                samples_column_row <-
                    intersect(colnames(data),
                            rownames(feature_specific_covariate))
                samples_row_column <-
                    intersect(rownames(data),
                            colnames(feature_specific_covariate))
            }
            
            if (length(samples_column_row) > 0 &
                length(samples_row_column) > 0) {
                logging::loginfo(
                    paste(
                        "Input format is feature_specific_covariate samples",
                        "as columns"
                    )
                )
                # transpose data frame so samples are rows
                feature_specific_covariate <-
                    as.data.frame(BiocGenerics::t(feature_specific_covariate))
                logging::logdebug("Transformed feature_specific_covariate 
                                so samples are rows")
            } else {
                logging::logerror(
                    paste(
                        "Unable to find samples in feature_specific_covariate.",
                        "Rows/columns do not match."
                    )
                )
                logging::logdebug("Data rows: %s",
                                paste(rownames(data), collapse = ","))
                logging::logdebug("Data columns: %s",
                                paste(colnames(data), collapse = ","))
                logging::logdebug(
                    "Feature specific covariate rows: %s",
                    paste(
                        rownames(feature_specific_covariate),
                        collapse = ","
                    )
                )
                logging::logdebug(
                    "Feature specific covariate columns: %s",
                    paste(
                        colnames(feature_specific_covariate),
                        collapse = ","
                    )
                )
                stop()
            }
        }
        
    }
    
    # replace unexpected characters in feature names
    colnames(data) <- make.names(colnames(data))
    if (!is.null(unscaled_abundance)) {
        colnames(unscaled_abundance) <-
            make.names(colnames(unscaled_abundance))
    }
    if (!is.null(feature_specific_covariate)) {
        colnames(feature_specific_covariate) <-
            make.names(colnames(feature_specific_covariate))
    }
    
    # get a set of the samples with both metadata and features
    intersect_samples <-
        intersect(rownames(data), rownames(metadata))
    logging::logdebug(
        "A total of %s samples were found in both the data and metadata",
        length(intersect_samples)
    )
    
    if (!is.null(feature_specific_covariate)) {
        intersect_samples <-
            intersect(intersect_samples,
                    rownames(feature_specific_covariate))
        logging::logdebug(
            "A total of %s samples were found in the data, metadata, 
            and feature specific covariates",
            length(intersect_samples)
        )
    }
    
    # check for samples without metadata
    extra_feature_samples <-
        setdiff(rownames(data), intersect_samples)
    if (length(extra_feature_samples) > 0)
        logging::loginfo(
            paste(
                "The following samples were found",
                "to have features but no metadata",
                "(or feature specific covariates if",
                "applicable).",
                "They will be removed. %s"
            ),
            paste(extra_feature_samples, collapse = ",")
        )
    
    # check for metadata samples without features
    extra_metadata_samples <-
        setdiff(rownames(metadata), intersect_samples)
    if (length(extra_metadata_samples) > 0)
        logging::loginfo(
            paste(
                "The following samples were found",
                "to have metadata but no features",
                "(or feature specific covariates if",
                "applicable).",
                "They will be removed. %s"
            ),
            paste(extra_metadata_samples, collapse = ",")
        )
    
    if (!is.null(feature_specific_covariate)) {
        extra_feature_specific_covariate_samples <-
            setdiff(rownames(feature_specific_covariate),
                    intersect_samples)
        if (length(extra_feature_specific_covariate_samples) > 0)
            logging::loginfo(
                paste(
                    "The following samples were found",
                    "to have feature specific covariates",
                    "but no features or no metadata.",
                    "They will be removed. %s"
                ),
                paste(extra_feature_specific_covariate_samples, collapse = ",")
            )
    }
    
    if (!is.null(unscaled_abundance)) {
        extra_unscaled_abundance_samples <-
            setdiff(rownames(unscaled_abundance), rownames(data))
        if (length(extra_unscaled_abundance_samples) > 0)
            logging::logdebug(
                paste(
                    "The following samples were found",
                    "to have unscaled abundances but no features.",
                    "They will be removed. %s"
                ),
                paste(extra_unscaled_abundance_samples, collapse = ",")
            )
    }
    
    if (!is.null(unscaled_abundance))  {
        if (!all(rownames(data) %in% rownames(unscaled_abundance))) {
            stop("some data samples do not have an unscaled abundance")
        } else if (length(colnames(unscaled_abundance)) > 1) {
            stop("there is more than 1 column in 
                the unscaled abundance data frame")
        } else if (colnames(unscaled_abundance) %in% colnames(data)) {
            logging::logdebug("Using unscaled abundance as spike-in feature")
        } else if (colnames(unscaled_abundance) == 'total') {
            logging::logdebug("Using unscaled abundance as total abundances")
        } else {
            stop("unscaled abundance column must be a feature name or 'total'")
        }
    }
    
    # now order both data and metadata with the same sample ordering
    logging::logdebug("Reordering data/metadata to use same sample ordering")
    data <- data[intersect_samples, , drop = FALSE]
    metadata <- metadata[intersect_samples, , drop = FALSE]
    
    if (!is.null(unscaled_abundance)) {
        unscaled_abundance <-
            unscaled_abundance[intersect_samples, , drop = FALSE]
    }
    
    if (!is.null(feature_specific_covariate)) {
        feature_specific_covariate <-
            feature_specific_covariate[intersect_samples, , drop = FALSE]
    }
    
    return(
        list(
            "data" = data,
            "metadata" = metadata,
            "feature_specific_covariate" = feature_specific_covariate,
            "unscaled_abundance" = unscaled_abundance
        )
    )
}

###########################################
# Compute the formula based on user input #
###########################################

maaslin_compute_formula <- function(data,
                                    metadata,
                                    fixed_effects = NULL,
                                    random_effects = NULL,
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL) {
    random_effects_formula <- NULL
    # use all metadata if no fixed effects are provided
    if (is.null(fixed_effects)) {
        fixed_effects <- colnames(metadata)
    } else {
        fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
        # remove any fixed effects not found in metadata names
        to_remove <- setdiff(fixed_effects, colnames(metadata))
        if (length(to_remove) > 0) {
            logging::logerror(
                paste(
                    "Variable name not found in metadata",
                    "so not applied to formula as fixed effect: %s"
                ),
                paste(to_remove, collapse = " , ")
            )
            stop()
        }
    }
    
    if (!is.null(random_effects)) {
        random_effects <-
            unlist(strsplit(random_effects, ",", fixed = TRUE))
        
        common_variables <- intersect(fixed_effects, random_effects)
        if (length(common_variables) > 0) {
            logging::logwarn(
                paste(
                    "Feature name included as fixed and random effect,",
                    "check that this is intended: %s"
                ),
                paste(common_variables, collapse = " , ")
            )
        }
        
        # remove any random effects not found in metadata
        to_remove <- setdiff(random_effects, colnames(metadata))
        if (length(to_remove) > 0) {
            logging::logerror(
                paste(
                    "Effect name not found in metadata",
                    "so not applied to formula as random effect: %s"
                ),
                paste(to_remove, collapse = " , ")
            )
            stop()
        }
        
        # create formula
        if (length(random_effects) > 0) {
            random_effects_formula_text <-
                paste("expr ~ (1 | ",
                    paste(
                        random_effects,
                        ")",
                        sep = '',
                        collapse = " + (1 | "
                    ),
                    sep = '')
            logging::loginfo("Formula for random effects: %s",
                            random_effects_formula_text)
            random_effects_formula <-
                tryCatch(
                    as.formula(random_effects_formula_text),
                    error = function(e)
                        stop(
                            sprintf(
                                "Invalid formula for random effects: %s",
                                random_effects_formula_text
                            )
                        )
                )
        }
    }
    
    if (!is.null(group_effects) | !is.null(ordered_effects)) {
        multi_effects <- c()
        if (!is.null(group_effects)) {
            multi_effects <-
                c(multi_effects,
                strsplit(group_effects, ",", fixed = TRUE))
        }
        if (!is.null(ordered_effects)) {
            multi_effects <-
                c(multi_effects,
                strsplit(ordered_effects, ",", fixed = TRUE))
        }
        
        common_variables <- intersect(fixed_effects, multi_effects)
        if (length(common_variables) > 0) {
            logging::logerror(
                paste(
                    "Feature name included as fixed and group/ordered effect,",
                    "this is not allowed: %s"
                ),
                paste(common_variables, collapse = " , ")
            )
            stop()
        }
        
        # remove any random effects not found in metadata
        to_remove <- setdiff(multi_effects, colnames(metadata))
        if (length(to_remove) > 0) {
            logging::logerror(paste0(
                "Effect name not found in metadata: ",
                paste0(to_remove, collapse = ", ")
            ))
            stop()
        }
    }
    
    if (length(fixed_effects) == 0 &
        length(group_effects) == 0 &
        length(ordered_effects) == 0 &
        is.null(feature_specific_covariate_name)) {
        logging::logerror("No fixed/group/ordered/
                        feature-specific effects provided.")
        stop()
    }
    
    # reduce metadata to only include fixed/group/random effects in formula
    effects_names <-
        unique(
            c(
                fixed_effects,
                random_effects,
                group_effects,
                ordered_effects,
                strata_effects
            )
        )
    metadata <- metadata[, effects_names, drop = FALSE]
    
    # create the fixed effects formula text
    formula_effects <- fixed_effects
    if (length(group_effects) > 0) {
        formula_effects <-
            union(formula_effects, paste0("group(", group_effects, ")"))
    }
    if (length(ordered_effects) > 0) {
        formula_effects <-
            union(formula_effects,
                paste0("ordered(", ordered_effects, ")"))
    }
    if (length(strata_effects) > 0) {
        formula_effects <-
            union(formula_effects,
                paste0("strata(", strata_effects, ")"))
    }
    if (!is.null(feature_specific_covariate_name)) {
        formula_effects <-
            union(formula_effects, feature_specific_covariate_name)
    }
    
    formula_text <-
        paste("expr ~ ", paste(formula_effects, collapse = " + "))
    logging::loginfo("Formula for fixed effects: %s", formula_text)
    formula <-
        tryCatch(
            as.formula(formula_text),
            error = function(e)
                stop(
                    sprintf(
                        "Invalid formula. Please provide 
                        a different formula: %s",
                        formula_text
                    )
                )
        )
    
    if (!(is.null(random_effects_formula))) {
        formula <-
            paste('. ~',
                paste(formula_effects, collapse = ' + '),
                '.',
                sep = ' + ')
        formula <- update(random_effects_formula, formula)
    }
    
    return(list(
        "formula" = formula,
        "random_effects_formula" = random_effects_formula
    ))
}

##############################
# Check a user input formula #
##############################

maaslin_check_formula <- function(data,
                                metadata,
                                input_formula = NULL,
                                feature_specific_covariate_name = NULL) {
    if (methods::is(input_formula, "formula")) {
        input_formula <-
            paste0(trimws(deparse(input_formula)), collapse = " ")
    }
    
    random_effects_formula <- NULL
    
    if (is.null(input_formula)) {
        logging::logerror(paste("No user formula provided"))
    }
    
    # Remove anything before the tilde if necessary
    input_formula <- sub(".*~\\s*", "", input_formula)
    
    if (!is.null(feature_specific_covariate_name)) {
        if (!grepl(feature_specific_covariate_name, input_formula)) {
            input_formula <-
                paste0("expr ~ ",
                    feature_specific_covariate_name,
                    ' + ',
                    input_formula)
        } else {
            input_formula <- paste0("expr ~ ", input_formula)
        }
    } else {
        input_formula <- paste0("expr ~ ", input_formula)
    }
    
    formula <-
        tryCatch(
            as.formula(input_formula),
            error = function(e)
                stop(sprintf("Invalid formula: %s",
                            input_formula))
        )
    
    formula_terms <- all.vars(formula)
    formula_terms <-
        formula_terms[!formula_terms %in% c("expr", 
                                            feature_specific_covariate_name)]
    
    to_remove <- setdiff(formula_terms, colnames(metadata))
    if (length(to_remove) > 0) {
        logging::logerror(
            paste("Effect name not found in metadata: %s"),
            paste(to_remove, collapse = ", ")
        )
        stop()
    }
    
    term_labels <- attr(terms(formula), "term.labels")
    
    if (sum(!grepl("strata\\(|\\|", term_labels)) == 0 &
        is.null(feature_specific_covariate_name)) {
        logging::logerror("No fixed, group, or 
                        ordered effects included in formula.")
        stop()
    }
    
    # create formula
    if (sum(grepl("\\|", term_labels)) > 0) {
        logging::loginfo("Formula for random effects: %s",
                        input_formula)
        random_effects_formula <- formula
    } else {
        # create the fixed effects formula text
        formula_text <- deparse(formula)
        logging::loginfo("Formula for fixed effects: %s", formula_text)
        random_effects_formula <- NULL
    }
    
    # reduce metadata to only include fixed/group/ordered/
    # strata/random effects in formula
    metadata <- metadata[, formula_terms, drop = FALSE]
    
    return(list(
        "formula" = formula,
        "random_effects_formula" = random_effects_formula
    ))
}

#################
# Normalization #
#################

maaslin_normalize <- function(data,
                            output,
                            zero_threshold = 0,
                            normalization = 'TSS',
                            unscaled_abundance = NULL) {
    features <- data
    
    normalization <- toupper(normalization)
    
    logging::loginfo("Running selected normalization method: %s", normalization)
    
    if (normalization == 'TSS') {
        features <- TSSnorm(features, zero_threshold)
    }
    if (normalization == 'CLR') {
        features <- CLRnorm(features, zero_threshold)
    }
    if (normalization == 'NONE') {
        features <- NONEnorm(features, zero_threshold)
    }
    
    if (!is.null(unscaled_abundance)) {
        features <-
            UNSCALEDnorm(features, unscaled_abundance, zero_threshold)
    }
    
    features_folder <- file.path(output, "features")
    if (!file.exists(features_folder)) {
        logging::loginfo("Creating output feature tables folder")
        dir.create(features_folder, recursive = TRUE)
    }
    
    data_norm_file <- file.path(features_folder, "data_norm.tsv")
    logging::loginfo("Writing normalized data to file %s", data_norm_file)
    write.table(
        data.frame("feature" = rownames(features), features),
        file = data_norm_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    
    return(features)
}


########################################################################
# Filter data based on min abundance, min prevalence, and min variance #
########################################################################

maaslin_filter <- function(normalized_data,
                        output,
                        min_abundance = 0,
                        min_prevalence = 0.0,
                        zero_threshold = 0,
                        min_variance = 0) {
    unfiltered_data <- normalized_data
    
    # require at least total samples * min prevalence values
    # for each feature to be greater than min abundance
    logging::loginfo("Filter data based on min abundance and min prevalence")
    total_samples <- nrow(unfiltered_data)
    logging::loginfo("Total samples in data: %d", total_samples)
    min_samples <- total_samples * min_prevalence
    logging::loginfo(
        paste(
            "Min samples required with min abundance",
            "for a feature not to be filtered: %f"
        ),
        min_samples
    )
    
    # Filter by abundance
    data_zeros <- unfiltered_data
    for (col_index in seq_along(data_zeros)) {
        data_zeros[, col_index][is.na(data_zeros[, col_index])] <-
            min(min_abundance, zero_threshold) - 1
    }
    
    ##########################################
    # Apply the non-zero abundance threshold #
    ##########################################
    
    prevalence_mask <- ifelse(data_zeros > zero_threshold, 1, 0)
    data_zeros <- data_zeros * prevalence_mask
    
    filtered_data <-
        unfiltered_data[, colSums(data_zeros > min_abundance) > min_samples,
                        drop = FALSE]
    total_filtered_features <-
        ncol(unfiltered_data) - ncol(filtered_data)
    logging::loginfo("Total filtered features: %d", total_filtered_features)
    filtered_feature_names <-
        setdiff(names(unfiltered_data), names(filtered_data))
    logging::loginfo(
        "Filtered feature names from abundance and prevalence filtering: %s",
        toString(filtered_feature_names)
    )
    
    #################################
    # Filter data based on variance #
    #################################
    
    vars <- apply(filtered_data, 2, var, na.rm = TRUE)
    variance_filtered_data <-
        filtered_data[, which(vars > min_variance), drop = FALSE]
    variance_filtered_features <-
        ncol(filtered_data) - ncol(variance_filtered_data)
    logging::loginfo(
        "Total features filtered by non-zero variance filtering: %d",
        variance_filtered_features)
    variance_filtered_feature_names <-
        setdiff(names(filtered_data), names(variance_filtered_data))
    logging::loginfo(
        "Filtered feature names from variance filtering: %s",
        toString(variance_filtered_feature_names)
    )
    filtered_data <- variance_filtered_data
    
    #######################
    # Write filtered data #
    #######################
    
    features_folder <- file.path(output, "features")
    if (!file.exists(features_folder)) {
        logging::loginfo("Creating output feature tables folder")
        dir.create(features_folder, recursive = TRUE)
    }
    
    filtered_file <- file.path(features_folder, "filtered_data.tsv")
    logging::loginfo("Writing filtered data to file %s", filtered_file)
    write.table(
        data.frame("feature" = rownames(filtered_data), filtered_data),
        file = filtered_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    
    return(filtered_data)
}

#######################################
# Filter and standardize the metadata #
#######################################

maaslin_process_metadata <- function(metadata,
                                    formula = NULL,
                                    fixed_effects = NULL,
                                    reference = NULL,
                                    feature_specific_covariate_name = NULL,
                                    standardize = TRUE) {
    if (is.null(reference)) {
        reference <- ","
    }
    split_reference <- unlist(strsplit(reference, "[,;]"))
    
    # Extract fixed effects from formula if null
    if (is.null(fixed_effects)) {
        term_labels <- attr(terms(formula), "term.labels")
        term_labels <-
            term_labels[!grepl("group|ordered|strata|\\|", term_labels)]
        tmp_formula <-
            formula(paste0("~ ", paste0(term_labels, collapse = " + ")))
        formula_terms <- all.vars(tmp_formula)
        if (is.null(feature_specific_covariate_name)) {
            fixed_effects <- formula_terms
        } else {
            fixed_effects <-
                formula_terms[formula_terms != feature_specific_covariate_name]
        }
    }
    
    # for each fixed effect, check that a 
    # reference level has been set if necessary:
    # number of levels > 2 and metadata isn't already an ordered factor
    for (i in fixed_effects) {
        # don't check for or require reference levels for numeric metadata
        if (is.numeric(metadata[, i])) {
            next
        }
        # respect ordering if a factor is 
        # explicitly passed in with no reference set
        if (is.factor(metadata[, i]) && !(i %in% split_reference)) {
            logging::loginfo(
                paste(
                    "Factor detected for categorial metadata '",
                    i,
                    "'. Using as-is.",
                    sep = ""
                )
            )
            next
        }
        
        # set metadata as a factor (ordered alphabetically)
        metadata[, i] <- as.factor(metadata[, i])
        mlevels <- levels(metadata[, i])
        
        # get reference level for variable being considered, 
        # returns NA if not found
        ref <- split_reference[match(i, split_reference) + 1]

        # if metadata has 2 levels, allow but don't require 
        # setting reference level, otherwise require it
        if ((length(mlevels) == 2)) {
            if (!is.na(ref)) {
                metadata[, i] <- relevel(metadata[, i], ref = ref)
            }
        } else if (length(mlevels) > 2) {
            if (!is.na(ref)) {
                metadata[, i] <- relevel(metadata[, i], ref = ref)
            } else {
                stop(
                    paste(
                        "Please provide the reference for the variable '",
                        i,
                        "' which includes more than 2 levels: ",
                        paste(as.character(mlevels), collapse = ", "),
                        ". ",
                        "Alternatively, set the variable as 
                        a factor beforehand.",
                        sep = ""
                    )
                )
            }
        } else {
            stop("Provided categorical metadata has 
                fewer than 2 unique, non-NA values.")
        }
    }
    
    ################################
    # Standardize metadata, if set #
    ################################
    
    if (standardize) {
        logging::loginfo("Applying z-score to standardize continuous metadata")
        metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
    } else {
        logging::loginfo("Bypass z-score application to metadata")
    }
    
    return(metadata)
}

##################
# Transformation #
##################

maaslin_transform <- function(filtered_data,
                            output,
                            transform = 'LOG') {
    features <- filtered_data
    
    logging::loginfo("Running selected transform method: %s", transform)
    
    if (transform == 'LOG') {
        features <- LOG(features)
    }
    if (transform == 'PLOG') {
        features <- PLOG(features)
    }
    
    features_folder <- file.path(output, "features")
    if (!file.exists(features_folder)) {
        logging::loginfo("Creating output feature tables folder")
        dir.create(features_folder, recursive = TRUE)
    }
    
    filtered_data_norm_transformed_file <-
        file.path(features_folder, "data_transformed.tsv")
    logging::loginfo(
        "Writing normalized, filtered, transformed data to file %s",
        filtered_data_norm_transformed_file
    )
    write.table(
        data.frame("feature" = rownames(features), features),
        file = filtered_data_norm_transformed_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    
    return(features)
}

##############
# Fit linear #
##############

maaslin_fit <- function(filtered_data,
                        transformed_data,
                        metadata,
                        formula,
                        random_effects_formula,
                        feature_specific_covariate = NULL,
                        feature_specific_covariate_name = NULL,
                        feature_specific_covariate_record = NULL,
                        zero_threshold = 0,
                        correction = 'BH',
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        median_comparison_abundance_threshold = 0.25,
                        median_comparison_prevalence_threshold = 0.25,
                        subtract_median = FALSE,
                        augment = TRUE,
                        evaluate_only = NULL,
                        cores = 1,
                        save_models = FALSE) {
    logging::loginfo("Running the linear model component")
    
    if (!is.null(feature_specific_covariate)) {
        tryCatch({
            feature_specific_covariate <-
                feature_specific_covariate[rownames(filtered_data), 
                                        colnames(filtered_data)]
        }, error = function(e) {
            stop(
                "feature_specific_covariate does not contain the 
                features and samples of the filtered data."
            )
        })
    }
    
    # Match variable ignoring case then set correctly as required for p.adjust
    correction <- correction_choices[match(toupper(correction),
                                        toupper(correction_choices))]
    
    if (is.null(evaluate_only) || evaluate_only == "abundance") {
        #######################
        # For non-zero models #
        #######################
        
        fit_data_abundance <-
            fit.model(
                features = transformed_data,
                metadata = metadata,
                model = 'linear',
                formula = formula,
                random_effects_formula = random_effects_formula,
                correction = correction,
                save_models = save_models,
                augment = augment,
                cores = cores,
                median_comparison = median_comparison_abundance,
                median_comparison_threshold = 
                    median_comparison_abundance_threshold,
                subtract_median = subtract_median,
                feature_specific_covariate = feature_specific_covariate,
                feature_specific_covariate_name = 
                    feature_specific_covariate_name,
                feature_specific_covariate_record = 
                    feature_specific_covariate_record
            )
        
        #################################################################
        # Count the total values for each feature (untransformed space) #
        #################################################################
        
        logging::loginfo("Counting total values for each feature")
        
        fit_data_abundance$results$N <-
            apply(
                fit_data_abundance$results,
                1,
                FUN = function(x)
                    length(transformed_data[, x[1]])
            )
        fit_data_abundance$results$N.not.zero <-
            apply(
                fit_data_abundance$results,
                1,
                FUN = function(x)
                    length(which(filtered_data[, x[1]] > zero_threshold))
            )
        
    }
    
    if (is.null(evaluate_only) || evaluate_only == "prevalence") {
        #####################
        # For binary models #
        #####################
        
        logging::loginfo("Running the logistic model component")
        
        prevalence_mask <- ifelse(!is.na(filtered_data), 1, 0)
        
        fit_data_prevalence <-
            fit.model(
                features = prevalence_mask,
                metadata = metadata,
                model = 'logistic',
                formula = formula,
                random_effects_formula = random_effects_formula,
                correction = correction,
                save_models = save_models,
                augment = augment,
                cores = cores,
                median_comparison = median_comparison_prevalence,
                median_comparison_threshold = 
                    median_comparison_prevalence_threshold,
                subtract_median = subtract_median,
                feature_specific_covariate = feature_specific_covariate,
                feature_specific_covariate_name = 
                    feature_specific_covariate_name,
                feature_specific_covariate_record = 
                    feature_specific_covariate_record
            )
        
        logging::loginfo("Counting total values for each feature")
        
        fit_data_prevalence$results$N <-
            apply(
                fit_data_prevalence$results,
                1,
                FUN = function(x)
                    length(transformed_data[, x[1]])
            )
        fit_data_prevalence$results$N.not.zero <-
            apply(
                fit_data_prevalence$results,
                1,
                FUN = function(x)
                    length(which(filtered_data[, x[1]] > zero_threshold))
            )
    }
    
    if (is.null(evaluate_only)) {
        results <-
            add_joint_signif(fit_data_abundance,
                            fit_data_prevalence,
                            'linear',
                            correction)
        fit_data_abundance$results <- results[[1]]
        fit_data_prevalence$results <- results[[2]]
    } else if (evaluate_only == 'abundance') {
        fit_data_abundance$results$pval_joint <-
            fit_data_abundance$results$pval
        fit_data_abundance$results$qval_joint <-
            fit_data_abundance$results$qval
        fit_data_abundance$results <- fit_data_abundance$results %>%
            dplyr::rename(pval_individual = .data$pval,
                        qval_individual = .data$qval)
    } else if (evaluate_only == 'prevalence') {
        fit_data_prevalence$results$pval_joint <-
            fit_data_prevalence$results$pval
        fit_data_prevalence$results$qval_joint <-
            fit_data_prevalence$results$qval
        fit_data_prevalence$results <-
            fit_data_prevalence$results %>%
            dplyr::rename(pval_individual = .data$pval,
                        qval_individual = .data$qval)
    }
    
    if (is.null(evaluate_only) || evaluate_only == "prevalence") {
        current_likely_error_subsetter <-
            !is.na(fit_data_prevalence$results$N.not.zero) &
            (
                fit_data_prevalence$results$N.not.zero < 50 &
                    fit_data_prevalence$results$N.not.zero / 
                    fit_data_prevalence$results$N < 0.05
            ) &
            ((
                !is.na(fit_data_prevalence$results$coef) &
                    abs(fit_data_prevalence$results$coef) > 15
            ) |
                (
                    !is.na(fit_data_prevalence$results$pval_individual) &
                        fit_data_prevalence$results$pval_individual < 10 ^ -10
                )
            )
        current_errors_for_likely_issues <-
            fit_data_prevalence$results$error[current_likely_error_subsetter]
        fit_data_prevalence$results$error[current_likely_error_subsetter] <-
            ifelse(
                !is.na(current_errors_for_likely_issues),
                current_errors_for_likely_issues,
                "A large coefficient (>15 in absolute value) or small 
                p-value (< 10^-10) was obtained from a feature present 
                in <5% of samples. Check this is intended."
            )
        
        current_likely_error_subsetter <-
            !is.na(fit_data_prevalence$results$N.not.zero) &
            (
                fit_data_prevalence$results$N - 
                    fit_data_prevalence$results$N.not.zero < 50 &
                    fit_data_prevalence$results$N.not.zero / 
                    fit_data_prevalence$results$N > 0.95
            ) &
            ((
                !is.na(fit_data_prevalence$results$coef) &
                    abs(fit_data_prevalence$results$coef) > 15
            ) |
                (
                    !is.na(fit_data_prevalence$results$pval_individual) &
                        fit_data_prevalence$results$pval_individual < 10 ^ -10
                )
            )
        current_errors_for_likely_issues <-
            fit_data_prevalence$results$error[current_likely_error_subsetter]
        fit_data_prevalence$results$error[current_likely_error_subsetter] <-
            ifelse(
                !is.na(current_errors_for_likely_issues),
                current_errors_for_likely_issues,
                "A large coefficient (>15 in absolute value) or small p-value 
                (< 10^-10) was obtained from a feature present in >95% of 
                samples. Check this is intended."
            )
        
        fit_data_prevalence$results <-
            fit_data_prevalence$results[
                order(fit_data_prevalence$results$qval_joint), ]
        fit_data_prevalence$results <-
            fit_data_prevalence$results[
                order(!is.na(fit_data_prevalence$results$error)), ] 
        # Move all that had errors to the end
    } else {
        fit_data_prevalence <- NULL
    }
    
    if (is.null(evaluate_only) || evaluate_only == "abundance") {
        fit_data_abundance$results <-
            fit_data_abundance$results[
                order(fit_data_abundance$results$qval_joint), ]
        fit_data_abundance$results <-
            fit_data_abundance$results[
                order(!is.na(fit_data_abundance$results$error)), ] 
        # Move all that had errors to the end
    } else {
        fit_data_abundance <- NULL
    }
    
    return(
        list(
            "fit_data_abundance" = fit_data_abundance,
            "fit_data_prevalence" = fit_data_prevalence
        )
    )
}

###########################
# Write tables of results #
###########################

maaslin_write_results <- function(output,
                                fit_data_abundance,
                                fit_data_prevalence,
                                random_effects_formula = NULL,
                                max_significance = 0.1,
                                save_models = FALSE) {
    # create an output folder if it does not exist
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    
    write_fits(
        output,
        fit_data_abundance,
        fit_data_prevalence,
        random_effects_formula,
        save_models
    )
    
    write_results(output,
                fit_data_abundance,
                fit_data_prevalence,
                max_significance)
}

maaslin_write_results_lefse_format <- function(output,
                                            fit_data_abundance,
                                            fit_data_prevalence) {
    # create an output folder if it does not exist
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    
    if (!is.null(fit_data_abundance$results)) {
        write_results_in_lefse_format(
            fit_data_abundance$results,
            file.path(output, 'lefse_style_results_abundance.res')
        )
    }
    if (!is.null(fit_data_prevalence$results)) {
        write_results_in_lefse_format(
            fit_data_prevalence$results,
            file.path(output, 'lefse_style_results_prevalence.res')
        )
    }
}

#######################################################
# Create visualizations for results passing threshold #
#######################################################

maaslin_plot_results <- function(output,
                                transformed_data,
                                unstandardized_metadata,
                                fit_data_abundance,
                                fit_data_prevalence,
                                normalization,
                                transform,
                                feature_specific_covariate = NULL,
                                feature_specific_covariate_name = NULL,
                                feature_specific_covariate_record = NULL,
                                median_comparison_abundance = TRUE,
                                median_comparison_prevalence = FALSE,
                                max_significance = 0.1,
                                plot_summary_plot = TRUE,
                                summary_plot_first_n = 25,
                                coef_plot_vars = NULL,
                                heatmap_vars = NULL,
                                plot_associations = TRUE,
                                max_pngs = 30) {
    # create an output folder and figures folder if it does not exist
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    if (plot_summary_plot || plot_associations) {
        figures_folder <- file.path(output, "figures")
        if (!file.exists(figures_folder)) {
            logging::loginfo("Creating output figures folder")
            dir.create(figures_folder)
        }
    }
    
    if (is.null(fit_data_abundance$results)) {
        merged_results <- fit_data_prevalence$results
    } else if (is.null(fit_data_prevalence$results)) {
        merged_results <- fit_data_abundance$results
    } else {
        merged_results <- rbind(fit_data_abundance$results,
                                fit_data_prevalence$results)
    }
    
    if (plot_summary_plot) {
        summary_plot_file <- file.path(figures_folder, "summary_plot.pdf")
        logging::loginfo("Writing summary plot of significant 
                        results to file: %s",
                        summary_plot_file)
        
        if (!is.null(coef_plot_vars) &
            length(coef_plot_vars) == 1) {
            coef_plot_vars <- trimws(unlist(strsplit(coef_plot_vars, ',')))
        }
        if (!is.null(heatmap_vars) & length(heatmap_vars) == 1) {
            heatmap_vars <- trimws(unlist(strsplit(heatmap_vars, ',')))
        }
        
        maaslin3_summary_plot(
            merged_results,
            summary_plot_file,
            figures_folder,
            first_n = summary_plot_first_n,
            max_significance = max_significance,
            coef_plot_vars = coef_plot_vars,
            heatmap_vars = heatmap_vars,
            median_comparison_abundance = median_comparison_abundance,
            median_comparison_prevalence = median_comparison_prevalence
        )
    }
    
    if (plot_associations) {
        logging::loginfo(
            paste(
                "Writing association plots",
                "(one for each significant association)",
                "to output folder: %s"
            ),
            figures_folder
        )
        
        return(
            maaslin3_association_plots(
                merged_results = merged_results,
                metadata = unstandardized_metadata,
                features = transformed_data,
                max_significance = max_significance,
                figures_folder = figures_folder,
                max_pngs = max_pngs,
                normalization = normalization,
                transform = transform,
                feature_specific_covariate = feature_specific_covariate,
                feature_specific_covariate_name = 
                    feature_specific_covariate_name,
                feature_specific_covariate_record = 
                    feature_specific_covariate_record
            )
        )
    }
}

maaslin_plot_results_from_output <- function(output,
                                            metadata,
                                            normalization,
                                            transform,
                                            feature_specific_covariate = 
                                                NULL,
                                            feature_specific_covariate_name = 
                                                NULL,
                                            feature_specific_covariate_record =
                                                NULL,
                                            median_comparison_abundance = 
                                                TRUE,
                                            median_comparison_prevalence = 
                                                FALSE,
                                            max_significance = 0.1,
                                            plot_summary_plot = TRUE,
                                            summary_plot_first_n = 25,
                                            coef_plot_vars = NULL,
                                            heatmap_vars = NULL,
                                            plot_associations = TRUE,
                                            max_pngs = 30) {
    
    # create an output folder and figures folder if it does not exist
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    if (plot_summary_plot || plot_associations) {
        figures_folder <- file.path(output, "figures")
        if (!file.exists(figures_folder)) {
            logging::loginfo("Creating output figures folder")
            dir.create(figures_folder)
        }
    }
    
    all_results_file <-
        paste0(gsub('/$', '', output), '/', 'all_results.tsv')
    if (!file.exists(all_results_file)) {
        stop(sprintf(
            'Please generate the results file first: %s',
            all_results_file
        ))
    }
    merged_results <- utils::read.csv(all_results_file, sep = '\t')
    merged_results$model[merged_results$model == 'abundance'] <-
        'linear'
    merged_results$model[merged_results$model == 'prevalence'] <-
        'logistic'
    
    if (plot_summary_plot) {
        summary_plot_file <- file.path(figures_folder, "summary_plot.pdf")
        logging::loginfo("Writing summary plot of 
                        significant results to file: %s",
                        summary_plot_file)
        
        if (!is.null(coef_plot_vars) &
            length(coef_plot_vars) == 1) {
            coef_plot_vars <- trimws(unlist(strsplit(coef_plot_vars, ',')))
        }
        if (!is.null(heatmap_vars) & length(heatmap_vars) == 1) {
            heatmap_vars <- trimws(unlist(strsplit(heatmap_vars, ',')))
        }
        
        maaslin3_summary_plot(
            merged_results,
            summary_plot_file,
            figures_folder,
            first_n = summary_plot_first_n,
            max_significance = max_significance,
            coef_plot_vars = coef_plot_vars,
            heatmap_vars = heatmap_vars,
            median_comparison_abundance = median_comparison_abundance,
            median_comparison_prevalence = median_comparison_prevalence
        )
    }
    
    if (plot_associations) {
        features_file <-
            paste0(gsub('/$', '', output),
                '/',
                'features/data_transformed.tsv')
        if (!file.exists(features_file)) {
            stop(sprintf(
                'Please generate the results file first: %s',
                features_file
            ))
        }
        transformed_data <-
            utils::read.csv(
                features_file,
                sep = '\t',
                row.names = 1,
                check.names = FALSE
            )
        
        logging::loginfo(
            paste(
                "Writing association plots",
                "(one for each significant association)",
                "to output folder: %s"
            ),
            figures_folder
        )
        
        # Need to redo this if not fitting the model
        if (!is.null(feature_specific_covariate)) {
            tryCatch({
                feature_specific_covariate <-
                    feature_specific_covariate[rownames(transformed_data), 
                                            colnames(transformed_data)]
            }, error = function(e) {
                stop(
                    "feature_specific_covariate does not contain the features 
                    and samples of the filtered data."
                )
            })
        }
        
        if (missing("normalization") |
            missing("transform") | missing("metadata")) {
            stop(
                "Missing normalization, transform, or metadata argument to 
                maaslin_plot_results_from_output"
            )
        }
        
        plots_out <- maaslin3_association_plots(
            merged_results = merged_results,
            metadata = metadata,
            features = transformed_data,
            max_significance = max_significance,
            figures_folder = figures_folder,
            max_pngs = max_pngs,
            normalization = normalization,
            transform = transform,
            feature_specific_covariate = feature_specific_covariate,
            feature_specific_covariate_name = feature_specific_covariate_name,
            feature_specific_covariate_record = 
                feature_specific_covariate_record
        )
    } else {
        plots_out <- NULL
    }
    
    if ('logging::writeToFile' %in% names(logging::getLogger()[['handlers']])) {
        logging::removeHandler('logging::writeToFile')
    }
    
    return(plots_out)
}

#######################################################
# Main maaslin3 function (defaults same command line) #
#######################################################

maaslin3 <- function(input_data,
                    input_metadata,
                    output,
                    formula = NULL,
                    fixed_effects = NULL,
                    reference = NULL,
                    random_effects = NULL,
                    group_effects = NULL,
                    ordered_effects = NULL,
                    strata_effects = NULL,
                    feature_specific_covariate = NULL,
                    feature_specific_covariate_name = NULL,
                    feature_specific_covariate_record = NULL,
                    min_abundance = 0,
                    min_prevalence = 0.0,
                    zero_threshold = 0,
                    min_variance = 0,
                    max_significance = 0.1,
                    normalization = 'TSS',
                    transform = 'LOG',
                    correction = 'BH',
                    standardize = TRUE,
                    unscaled_abundance = NULL,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    median_comparison_abundance_threshold = 0.25,
                    median_comparison_prevalence_threshold = 0.25,
                    subtract_median = FALSE,
                    augment = TRUE,
                    evaluate_only = NULL,
                    plot_summary_plot = TRUE,
                    summary_plot_first_n = 25,
                    coef_plot_vars = NULL,
                    heatmap_vars = NULL,
                    plot_associations = TRUE,
                    max_pngs = 30,
                    cores = 1,
                    save_models = FALSE) {
    # Allow for lower case variables
    normalization <- toupper(normalization)
    transform <- toupper(transform)
    
    # Match variable ignoring case then set correctly as required for p.adjust
    correction <- correction_choices[match(toupper(correction),
                                        toupper(correction_choices))]
    
    # If formula is a formula object, convert it back to a string
    if (methods::is(formula, "formula")) {
        formula <- paste0(trimws(deparse(formula)), collapse = " ")
    }
    
    # Create log file, log arguments, and check arguments
    maaslin_log_arguments(
        input_data,
        input_metadata,
        output,
        formula,
        fixed_effects,
        reference,
        random_effects,
        group_effects,
        ordered_effects,
        strata_effects,
        feature_specific_covariate,
        feature_specific_covariate_name,
        feature_specific_covariate_record,
        min_abundance,
        min_prevalence,
        zero_threshold,
        min_variance,
        max_significance,
        normalization,
        transform,
        correction,
        standardize,
        unscaled_abundance,
        median_comparison_abundance,
        median_comparison_prevalence,
        median_comparison_abundance_threshold,
        median_comparison_prevalence_threshold,
        subtract_median,
        augment,
        evaluate_only,
        plot_summary_plot,
        summary_plot_first_n,
        coef_plot_vars,
        heatmap_vars,
        plot_associations,
        max_pngs,
        cores,
        save_models
    )
    
    read_data_list <- maaslin_read_data(input_data,
                                        input_metadata,
                                        feature_specific_covariate,
                                        unscaled_abundance)
    
    read_data_list <- maaslin_reorder_data(
        read_data_list$data,
        read_data_list$metadata,
        read_data_list$feature_specific_covariate,
        read_data_list$unscaled_abundance
    )
    
    data <- read_data_list$data
    metadata <- read_data_list$metadata
    unscaled_abundance <- read_data_list$unscaled_abundance
    feature_specific_covariate <-
        read_data_list$feature_specific_covariate
    
    if (is.null(formula)) {
        formulas <- maaslin_compute_formula(
            data,
            metadata,
            fixed_effects,
            random_effects,
            group_effects,
            ordered_effects,
            strata_effects,
            feature_specific_covariate_name
        )
    } else {
        formulas <- maaslin_check_formula(data,
                                        metadata,
                                        formula,
                                        feature_specific_covariate_name)
    }
    
    formula <- formulas$formula
    random_effects_formula <- formulas$random_effects_formula
    
    normalized_data <- maaslin_normalize(data,
                                        output,
                                        zero_threshold,
                                        normalization,
                                        unscaled_abundance)
    
    filtered_data <- maaslin_filter(
        normalized_data,
        output,
        min_abundance,
        min_prevalence,
        zero_threshold,
        min_variance
    )
    
    transformed_data <- maaslin_transform(filtered_data,
                                        output,
                                        transform)
    
    standardized_metadata <- maaslin_process_metadata(
        metadata,
        formula,
        fixed_effects,
        reference,
        feature_specific_covariate_name,
        standardize
    )
    
    maaslin_results <- maaslin_fit(
        filtered_data,
        transformed_data,
        standardized_metadata,
        formula,
        random_effects_formula,
        feature_specific_covariate,
        feature_specific_covariate_name,
        feature_specific_covariate_record,
        zero_threshold,
        correction,
        median_comparison_abundance,
        median_comparison_prevalence,
        median_comparison_abundance_threshold,
        median_comparison_prevalence_threshold,
        subtract_median,
        augment,
        evaluate_only,
        cores,
        save_models
    )
    
    maaslin_write_results(
        output,
        maaslin_results$fit_data_abundance,
        maaslin_results$fit_data_prevalence,
        random_effects_formula,
        max_significance,
        save_models
    )
    
    if (plot_summary_plot | plot_associations) {
        maaslin_plot_results(
            output,
            transformed_data,
            metadata,
            maaslin_results$fit_data_abundance,
            maaslin_results$fit_data_prevalence,
            normalization,
            transform,
            feature_specific_covariate,
            feature_specific_covariate_name,
            feature_specific_covariate_record,
            median_comparison_abundance,
            median_comparison_prevalence,
            max_significance,
            plot_summary_plot,
            summary_plot_first_n,
            coef_plot_vars,
            heatmap_vars,
            plot_associations,
            max_pngs
        )
    }
    
    if ('logging::writeToFile' %in% names(logging::getLogger()[['handlers']])) {
        logging::removeHandler('logging::writeToFile')
    }
    
    return(
        list(
            "data" = data,
            "normalized_data" = normalized_data,
            "filtered_data" = filtered_data,
            "transformed_data" = transformed_data,
            "metadata" = metadata,
            "standardized_metadata" = standardized_metadata,
            "formula" = formulas,
            "fit_data_abundance" = maaslin_results$fit_data_abundance,
            "fit_data_prevalence" = maaslin_results$fit_data_prevalence
        )
    )
}

###########################################################################
# If running on the command line, get arguments and call maaslin function #
###########################################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {

    # get command line options and positional arguments
    parsed_arguments <- optparse::parse_args(options,
                                            positional_arguments = TRUE)
    current_args <- parsed_arguments[["options"]]
    positional_args <- parsed_arguments[["args"]]
    # check three positional arguments are provided
    if (length(positional_args) != 3) {
        optparse::print_help(options)
        stop(
            sprintf(
                "Please provide the required positional arguments 
                <data.tsv> <metadata.tsv> <output_folder>"
            )
        )
    }
    
    # call maaslin with the command line options
    fit_data <-
        maaslin3(
            input_data = positional_args[1],
            input_metadata = positional_args[2],
            output = positional_args[3],
            min_abundance = current_args$min_abundance,
            zero_threshold = current_args$zero_threshold,
            min_prevalence = current_args$min_prevalence,
            min_variance = current_args$min_variance,
            max_significance = current_args$max_significance,
            normalization = current_args$normalization,
            transform = current_args$transform,
            random_effects = current_args$random_effects,
            fixed_effects = current_args$fixed_effects,
            group_effects = current_args$group_effects,
            ordered_effects = current_args$ordered_effects,
            strata_effects = current_args$strata_effects,
            feature_specific_covariate = 
                current_args$feature_specific_covariate,
            feature_specific_covariate_name = 
                current_args$feature_specific_covariate_name,
            feature_specific_covariate_record = 
                current_args$feature_specific_covariate_record,
            median_comparison_abundance = 
                current_args$median_comparison_abundance,
            median_comparison_prevalence = 
                current_args$median_comparison_prevalence,
            median_comparison_abundance_threshold = 
                current_args$median_comparison_abundance_threshold,
            median_comparison_prevalence_threshold = 
                current_args$median_comparison_prevalence_threshold,
            subtract_median = current_args$subtract_median,
            formula = current_args$formula,
            correction = current_args$correction,
            standardize = current_args$standardize,
            cores = current_args$cores,
            plot_summary_plot = current_args$plot_summary_plot,
            summary_plot_first_n = current_args$summary_plot_first_n,
            coef_plot_vars = current_args$coef_plot_vars,
            heatmap_vars = current_args$heatmap_vars,
            plot_associations = current_args$plot_associations,
            max_pngs = current_args$max_pngs,
            save_models = current_args$save_models,
            augment = current_args$augment,
            evaluate_only = current_args$evaluate_only,
            reference = current_args$reference,
            unscaled_abundance = current_args$unscaled_abundance
        )
}
