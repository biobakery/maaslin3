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

# load in the required libraries, report an error if they are not installed

for (lib in c('optparse', 'logging', 'data.table', 'dplyr')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

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
    
    for (R_file in dir(script_dir, pattern = "*.R$"))
    {
        if (!(R_file == script_name))
            source(file.path(script_dir, R_file))
    }
}

#### Set the default options ####

normalization_choices <- c("TSS", "CLR", "CSS", "NONE", "TMM")
transform_choices <- c("LOG", "LOGIT", "AST", "NONE")
valid_choice_transform_norm <- hash::hash()
valid_choice_transform_norm[[transform_choices[2]]] <-
  normalization_choices[c(4)]
valid_choice_transform_norm[[transform_choices[3]]] <-
  normalization_choices[c(1, 4)]

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
args$random_effects <- NULL
args$group_effects <- NULL
args$ordered_effects <- NULL
args$fixed_effects <- NULL
args$formula <- NULL
args$standardize <- TRUE
args$median_comparison_abundance <- TRUE
args$median_comparison_prevalence <- FALSE
args$median_comparison_abundance_threshold <- 0.25
args$median_comparison_prevalence_threshold <- 0.25
args$augment <- TRUE
args$unscaled_abundance <- NULL
args$plot_heatmap <- TRUE
args$heatmap_first_n <- 50
args$plot_scatter <- TRUE
args$max_pngs <- 10
args$save_scatter <- FALSE
args$cores <- 1
args$save_models <- FALSE
args$reference <- NULL

#### end ####

#### Add command line arguments ####

options <-
    optparse::OptionParser(usage = paste(
        "%prog [options]",
        " <data.tsv> ",
        "<metadata.tsv> ",
        "<output_folder>" 
        )
    )
options <-
    optparse::add_option(
        options,
        c("-a", "--min_abundance"),
        type = "double",
        dest = "min_abundance",
        default = args$min_abundance,
        help = paste0("The minimum abundance for each feature (before normalization and transformation)",
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
    help = paste0("The minimum abundance to be considered non-zero",
                  "[ Default: %default ]"
    )
  )
options <-
    optparse::add_option(
        options,
        c("-p", "--min_prevalence"),
        type = "double",
        dest = "min_prevalence",
        default = args$min_prevalence,
        help = paste0("The minimum percent of samples for which",
            "a feature is detected at minimum abundance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-b", "--min_variance"),
        type = "double",
        dest = "min_variance",
        default = args$min_variance,
        help = paste0("Keep features with variances",
            "greater than value",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-s", "--max_significance"),
        type = "double",
        dest = "max_significance",
        default = args$max_significance,
        help = paste0("The q-value threshold for significance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-n", "--normalization"),
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
        c("-t", "--transform"),
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
        c("-r", "--random_effects"),
        type = "character",
        dest = "random_effects",
        default = args$random_effects,
        help = paste("The random effects for the model, ",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
  optparse::add_option(
    options,
    c("-f", "--fixed_effects"),
    type = "character",
    dest = "fixed_effects",
    default = args$fixed_effects,
    help = paste("The fixed effects for the model,",
                 "comma-delimited for multiple effects",
                 "[ Default: all ]"
    )
  )
options <-
  optparse::add_option(
    options,
    c("--group_effects"),
    type = "character",
    dest = "group_effects",
    default = args$group_effects,
    help = paste("The group effects for the model, ",
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
    help = paste("The ordered effects for the model, ",
                 "comma-delimited for multiple effects",
                 "[ Default: none ]"
    )
  )
options <-
    optparse::add_option(
        options,
        c("--formula"),
        type = "character",
        dest = "formula",
        default = args$formula,
        help = paste("The formula for the model,",
                     "[ Default: all variables fixed ]"
      )
  )
options <-
    optparse::add_option(
        options,
        c("-c", "--correction"),
        type = "character",
        dest = "correction",
        default = args$correction,
        help = paste("The correction method for computing",
            "the q-value [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-z", "--standardize"),
        type = "logical",
        dest = "standardize",
        default = args$standardize,
        help = paste("Apply z-score so continuous metadata are on",
            "the same scale [ Default: %default ]"
        )
    )
options <-
  optparse::add_option(
    options,
    c("--median_comparison_abundance"),
    type = "logical",
    dest = "median_comparison_abundance",
    default = args$median_comparison_abundance,
    help = paste("Test abundance coefficients against the median", 
                 "association [ Default: %default ]")
  )
options <-
  optparse::add_option(
    options,
    c("--median_comparison_prevalence"),
    type = "logical",
    dest = "median_comparison_prevalence",
    default = args$median_comparison_prevalence,
    help = paste("Test prevalence coefficients against the median", 
                 "association [ Default: %default ]")
  )
options <-
  optparse::add_option(
    options,
    c("--median_comparison_abundance_threshold"),
    type = "logical",
    dest = "median_comparison_abundance_threshold",
    default = args$median_comparison_abundance_threshold,
    help = paste("Radius within which the median adjustment", 
                 "gives a p-value of 1 [ Default: %default ]")
  )
options <-
  optparse::add_option(
    options,
    c("--median_comparison_prevalence_threshold"),
    type = "logical",
    dest = "median_comparison_prevalence_threshold",
    default = args$median_comparison_prevalence_threshold,
    help = paste("Radius within which the median adjustment", 
                 "gives a p-value of 1 [ Default: %default ]")
  )
options <-
  optparse::add_option(
    options,
    c("--augment"),
    type = "logical",
    dest = "augment",
    default = args$augment,
    help = paste("Add weighted extra 0s and 1s to avoid linear", 
                 "separability [ Default: %default ]")
  )
options <-
    optparse::add_option(
        options,
        c("-l", "--plot_heatmap"),
        type = "logical",
        dest = "plot_heatmap",
        default = args$plot_heatmap,
        help = paste("Generate a heatmap for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-i", "--heatmap_first_n"),
        type = "double",
        dest = "heatmap_first_n",
        default = args$heatmap_first_n,
        help = paste("In heatmap, plot top N features with significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-o", "--plot_scatter"),
        type = "logical",
        dest = "plot_scatter",
        default = args$plot_scatter,
        help = paste("Generate scatter plots for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-g", "--max_pngs"),
        type = "double",
        dest = "max_pngs",
        default = args$max_pngs,
        help = paste("The maximum number of scatterplots for signficant",
                     "associations to save as png files [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-O", "--save_scatter"),
        type = "logical",
        dest = "save_scatter",
        default = args$save_scatter,
        help = paste("Save all scatter plot ggplot objects",
                     "to an RData file [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-e", "--cores"),
        type = "double",
        dest = "cores",
        default = args$cores,
        help = paste("The number of R processes to ",
            "run in parallel [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-j", "--save_models"),
        type = "logical",
        dest = "save_models",
        default = args$save_models,
        help = paste("Return the full model outputs ",
                     "and save to an RData file [ Default: %default ]"
        )
    )
options <-
  optparse::add_option(
    options,
    c("--unscaled_abundance"),
    type = "character",
    dest = "unscaled_abundance",
    default = args$unscaled_abundance,
    help = paste("The table to use as an unscaled",
                 "abundance reference (the single column",
                 "name must be the same as one of the",
                 "features or 'total')"
    )
  )
options <-
    optparse::add_option(
        options,
        c("-d", "--reference"),
        type = "character",
        dest = "reference",
        default = args$reference,
        help = paste("The factor to use as a reference for",
            "a variable with more than two levels",
            "provided as a string of 'variable,reference'",
            "semi-colon delimited for multiple variables [ Default: NA ]"
        )
    )

option_not_valid_error <- function(message, valid_options) {
  logging::logerror(paste(message, ": %s"), toString(valid_options))
  stop("Option not valid", call. = FALSE)
}

#### end ####

##################################################
# Fill in parameter list with missing parameters #
##################################################

maaslin_parse_param_list <- function(param_list) {
  if (sum(c("input_data", "input_metadata", "output") %in% names(param_list)) != 3) {
    stop("param_list must include input_data, input_metadata, and output")
  }
  default_params = list(min_abundance = args$min_abundance,
                        zero_threshold = args$zero_threshold,
                        min_prevalence = args$min_prevalence,
                        min_variance = args$min_variance,
                        normalization = args$normalization,
                        transform = args$transform,
                        max_significance = args$max_significance,
                        random_effects = args$random_effects,
                        fixed_effects = args$fixed_effects,
                        group_effects = args$group_effects,
                        ordered_effects = args$ordered_effects,
                        formula = args$formula,
                        correction = args$correction,
                        standardize = args$standardize,
                        median_comparison_abundance = args$median_comparison_abundance,
                        median_comparison_prevalence = args$median_comparison_prevalence,
                        median_comparison_abundance_threshold = args$median_comparison_abundance_threshold,
                        median_comparison_prevalence_threshold = args$median_comparison_prevalence_threshold,
                        augment = args$augment,
                        cores = args$cores,
                        plot_heatmap = args$plot_heatmap,
                        heatmap_first_n = args$heatmap_first_n,
                        plot_scatter = args$plot_scatter,
                        max_pngs = args$max_pngs,
                        save_scatter = args$save_scatter,
                        save_models = args$save_models,
                        unscaled_abundance = args$unscaled_abundance,
                        reference = args$reference)
  
  input_param_list_names <- names(param_list)
  missing_params <- setdiff(names(default_params), input_param_list_names)
  extra_params <- setdiff(input_param_list_names, c("input_data", "input_metadata", "output", names(default_params)))
  
  if (length(extra_params) > 0) {
    stop(paste0('Extra parameters were included in the parameter list: ', 
                paste0(extra_params, collapse = ', ')))
  }
  
  for (param in missing_params) {
    param_list[[param]] <- default_params[[param]]
  }
  
  if (!is.null(param_list[['unscaled_abundance']]) & 
      param_list[["median_comparison_abundance"]] & 
      !('median_comparison_abundance' %in% input_param_list_names)) {
    stop("`median_comparison_abundance` usually should not be TRUE (default) with unscaled abundances. 
         Set `median_comparison_abundance` to TRUE in the parameters to bypass this error.")
  }
  
  # Allow for lower case variables
  param_list[["normalization"]] <- toupper(param_list[["normalization"]])
  param_list[["transform"]] <- toupper(param_list[["transform"]])

  # Match variable ignoring case then set correctly as required for p.adjust
  param_list[["correction"]] <- correction_choices[match(toupper(param_list[["correction"]]), 
                                                         toupper(correction_choices))]
  
  return(param_list)
}

####################################
# Check valid options are selected #
####################################

maaslin_check_arguments <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  
  # Check valid normalization option selected
  logging::loginfo("Verifying options selected are valid")
  if (!param_list[["normalization"]] %in% normalization_choices) {
    option_not_valid_error(
      paste(
        "Please select a normalization",
        "from the list of available options"),
      toString(normalization_choices)
    )
  }
  
  if (!is.null(param_list[['unscaled_abundance']]) & param_list[["normalization"]] != 'TSS') {
    stop(
      paste("Normalization must be TSS if using unscaled abundance")
    )
  }
  
  # check valid transform option selected
  if (!param_list[["transform"]] %in% transform_choices) {
    option_not_valid_error(
      "Please select a transform from the list of available options",
      toString(transform_choices)
    )
  }
  
  # check valid correction method selected
  if (!param_list[["correction"]] %in% correction_choices) {
    option_not_valid_error(
      paste("Please select a correction method",
            "from the list of available options"),
      toString(correction_choices)
    )
  }
  
  # check a valid choice combination is selected
  for (limited_transform in hash::keys(valid_choice_transform_norm)) {
    if (param_list[["transform"]] == limited_transform) {
      if (!param_list[["normalization"]] %in% 
          valid_choice_transform_norm[[limited_transform]]) {
        option_not_valid_error(
          paste0("This transform can only be used",
                 " with a subset of normalizations. ",
                 "Please select from the following list"
          ),
          toString(
            valid_choice_transform_norm[[limited_transform]])
        )
      }
    }
  }
  
  # check that plots are generated if to be saved
  if (!param_list[["plot_scatter"]] && param_list[["save_scatter"]]) {
    logging::logerror("Scatter plots cannot be saved if they are not plotted")
    stop("Option not valid", call. = FALSE)
  }
  
  return(param_list)
}

#####################################
# Create log file and log arguments #
#####################################

maaslin_log_arguments <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  output <- param_list[["output"]]
  
  # create an output folder
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  # create log file (write info to stdout and debug level to log file)
  # set level to finest so all log levels are reviewed
  log_file <- file.path(output, "maaslin3.log")
  
  # remove log file if already exists (to avoid append)
  if (file.exists(log_file)) {
    print(paste("Warning: Deleting existing log file:", log_file))
    unlink(log_file)
  }
  
  logging::basicConfig(level = 'FINEST')
  logging::addHandler(logging::writeToFile, 
                      file = log_file, level = "DEBUG")
  logging::setLevel(20, logging::getHandler('basic.stdout'))
  
  logging::loginfo("Writing function arguments to log file")
  logging::logdebug("Function arguments")
  if (is.character(param_list[["input_data"]])) {
    logging::logdebug("Input data file: %s", param_list[["input_data"]])
  }
  if (is.character(param_list[["input_metadata"]])) {
    logging::logdebug("Input metadata file: %s", param_list[["input_metadata"]])
  }
  logging::logdebug("Output folder: %s", param_list[["output"]])
  logging::logdebug("Min Abundance: %f", param_list[["min_abundance"]])
  logging::logdebug("Zero Threshold: %f", param_list[["zero_threshold"]])
  logging::logdebug("Min Prevalence: %f", param_list[["min_prevalence"]])
  logging::logdebug("Normalization: %s", param_list[["normalization"]])
  logging::logdebug("Transform: %s", param_list[["transform"]])
  logging::logdebug("Max significance: %f", param_list[["max_significance"]])
  logging::logdebug("Random effects: %s", param_list[["random_effects"]])
  logging::logdebug("Fixed effects: %s", param_list[["fixed_effects"]])
  logging::logdebug("Group effects: %s", param_list[["group_effects"]])
  logging::logdebug("Ordered effects: %s", param_list[["ordered_effects"]])
  logging::logdebug("Formula: %s", param_list[["formula"]])
  logging::logdebug("Correction method: %s", param_list[["correction"]])
  logging::logdebug("Standardize: %s", param_list[["standardize"]])
  logging::logdebug("Augment: %s", param_list[["augment"]])
  logging::logdebug("Cores: %d", param_list[["cores"]])
  logging::logdebug("Abundance median comparison: %s", param_list[["median_comparison_abundance"]])
  logging::logdebug("Prevalence median comparison: %s", param_list[["median_comparison_prevalence"]])
  logging::logdebug("Abundance median comparison threshold: %s", param_list[["median_comparison_abundance_threshold"]])
  logging::logdebug("Prevalence median comparison threshold: %s", param_list[["median_comparison_prevalence_threshold"]])
  if (is.character(param_list[["unscaled_abundance"]])) {
    logging::logdebug("Unscaled abundance: %s", param_list[["unscaled_abundance"]])
  }
  
  maaslin_check_arguments(param_list)
  
  return(param_list)
}

#################################
# Read in the data and metadata #
#################################

maaslin_read_data <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  input_data <- param_list[["input_data"]]
  input_metadata <- param_list[["input_metadata"]]
  unscaled_abundance <- param_list[["unscaled_abundance"]]
  
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(input_data) && file.exists(input_data)) {
    data <-
      data.frame(data.table::fread(
        input_data, header = TRUE, sep = "\t"),
        row.names = 1)
    if (nrow(data) == 1) {
      # read again to get row name
      data <- read.table(input_data, header = TRUE, row.names = 1)
    }
  } else if (is.data.frame(input_data)) {
    if (!tibble::has_rownames(input_data)) {
      stop("If supplying input_data as a data frame, it must have appropriate rownames!")
    }
    data <- as.data.frame(input_data) # in case it's a tibble or something
  } else if (is.matrix(input_data)) {
    logging::logwarn("Input is a matrix, passing through as.data.frame() .")
    data <- as.data.frame(input_data)
  } else {
    stop("input_data is neither a file nor a data frame!")
  }
  
  if (is.character(input_metadata) && file.exists(input_metadata)) {
    metadata <-
      data.frame(data.table::fread(
        input_metadata, header = TRUE, sep = "\t"),
        row.names = 1)
    if (nrow(metadata) == 1) {
      metadata <- read.table(input_metadata,
                             header = TRUE,
                             row.names = 1)
    }
  } else if (is.data.frame(input_metadata)) {
    if (!tibble::has_rownames(input_metadata)) {
      stop("If supplying input_metadata as a data frame, it must have appropriate rownames!")
    }
    metadata <- as.data.frame(input_metadata) # in case it's a tibble or something
  } else {
    stop("input_metadata is neither a file nor a data frame!")
  }
  
  if (is.character(unscaled_abundance) && file.exists(unscaled_abundance)) {
    unscaled_abundance <-
      data.frame(data.table::fread(
        unscaled_abundance, header = TRUE, sep = "\t"),
        row.names = 1)
  } else if (is.data.frame(unscaled_abundance)) {
    if (!tibble::has_rownames(unscaled_abundance)) {
      stop("If supplying unscaled_abundance as a data frame, it must have appropriate rownames!")
    }
    unscaled_abundance <- as.data.frame(unscaled_abundance) # in case it's a tibble or something
  } else if (!is.null(unscaled_abundance)) {
    stop("unscaled_abundance is not a file or data frame!")
  }
  
  param_list[["unscaled_abundance"]] <- unscaled_abundance
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata))
}

###############################################################
# Determine orientation of data in input and reorder to match #
###############################################################

maaslin_reorder_data <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  unscaled_abundance <- param_list[["unscaled_abundance"]]
  
  logging::loginfo("Determining format of input files")
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
    logging::loginfo(
      paste(
        "Input format is data samples",
        "as rows and metadata samples as rows"))
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata))
    
    if (length(samples_column_row) == 0) {
      # modify possibly included special chars in sample names in metadata
      rownames(metadata) <- make.names(rownames(metadata))
      
      samples_column_row <- intersect(colnames(data), rownames(metadata))
    }
    
    if (length(samples_column_row) > 0) {
      logging::loginfo(
        paste(
          "Input format is data samples",
          "as columns and metadata samples as rows"))
      # transpose data frame so samples are rows
      data <- as.data.frame(t(data))
      logging::logdebug(
        "Transformed data so samples are rows")
    } else {
      samples_column_column <- 
        intersect(colnames(data), colnames(metadata))
      if (length(samples_column_column) > 0) {
        logging::loginfo(
          paste(
            "Input format is data samples",
            "as columns and metadata samples as columns"))
        data <- as.data.frame(t(data))
        metadata <- type.convert(as.data.frame(t(metadata)))
        logging::logdebug(
          "Transformed data and metadata so samples are rows")
      } else {
        samples_row_column <- 
          intersect(rownames(data), colnames(metadata))
        
        if (length(samples_row_column) == 0) {
          # modify possibly included special chars in sample names in data
          rownames(data) <- make.names(rownames(data))
          
          samples_row_column <- intersect(rownames(data), colnames(metadata))
        }
        
        if (length(samples_row_column) > 0) {
          logging::loginfo(
            paste(
              "Input format is data samples",
              "as rows and metadata samples as columns"))
          metadata <- type.convert(as.data.frame(t(metadata)))
          logging::logdebug(
            "Transformed metadata so samples are rows")
        } else {
          logging::logerror(
            paste("Unable to find samples in data and",
                  "metadata files.",
                  "Rows/columns do not match."))
          logging::logdebug(
            "Data rows: %s", 
            paste(rownames(data), collapse = ","))
          logging::logdebug(
            "Data columns: %s", 
            paste(colnames(data), collapse = ","))
          logging::logdebug(
            "Metadata rows: %s", 
            paste(rownames(metadata), collapse = ","))
          logging::logdebug(
            "Metadata columns: %s",
            paste(colnames(data), collapse = ","))
          stop()
        }
      }
    }
  }
  
  # replace unexpected characters in feature names
  colnames(data) <- make.names(colnames(data))
  if (!is.null(unscaled_abundance)) {
    colnames(unscaled_abundance) <- make.names(colnames(unscaled_abundance))
  }

  # check for samples without metadata
  extra_feature_samples <-
    setdiff(rownames(data), rownames(metadata))
  if (length(extra_feature_samples) > 0)
    logging::logdebug(
      paste("The following samples were found",
            "to have features but no metadata.",
            "They will be removed. %s"),
      paste(extra_feature_samples, collapse = ",")
    )
  
  # check for metadata samples without features
  extra_metadata_samples <-
    setdiff(rownames(metadata), rownames(data))
  if (length(extra_metadata_samples) > 0)
    logging::logdebug(
      paste("The following samples were found",
            "to have metadata but no features.",
            "They will be removed. %s"),
      paste(extra_metadata_samples, collapse = ",")
    )
  
  if (!is.null(unscaled_abundance)) {
    extra_unscaled_abundance_samples <-
      setdiff(rownames(unscaled_abundance), rownames(data))
    if (length(extra_unscaled_abundance_samples) > 0)
      logging::logdebug(
        paste("The following samples were found",
              "to have unscaled abundances but no features.",
              "They will be removed. %s"),
        paste(extra_unscaled_abundance_samples, collapse = ",")
      )
  }
  
  # get a set of the samples with both metadata and features
  intersect_samples <- intersect(rownames(data), rownames(metadata))
  logging::logdebug(
    "A total of %s samples were found in both the data and metadata",
    length(intersect_samples)
  )
  
  if (!is.null(unscaled_abundance))  {
    if (!all(rownames(data) %in% rownames(unscaled_abundance))) {
      stop("some data samples do not have an unscaled abundance")
    } else if (length(colnames(unscaled_abundance)) > 1) {
      stop("there is more than 1 column in the unscaled abundance data frame")
    } else if (colnames(unscaled_abundance) %in% colnames(data)) {
      logging::logdebug(
        "Using unscaled abundance as spike-in feature"
      )
    } else if (colnames(unscaled_abundance) == 'total') {
      logging::logdebug(
        "Using unscaled abundance as total abundances"
      )
    } else {
      stop("unscaled abundance column must be a feature name or 'total'")
    }
  }
  
  # now order both data and metadata with the same sample ordering
  logging::logdebug(
    "Reordering data/metadata to use same sample ordering")
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]
  if (!is.null(unscaled_abundance)) {
    unscaled_abundance <- unscaled_abundance[intersect_samples, , drop = FALSE]
  }
  param_list[["unscaled_abundance"]] <- unscaled_abundance
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata))
}

###########################################
# Compute the formula based on user input #
###########################################

maaslin_compute_formula <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  
  fixed_effects <- param_list[["fixed_effects"]]
  group_effects <- param_list[["group_effects"]]
  ordered_effects <- param_list[["ordered_effects"]]
  random_effects <- param_list[["random_effects"]]
  
  if (!is.null(param_list[["formula"]])) {
    if (!is.null(param_list[["fixed_effects"]]) | 
        !is.null(param_list[["random_effects"]]) |
        !is.null(param_list[["group_effects"]]) | 
        !is.null(param_list[["ordered_effects"]])) {
      logging::logwarn(
        paste("fixed_effects, random_effects, group_effects, or ordered_effects provided in addition to formula,", 
              "using just the fixed_effects, random_effects, group_effects, and ordered_effects"))
    } else {
      logging::logwarn(
        paste("maaslin_compute_formula called even though a formula is provided",
              "without fixed_effects or random_effects,",
              "creating new formula based on metadata"))
    }
  }
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(fixed_effects)) {
    fixed_effects <- colnames(metadata)
  } else {
    fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
    # remove any fixed effects not found in metadata names
    to_remove <- setdiff(fixed_effects, colnames(metadata))
    if (length(to_remove) > 0)
      logging::logwarn(
        paste("Feature name not found in metadata",
              "so not applied to formula as fixed effect: %s"),
        paste(to_remove, collapse = " , ")
      )
    
    fixed_effects <- setdiff(fixed_effects, to_remove)
  }
  
  if (!is.null(random_effects)) {
    random_effects <-
      unlist(strsplit(random_effects, ",", fixed = TRUE))

    common_variables <- intersect(fixed_effects, random_effects)
    if (length(common_variables) > 0) {
      logging::logwarn(
        paste("Feature name included as fixed and random effect,",
              "check that this is intended: %s"),
        paste(common_variables, collapse = " , ")
      )
    }
        
    # remove any random effects not found in metadata
    to_remove <- setdiff(random_effects, colnames(metadata))
    if (length(to_remove) > 0) {
      logging::logwarn(
        paste("Feature name not found in metadata",
              "so not applied to formula as random effect: %s"),
        paste(to_remove, collapse = " , ")
      )
      random_effects <- setdiff(random_effects, to_remove)
    }
    
    # create formula
    if (length(random_effects) > 0) {
      random_effects_formula_text <-
        paste(
          "expr ~ (1 | ",
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
              paste(
                "Invalid formula for random effects: ",
                random_effects_formula_text
              )
            )
        )
    }
  }
  
  if (!is.null(group_effects) | !is.null(ordered_effects)) {
    multi_effects <- c()
    if (!is.null(group_effects)) {
      multi_effects <- c(multi_effects, strsplit(group_effects, ",", fixed = TRUE))
    }
    if (!is.null(ordered_effects)) {
      multi_effects <- c(multi_effects, strsplit(ordered_effects, ",", fixed = TRUE))
    }

    common_variables <- intersect(fixed_effects, multi_effects)
    if (length(common_variables) > 0) {
      logging::logerror(
        paste("Feature name included as fixed and group/ordered effect,",
              "this is not allowed: %s"),
        paste(common_variables, collapse = " , ")
      )
      stop()
    }
    
    # remove any random effects not found in metadata
    to_remove <- setdiff(multi_effects, colnames(metadata))
    if (length(to_remove) > 0) {
      logging::logerror(paste0("Feature name not found in metadata: ", 
                               paste0(to_remove, collapse = ", ")))
      stop()
    }
  }
  
  if (length(fixed_effects) == 0 & length(group_effects) == 0 & length(ordered_effects) == 0) {
    logging::logerror("No fixed/group/ordered effects provided.")
    stop()
  }
  
  # reduce metadata to only include fixed/group/random effects in formula
  effects_names <- unique(c(fixed_effects, random_effects, group_effects, ordered_effects))
  metadata <- metadata[, effects_names, drop = FALSE]
  
  # create the fixed effects formula text
  formula_effects <- fixed_effects
  if (length(group_effects) > 0) {
    formula_effects <- union(formula_effects, paste0("group(", group_effects, ")"))
  }
  if (length(ordered_effects) > 0) {
    formula_effects <- union(formula_effects, paste0("ordered(", ordered_effects, ")"))
  }
  
  formula_text <-
    paste("expr ~ ", paste(formula_effects, collapse = " + "))
  logging::loginfo("Formula for fixed effects: %s", formula_text)
  formula <-
    tryCatch(
      as.formula(formula_text),
      error = function(e)
        stop(
          paste(
            "Invalid formula.",
            "Please provide a different formula: ",
            formula_text
          )
        )
    )
  
  if (!(is.null(random_effects_formula))) {
    formula <-
      paste(
        '. ~', 
        paste(formula_effects, collapse = ' + '), 
        '.', 
        sep = ' + ')
    formula <- update(random_effects_formula, formula)
  }
  
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata, 
              "formula" = list("formula" = formula, 
                               "random_effects_formula" = random_effects_formula)))
}

##############################
# Check a user input formula #
##############################

maaslin_check_formula <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  
  input_formula <- param_list[["formula"]]
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(input_formula)) {
    logging::logwarn(
      paste("No user formula provided,",
            "building one instead"))
    return(maaslin_compute_formula(params_and_data))
  }
  
  if (!is.null(param_list[["fixed_effects"]]) | 
      !is.null(param_list[["random_effects"]]) | 
      !is.null(param_list[["group_effects"]]) | 
      !is.null(param_list[["ordered_effects"]])) {
    logging::logwarn(
      paste("fixed_effects, random_effects, group_effects, ordered_effects provided in addition to formula,", 
            "using only formula"))
  }
  
  # Remove anything before the tilde if necessary
  input_formula <- sub(".*~\\s*", "", input_formula)
  input_formula <- paste0("expr ~ ", input_formula)
  
  formula <-
    tryCatch(
      as.formula(input_formula),
      error = function(e)
        stop(
          paste(
            "Invalid formula: ",
            input_formula
          )
        )
    )
  
  formula_terms <- all.vars(formula)
  formula_terms <- formula_terms[formula_terms != "expr"]
  
  to_remove <- setdiff(formula_terms, colnames(metadata))
  if (length(to_remove) > 0) {
    logging::logerror(
      paste("Feature name not found in metadata: %s"),
      paste(to_remove, collapse = " , ")
    )
    stop()
  }
  
  term_labels <- attr(terms(formula), "term.labels")
  
  if (sum(!grepl("\\|", term_labels)) == 0) {
    logging::logerror("No fixed, group, or ordered effects included in formula.")
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
  
  # reduce metadata to only include fixed/group/ordered/random effects in formula
  metadata <- metadata[, formula_terms, drop = FALSE]
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata, 
              "formula" = list("formula" = formula, 
                               "random_effects_formula" = random_effects_formula)))
}

##############################################################################################
# Filter data based on min abundance, min prevalence, and min variance; standardize metadata #
##############################################################################################

maaslin_filter_and_standardize <- function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  data <- params_and_data_and_formula[["data"]]
  metadata <- params_and_data_and_formula[["metadata"]]
  
  #########################################################
  # Filter data based on min abundance and min prevalence #
  #########################################################
  
  reference <- param_list[["reference"]]
  if (is.null(reference)) {
    reference <- ","
  }
  split_reference <- unlist(strsplit(reference, "[,;]"))
  
  fixed_effects <- param_list[["fixed_effects"]]
  # for each fixed effect, check that a reference level has been set if necessary: 
  # number of levels > 2 and metadata isn't already an ordered factor
  for (i in fixed_effects) {
    # don't check for or require reference levels for numeric metadata
    if (is.numeric(metadata[,i])) {
      next
    }
    # respect ordering if a factor is explicitly passed in with no reference set
    if (is.factor(metadata[,i]) && !(i %in% split_reference)) {
      logging::loginfo(paste("Factor detected for categorial metadata '", 
                             i, "'. Provide a reference argument or manually set factor ordering to change reference level.", sep=""))
      next
    }
    
    # set metadata as a factor (ordered alphabetically)
    metadata[,i] <- as.factor(metadata[,i])
    mlevels <- levels(metadata[,i])
    
    # get reference level for variable being considered, returns NA if not found
    ref <- split_reference[match(i, split_reference)+1]
    
    # if metadata has 2 levels, allow but don't require setting reference level, otherwise require it
    if ((length(mlevels) == 2)) {
      if(!is.na(ref)) {
        metadata[,i] = relevel(metadata[,i], ref = ref)
      }
    } else if (length(mlevels) > 2) {
      if (!is.na(ref)) {
        metadata[,i] = relevel(metadata[,i], ref = ref)
      } else {
        stop(paste("Please provide the reference for the variable '",
                   i, "' which includes more than 2 levels: ",
                   paste(as.character(mlevels), collapse=", "), ". ",
                   "Alternatively, set the variable as a factor beforehand.", sep=""))   
      } 
    } else {
      stop("Provided categorical metadata has fewer than 2 unique, non-NA values.")
    }
  }
  
  unfiltered_data <- data
  unfiltered_metadata <- metadata
  
  # require at least total samples * min prevalence values 
  # for each feature to be greater than min abundance
  logging::loginfo(
    "Filter data based on min abundance and min prevalence")
  total_samples <- nrow(unfiltered_data)
  logging::loginfo("Total samples in data: %d", total_samples)
  min_samples <- total_samples * param_list[["min_prevalence"]]
  logging::loginfo(
    paste("Min samples required with min abundance",
          "for a feature not to be filtered: %f"),
    min_samples
  )
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data
  data_zeros[is.na(data_zeros)] <- 0
  filtered_data <-
    unfiltered_data[, 
                    colSums(data_zeros > param_list[["min_abundance"]]) > min_samples,
                    drop = FALSE]
  total_filtered_features <-
    ncol(unfiltered_data) - ncol(filtered_data)
  logging::loginfo("Total filtered features: %d", total_filtered_features)
  filtered_feature_names <-
    setdiff(names(unfiltered_data), names(filtered_data))
  logging::loginfo("Filtered feature names from abundance and prevalence filtering: %s",
                   toString(filtered_feature_names))
  
  ##################################################################################
  # Apply the non-zero abundance threshold to split the data into 0s and non-zeros #
  ##################################################################################
  
  prevalence_mask <- ifelse(filtered_data > param_list[["zero_threshold"]], 1, 0)
  filtered_data <- filtered_data * prevalence_mask
  
  #################################
  # Filter data based on variance #
  #################################
  
  vars <- apply(filtered_data, 2, var)
  variance_filtered_data <- filtered_data[, which(vars > param_list[["min_variance"]]), drop = FALSE]
  variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
  logging::loginfo("Total filtered features with variance filtering: %d", variance_filtered_features)
  variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
  logging::loginfo("Filtered feature names from variance filtering: %s",
                   toString(variance_filtered_feature_names))
  filtered_data <- variance_filtered_data

  #######################
  # Write filtered data #
  #######################
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_file = file.path(features_folder, "filtered_data.tsv")
  logging::loginfo("Writing filtered data to file %s", filtered_file)
  write.table(
    data.frame("feature" = rownames(filtered_data), filtered_data), 
    file = filtered_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  ################################
  # Standardize metadata, if set #
  ################################
  
  if (param_list[["standardize"]]) {
    logging::loginfo(
      "Applying z-score to standardize continuous metadata")
    metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
  } else {
    logging::loginfo("Bypass z-score application to metadata")
  }
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]],
              "data" = data, 
              "metadata" = metadata,
              "filtered_data" = filtered_data, 
              "unfiltered_metadata" = unfiltered_metadata, 
              "formula" = params_and_data_and_formula[["formula"]]))
}

#################
# Normalization #
#################

maaslin_normalize = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  features <- params_and_data_and_formula[["filtered_data"]]
  normalization <- param_list[["normalization"]]

  logging::loginfo(
    "Running selected normalization method: %s", normalization)
  
  if (normalization == 'TSS') {features <- TSSnorm(features)}
  if (normalization == 'CLR') {features <- CLRnorm(features)}
  if (normalization == 'CSS') {features <- CSSnorm(features)}
  if (normalization == 'TMM') {features <- TMMnorm(features)}
  if (normalization == 'NONE') {features <- NONEnorm(features)}
  
  if (!is.null(param_list[['unscaled_abundance']])) {
    features <- UNSCALEDnorm(features, param_list[['unscaled_abundance']])
  }
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_data_norm_file = file.path(features_folder, "filtered_data_norm.tsv")
  logging::loginfo("Writing filtered, normalized data to file %s", filtered_data_norm_file)
  write.table(
    data.frame("feature" = rownames(features), features), 
    file = filtered_data_norm_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]][,colnames(features), drop = F], 
              "filtered_data_norm" = features,
              "formula" = params_and_data_and_formula[["formula"]]))
}

##################
# Transformation #
##################

maaslin_transform = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  features <- params_and_data_and_formula[["filtered_data_norm"]]
  transformation <- param_list[["transform"]]
  
  logging::loginfo("Running selected transform method: %s", transformation)
  
  if (transformation == 'LOG') {features <- apply(features, 2, LOG)}
  if (transformation == 'LOGIT') {features <- apply(features, 2, LOGIT)}
  if (transformation == 'AST') {features <- apply(features, 2, AST)}
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_data_norm_transformed_file = file.path(features_folder, "filtered_data_norm_transformed.tsv")
  logging::loginfo("Writing filtered, normalized, transformed data to file %s", filtered_data_norm_transformed_file)
  write.table(
    data.frame("feature" = rownames(features), features), 
    file = filtered_data_norm_transformed_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]], 
              "filtered_data_norm_transformed" = features,
              "formula" = params_and_data_and_formula[["formula"]]))
}

##########
# Fit LM #
##########

maaslin_fit = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])

  logging::loginfo("Running the linear model (LM) component")
  
  prevalence_mask <- ifelse(params_and_data_and_formula[["filtered_data"]] > 0, 1, 0)
  
  #######################
  # For non-zero models #
  #######################
  
  fit_data_non_zero <-
    fit.model(
      features = params_and_data_and_formula[["filtered_data_norm_transformed"]],
      metadata = params_and_data_and_formula[["metadata"]],
      model = 'LM',
      formula = params_and_data_and_formula[["formula"]][["formula"]],
      random_effects_formula = params_and_data_and_formula[["formula"]][["random_effects_formula"]],
      correction = param_list[["correction"]],
      save_models = param_list[["save_models"]],
      augment = param_list[["augment"]],
      cores = param_list[["cores"]],
      median_comparison = param_list[["median_comparison_abundance"]],
      median_comparison_threshold = param_list[["median_comparison_abundance_threshold"]]
    )
  
  #################################################################
  # Count the total values for each feature (untransformed space) #
  #################################################################
  
  logging::loginfo("Counting total values for each feature")
  
  fit_data_non_zero$results$N <-
    apply(
      fit_data_non_zero$results,
      1,
      FUN = function(x)
        length(params_and_data_and_formula[["filtered_data_norm_transformed"]][, x[1]])
    )
  fit_data_non_zero$results$N.not.zero <-
    apply(
      fit_data_non_zero$results,
      1,
      FUN = function(x)
        length(which(params_and_data_and_formula[["filtered_data"]][, x[1]] != 0))
    )
  
  #####################
  # For binary models #
  #####################
  
  logging::loginfo("Running the logistic model component")
  
  fit_data_binary <-
    fit.model(
      features = prevalence_mask,
      metadata = params_and_data_and_formula[["metadata"]],
      model = 'logistic',
      formula = params_and_data_and_formula[["formula"]][["formula"]],
      random_effects_formula = params_and_data_and_formula[["formula"]][["random_effects_formula"]],
      correction = param_list[["correction"]],
      save_models = param_list[["save_models"]],
      augment = param_list[["augment"]],
      cores = param_list[["cores"]],
      median_comparison = param_list[["median_comparison_prevalence"]],
      median_comparison_threshold = param_list[["median_comparison_prevalence_threshold"]]
    )

  logging::loginfo("Counting total values for each feature")
  
  fit_data_binary$results$N <-
    apply(
      fit_data_binary$results,
      1,
      FUN = function(x)
        length(params_and_data_and_formula[["filtered_data_norm_transformed"]][, x[1]])
    )
  fit_data_binary$results$N.not.zero <-
    apply(
      fit_data_binary$results,
      1,
      FUN = function(x)
        length(which(params_and_data_and_formula[["filtered_data"]][, x[1]] != 0))
    )
  
  results <- add_joint_signif(fit_data_non_zero, fit_data_binary, 'LM', param_list[["correction"]])
  fit_data_non_zero$results <- results[[1]]
  fit_data_binary$results <- results[[2]]
  
  current_errors_for_likely_issues <- fit_data_binary$results$error[!is.na(fit_data_binary$results$N.not.zero) & 
                                                                      (fit_data_binary$results$N.not.zero < 50 &
                                                                         fit_data_binary$results$N.not.zero / fit_data_binary$results$N < 0.05) &
                                                                      ((!is.na(fit_data_binary$results$coef) & 
                                                                          abs(fit_data_binary$results$coef) > 15) |
                                                                         (!is.na(fit_data_binary$results$pval_individual) & 
                                                                         fit_data_binary$results$pval_individual < 10^-10))]
  fit_data_binary$results$error[!is.na(fit_data_binary$results$N.not.zero) & 
                                  (fit_data_binary$results$N.not.zero < 50 &
                                  fit_data_binary$results$N.not.zero / fit_data_binary$results$N < 0.05) &
                                  ((!is.na(fit_data_binary$results$coef) & 
                                      abs(fit_data_binary$results$coef) > 15) |
                                     (!is.na(fit_data_binary$results$pval_individual) & 
                                        fit_data_binary$results$pval_individual < 10^-10))] <- 
    ifelse(!is.na(current_errors_for_likely_issues),
           current_errors_for_likely_issues,
           "A large coefficient (>15 in absolute value) or small p-value (< 10^-10) was obtained from a feature present in <5% of samples. Check this is intended.")
  
  fit_data_non_zero$results <- fit_data_non_zero$results[order(fit_data_non_zero$results$qval_joint),]
  fit_data_non_zero$results <- fit_data_non_zero$results[order(!is.na(fit_data_non_zero$results$error)),] # Move all that had errors to the end
  
  fit_data_binary$results <- fit_data_binary$results[order(fit_data_binary$results$qval_joint),]
  fit_data_binary$results <- fit_data_binary$results[order(!is.na(fit_data_binary$results$error)),] # Move all that had errors to the end
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]], 
              "filtered_data_norm_transformed" = params_and_data_and_formula[["filtered_data_norm_transformed"]],
              "formula" = params_and_data_and_formula[["formula"]],
              "fit_data_non_zero" = fit_data_non_zero,
              "fit_data_binary" = fit_data_binary
              ))
}

###########################
# Write tables of results #
###########################

maaslin_write_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  
  # create an output folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  write_fits(params_data_formula_fit)
  write_results(params_data_formula_fit)
  
  return(params_data_formula_fit)
}

maaslin_write_results_lefse_format <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  
  # create an output folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  write_results_in_lefse_format(params_data_formula_fit$fit_data_non_zero$results, 
                                file.path(output, 'lefse_style_results_abundance.res'))
  write_results_in_lefse_format(params_data_formula_fit$fit_data_binary$results, 
                                file.path(output, 'lefse_style_results_prevalence.res'))
  
  return(params_data_formula_fit)
}

#######################################################
# Create visualizations for results passing threshold #
#######################################################

# TODO
maaslin_plot_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  
  # create an output folder and figures folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  if (param_list[["plot_heatmap"]] || param_list[["plot_scatter"]]) {
    figures_folder <- file.path(output, "figures")
    if (!file.exists(figures_folder)) {
      print("Creating output figures folder")
      dir.create(figures_folder)
    }
  }
  
  fit_out_lm <- params_data_formula_fit$fit_data_non_zero$results
  fit_out_lm <- fit_out_lm[c("feature", "value", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
  fit_out_lm$association <- "abundance"
  
  fit_out_binary <- params_data_formula_fit$fit_data_binary$results
  fit_out_binary <- fit_out_binary[c("feature", "value", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
  fit_out_binary$association <- "prevalence"
  
  fit_out <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))
  
  if (param_list[["plot_heatmap"]]) {
    heatmap_file <- file.path(output, "heatmap.pdf")
    logging::loginfo(
      "Writing heatmap of significant results to file: %s",
      heatmap_file)
    save_heatmap(fit_out, heatmap_file, figures_folder,
                 first_n = param_list[["heatmap_first_n"]],
                 max_significance = param_list[['max_significance']])
  }
  
  # if (param_list[["plot_scatter"]]) {
  #   logging::loginfo(
  #     paste("Writing association plots",
  #           "(one for each significant association)",
  #           "to output folder: %s"),
  #     output
  #   )
  #   significant_results_file <-
  #     file.path(output, "significant_results.tsv")
  #   saved_plots <- maaslin2_association_plots(
  #     params_data_formula_fit[["unfiltered_metadata"]],
  #     params_data_formula_fit[["filtered_data"]],
  #     significant_results_file,
  #     output,
  #     figures_folder,
  #     param_list[["max_pngs"]],
  #     param_list[["save_scatter"]])
  #   if (param_list[["save_scatter"]]) {
  #     scatter_file <- file.path(figures_folder, "scatter_plots.rds")
  #     # remove plots file if already exists
  #     if (file.exists(scatter_file)) {
  #       logging::logwarn(
  #         "Deleting existing scatter plot objects file: %s", scatter_file)
  #       unlink(scatter_file)
  #     }
  #     logging::loginfo("Writing scatter plot objects to file %s", scatter_file)
  #     saveRDS(saved_plots, file = scatter_file)   
  #   }
  # }
}

#######################################################
# Main maaslin3 function (defaults same command line) #
#######################################################

maaslin3 <- function(param_list = list()) {
  # Create log file, log arguments, and check arguments
  params_tmp <- maaslin_log_arguments(param_list) %>%
    maaslin_read_data() %>%
    maaslin_reorder_data()
  
  if (is.null(params_tmp[['param_list']][['formula']])) {
    params_tmp_2 <- params_tmp %>%
      maaslin_compute_formula()
  } else {
    params_tmp_2 <- params_tmp %>%
      maaslin_check_formula()
  }
  
  params_and_data_and_formula <- params_tmp_2 %>%
    maaslin_filter_and_standardize() %>%
    maaslin_normalize() %>%
    maaslin_transform()
  
  params_data_formula_fit <- params_and_data_and_formula %>%
    maaslin_fit()

  maaslin_write_results(params_data_formula_fit)
  
  # TODO: plot if called, write results if called
  maaslin_plot_results(params_data_formula_fit)
  
  if ('logging::writeToFile' %in% names(logging::getLogger()[['handlers']])) {
    logging::removeHandler('logging::writeToFile')
  }

  return(params_data_formula_fit)
}

###########################################################################
# If running on the command line, get arguments and call maaslin function #
###########################################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  # check three positional arguments are provided
  if (length(positional_args) != 3) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments",
            "<data.tsv> <metadata.tsv> <output_folder>")
    )
  }
  
  # call maaslin with the command line options
  fit_data <-
    Maaslin3(list(
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
      median_comparison_abundance = current_args$median_comparison_abundance,
      median_comparison_prevalence = current_args$median_comparison_prevalence,
      median_comparison_abundance_threshold = current_args$median_comparison_abundance_threshold,
      median_comparison_prevalence_threshold = current_args$median_comparison_prevalence_threshold,
      formula = current_args$formula,
      correction = current_args$correction,
      standardize = current_args$standardize,
      cores = current_args$cores,
      plot_heatmap = current_args$plot_heatmap,
      heatmap_first_n = current_args$heatmap_first_n,
      plot_scatter = current_args$plot_scatter,
      max_pngs = current_args$max_pngs,
      save_scatter = current_args$save_scatter,
      save_models = current_args$save_models,
      augment = current_args$augment,
      reference = current_args$reference,
      unscaled_abundance = current_args$unscaled_abundance
    ))
}


