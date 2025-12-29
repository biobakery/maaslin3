###############################################################################
# MaAsLin3 utility scripts

# Copyright (c) 2025 Harvard School of Public Health

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

######################
## TSS Normalization #
######################

TSSnorm <- function(features, zero_threshold) {
    # Convert to Matrix from Data Frame
    features_norm <- as.matrix(features)
    X_mask <- ifelse(features > zero_threshold, 1, 0)
    dd <- colnames(features_norm)
    
    ##############
    # From vegan #
    ##############
    
    x <- as.matrix(features_norm)
    if (any(x < 0, na.rm = TRUE)) {
        k <- min(x, na.rm = TRUE)
        warning("input data contains negative entries: result may be non-sense")
    } else {
        k <- .Machine$double.eps
    }
    
    MARGIN <- 1
    
    tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = TRUE))
    x <- sweep(x, MARGIN, tmp, "/")
    if (any(is.nan(x)))
        warning("result contains NaN, perhaps due to impossible mathematical\n
            operation\n")
    
    x <- ifelse(X_mask, x, NA)
    
    # Convert back to data frame
    features_TSS <- as.data.frame(x)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_TSS) <- dd
    
    # Return
    return(features_TSS)
}

######################
## CLR Normalization #
######################

# Only CLRs on the non-zero part of the vector
CLRnorm <- function(features, zero_threshold) {
    # Convert to Matrix from Data Frame
    features_norm <- as.matrix(features)
    dd <- colnames(features_norm)
    
    #####################
    # from chemometrics #
    #####################
    
    # CLR Normalizing the Data
    X <- features_norm
    X_mask <- ifelse(X > zero_threshold, 1, 0)
    Xgeom <-
        exp(1) ^ apply(X, 1, function(x) {
            mean(log(x[x > zero_threshold]))
        })
    features_CLR <- log(X / Xgeom)
    features_CLR <- ifelse(X_mask, features_CLR, NA)
    
    # Convert back to data frame
    features_CLR <- as.data.frame(features_CLR)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_CLR) <- dd
    
    # Return
    return(features_CLR)
}

######################
# NONE Normalization #
######################

NONEnorm <- function(features, zero_threshold) {
    X <- as.matrix(features)
    X_mask <- ifelse(X > zero_threshold, 1, 0)
    features_NONE <-
        data.frame(ifelse(X_mask, X, NA), check.names = FALSE)
    return(features_NONE)
}

##########################
# Unscaled Normalization #
##########################

UNSCALEDnorm <- function(features, abs_abundances, zero_threshold) {
    # Convert to Matrix from Data Frame
    if (colnames(abs_abundances) == 'total') {
        X_mask <- ifelse(features > zero_threshold, 1, 0)
        abs_mult_fact <- abs_abundances[rownames(features), 1]
    } else {
        abs_feature <- colnames(abs_abundances)
        X_mask <-
            ifelse(features[, colnames(features) != abs_feature] > 
                    zero_threshold, 1, 0)
        abs_mult_fact <-
            abs_abundances[rownames(features), 1] / features[, abs_feature]
        features <- features[, colnames(features) != abs_feature]
    }
    
    if (any(is.na(abs_mult_fact)) |
        any(is.infinite(abs_mult_fact))) {
        stop(
            '1+ of the unscaled abundance multipliers are NA or infinite.
        Check that the spike-in is present in all samples and the
        `unscaled_abundance` table is entirely non-NA.'
        )
    }
    
    features_norm <- as.matrix(features)
    dd <- colnames(features_norm)
    x <- as.matrix(features_norm)
    
    MARGIN <- 1
    
    x <- sweep(x, MARGIN, abs_mult_fact, "*")
    if (any(is.nan(x)))
        warning("result contains NaN, perhaps due to impossible mathematical\n
            operation\n")
    
    if (any(is.infinite(x)))
        warning("result contains Inf, perhaps due to impossible mathematical\n
            operation\n")
    
    x <- ifelse(X_mask, x, NA)
    
    # Convert back to data frame
    features_ABS <- as.data.frame(x)
    
    # Rename the features - Same Format as Before
    colnames(features_ABS) <- dd
    
    # Return
    return(features_ABS)
}

######################
# Log Transformation #
######################

LOG <- function(x) {
    if (any(x[!is.na(x)] <= 0)) {
        stop("Log transformation is only valid for values above 0")
    }
    return(log2(x))
}

#############################
# Pseudo-log Transformation #
#############################

PLOG <- function(x) {
    if (any(x[!is.na(x)] < 0)) {
        stop("Pseudo-log transformation is only valid for values >= 0")
    }
    min_observed <- min(x[!is.na(x) & x > 0])
    x[!is.na(x) & x == 0] <- min_observed / 2
    return(log2(x))
}

#################
# Write outputs #
#################

write_fits <- function(output,
                    fit_data_abundance,
                    fit_data_prevalence,
                    random_effects_formula = NULL,
                    save_models = FALSE) {
    fits_folder <- file.path(output, "fits")
    if (!file.exists(fits_folder)) {
        logging::loginfo("Creating fits folder")
        dir.create(fits_folder)
    }

    model_types <- c('linear', 'logistic')
    
    lapply(model_types, function(model_type) {
        fit_data <- if (model_type == 'linear') {
            fit_data_abundance
        } else {
            fit_data_prevalence
        }
        
        if (is.null(fit_data)) return(NULL)
        
        # Write out the raw model fits
        if (save_models) {
            model_file <- file.path(fits_folder, 
                                    paste0("models_", model_type, ".rds"))
            if (file.exists(model_file)) {
                logging::logwarn("Deleting existing model objects file: %s", 
                                    model_file)
                unlink(model_file)
            }
            logging::loginfo("Writing model objects to file %s", model_file)
            saveRDS(fit_data$fits, file = model_file)
        }
        
        # Write residuals to file
        residuals_file <- file.path(fits_folder, 
                                    paste0("residuals_", model_type, ".rds"))
        if (file.exists(residuals_file)) {
            logging::logwarn("Deleting existing residuals file: %s", 
                            residuals_file)
            unlink(residuals_file)
        }
        logging::loginfo("Writing residuals to file %s", residuals_file)
        saveRDS(fit_data$residuals, file = residuals_file)
        
        # Write fitted values to file
        fitted_file <- file.path(fits_folder, 
                                paste0("fitted_", model_type, ".rds"))
        if (file.exists(fitted_file)) {
            logging::logwarn("Deleting existing fitted file: %s", fitted_file)
            unlink(fitted_file)
        }
        logging::loginfo("Writing fitted values to file %s", fitted_file)
        saveRDS(fit_data$fitted, file = fitted_file)
        
        # Write extracted random effects to file (if specified)
        if (!is.null(random_effects_formula)) {
            ranef_file <- file.path(fits_folder, 
                                    paste0("ranef_", model_type, ".rds"))
            if (file.exists(ranef_file)) {
                logging::logwarn("Deleting existing ranef file: %s", ranef_file)
                unlink(ranef_file)
            }
            logging::loginfo("Writing extracted random effects to file %s", 
                            ranef_file)
            saveRDS(fit_data$ranef, file = ranef_file)
        }
    })
}

write_results <- function(output,
                        fit_data_abundance,
                        fit_data_prevalence,
                        max_significance = 0.1) {
    if (is.null(fit_data_abundance)) {
        fit_data <- fit_data_prevalence$results
    } else if (is.null(fit_data_prevalence)) {
        fit_data <- fit_data_abundance$results
    } else {
        fit_data <- rbind(fit_data_abundance$results,
                        fit_data_prevalence$results)
    }
    
    fit_data$model <- data.table::fcase(
            fit_data$model == 'linear', 'abundance',
            fit_data$model == 'logistic',  'prevalence'
    )
    
    small_random_effects_warning <- 
        paste0("<4 average observations per random effect group often inflates",
            " coefficients and deflates p-values: consider setting ",
            "small_random_effects=TRUE and see tutorial")
    
    fit_data <- fit_data[order(fit_data$qval_individual), ]
    # Move all that had errors to the end
    fit_data <- fit_data[order(!is.na(fit_data$error) & 
        fit_data$error != small_random_effects_warning), ] 
    
    #############################
    # Write all results to file #
    #############################
    
    results_file <- file.path(output, "all_results.tsv")
    
    logging::loginfo(
        paste(
            "Writing all the results to file (ordered 
            by increasing individual q-values): %s"
        ),
        results_file
    )
    
    write.table(
        fit_data,
        file = results_file,
        sep = "\t",
        quote = TRUE,
        row.names = FALSE
    )
    
    ###########################################
    # Write results passing threshold to file #
    ###########################################
    
    significant_results <-
        fit_data[(fit_data$qval_joint <= max_significance &
                    !is.na(fit_data$qval_joint)) | 
                    (fit_data$qval_individual <= max_significance &
                    !is.na(fit_data$qval_individual)) &
                    (is.na(fit_data$error) | 
                    fit_data$error == small_random_effects_warning),]
    significant_results_file <-
        file.path(output, "significant_results.tsv")
    
    logging::loginfo(
        paste(
            "Writing the significant results without errors",
            "(those which have joint q-values less than",
            "or equal to the threshold",
            "of %f ) to file (ordered by increasing individual q-values): %s"
        ),
        max_significance,
        significant_results_file
    )
    
    write.table(
        significant_results,
        file = significant_results_file,
        sep = "\t",
        quote = TRUE,
        row.names = FALSE
    )
}

write_results_in_lefse_format <-
    function(results, output_file_name) {
        lines_vec <- vapply(seq(nrow(results)), function(i) {
            if (is.na(results[i, ]$error) & 
                !is.na(results[i, ]$qval_individual)) {
                if (results[i, ]$qval_individual < 0.1) {
                    return(paste0(
                        c(
                            results[i, ]$feature,
                            results[i, ]$coef,
                            paste0(results[i, ]$metadata, '_', 
                                    results[i, ]$value),
                            results[i, ]$coef,
                            results[i, ]$pval_individual
                        ),
                        collapse = '\t'
                    ))
                } else {
                    return(paste0(c(results[i, ]$feature, 
                            results[i, ]$coef, "", "", "-"), collapse = '\t'))
                }
            } else {
                # return NA for rows that don't meet the condition
                return(NA_character_)
            }
        }, FUN.VALUE = character(1))
        lines_vec <- lines_vec[lines_vec != 'FALSE']
        
        writeLines(sort(lines_vec), con = output_file_name)
    }

####################
# Contrast testing #
####################

# add_qvals
# create_combined_pval
# flag_abundance_turned_prevalence
# add_joint_signif
# append_joint
# Check for highly significant likely model misfits
# Warn about prevalence associations induced by abundances changes
# Small random effects warning

maaslin_contrast_test <- function(
                        maaslin3_fit,
                        contrast_mat,
                        rhs = NULL,
                        max_significance = 0.1,
                        correction = 'BH',
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        subtract_median = FALSE,
                        small_random_effects = FALSE,
                        evaluate_only = NULL) {
    
    match.arg(correction, correction_choices)
    
    # Match variable ignoring case then set correctly as required for p.adjust
    correction <- correction_choices[match(toupper(correction),
                                        toupper(correction_choices))]
    
    # Run linear model component
    if (is.null(evaluate_only) || evaluate_only == "abundance") {
        if (is.null(maaslin3_fit$fit_data_abundance$fits)) {
            error_message <- paste0(
                "maaslin3_fit$fit_data_abundance$fits cannot be null,",
                " check if evaluate_only should be 'prevalence'")
            stop(error_message)
        }
        
        #######################
        # For non-zero models #
        #######################
        
        fit_data_abundance <- list()
        fit_data_abundance$results <- maaslin_contrast_test_runner(
            maaslin3_fit$fit_data_abundance$fits,
            contrast_mat,
            rhs = rhs,
            median_comparison = median_comparison_abundance)

        #################################################################
        # Count the total values for each feature (untransformed space) #
        #################################################################
        
        fit_data_abundance$results$N <- as.numeric(fmapvalues(
            fit_data_abundance$results$feature,
            maaslin3_fit$fit_data_abundance$results$feature,
            maaslin3_fit$fit_data_abundance$results$N))
        
        fit_data_abundance$results$N_not_zero <- as.numeric(fmapvalues(
            fit_data_abundance$results$feature,
            maaslin3_fit$fit_data_abundance$results$feature,
            maaslin3_fit$fit_data_abundance$results$N_not_zero))
    }
    
    # Run logistic model component
    if (is.null(evaluate_only) || evaluate_only == "prevalence") {
        if (is.null(maaslin3_fit$fit_data_prevalence$fits)) {
            error_message <- paste0(
                "maaslin3_fit$fit_data_prevalence$fits cannot be null,",
                " check if evaluate_only should be 'abundance'")
            stop(error_message)
        }
        
        #####################
        # For binary models #
        #####################
        
        fit_data_prevalence <- list()
        fit_data_prevalence$results <- maaslin_contrast_test_runner(
            maaslin3_fit$fit_data_prevalence$fits,
            contrast_mat,
            rhs = rhs,
            median_comparison = median_comparison_prevalence)
        
        fit_data_prevalence$results$N <- as.numeric(fmapvalues(
            fit_data_prevalence$results$feature,
            maaslin3_fit$fit_data_prevalence$results$feature,
            maaslin3_fit$fit_data_prevalence$results$N))
        
        fit_data_prevalence$results$N_not_zero <- as.numeric(fmapvalues(
            fit_data_prevalence$results$feature,
            maaslin3_fit$fit_data_prevalence$results$feature,
            maaslin3_fit$fit_data_prevalence$results$N_not_zero))
    }
    
    # Check for highly significant likely model misfits
    if (is.null(evaluate_only) || evaluate_only == "prevalence") {
        current_likely_error_subsetter <-
            !is.na(fit_data_prevalence$results$N_not_zero) &
            (
                fit_data_prevalence$results$N_not_zero < 50 &
                    fit_data_prevalence$results$N_not_zero /
                    fit_data_prevalence$results$N < 0.05
            ) &
            ((
                !is.na(fit_data_prevalence$results$coef) &
                    abs(fit_data_prevalence$results$coef) > 15
            ) |
                (
                    !is.na(fit_data_prevalence$results$pval) &
                        fit_data_prevalence$results$pval < 10 ^ -10
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
            !is.na(fit_data_prevalence$results$N_not_zero) &
            (
                fit_data_prevalence$results$N -
                    fit_data_prevalence$results$N_not_zero < 50 &
                    fit_data_prevalence$results$N_not_zero /
                    fit_data_prevalence$results$N > 0.95
            ) &
            ((
                !is.na(fit_data_prevalence$results$coef) &
                    abs(fit_data_prevalence$results$coef) > 15
            ) |
                (
                    !is.na(fit_data_prevalence$results$pval) &
                        fit_data_prevalence$results$pval < 10 ^ -10
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
    } else {
        fit_data_prevalence <- NULL
    }
    
    results <- add_qvals(fit_data_abundance,fit_data_prevalence, correction)
    if (!is.null(fit_data_abundance)) {
        fit_data_abundance$results <- results[[1]]
    }
    if (!is.null(fit_data_prevalence)) {
        fit_data_prevalence$results <- results[[2]]
    }
    if (!is.null(fit_data_abundance)) {
        fit_data_abundance$results <- fit_data_abundance$results |>
            collapse::fmutate(null_hypothesis = rhs)
        
        if (subtract_median & median_comparison_abundance) {
            fit_data_abundance$results$coef <- 
                fit_data_abundance$results$coef - 
                fit_data_abundance$results$null_hypothesis
            fit_data_abundance$results$null_hypothesis <- 0
        }
    }
    if (!is.null(fit_data_prevalence)) {
        fit_data_prevalence$results <- fit_data_prevalence$results |>
            collapse::fmutate(null_hypothesis = rhs)
        
        if (subtract_median & median_comparison_prevalence) {
            fit_data_prevalence$results$coef <- 
                fit_data_prevalence$results$coef - 
                fit_data_prevalence$results$null_hypothesis
            fit_data_prevalence$results$null_hypothesis <- 0
        }
    }
    
    # Add in joint p/q-values
    if (is.null(evaluate_only)) {
        if (!is.null(fit_data_abundance)) {
            
            dtest = fit_data_abundance$results$test
            
            fit_data_abundance$results <- fit_data_abundance$results |> 
                collapse::fmutate(metadata = dtest,
                                  value    = dtest,
                                  name     = dtest)
            
        }
        if (!is.null(fit_data_prevalence)) {
            
            dtest = fit_data_prevalence$results$test
            
            fit_data_prevalence$results <- fit_data_prevalence$results |> 
                collapse::fmutate(metadata = dtest,
                                  value    = dtest,
                                  name     = dtest)
        }
        results <-
            add_joint_signif(fit_data_abundance,
                            fit_data_prevalence,
                            NULL,
                            max_significance,
                            correction)
        if (!is.null(fit_data_abundance)) {
            fit_data_abundance$results <- results[[1]]
            
            collapse::fselect(fit_data_abundance$results,
                              c("metadata", "value", "name")) <- NULL
        }
        if (!is.null(fit_data_prevalence)) {
            fit_data_prevalence$results <- results[[2]]
            collapse::fselect(fit_data_prevalence$results,
                              c("metadata", "value", "name")) <- NULL
        }
    } else if (evaluate_only == 'abundance') {
        if (!is.null(fit_data_abundance)) {
            fit_data_abundance$results$pval_joint <-
                fit_data_abundance$results$pval
            
            fit_data_abundance$results$qval_joint <-
                fit_data_abundance$results$qval
            
            collapse::setrename(fit_data_abundance$results,
                                 "pval_individual" = "pval",
                                 "qval_individual" = "qval")
        }
    } else if (evaluate_only == 'prevalence') {
        fit_data_prevalence$results$pval_joint <-
            fit_data_prevalence$results$pval
        
        fit_data_prevalence$results$qval_joint <-
            fit_data_prevalence$results$qval
        
        collapse::setrename(fit_data_prevalence$results,
                             "pval_individual" = "pval",
                             "qval_individual" = "qval")
    }
    
    if (!is.null(fit_data_prevalence)) {
        used_random_effects <- 
            any(unlist(lapply(maaslin3_fit$fit_data_prevalence$fits, 
                FUN = function(x){is(x,"glmerMod")})))
        
        if (used_random_effects) {
            if (any(grepl('<4 average observations per random effect', 
                maaslin3_fit$fit_data_prevalence$results$error))) {
                fit_data_prevalence$results$error <- 
                    ifelse(is.na(fit_data_prevalence$results$error),
                            paste0("<4 average observations per random ",
                                    "effect group often inflates ",
                                    "coefficients and deflates p-values: ",
                                    "consider setting ",
                                    "small_random_effects=TRUE and see ",
                                    "tutorial"), 
                            fit_data_prevalence$results$error)
            }
        }
    }
    
    # Set column order
    col_order <- c("feature", "test", "coef", 
                    "null_hypothesis", "stderr",
                    "pval_individual", "qval_individual", "pval_joint",
                    "qval_joint", "error", "model", "N", "N_not_zero")
    
    # Reorder outputs
    if (is.null(evaluate_only) || evaluate_only == "abundance") {
        fit_data_abundance$results <-
            fit_data_abundance$results[
                order(fit_data_abundance$results$qval_individual), ]
        # Move all that had errors to the end
        fit_data_abundance$results <-
            fit_data_abundance$results[
                order(!is.na(fit_data_abundance$results$error)), ]
        
        # Reorder columns
        fit_data_abundance$results <- fit_data_abundance$results[,col_order]
    } else {
        fit_data_abundance <- NULL
    }
    
    if (is.null(evaluate_only) || evaluate_only == "prevalence") {
        fit_data_prevalence$results <-
            fit_data_prevalence$results[
                order(fit_data_prevalence$results$qval_individual), ]
        # Move all that had errors to the end
        fit_data_prevalence$results <-
            fit_data_prevalence$results[
                order(!is.na(fit_data_prevalence$results$error)), ]
        
        # Reorder columns
        fit_data_prevalence$results <- fit_data_prevalence$results[,col_order]
    } else {
        fit_data_prevalence <- NULL
    }
    
    return(
        list(
            "fit_data_abundance" = fit_data_abundance,
            "fit_data_prevalence" = fit_data_prevalence
        )
    )
}

maaslin_contrast_test_runner <- function(fits,
                                contrast_mat,
                                rhs = NULL,
                                median_comparison = NULL) {
    
    if (is.null(fits)) {
        stop("fits must not be null")
    }
    
    if (is.null(colnames(contrast_mat))) {
        stop("contrast_mat columns must be named to match model variables")
    }
    
    # Identify model type
    if (any(unlist(lapply(fits, FUN = 
                                function(x){is(x,"lm") & !is(x,"glm")})))) {
        model <- "abundance"
        uses_random_effects <- FALSE
    } else if (any(unlist(lapply(fits, FUN = function(x){is(x,"glm")})))) {
        model <- "prevalence"
        uses_random_effects <- FALSE
    } else if (any(unlist(lapply(fits, FUN = function(x){is(x,"lmerMod")})))) {
        model <- "abundance"
        uses_random_effects <- TRUE
    } else if (any(unlist(lapply(fits, FUN = function(x){is(x,"glmerMod")})))) {
        model <- "prevalence"
        uses_random_effects <- TRUE
    } else if (any(unlist(lapply(fits, FUN = function(x){is(x,"coxph")})))) {
        model <- "prevalence"
        uses_random_effects <- FALSE
    }

    if (is.null(median_comparison)) {
        if (model == 'abundance') {
            median_comparison <- TRUE
        } else {
            median_comparison <- FALSE
        }
    }
    
    if (median_comparison & !is.null(rhs)) {
        stop("rhs must be NULL with median_comparison")
    }
    
    if (is.null(rhs)) {
        rhs <- rep(0, nrow(contrast_mat))
    }
    
    if (nrow(contrast_mat) != length(rhs)) {
        stop("Number of rows in contrast matrix must equal the length of rhs")
    }
    
    # Run contrast test on each  model
    test_out_joined <- do.call(rbind, 
        lapply(seq_along(fits), function(fit_index) {
                fit <- fits[[fit_index]]
        if (!is(fit, 'lmerMod') & !is(fit, 'glmerMod') & !is(fit, "coxph")) {
            if (all(is.na(fit))) {
                return(matrix(rep(c(NA, NA, NA, "Fit is NA"), 
                            nrow(contrast_mat)),
                            nrow = nrow(contrast_mat), byrow = TRUE))
            }
        }
        
        # Fill in contrast matrix gaps if necessary
        if (uses_random_effects) {
            included_coefs <- names(lme4::fixef(fit))
        } else {
            included_coefs <- names(coef(fit, complete = FALSE))
        }
        contrast_mat_cols <- colnames(contrast_mat)
        contrast_mat_tmp <- contrast_mat
        contrast_mat_tmp <- cbind(contrast_mat_tmp, 
                                matrix(0, 
                                    nrow = nrow(contrast_mat_tmp),
                                    ncol = length(
                                    setdiff(included_coefs,
                                    contrast_mat_cols))))
        colnames(contrast_mat_tmp) <- c(contrast_mat_cols,
                                        setdiff(included_coefs,
                                                contrast_mat_cols))
        if (any(contrast_mat_tmp[,setdiff(contrast_mat_cols, 
                                        included_coefs)] != 0)) {
            return(matrix(rep(c(NA, NA, NA, 
                "Predictors not in the model had non-zero contrast values"), 
                    nrow(contrast_mat)),
                    nrow = nrow(contrast_mat), byrow = TRUE))
        }
        
        contrast_mat_tmp <- contrast_mat_tmp[,included_coefs, drop = FALSE]
        
        # Run contrast test
        if (!uses_random_effects | model == 'prevalence') {
            test_out_joined <- vapply(X = seq(nrow(contrast_mat)), 
                FUN = function(row_num) {
                    contrast_vec <- t(matrix(contrast_mat_tmp[row_num,]))
                    
                    error_message <- NA
                    calling_env <- environment()
                    test_out <- tryCatch({
                        if (!uses_random_effects) {
                            summary_out <- summary(multcomp::glht(
                                fit,
                                linfct = contrast_vec,
                                rhs = rhs[row_num],
                                coef. = function(x) {
                                    coef(x, complete = FALSE)
                                }
                            )
                            )$test
                        } else {
                            summary_out <- summary(multcomp::glht(
                                fit,
                                linfct = contrast_vec,
                                rhs = rhs[row_num],
                            )
                            )$test
                        }
                        
                        c(summary_out$coefficients,
                            summary_out$sigma,
                            summary_out$pvalues)
                    }, warning = function(w) {
                        message(sprintf("Feature %s : %s", 
                                        names(fit), w))
                        
                        assign("error_message",
                                conditionMessage(w),
                                envir = calling_env)
                        invokeRestart("muffleWarning")
                    },
                    error = function(err) {
                        assign("error_message", err$message, 
                            envir = calling_env)
                        error_obj <- c(NA, NA, NA)
                        return(error_obj)
                    })
                    return(as.character(c(test_out, error_message)))
                }, character(4))
        } else {
            test_out_joined <- vapply(X = seq(nrow(contrast_mat)), 
                FUN = function(row_num) {
                    contrast_vec <- t(matrix(contrast_mat_tmp[row_num,]))
                    
                    error_message <- NA
                    calling_env <- environment()
                    test_out <- tryCatch({
                        pval <- lmerTest::contest(fit,
                            matrix(contrast_vec, TRUE), 
                            rhs = rhs[row_num])[['Pr(>F)']]
                        coef <- contrast_vec %*% lme4::fixef(fit)
                        sigma <- sqrt((contrast_vec %*% vcov(fit) %*% 
                                        t(contrast_vec))[1, 1])
                        c(coef, sigma, pval)
                    }, warning = function(w) {
                        message(sprintf("Feature %s : %s", names(fit), w))
                        
                        assign("error_message",
                                conditionMessage(w),
                                envir = calling_env)
                        invokeRestart("muffleWarning")
                    },
                    error = function(err) {
                        assign("error_message", err$message,
                            envir = calling_env)
                        error_obj <- c(NA, NA, NA)
                        return(error_obj)
                    })
                    return(as.character(c(test_out, error_message)))
                }, character(4))
        }
        return(t(test_out_joined))
    }))
    
    test_out_joined <- data.frame(test_out_joined, check.names = FALSE)
    colnames(test_out_joined) <- c("coefs_new", "sigmas_new", 
        "pvals_new", "errors")
    test_out_joined[, seq(3)] <- apply(test_out_joined[, seq(3)], 2, as.numeric)

    if (all(is.null(rownames(contrast_mat)))) {
        test_names <- seq(nrow(contrast_mat))
    } else {
        test_names <- rownames(contrast_mat)
    }
    
    # Prepare results for return
    paras <- data.frame(
        feature = rep(names(fits), each = nrow(contrast_mat)),
        test = rep(test_names, length(fits)),
        coef = test_out_joined$coefs_new,
        rhs = rep(rhs, length(fits)),
        stderr = test_out_joined$sigmas_new,
        pval = test_out_joined$pvals_new,
        error = test_out_joined$errors,
        model = model
    )
    
    # Perform median comparison as though these were the original fits
    if (median_comparison) {
        paras <- do.call(rbind, lapply(unique(paras$test), 
                                        function(selected_test) {
            paras_sub <- paras[paras$test == selected_test, , drop = FALSE]
            paras_sub_out <- paras_sub[is.na(paras_sub$coef) | 
                                        is.na(paras_sub$stderr), , drop = FALSE]
            paras_sub <- paras_sub[!is.na(paras_sub$coef) & 
                                    !is.na(paras_sub$stderr), , drop = FALSE]
            if (nrow(paras_sub) == 0) {
                return(paras_sub_out)
            }
            
            # Make sure other p-values are set to NA to not confuse
            # median comparison with non-comparison
            paras_sub_out$pval <- rep(NA, 
                length(paras_sub_out$pval))
            
            use_this_coef <- !is.na(paras_sub$pval) & 
                paras_sub$pval < 0.95
            
            n_coefs <- nrow(paras_sub)
            sigmas <- paras_sub$stderr
            coefs <- paras_sub$coef
            sigma_sq_med <- var(coefs[use_this_coef], na.rm=TRUE)
            
            cur_median <- median(coefs[use_this_coef])
            
            # Variance from asymptotic distribution
            sd_median <- sqrt(0.25 * 2 * base::pi * 
                            sigma_sq_med / sum(use_this_coef))
            
            # MC for covariance
            nsims <- 10000
            
            sim_results <- replicate(nsims, {
                sim_coefs <- rnorm(n_coefs, coefs, sigmas)
                sim_median <- median(sim_coefs[use_this_coef])
                c(sim_median, sim_coefs)
            })
            
            sim_medians <- sim_results[1, ]
            all_sims <- sim_results[-1, , drop = FALSE]
            cov_adjust <- apply(all_sims, 1, function(x){cov(x, sim_medians)})
            
            # Necessary offsets for contrast testing
            offsets_to_test <- abs(cur_median - coefs) * 
                sqrt((sigmas^2) / (sigmas^2+ sd_median^2 - 2 * cov_adjust)) + 
                coefs
            
            pvals_new <- vapply(seq_len(nrow(paras_sub)), function(row_index) {
                feature <- paras_sub$feature[row_index]
                fit <- fits[[feature]]
                
                # Fill in contrast matrix gaps if necessary
                if (uses_random_effects) {
                    included_coefs <- names(lme4::fixef(fit))
                } else {
                    included_coefs <- names(coef(fit, complete = FALSE))
                }
                contrast_mat_cols <- colnames(contrast_mat)
                contrast_mat_tmp <- contrast_mat
                contrast_mat_tmp <- cbind(contrast_mat_tmp, 
                                        matrix(0, 
                                        nrow = nrow(contrast_mat_tmp),
                                        ncol = length(setdiff(
                                        included_coefs, contrast_mat_cols))))
                colnames(contrast_mat_tmp) <- 
                    c(contrast_mat_cols, 
                        setdiff(included_coefs, contrast_mat_cols))
                
                contrast_mat_tmp <- 
                    contrast_mat_tmp[, included_coefs, drop = FALSE]
                
                # Run contrast test
                contrast_vec <- t(matrix(contrast_mat_tmp[selected_test,]))
                
                error_message <- NA
                calling_env <- environment()
                test_out <- tryCatch({
                    if (!uses_random_effects | model == 'prevalence') {
                        summary_out <- summary(multcomp::glht(
                            fit,
                            linfct = contrast_vec,
                            rhs = offsets_to_test[row_index],
                            coef. = function(x) { coef(x, complete = FALSE) }
                        ))$test
                    } else {
                        summary_out <- summary(multcomp::glht(
                            fit,
                            linfct = contrast_vec,
                            rhs = offsets_to_test[row_index]
                        ))$test
                    }
                    
                    c(summary_out$pvalues, summary_out$coefficients, 
                        summary_out$sigma)
                }, warning = function(w) {
                    message(sprintf("Feature %s : %s", names(fit), w))
                    
                    assign("error_message", conditionMessage(w), 
                        envir = calling_env)
                    invokeRestart("muffleWarning")
                }, error = function(err) {
                    assign("error_message", err$message, envir = calling_env)
                    error_obj <- c(NA, NA, NA)
                    return(error_obj)
                })
                
                return(test_out[1])
            }, numeric(1))
            
            paras_sub$error <- 
                ifelse(!is.na(paras_sub$pval) & 
                        is.na(pvals_new),
                    "P-value became NA during median comparison",
                    paras_sub$error
                )
            
            paras_sub$pval <- pvals_new
            paras_sub$rhs <- cur_median
            
            return(rbind(paras_sub, paras_sub_out))
        }))
    }
    
    if (!is.null(rownames(paras)) && 
        !any(is.na(as.numeric(rownames(paras))))) {
        paras <- paras[order(as.numeric(rownames(paras))),]
    }
    
    return(paras)
}

##############################
# DNA pre-processing for MTX #
##############################

# rna_table must be samples (rows) by features (cols)
preprocess_dna_mtx <- function(dna_table, rna_table) {
    samples_row_row <-
        intersect(rownames(dna_table), rownames(rna_table))
    if (length(samples_row_row) == 0) {
        samples_column_row <-
            intersect(colnames(dna_table), rownames(rna_table))
        
        if (length(samples_column_row) > 0) {
            dna_table <- as.data.frame(t(dna_table))
        } else {
            stop(paste0(
                paste0("Rows/columns do not match."),
                paste0("DNA rows: ",
                    paste(
                        rownames(dna_table), collapse = ","
                    )),
                paste0("DNA columns: ",
                    paste(
                        colnames(dna_table), collapse = ","
                    )),
                paste0("RNA rows: ",
                    paste(
                        rownames(rna_table), collapse = ","
                    )),
                paste0("RNA columns: ",
                    paste(
                        colnames(rna_table), collapse = ","
                    )),
                collapse = '\n'
            ))
        }
    }
    intersect_samples <-
        intersect(rownames(dna_table), rownames(rna_table))
    
    dna_table <- dna_table[intersect_samples, , drop = FALSE]
    rna_table <- rna_table[intersect_samples, , drop = FALSE]
    
    # At this point, DNA and RNA tables are 
    # samples x features with same samples
    
    dna_table <- TSSnorm(dna_table, 0)
    dna_table <- apply(dna_table, 2, function(x) {
        x[is.na(x)] <- 0
        return(x)
    })
    dna_table <- as.data.frame(dna_table)

    rna_table <- TSSnorm(rna_table, 0)
    rna_table <- apply(rna_table, 2, function(x) {
        x[is.na(x)] <- 0
        return(x)
    })
    rna_table <- as.data.frame(rna_table)
    
    # Transforming DNA table
    impute_val <- log2(min(dna_table[dna_table > 0]) / 2)

    intersect_features <-
        intersect(colnames(dna_table), colnames(rna_table))
    
    dna_table <- dna_table[, intersect_features, drop = FALSE]
    rna_table <- rna_table[, intersect_features, drop = FALSE]
    
    if (!all(colnames(dna_table) == colnames(rna_table)) |
        !all(rownames(dna_table) == rownames(rna_table))) {
        stop("Something went wrong in preprocessing")
    }
    
    # Transform DNA 0s as necessary
    dna_table <- log2(dna_table)
    dna_table[dna_table == -Inf & rna_table > 0] <- impute_val
    dna_table[dna_table == -Inf & rna_table == 0] <- NA
    
    return(list("dna_table" = dna_table,
                "rna_table" = rna_table))
}

###############################
# Taxa pre-processing for MTX #
###############################

# rna_table must be samples (rows) by features (cols)
preprocess_taxa_mtx <- function(taxa_table, rna_table, rna_per_taxon) {
    if (any(colnames(rna_per_taxon) != c("RNA", 'taxon'))) {
        stop("colnames of rna_per_taxon should be RNA and taxon")
    }
    
    rna_vec <- c(unlist(rna_per_taxon[,'RNA']))
    taxon_vec <- c(unlist(rna_per_taxon[,'taxon']))

    samples_row_row <-
        intersect(rownames(taxa_table), rownames(rna_table))
    if (length(samples_row_row) == 0) {
        samples_column_row <-
            intersect(colnames(taxa_table), rownames(rna_table))
        
        if (length(samples_column_row) > 0) {
            taxa_table <- as.data.frame(t(taxa_table))
        } else {
            stop(paste0(
                paste0("Rows/columns do not match."),
                paste0("Taxa rows: ",
                    paste(rownames(taxa_table), collapse = ",")),
                paste0("Taxa columns: ",
                    paste(colnames(taxa_table), collapse = ",")),
                paste0("RNA rows: ",
                    paste(rownames(rna_table), collapse = ",")),
                paste0("RNA columns: ",
                    paste(colnames(rna_table), collapse = ",")),
                collapse = '\n'
            ))
        }
    }
    intersect_samples <-
        intersect(rownames(taxa_table), rownames(rna_table))
    
    taxa_table <- taxa_table[intersect_samples, , drop = FALSE]
    rna_table <- rna_table[intersect_samples, , drop = FALSE]
    
    if (any(!colnames(rna_table) %in% rna_vec)) {
        stop("rna_per_taxon RNA must contain all columns of rna_table")
    }
    if (any(!colnames(taxa_table) %in% taxon_vec)) {
        stop("rna_per_taxon taxon must contain all columns of taxa_table")
    }
    
    # At this point, taxa and RNA tables are 
    # samples x features with same samples
    
    taxa_table <- TSSnorm(taxa_table, 0)
    taxa_table <- apply(taxa_table, 2, function(x) {
        x[is.na(x)] <- 0
        return(x)
    })
    taxa_table <- as.data.frame(taxa_table)
    
    rna_table <- TSSnorm(rna_table, 0)
    rna_table <- apply(rna_table, 2, function(x) {
        x[is.na(x)] <- 0
        return(x)
    })
    rna_table <- as.data.frame(rna_table)
    
    # Create a dna table by choosing the taxa to match the rna table
    dna_table <- taxa_table[, 
        fmapvalues(colnames(rna_table), rna_vec, taxon_vec)]
    colnames(dna_table) <- colnames(rna_table)

    # Transforming DNA table
    impute_val <- log2(min(dna_table[dna_table > 0]) / 2)
    
    if (!all(colnames(dna_table) == colnames(rna_table)) |
        !all(rownames(dna_table) == rownames(rna_table))) {
        stop("Something went wrong in preprocessing")
    }
    
    # Transform DNA 0s as necessary
    dna_table <- log2(dna_table)
    
    max_rna_table <- vapply(split(rna_per_taxon$RNA, rna_per_taxon$taxon), 
        function(rna_cols) {
            rna_cols <- rna_cols[rna_cols %in% colnames(rna_table)]
            apply(rna_table[, rna_cols, drop = FALSE], 1, max, na.rm = TRUE)
        }, FUN.VALUE = numeric(nrow(rna_table)))
    max_rna_table <- max_rna_table[, 
        fmapvalues(colnames(rna_table), rna_vec, taxon_vec)]
    colnames(max_rna_table) <- colnames(rna_table)

    dna_table[dna_table == -Inf & max_rna_table > 0] <- impute_val
    if (any(dna_table == -Inf & rna_table == 0)) {
        dna_table[dna_table == -Inf & rna_table == 0] <- NA
    }
    
    return(list("dna_table" = dna_table,
                "rna_table" = rna_table))
}

##################################
# Read from SummarizedExperiment #
##################################

maaslin_read_summarized_experiment_data <- 
    function(summarized_experiment, assay.type = 1) {
    if (!inherits(summarized_experiment, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object")
    }
    
    data <- as.data.frame(t(SummarizedExperiment::assay(
        summarized_experiment, assay.type)))
    metadata <- as.data.frame(
        SummarizedExperiment::colData(summarized_experiment))
    return(list(
        "data" = data,
        "metadata" = metadata
    ))
}

fmapvalues <- function(x, from, to) {
    # To avoid dependency on plyr for mapvalues()
   
    stopifnot(length(from) == length(to)) 
    
    from_matches = collapse::fmatch(x, from)
    
    mch_id = collapse::whichNA(from_matches, invert = TRUE)
    
    to_set = to[from_matches[mch_id]]
    
    return(collapse::copyv(x, mch_id, to_set)) 
}
