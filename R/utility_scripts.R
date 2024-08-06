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
        data.frame(ifelse(X_mask > zero_threshold, X, NA))
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

############################
# Write out the model fits #
############################

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
    
    
    for (model_type in c('LM', 'logistic')) {
        if (model_type == 'LM') {
            fit_data <- fit_data_abundance
        } else {
            fit_data <- fit_data_prevalence
        }
        
        if (is.null(fit_data)) {
            next
        }
        
        ################################
        # Write out the raw model fits #
        ################################
        
        if (save_models) {
            model_file <-
                file.path(fits_folder,
                        paste0("models_", model_type, ".rds"))
            # remove models file if already exists (since models append)
            if (file.exists(model_file)) {
                logging::logwarn("Deleting existing model objects file: %s",
                                model_file)
                unlink(model_file)
            }
            logging::loginfo("Writing model objects to file %s", model_file)
            saveRDS(fit_data$fits, file = model_file)
        }
        
        ###########################
        # Write residuals to file #
        ###########################
        
        residuals_file <-
            file.path(fits_folder, paste0("residuals_", model_type, ".rds"))
        # remove residuals file if already exists (since residuals append)
        if (file.exists(residuals_file)) {
            logging::logwarn("Deleting existing residuals file: %s",
                            residuals_file)
            unlink(residuals_file)
        }
        logging::loginfo("Writing residuals to file %s", residuals_file)
        saveRDS(fit_data$residuals, file = residuals_file)
        
        ###############################
        # Write fitted values to file #
        ###############################
        
        fitted_file <-
            file.path(fits_folder, paste0("fitted_", model_type, ".rds"))
        # remove fitted file if already exists (since fitted append)
        if (file.exists(fitted_file)) {
            logging::logwarn("Deleting existing fitted file: %s", fitted_file)
            unlink(fitted_file)
        }
        logging::loginfo("Writing fitted values to file %s", fitted_file)
        saveRDS(fit_data$fitted, file = fitted_file)
        
        #########################################################
        # Write extracted random effects to file (if specified) #
        #########################################################
        
        if (!is.null(random_effects_formula)) {
            ranef_file <-
                file.path(fits_folder, paste0("ranef_", model_type, ".rds"))
            # remove ranef file if already exists (since ranef append)
            if (file.exists(ranef_file)) {
                logging::logwarn("Deleting existing ranef file: %s", ranef_file)
                unlink(ranef_file)
            }
            logging::loginfo("Writing extracted random effects to file %s",
                            ranef_file)
            saveRDS(fit_data$ranef, file = ranef_file)
        }
    }
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
    
    fit_data$model <-
        dplyr::case_when(
            fit_data$model == 'LM' ~ 'abundance',
            fit_data$model == 'logistic' ~ 'prevalence',
            TRUE ~ NA
        )
    
    fit_data <- fit_data[order(fit_data$qval_joint), ]
    # Move all that had errors to the end
    fit_data <- fit_data[order(!is.na(fit_data$error)), ] 
    
    #############################
    # Write all results to file #
    #############################
    
    results_file <- file.path(output, "all_results.tsv")
    
    logging::loginfo(
        paste(
            "Writing all the results to file (ordered 
            by increasing joint q-values): %s"
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
        fit_data[fit_data$qval_joint <= max_significance &
                    is.na(fit_data$error),]
    significant_results$error <- NULL
    significant_results_file <-
        file.path(output, "significant_results.tsv")
    
    logging::loginfo(
        paste(
            "Writing the significant results without errors",
            "(those which are less than or equal to the threshold",
            "of %f ) to file (ordered by increasing joint q-values): %s"
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
        lines_vec <- vector(length = nrow(results))
        for (i in seq(nrow(results))) {
            if (is.na(results[i, ]$error) &
                !is.na(results[i, ]$qval_individual)) {
                if (results[i, ]$qval_individual < 0.1) {
                    lines_vec[i] <- paste0(
                        c(
                            results[i, ]$feature,
                            results[i, ]$coef,
                            paste0(results[i, ]$metadata, '_', 
                                results[i, ]$value),
                            results[i, ]$coef,
                            results[i, ]$pval_individual
                        ),
                        collapse = '\t'
                    )
                } else {
                    lines_vec[i] <- paste0(c(results[i, ]$feature,
                                            results[i, ]$coef,
                                            "",
                                            "",
                                            "-"),
                                        collapse = '\t')
                }
            }
        }
        lines_vec <- lines_vec[lines_vec != 'FALSE']
        
        writeLines(sort(lines_vec), con = output_file_name)
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
        
        if (length(samples_column_row) == 0) {
            # modify possibly included special chars in sample names in metadata
            rownames(rna_table) <- make.names(rownames(rna_table))
            
            samples_column_row <-
                intersect(colnames(dna_table), rownames(rna_table))
        }
        
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
    
    # replace unexpected characters in feature names
    colnames(dna_table) <- make.names(colnames(dna_table))
    
    intersect_samples <-
        intersect(rownames(dna_table), rownames(rna_table))
    
    # check for samples without RNA abundances
    extra_dna_samples <-
        setdiff(rownames(dna_table), intersect_samples)
    
    # check for samples without DNA abundances
    extra_rna_samples <-
        setdiff(rownames(rna_table), intersect_samples)
    
    dna_table <- dna_table[intersect_samples, , drop = FALSE]
    rna_table <- rna_table[intersect_samples, , drop = FALSE]
    
    intersect_features <-
        intersect(colnames(dna_table), colnames(rna_table))
    
    # check for features without RNA abundances
    extra_dna_samples <-
        setdiff(colnames(dna_table), intersect_features)
    
    # check for features without DNA abundances
    extra_rna_samples <-
        setdiff(colnames(rna_table), intersect_features)
    
    dna_table <- dna_table[, intersect_features, drop = FALSE]
    rna_table <- rna_table[, intersect_features, drop = FALSE]
    
    if (!all(colnames(dna_table) == colnames(rna_table)) |
        !all(rownames(dna_table) == rownames(rna_table))) {
        stop("Something went wrong in preprocessing")
    }
    
    # At this point, DNA and RNA tables are 
    # samples x features with same features and samples
    
    dna_table <- TSSnorm(dna_table, 0)
    for (col_index in seq(ncol(dna_table))) {
        dna_table[, col_index][is.na(dna_table[, col_index])] <- 0
    }
    
    rna_table <- TSSnorm(rna_table, 0)
    for (col_index in seq(ncol(rna_table))) {
        rna_table[, col_index][is.na(rna_table[, col_index])] <- 0
    }
    
    # Transforming DNA table
    impute_val <- log2(min(dna_table[dna_table > 0]) / 2)
    dna_table <- log2(dna_table)
    dna_table[dna_table == -Inf & rna_table > 0] <- impute_val
    dna_table[dna_table == -Inf & rna_table == 0] <- NA
    
    return(list("dna_table" = dna_table,
                "rna_table" = rna_table))
}
