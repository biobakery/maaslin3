######################
## TSS Normalization #
######################

TSSnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  X_mask <- ifelse(features > 0, 1, 0)
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
CLRnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  #####################
  # from chemometrics #
  #####################
  
  # CLR Normalizing the Data
  X <- features_norm
  X_mask <- ifelse(X > 0, 1, 0)
  Xgeom <- exp(1)^apply(X, 1, function(x) {mean(log(x[x > 0]))})
  features_CLR <- log(X/Xgeom)
  features_CLR <- ifelse(X_mask, features_CLR, NA)

  # Convert back to data frame
  features_CLR <- as.data.frame(features_CLR)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CLR) <- dd
  
  # Return
  return(features_CLR)
}

######################
## CSS Normalization #
######################

######################
# From metagenomeSeq #
######################

calcNormFactors <- function(x, p = cumNormStat(x)) {
  xx <- x
  xx[xx == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  normFactors <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  names(normFactors) <- colnames(x)
  as.data.frame(normFactors)
}

cumNormStat <- function (counts, qFlag = TRUE, pFlag = FALSE, rel = 0.1, ...) {
  mat = counts
  if (any(colSums(mat) == 0)) 
    stop("Warning empty sample")
  smat = sapply(1:ncol(mat), function(i) {
    sort(mat[, i], decreasing = FALSE)
  })
  ref = rowMeans(smat)
  yy = mat
  yy[yy == 0] = NA
  ncols = ncol(mat)
  refS = sort(ref)
  k = which(refS > 0)[1]
  lo = (length(refS) - k + 1)
  if (qFlag == TRUE) {
    diffr = sapply(1:ncols, function(i) {
      refS[k:length(refS)] - quantile(yy[, i], p = seq(0, 
                                                       1, length.out = lo), na.rm = TRUE)
    })
  }
  if (qFlag == FALSE) {
    diffr = sapply(1:ncols, function(i) {
      refS[k:length(refS)] - approx(sort(yy[, i], decreasing = FALSE), 
                                    n = lo)$y
    })
  }
  diffr2 = matrixStats::rowMedians(abs(diffr), na.rm = TRUE)
  if (pFlag == TRUE) {
    plot(abs(diff(diffr2[diffr2 > 0]))/diffr2[diffr2 > 0][-1], 
         type = "h", ylab = "Relative difference for reference", 
         xaxt = "n", ...)
    graphics::abline(h = rel)
    graphics::axis(1, at = seq(0, length(diffr2), length.out = 5), 
         labels = seq(0, 1, length.out = 5))
  }
  x = which(abs(diff(diffr2))/diffr2[-1] > rel)[1]/length(diffr2)
  if (x <= 0.5) {
    message("Default value being used.")
    x = 0.5
  }
  return(x)
}

cumNormMat <- function (x, p = cumNormStat(x), sl = 1000) {
  xx <- x
  xx[xx == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  newMat <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  nmat <- sweep(x, 2, newMat/sl, "/")
  return(nmat)
}

MRcounts <- function (counts, norm_factors, sl = 1000)  {
  if (any(is.na(norm_factors))) {
    x = cumNormMat(as.matrix(counts), sl = sl)
  } else {
    x = sweep(as.matrix(counts), 2, norm_factors/sl, "/")
  }
  return(x)
}

CSSnorm = function(features) {
  features_norm = as.matrix(features)
  X_mask <- ifelse(features_norm > 0, 1, 0)
  dd <- colnames(features_norm)
  
  counts = t(features_norm)
  counts = counts[,colSums(counts, na.rm = T) > 0]
  norm_factors <- calcNormFactors(counts)$normFactors
  features_CSS <- as.data.frame(t(MRcounts(counts, norm_factors)))
  
  features_CSS[setdiff(rownames(features_norm),rownames(features_CSS)),] <- NA
  features_CSS <- features_CSS[match(rownames(features_norm), rownames(features_CSS)), ]
  
  colnames(features_CSS) <- dd
  features_CSS <- data.frame(ifelse(X_mask > 0, as.matrix(features_CSS), NA))
  
  return(features_CSS)
}

######################
## TMM Normalization #
######################

##############
# From edgeR #
##############

.calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	Created 16 Aug 2010. Last modified 12 Sep 2020.
{
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
  if(min(f)==0) warning("One or more quantiles are zero")
  f / lib.size
}

.calcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM between two libraries
  #	Mark Robinson
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  
  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref
  
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  if(max(abs(logR)) < 1e-6) return(1)
  
  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)
  
  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

TMMnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  X_mask <- ifelse(features_norm > 0, 1, 0)
  dd <- colnames(features_norm)
  
  # TMM Normalizing the Data
  X <- t(features_norm)
  X = X[,colSums(X, na.rm = T) > 0]
  x <- as.matrix(X)
  if (any(is.na(x)))
    stop("NA counts not permitted")
  nsamples <- ncol(x)
  lib.size <- colSums(x)
  method <- "TMM"
  allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L
  if (any(allzero))
    x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || nsamples == 1)
    method = "none"
  
  f <- switch(method, TMM = {
    f75 <- suppressWarnings(.calcFactorQuantile(data = x, 
                                                lib.size = lib.size, p = 0.75))
    if (median(f75) < 1e-20) {
      refColumn <- which.max(colSums(sqrt(x)))
    } else {
      refColumn <- which.min(abs(f75 - mean(f75)))
    }
    f <- rep_len(NA_real_, nsamples)
    for (i in 1:nsamples) {
      f[i] <- .calcFactorTMM(obs = x[,i], ref = x[, refColumn], libsize.obs = lib.size[i], 
                             libsize.ref = lib.size[refColumn], logratioTrim = 0.3, 
                             sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
    }
    f
  }, 
  none = rep_len(1, nsamples))
  f <- f/exp(mean(log(f)))
  names(f) <- colnames(x)
  libSize <- f
  
  eff.lib.size = colSums(X) * libSize
  
  ref.lib.size = mean(eff.lib.size)
  #Use the mean of the effective library sizes as a reference library size
  X.output = sweep(X, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
  #Normalized read counts
  
  # Convert back to data frame
  features_TMM <- as.data.frame(t(X.output))
  
  features_TMM[setdiff(rownames(features_norm), rownames(features_TMM)),] <- NA
  features_TMM <- features_TMM[match(rownames(features_norm), rownames(features_TMM)), ]
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TMM) <- dd
  
  features_TMM <- data.frame(ifelse(X_mask > 0, as.matrix(features_TMM), NA))
  
  # Return as list
  return(features_TMM)
}

######################
# NONE Normalization #
######################

NONEnorm = function(features) {
  X <- as.matrix(features)
  X_mask <- ifelse(X > 0, 1, 0)
  features_NONE <- data.frame(ifelse(X_mask > 0, X, NA))
  return(features_NONE)
}

##########################
# Unscaled Normalization #
##########################

UNSCALEDnorm = function(features, abs_abundances) {
  # Convert to Matrix from Data Frame
  if (colnames(abs_abundances) == 'total') {
    X_mask <- ifelse(features > 0, 1, 0)
    abs_mult_fact <- abs_abundances[rownames(features),1]
  } else {
    abs_feature <- colnames(abs_abundances)
    X_mask <- ifelse(features[,colnames(features) != abs_feature] > 0, 1, 0)
    abs_mult_fact <- abs_abundances[rownames(features),1] / features[,abs_feature]
    features <- features[,colnames(features) != abs_feature]
  }
  
  if (any(is.na(abs_mult_fact)) | any(is.infinite(abs_mult_fact))) {
    stop('1+ of the unscaled abundance multipliers are NA or infinite. 
         Check that the spike-in is present in all samples and the 
         `unscaled_abundance` table is entirely non-NA.')
  }
  
  features_norm = as.matrix(features)
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

#######################################
# Arc-Sine Square Root Transformation #
#######################################

AST <- function(x) {
  input_na_count <- sum(is.na(x))
  y <- sign(x) * asin(sqrt(abs(x)))
  if(input_na_count < sum(is.na(y))) {
    logging::logerror(
      paste0("AST transform is only valid for values between -1 and 1. ",
             "Please select an appropriate normalization option or ",
             "normalize your data prior to running."))
    stop()
  }
  return(y)
}

########################
# Logit Transformation #
########################

# Zero-inflated Logit Transformation (Does not work well for microbiome data)
LOGIT <- function(p) {
  
  ########################
  # From the car package #
  ########################
  
  range.p <- range(p, na.rm = TRUE)
  if (range.p[2] > 1) {
    percents <- TRUE
    logging::loginfo("Note: largest value of p > 1 so values of p interpreted as percents")
  } else {
    percents <- FALSE
  }
  if (percents) {
    if (range.p[1] < 0 || range.p[2] > 100) 
      stop("p must be in the range 0 to 100")
    p <- p/100
    range.p <- range.p/100
  } else if (range.p[1] < 0 || range.p[2] > 1)  {
    stop("p must be in the range 0 to 1")
  }
  a <- 1
  y <- log((0.5 + a * (p - 0.5))/(1 - (0.5 + a * (p - 0.5))))
  if(any(y[!is.na(y)] == -Inf)) {
    stop("Logit transformation is only valid for values above 0")
  }
  return(y)
}

######################
# Log Transformation #
######################

LOG <- function(x) {
  if(any(x[!is.na(x)] <= 0)) {
    stop("Log transformation is only valid for values above 0")
  }
  return(log2(x))
}

############################
# Write out the model fits #
############################

write_fits <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]

  fits_folder <- file.path(output, "fits")
  if (!file.exists(fits_folder)) {
    print("Creating output fits folder")
    dir.create(fits_folder)
  }
  
  for (model_type in c('LM', 'logistic')) {
    if (model_type == 'LM') {
      fit_data <- params_data_formula_fit[["fit_data_abundance"]]
    } else {
      fit_data <- params_data_formula_fit[["fit_data_prevalence"]]
    }

    ################################
    # Write out the raw model fits #
    ################################
    
    if (param_list[["save_models"]]) {
      model_file = file.path(fits_folder, paste0("models_", model_type, ".rds"))
      # remove models file if already exists (since models append)
      if (file.exists(model_file)) {
        logging::logwarn(
          "Deleting existing model objects file: %s", model_file)
        unlink(model_file)
      }
      logging::loginfo("Writing model objects to file %s", model_file)
      saveRDS(fit_data$fits, file = model_file)   
    }
    
    ###########################
    # Write residuals to file #
    ###########################
    
    residuals_file = file.path(fits_folder, paste0("residuals_", model_type, ".rds"))
    # remove residuals file if already exists (since residuals append)
    if (file.exists(residuals_file)) {
      logging::logwarn(
        "Deleting existing residuals file: %s", residuals_file)
      unlink(residuals_file)
    }
    logging::loginfo("Writing residuals to file %s", residuals_file)
    saveRDS(fit_data$residuals, file = residuals_file)
    
    ###############################
    # Write fitted values to file #
    ###############################
    
    fitted_file = file.path(fits_folder, paste0("fitted_", model_type, ".rds"))
    # remove fitted file if already exists (since fitted append)
    if (file.exists(fitted_file)) {
      logging::logwarn(
        "Deleting existing fitted file: %s", fitted_file)
      unlink(fitted_file)
    }
    logging::loginfo("Writing fitted values to file %s", fitted_file)
    saveRDS(fit_data$fitted, file = fitted_file)
    
    #########################################################
    # Write extracted random effects to file (if specified) #
    #########################################################
    
    if (!is.null(param_list[["random_effects"]])) {
      ranef_file = file.path(fits_folder, paste0("ranef_", model_type, ".rds"))
      # remove ranef file if already exists (since ranef append)
      if (file.exists(ranef_file)) {
        logging::logwarn(
          "Deleting existing ranef file: %s", ranef_file)
        unlink(ranef_file)
      }
      logging::loginfo("Writing extracted random effects to file %s", ranef_file)
      saveRDS(fit_data$ranef, file = ranef_file)
    }
  }
}

write_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  max_significance <- param_list[["max_significance"]]
  fit_data <- rbind(params_data_formula_fit[["fit_data_abundance"]]$results,
                    params_data_formula_fit[["fit_data_prevalence"]]$results)
  
  fit_data$model <- dplyr::case_when(fit_data$model == 'LM' ~ 'abundance',
                                     fit_data$model == 'logistic' ~ 'prevalence',
                                     TRUE ~ NA)
  
  fit_data <- fit_data[order(fit_data$qval_joint),]
  fit_data <- fit_data[order(!is.na(fit_data$error)),] # Move all that had errors to the end

  #############################
  # Write all results to file #
  #############################
  
  results_file <- file.path(output, "all_results.tsv")
  
  logging::loginfo(
    paste("Writing all the results to file (ordered by increasing joint q-values): %s"),
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
  
  significant_results <- fit_data[fit_data$qval_joint <= max_significance & is.na(fit_data$error), ]
  significant_results$error <- NULL
  significant_results_file <- file.path(output, "significant_results.tsv")
  
  logging::loginfo(
    paste("Writing the significant results without errors",
          "(those which are less than or equal to the threshold",
          "of %f ) to file (ordered by increasing joint q-values): %s"),
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

write_results_in_lefse_format <- function(results, output_file_name) {
  lines_vec <- vector(length = nrow(results))
  for (i in 1:nrow(results)) {
    if (is.na(results[i,]$error) & !is.na(results[i,]$qval_individual)) {
      if(results[i,]$qval_individual < 0.1) {
        lines_vec[i] <- paste0(c(results[i,]$feature, 
                                 results[i,]$coef, 
                                 paste0(results[i,]$metadata, '_', results[i,]$value), 
                                 results[i,]$coef, 
                                 results[i,]$pval_individual), 
                               collapse = '\t')
      } else {
        lines_vec[i] <- paste0(c(results[i,]$feature, 
                                 results[i,]$coef, 
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
  samples_row_row <- intersect(rownames(dna_table), rownames(rna_table))
  if (length(samples_row_row) == 0) {
    samples_column_row <- intersect(colnames(dna_table), rownames(rna_table))
    
    if (length(samples_column_row) == 0) {
      # modify possibly included special chars in sample names in metadata
      rownames(rna_table) <- make.names(rownames(rna_table))
      
      samples_column_row <- intersect(colnames(dna_table), rownames(rna_table))
    }
    
    if (length(samples_column_row) > 0) {
      dna_table <- as.data.frame(t(dna_table))
      print("Transformed data so samples are rows")
    } else {
          print("Rows/columns do not match.")
          print(
            paste0("DNA rows: ", 
            paste(rownames(dna_table), collapse = ",")))
          print(
            paste0("DNA columns: ", 
            paste(colnames(dna_table), collapse = ",")))
          print(
            paste0("RNA rows: ", 
            paste(rownames(rna_table), collapse = ",")))
          print(
            paste0("RNA columns: ",
            paste(colnames(rna_table), collapse = ",")))
          stop()
    }
  }
  
  # replace unexpected characters in feature names
  colnames(dna_table) <- make.names(colnames(dna_table))

  intersect_samples <- intersect(rownames(dna_table), rownames(rna_table))
  print(paste0("A total of ", length(intersect_samples), 
               " samples were found in both the data and metadata"))
  
  # check for samples without RNA abundances
  extra_dna_samples <-
    setdiff(rownames(dna_table), intersect_samples)
  if (length(extra_dna_samples) > 0)
    print(
      paste("The following samples were found",
            "to have DNA but no RNA",
            "They will be removed. ", paste(extra_dna_samples, collapse = ","))
    )
  
  # check for samples without DNA abundances
  extra_rna_samples <-
    setdiff(rownames(rna_table), intersect_samples)
  if (length(extra_rna_samples) > 0)
    print(
      paste("The following samples were found",
            "to have RNA but no DNA.",
            "They will be removed. ", paste(extra_rna_samples, collapse = ","))
    )
  
  print("Reordering DNA/RNA to use same sample ordering")
  dna_table <- dna_table[intersect_samples, , drop = FALSE]
  rna_table <- rna_table[intersect_samples, , drop = FALSE]

  intersect_features <- intersect(colnames(dna_table), colnames(rna_table))
  print(paste0("A total of ", length(intersect_features), 
               " features were found in both the data and metadata"))
  
  # check for features without RNA abundances
  extra_dna_samples <-
    setdiff(colnames(dna_table), intersect_features)
  if (length(extra_dna_samples) > 0)
    print(
      paste("The following samples were found",
            "to have DNA but no RNA",
            "They will be removed. ", paste(extra_dna_samples, collapse = ","))
    )
  
  # check for features without DNA abundances
  extra_rna_samples <-
    setdiff(colnames(rna_table), intersect_features)
  if (length(extra_rna_samples) > 0)
    print(
      paste("The following samples were found",
            "to have RNA but no DNA.",
            "They will be removed. ", paste(extra_rna_samples, collapse = ","))
    )
  
  print("Reordering DNA/RNA to use same feature ordering")
  dna_table <- dna_table[, intersect_features, drop = FALSE]
  rna_table <- rna_table[, intersect_features, drop = FALSE]

  if (!all(colnames(dna_table) == colnames(rna_table)) | 
      !all(rownames(dna_table) == rownames(rna_table))) {
    stop("Something went wrong in preprocessing")
  }
  
  # At this point, DNA and RNA tables are samples x features with same features and samples
  
  dna_table <- TSSnorm(dna_table)
  for (col_index in 1:ncol(dna_table)) {
    dna_table[,col_index][is.na(dna_table[,col_index])] <- 0
  }

  rna_table <- TSSnorm(rna_table)
  for (col_index in 1:ncol(rna_table)) {
    rna_table[,col_index][is.na(rna_table[,col_index])] <- 0
  }
  
  # Transforming DNA table
  print("Transforming DNA table...this can take a while")
  for (row_name in rownames(dna_table)) {
    dna_table[row_name,] <- 
      ifelse(dna_table[row_name,] > 0,
             log2(dna_table[row_name,]), # If present, just log2 transform it
             ifelse(rna_table[row_name,] > 0, # Otherwise, if RNA is present, use pseudocount
                    log2(min(dna_table[row_name,][dna_table[row_name,] > 0]) / 2), # sample minimum / 2
                    NA) # If both are missing, set NA to exclude in analysis
            )
  }

  return(list("dna_table" = dna_table, 
              "rna_table" = rna_table))
}










