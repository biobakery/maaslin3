#!/usr/bin/env Rscript
# Load Required Packages
for (lib in c(
  'dplyr',
  'pbapply',
  'lmerTest',
  'parallel',
  'lme4',
  'plyr',
  'TcGSA'
)) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# Function to augment data for logistic fitting
augment_data <- function(formula, random_effects_formula, dat_sub){
  if (is.null(random_effects_formula)) { # No random effects
    formula <- formula(formula)
    
    # Get model matrix
    mm <- data.frame(model.matrix(formula, dat_sub), check.names = F)
    org_data_length <- nrow(mm)
    org_col_names <- colnames(mm)
    
    # Add in feature abundance
    mm$expr <- as.numeric(mapvalues(rownames(mm), rownames(dat_sub), dat_sub$expr, warn_missing = F))
    
    # Generate the new rows
    new_rows <- matrix(nrow = 0, ncol = ncol(mm))
    for (colname in colnames(mm)[!(colnames(mm) %in% c("(Intercept)", "expr"))]) {
      col = which(colnames(mm) == colname)
      new_rows_addition <- matrix(colMeans(mm), nrow = 4, ncol = ncol(mm), byrow = T)
      new_rows_addition[,col] <- new_rows_addition[,col] + c(-1, 1, -1, 1) * sd(mm[,col])
      new_rows_addition[,ncol(new_rows_addition)] <- c(1, 1, 0, 0)
      new_rows <- rbind(new_rows, new_rows_addition)
    }
    
    # Append these new rows to the model matrix
    mm_names <- colnames(mm)
    mm <- rbind(as.matrix(mm), new_rows)
    colnames(mm) <- mm_names
    mm <- data.frame(mm, check.names = F)
    if ("(Intercept)" %in% colnames(mm)) { mm$`(Intercept)` <- NULL }
    
    # Calculate the weights
    p <- ncol(mm) - 1
    weight_scheme <- c(rep(1, org_data_length), rep((p + 1) / (4 * p), nrow(mm) - org_data_length))
    
    # Get ready to return augmented data
    mm_input = mm
    new_formula <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"))
    return(list(mm_input = mm_input, weight_scheme = weight_scheme, new_formula = formula(new_formula)))
  } else { # With random effects
    input_formula <- formula(formula)
    
    # Generate the model matrix for fixed effects
    lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
    mm <- data.frame(lmer_parsed$X, check.names = F)
    org_data_dim <- dim(mm)
    org_col_names <- colnames(mm)
    
    # Add on feature abundance
    mm$expr <- as.numeric(mapvalues(rownames(mm), rownames(dat_sub), dat_sub$expr, warn_missing = F))
    
    # Get random effects
    random_effects <- names(lmer_parsed$reTrms$cnms)
    
    # Function to duplicate a random effects variable matrix for the 4 data points added
    duplicate_rows <- function(mat) {
      mat <- as.matrix(mat)
      num_rows <- nrow(mat)
      duplicated_mat <- matrix(, ncol = ncol(mat), nrow = num_rows * 4)
      for (i in 1:num_rows) {
        duplicated_mat[((i - 1) * 4 + 1):(i * 4), ] <- matrix(rep(mat[i, ], each = 4), ncol = ncol(mat))
      }
      return(duplicated_mat)
    }
    
    # Add the new rows
    new_rows <- matrix(nrow = 0, ncol = ncol(mm) + length(random_effects))
    for (colname in colnames(mm)[!(colnames(mm) %in% c("(Intercept)", "expr"))]) {
      if (length(random_effects) > 1) { # Get the set of possible random effects combinations
        new_re <- expand.grid(lapply(dat_sub[,random_effects], unique))
      } else {
        new_re <- matrix(unique(dat_sub[,random_effects]))
      }
      col = which(colnames(mm) == colname)
      
      # Add the 4 new points per metadata and per random effect combination
      new_rows_addition <- matrix(colMeans(mm), nrow = 4 * nrow(new_re), ncol = ncol(mm), byrow = T)
      new_rows_addition[,col] <- new_rows_addition[,col] + rep(c(-1, 1, -1, 1) * sd(mm[,col]), nrow(new_re))
      new_rows_addition[,ncol(new_rows_addition)] <- rep(c(1, 1, 0, 0), nrow(new_re))
      new_rows_addition <- cbind(new_rows_addition, duplicate_rows(new_re))
      new_rows <- rbind(new_rows, new_rows_addition)
    }
    new_rows <- data.frame(new_rows)
    colnames(new_rows) <- c(colnames(mm), random_effects)
    
    # Add on other necessary columns for the random effect (e.g., when using random slopes)
    extra_cols <- colnames(dat_sub)[!(colnames(dat_sub) %in% colnames(new_rows))]
    if (length(extra_cols) > 0) {
      # Identify missed columns and build a model matrix
      mm_2 <- model.matrix(formula(paste0("expr ~ ", paste0(extra_cols, collapse = " + "))), dat_sub)
      mm_2 <- mm_2[,-1] # Remove intercept
      new_colnames <- colnames(mm_2)[!(colnames(mm_2) %in% colnames(new_rows))]
      
      # Use the model matrix to add missed columns
      if (length(new_colnames) > 0) {
        new_rows[,new_colnames] <- rep(colMeans(as.matrix(mm_2[,new_colnames])), each = nrow(new_rows))
        prev_colnames <- colnames(mm)
        mm <- cbind(mm, mm_2[,new_colnames])
        colnames(mm) <- c(prev_colnames, new_colnames)
      }
    }
    
    # Join the new rows with the original data
    mm[,random_effects] <- dat_sub[,random_effects, drop=F][rownames(mm),]
    new_rows <- new_rows[,order(mapvalues(colnames(new_rows), colnames(mm), 1:ncol(mm)))]
    mm <- rbind(mm, new_rows)
    mm <- data.frame(mm, check.names = F)
    
    # Turn numeric columns to numeric
    for (colname in colnames(mm)) {
      mm[,colname] <- tryCatch({
        as.numeric(mm[,colname])
      }, warning = function(w) { 
        return(mm[,colname])
      })
    }
    
    if ("(Intercept)" %in% colnames(mm)) { mm$`(Intercept)` <- NULL }
    
    # Establish weighting scheme
    p <- ncol(mm) - 1
    weight_scheme <- c(rep(1, org_data_dim[1]), 
                       rep(org_data_dim[2] / (nrow(mm) - org_data_dim[1]), nrow(mm) - org_data_dim[1]))
    
    # Create new formula based on model matrix
    new_formula <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                          paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
    
    return(list(mm_input = mm, weight_scheme = weight_scheme, new_formula = formula(new_formula)))
  }
}

safe_deparse <- function(formula) {
  paste0(trimws(deparse(formula)), collapse = " ")
}

extract_special_predictor <- function(formula, predictor_type) {
  groups <- regmatches(safe_deparse(formula), gregexpr(paste0(predictor_type, "\\((.*?)\\)"), safe_deparse(formula)))[[1]]
  groups <- gsub("\\)$", "", gsub(paste0("^", predictor_type, "\\("), "", groups))
  formula_tmp <- trimws(gsub(paste0("(\\s*\\+\\s*)?", predictor_type, "\\(.*?\\)"), "", safe_deparse(formula)))
  if (substr(formula_tmp, nchar(formula_tmp), nchar(formula_tmp)) == '~') {
    formula_tmp <- paste0(formula_tmp, '1')
  }
  formula <- formula(formula_tmp)
  formula <- formula(gsub("~ \\+", "~", safe_deparse(formula)))
  return(list(formula, predictor_type = groups))
}

# Get all fixed effects
get_fixed_effects <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  names_to_include <- c()
  if (is.null(random_effects_formula)) { # Fixed and group effects only
    names_to_include <- colnames(model.matrix(formula(gsub("^expr ", "", safe_deparse(formula))), dat_sub))
    names_to_include <- names_to_include[names_to_include != "(Intercept)"]
  } else { # Random effects
    pattern <- paste0("\\b", paste(gsub("\\|", "\\\\|", findbars(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = "|"), "\\b")
    
    fixed_effects_only <- gsub(pattern, "", 
                               paste0(trimws(safe_deparse(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = " "), 
                               ignore.case = TRUE)
    fixed_effects_only <- trimws(gsub("\\(\\)", "", fixed_effects_only))
    while(grepl("\\+$", fixed_effects_only)) {
      fixed_effects_only <- trimws(gsub("\\+$", "", fixed_effects_only))
    }
    
    names_to_include <- colnames(model.matrix(formula(fixed_effects_only), dat_sub))
    names_to_include <- names_to_include[names_to_include != "(Intercept)"]
  }
  omp_levels <- c()
  for (omp in omps) {
    omp_levels <- c(omp_levels, paste0(omp, levels(dat_sub[[omp]])[-1]))
  }
  
  names_to_include <- c(names_to_include, groups, gomps, omp_levels)
  return(names_to_include)
}

get_character_cols <- function(dat_sub) {
  all_factors <- c()
  for (col in colnames(dat_sub)) {
    if (!is.numeric(dat_sub[,col])) {
      factor_levels <- levels(factor(dat_sub[,col]))
      # All factor levels except basline
      all_factors <- c(all_factors, paste0(col, unique(factor_levels[-1])))
    }
  }
  return(all_factors)
}

# For group effects
generate_new_formula_groups <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps, group) {
  if (is.null(random_effects_formula)) {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group], gomps, omps), collapse = " + ")))
    mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
    org_col_names <- colnames(mm)
    colnames_to_insert <- org_col_names[org_col_names != "(Intercept)"]
    if (length(colnames_to_insert) > 0) {
      formula_new <- paste0("expr ~ ", paste0("`", paste0(colnames_to_insert, collapse = "` + `"), "`"))
    } else {
      formula_new <- paste0("expr ~ 1")
    }
  } else {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group], gomps, omps), collapse = " + ")))
    lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
    mm <- data.frame(lmer_parsed$X, check.names = F)
    org_col_names <- colnames(mm)
    # Create new formula based on model matrix
    formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                          paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
  }
  return(formula(formula_new))
}

# For gomps
generate_new_formula_gomps <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps, gomp) {
  if (is.null(random_effects_formula)) {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + ")))
    mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
    org_col_names <- colnames(mm)
    colnames_to_insert <- org_col_names[org_col_names != "(Intercept)"]
    if (length(colnames_to_insert) > 0) {
      formula_new <- paste0("expr ~ ", paste0("`", paste0(colnames_to_insert, collapse = "` + `"), "`"))
    } else {
      formula_new <- paste0("expr ~ 1")
    }
  } else {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + ")))
    lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
    mm <- data.frame(lmer_parsed$X, check.names = F)
    org_col_names <- colnames(mm)
    # Create new formula based on model matrix
    formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                          paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
  }
  return(formula(formula_new))
}

# For omps
generate_new_formula_omps <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps, omp) {
  if (is.null(random_effects_formula)) {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + ")))
    mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
    org_col_names <- colnames(mm)
    org_col_names <- org_col_names[!org_col_names %in% omp]
    colnames_to_insert <- org_col_names[org_col_names != "(Intercept)"]
    if (length(colnames_to_insert) > 0) {
      formula_new <- paste0("expr ~ ", paste0("`", paste0(colnames_to_insert, collapse = "` + `"), "`"))
    } else {
      formula_new <- paste0("expr ~ 1")
    }
  } else {
    input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + ")))
    lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
    mm <- data.frame(lmer_parsed$X, check.names = F)
    org_col_names <- colnames(mm)
    org_col_names <- org_col_names[!org_col_names %in% omp]
    # Create new formula based on model matrix
    formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                          paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
  }
  return(formula(formula_new))
}

# Get joint significance for zeros and non-zeros
add_joint_signif <- function(fit_data_non_zero, fit_data_binary, analysis_method, correction) {
  # Subset to shared columns
  fit_data_binary_signif <- fit_data_binary$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_binary_signif) <- c("feature", "metadata", "value", "name", "logistic", "logistic_error")
  fit_data_non_zero_signif <- fit_data_non_zero$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_non_zero_signif) <- c("feature", "metadata", "value", "name", analysis_method, "LM_error")
  
  # Join and check linear and logistic pieces
  merged_signif <- full_join(unique(fit_data_binary_signif), unique(fit_data_non_zero_signif), 
                             by=c("feature", "metadata", "value", "name"))
  if (nrow(merged_signif) != nrow(unique(fit_data_binary_signif)) | nrow(merged_signif) != nrow(unique(fit_data_non_zero_signif))) {
    print(nrow(unique(merged_signif)))
    print(nrow(unique(fit_data_binary_signif)))
    print(nrow(unique(fit_data_non_zero_signif)))
    print(anti_join(unique(fit_data_binary_signif), unique(fit_data_non_zero_signif), by=c("feature", "metadata", "value", "name")))
    print(anti_join(unique(fit_data_non_zero_signif), unique(fit_data_binary_signif), by=c("feature", "metadata", "value", "name")))
    stop("Merged significance tables have different associations")
  }
  
  # Create a combined p-value
  merged_signif$pval_joint <- pbeta(pmin(merged_signif[,analysis_method], 
                                         merged_signif[,"logistic"]), 
                                    1, 2)
  
  # If NA or model errored, use the p-value of the non-NA
  merged_signif$pval_joint <- ifelse(is.na(merged_signif[,analysis_method]) | !is.na(merged_signif$LM_error), merged_signif[,"logistic"], merged_signif$pval_joint)
  merged_signif$pval_joint <- ifelse(is.na(merged_signif[,"logistic"]) | !is.na(merged_signif$logistic_error), merged_signif[,analysis_method], merged_signif$pval_joint)
  merged_signif$pval_joint <- ifelse((is.na(merged_signif[,"logistic"]) | !is.na(merged_signif$logistic_error)) & 
                                       (is.na(merged_signif[,analysis_method]) | !is.na(merged_signif$LM_error)), 
                                     NA, merged_signif$pval_joint)
  merged_signif$qval_joint <- as.numeric(p.adjust(merged_signif$pval_joint, method = correction))
  
  return(list(append_joint(fit_data_non_zero, merged_signif), append_joint(fit_data_binary, merged_signif)))
}

# Take logistic or LM component and add on the merged significance pieces
append_joint <- function(outputs, merged_signif) {
  merged_signif <- merged_signif[,c("feature", "metadata", "value", "name", "pval_joint", "qval_joint")]
  tmp_colnames <- colnames(outputs$results)
  tmp_colnames <- case_when(tmp_colnames == "pval" ~ "pval_single",
                            tmp_colnames == "qval" ~ "qval_single",
                            TRUE ~ tmp_colnames)
  colnames(outputs$results) <- tmp_colnames
  
  merged_signif <- merge(outputs$results, merged_signif, 
                         by=c("feature", "metadata", "value", "name"))
  
  merged_signif <- merged_signif[order(merged_signif$qval_joint),]
  
  return(merged_signif)
}

bind_and_reorder <- function(growing_resid_mat, outputs_tmp_residuals, params_and_data_and_formula_data) {
  tmp_out <- rbind(growing_resid_mat, outputs_tmp_residuals)
  tmp_out <- 
    tmp_out[order(as.numeric(mapvalues(rownames(tmp_out),
                                       colnames(params_and_data_and_formula_data),
                                       1:ncol(params_and_data_and_formula_data), warn_missing = F))), , drop=F]
  return(tmp_out)
}

##################
# Gomp functions #
##################

gomp_lm <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + "))
      lm_small <- lm(formula_less_one, dat_sub)
      beta_small <- coef(lm_small)
      X_small <- model.matrix(formula_less_one, dat_sub)
      nll_small <- nlog_likelihood_linear_regression(X_small, beta_small, dat_sub$expr, w=rep(1, nrow(X_small)))
      
      # Get the full model matrix and identify the gomp's columns
      X_big <- model.matrix(formula_new, dat_sub)
      gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
      gomp_columns <- gomp_columns[gomp_columns %in% colnames(X_big)]
      
      # Thermometer encode the factor levels
      X_big[,gomp_columns] <- fill_left_of_ones(X_big[,gomp_columns, drop=F])
      gomp_cols <- which(!colnames(X_big) %in% names(beta_small))
      
      nll_big <- get_nll_big(direction = 'inc', 
                             lik = 'linear', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = dat_sub, 
                             w=rep(1, nrow(X_big)))
      
      # Flipped order because negative log likelihoods
      pval_inc <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      
      # Same thing but for decreasing predictors
      nll_big <- get_nll_big(direction = 'dec', 
                             lik = 'linear', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = dat_sub, 
                             w=rep(1, nrow(X_big)))
      
      # Flipped order because negative log likelihoods
      pval_dec <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      pval[i,] <- c(pval_inc, pval_dec)
      pval[i,]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  return(outputs)
}

# Based on https://rpubs.com/bbolker/glmerconstr
gomp_lmer <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + "))
      
      for (direction in c('inc', 'dec')) {
        formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
        # Fit with  missing parameter first
        f0 <- glFormula(formula_less_one, dat_sub, family = gaussian, REML=FALSE)
        devfun <- do.call(mkGlmerDevfun, c(f0, list(family=gaussian())))
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, f0$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        beta_small_names <- colnames(rho$pp$X)
        theta <- rho$pp$theta
        opt_small <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf,nbeta)),
                               upper=c(rep(Inf,length(theta)),rep(Inf,nbeta)))
        beta_small <- opt_small$par
        fval_small <- opt_small$fval
        
        # Fit with constrained coefficients
        # Get model matrix
        lmer_parsed <- lme4::lFormula(formula_new, dat_sub)
        mm <- data.frame(lmer_parsed$X, check.names = F)
        # Copy other columns
        mm[,colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm)]] <- 
          dat_sub[, colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm)]]
        gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
        gomp_columns <- gomp_columns[gomp_columns %in% colnames(mm)]
        # For ensuring order
        mm[,gomp_columns] <- fill_left_of_ones(mm[,gomp_columns, drop=F])
        
        formula_new <- generate_new_formula_gomps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
        
        fmod <- glFormula(formula_new, mm, family = gaussian, REML=FALSE)
        devfun <- do.call(mkGlmerDevfun, c(fmod, list(family=gaussian())))
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, fmod$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        gomp_cols <- which(!colnames(rho$pp$X) %in% beta_small_names)
        theta <- rho$pp$theta
        
        if (direction == 'inc') {
          lb <- rep(-Inf, nbeta)
          lb[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,lb),
                               upper=c(rep(Inf,length(theta)), rep(Inf,nbeta)))
        } else {
          ub <- rep(Inf, nbeta)
          ub[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta,rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf, nbeta)),
                               upper=c(rep(Inf,length(theta)), ub))
        }
        
        compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                 opt_big$par[-c(1:length(theta))])
        
        # Test statistic
        fval <- opt_small$fval - opt_big$fval
        
        # Correct convergence issues
        fval <- correct_fval(fval, compare_results, gomp_cols)
        
        # Flipped order because negative log likelihoods
        pval_cur <- pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F)
        
        # For ~low p-values, check that the low p-value isn't because of a fitting issue
        if (!is.na(pval_cur) && pval_cur < 0.005) {
          refit <- lmer(formula_less_one, dat_sub, REML = F)
          opt_small$fval <- -2 * c(logLik(refit))
          opt_small$par[-c(1:length(theta))] <- fixef(refit)
          compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                   opt_big$par[-c(1:length(theta))])
          
          fval <- opt_small$fval - opt_big$fval
          fval <- correct_fval(fval, compare_results, gomp_cols)
          pval_cur <- max(pval_cur, pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F))
        }
        
        if (direction == 'inc') {
          pval[i, 1] <- pval_cur
        } else {
          pval[i, 2] <- pval_cur
        }
      }
      pval[i,]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  return(outputs)
}

gomp_glm <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + "))
      lm_small <- glm(formula_less_one, dat_sub, family = 'binomial')
      beta_small <- coef(lm_small)
      X_small <- model.matrix(formula_less_one, dat_sub)
      nll_small <- nlog_likelihood_logistic_regression(X_small, beta_small, dat_sub$expr, w=rep(1, nrow(X_small)))
      
      # Get the full model matrix and identify the gomp's columns
      X_big <- model.matrix(formula_new, dat_sub)
      gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
      gomp_columns <- gomp_columns[gomp_columns %in% colnames(X_big)]
      
      # Thermometer encode the factor levels
      X_big[,gomp_columns] <- fill_left_of_ones(X_big[,gomp_columns, drop=F])
      gomp_cols <- which(!colnames(X_big) %in% names(beta_small))
      
      nll_big <- get_nll_big(direction = 'inc', 
                             lik = 'logistic', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = dat_sub, 
                             w=rep(1, nrow(X_big)))
      
      # Flipped order because negative log likelihoods
      pval_inc <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      
      # Same thing but for decreasing predictors
      nll_big <- get_nll_big(direction = 'dec', 
                             lik = 'logistic', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = dat_sub, 
                             w=rep(1, nrow(X_big)))
      
      # Flipped order because negative log likelihoods
      pval_dec <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      pval[i,] <- c(pval_inc, pval_dec)
      pval[i, ]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  return(outputs)
}

gomp_glm_augment <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
  mm_input <- augmented_data[["mm_input"]]
  weight_scheme <- augmented_data[["weight_scheme"]]
  formula_new <- augmented_data[["new_formula"]]
  weight_scheme <- weight_scheme[rowSums(is.na(mm_input)) == 0]
  mm_input <- mm_input[rowSums(is.na(mm_input)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- generate_new_formula_gomps(formula, random_effects_formula, dat_sub, groups, gomps, omps, gomp)
      
      lm_small <- glm(formula_less_one, mm_input, family = 'binomial')
      beta_small <- coef(lm_small)
      X_small <- model.matrix(formula_less_one, mm_input)
      nll_small <- nlog_likelihood_logistic_regression(X_small, beta_small, mm_input$expr, w=weight_scheme)
      
      # Get the full model matrix and identify the gomp's columns
      X_big <- model.matrix(formula_new, mm_input)
      gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
      gomp_columns <- gomp_columns[gomp_columns %in% colnames(X_big)]
      
      # Thermometer encode the factor levels
      X_big[,gomp_columns] <- fill_left_of_ones(X_big[,gomp_columns, drop=F])
      gomp_cols <- which(!colnames(X_big) %in% names(beta_small))
      
      nll_big <- get_nll_big(direction = 'inc', 
                             lik = 'logistic', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = mm_input, 
                             w=weight_scheme)
      
      # Flipped order because negative log likelihoods
      pval_inc <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      
      # Same thing but for decreasing predictors
      nll_big <- get_nll_big(direction = 'dec', 
                             lik = 'logistic', 
                             X_big = X_big, 
                             gomp_cols = gomp_cols, 
                             dat_sub = mm_input, 
                             w=weight_scheme)
      
      # Flipped order because negative log likelihoods
      pval_dec <- pchisqmix(2 * (nll_small - nll_big), s=0, q=length(gomp_cols), lower.tail = F)
      pval[i,] <- c(pval_inc, pval_dec)
      pval[i,]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  return(outputs)
}

gomp_glmer <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- formula(paste0(c(safe_deparse(formula), groups, gomps[gomps != gomp], omps), collapse = " + "))
      
      for (direction in c('inc', 'dec')) {
        formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
        # Fit with  missing parameter first
        f0 <- glFormula(formula_less_one, dat_sub, family = binomial)
        devfun <- do.call(mkGlmerDevfun, f0)
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, f0$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        beta_small_names <- colnames(rho$pp$X)
        theta <- rho$pp$theta
        opt_small <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf,nbeta)),
                               upper=c(rep(Inf,length(theta)),rep(Inf,nbeta)))
        beta_small <- opt_small$par
        fval_small <- opt_small$fval
        
        # Fit with constrained coefficients
        # Get model matrix
        lmer_parsed <- lme4::lFormula(formula_new, dat_sub)
        mm <- data.frame(lmer_parsed$X, check.names = F)
        # Copy other columns
        mm[,colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm)]] <- 
          dat_sub[, colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm)]]
        gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
        gomp_columns <- gomp_columns[gomp_columns %in% colnames(mm)]
        # For ensuring order
        mm[,gomp_columns] <- fill_left_of_ones(mm[,gomp_columns, drop=F])
        
        formula_new <- generate_new_formula_gomps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
        
        fmod <- glFormula(formula_new, mm, family = binomial)
        devfun <- do.call(mkGlmerDevfun, fmod)
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, fmod$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        gomp_cols <- which(!colnames(rho$pp$X) %in% beta_small_names)
        theta <- rho$pp$theta
        
        if (direction == 'inc') {
          lb <- rep(-Inf, nbeta)
          lb[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,lb),
                               upper=c(rep(Inf,length(theta)), rep(Inf,nbeta)))
        } else {
          ub <- rep(Inf, nbeta)
          ub[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta,rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf, nbeta)),
                               upper=c(rep(Inf,length(theta)), ub))
        }
        
        compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                 opt_big$par[-c(1:length(theta))])
        
        # Test statistic
        fval <- opt_small$fval - opt_big$fval
        
        # Correct convergence issues
        fval <- correct_fval(fval, compare_results, gomp_cols)
        
        # Flipped order because negative log likelihoods
        pval_cur <- pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F)
        
        # For ~low p-values, check that the low p-value isn't because of a fitting issue
        if (!is.na(pval_cur) && pval_cur < 0.005) {
          refit <- glmer(formula_less_one, dat_sub, family = 'binomial')
          opt_small$fval <- -2 * c(logLik(refit))
          opt_small$par[-c(1:length(theta))] <- fixef(refit)
          compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                   opt_big$par[-c(1:length(theta))])
          
          fval <- opt_small$fval - opt_big$fval
          fval <- correct_fval(fval, compare_results, gomp_cols)
          pval_cur <- max(pval_cur, pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F))
        }
        
        if (direction == 'inc') {
          pval[i, 1] <- pval_cur
        } else {
          pval[i, 2] <- pval_cur
        }
      }
      pval[i,]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  return(outputs)
}

gomp_glmer_augment <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
  mm_input <- augmented_data[["mm_input"]]
  weight_scheme <<- augmented_data[["weight_scheme"]] # Needs to be global because of GLMER environment?!
  formula_new <- augmented_data[["new_formula"]]
  weight_scheme <<- weight_scheme[rowSums(is.na(mm_input)) == 0]
  mm_input <- mm_input[rowSums(is.na(mm_input)) == 0,]
  
  pval <- matrix(nrow = length(gomps), ncol = 2)
  for (i in seq_along(gomps)) {
    pval[i, ] <- tryCatch({
      gomp = gomps[i]
      # Fit on the data missing the current gomp
      formula_less_one <- generate_new_formula_gomps(formula, random_effects_formula, dat_sub, groups, gomps, omps, gomp)
      
      for (direction in c('inc', 'dec')) {
        formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
        
        # Fit with  missing parameter first
        f0 <- glFormula(formula_less_one, mm_input, family = binomial, weights = weight_scheme)
        devfun <- do.call(mkGlmerDevfun, f0)
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, f0$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        beta_small_names <- colnames(rho$pp$X)
        theta <- rho$pp$theta
        opt_small <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf,nbeta)),
                               upper=c(rep(Inf,length(theta)),rep(Inf,nbeta)))
        beta_small <- opt_small$par
        fval_small <- opt_small$fval
        
        # Fit with constrained coefficients
        # Get model matrix
        gomp_columns <- paste0(gomp, levels(dat_sub[[gomp]]))
        gomp_columns <- gomp_columns[gomp_columns %in% colnames(mm_input)]
        # For ensuring order
        mm_input_tmp <- mm_input
        mm_input_tmp[,gomp_columns] <- fill_left_of_ones(mm_input_tmp[,gomp_columns, drop=F])
        
        formula_new <- generate_new_formula_gomps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
        
        fmod <- glFormula(formula_new, mm_input_tmp, family = binomial, weights = weight_scheme)
        devfun <- do.call(mkGlmerDevfun, fmod)
        opt1 <- optimizeGlmer(devfun)
        theta.lwr <- environment(devfun)$lower
        devfun <- updateGlmerDevfun(devfun, fmod$reTrms)
        rho <- environment(devfun)
        nbeta <- ncol(rho$pp$X)
        gomp_cols <- which(!colnames(rho$pp$X) %in% beta_small_names)
        theta <- rho$pp$theta
        
        if (direction == 'inc') {
          lb <- rep(-Inf, nbeta)
          lb[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,lb),
                               upper=c(rep(Inf,length(theta)), rep(Inf,nbeta)))
        } else {
          ub <- rep(Inf, nbeta)
          ub[gomp_cols] <- 0
          opt_big <- nloptwrap(par=c(rho$pp$theta,rep(0,nbeta)),
                               fn=devfun,
                               lower=c(theta.lwr,rep(-Inf, nbeta)),
                               upper=c(rep(Inf,length(theta)), ub))
        }
        
        compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                 opt_big$par[-c(1:length(theta))])
        
        # Test statistic
        fval <- opt_small$fval - opt_big$fval
        
        # Correct convergence issues
        fval <- correct_fval(fval, compare_results, gomp_cols)
        
        # Flipped order because negative log likelihoods
        pval_cur <- pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F)
        
        # For ~low p-values, check that the low p-value isn't because of a fitting issue
        if (!is.na(pval_cur) && pval_cur < 0.005) {
          refit <- glmer(formula_less_one, mm_input, family = 'binomial', weights = weight_scheme)
          opt_small$fval <- -2 * c(logLik(refit))
          opt_small$par[-c(1:length(theta))] <- fixef(refit)
          compare_results <- rbind(insert_zeros_at_indices(opt_small$par[-c(1:length(theta))], gomp_cols),
                                   opt_big$par[-c(1:length(theta))])
          
          fval <- opt_small$fval - opt_big$fval
          fval <- correct_fval(fval, compare_results, gomp_cols)
          pval_cur <- max(pval_cur, pchisqmix(fval, s=0, q=length(gomp_cols), lower.tail = F))
        }
        
        if (direction == 'inc') {
          pval[i, 1] <- pval_cur
        } else {
          pval[i, 2] <- pval_cur
        }
      }
      pval[i, ]
    }, error = function(err) { 
      pval[i, ] <- rep(NA, 2)
      return(rep(NA, 2))
    })
  }
  
  pval_final <- apply(pval, 1, function(x) {min(p.adjust(x, method = 'BY'))})
  directions <- apply(pval, 1, function(x) {ifelse(x[1] <= x[2], 1, -1)})
  
  outputs <- data.frame(matrix(ncol = 4, nrow = length(gomps)))
  outputs[,1] <- directions
  outputs[,3] <- pval_final
  outputs[,4] <- gomps
  
  remove(weight_scheme, pos = ".GlobalEnv") # Stop from being global
  
  return(outputs)
}

nlog_likelihood_linear_regression <- function(X, beta, y, w) {
  residuals <- y - X %*% beta
  sigma_sq <- c(var(residuals))
  log_likelihood <- -0.5 * (sum(log(2 * pi * sigma_sq / w)) + sum(residuals^2 / (sigma_sq / w) ))
  return(-log_likelihood)
}

nlog_likelihood_logistic_regression <- function(X, beta, y, w) {
  lin_pred <- X %*% beta
  probabilities <- plogis(lin_pred)
  log_likelihood <- sum(w * ifelse(y == 1, log(probabilities), log(1 - probabilities)))
  return(-log_likelihood)
}

get_nll_big <- function(direction, lik, X_big, gomp_cols, dat_sub, w) {
  # Build constraint matrix
  constr_mat <- matrix(0, nrow = length(gomp_cols), ncol = ncol(X_big))
  coordinates <- cbind(1:length(gomp_cols), gomp_cols)
  for (i in 1:nrow(coordinates)) {
    coord <- coordinates[i,]
    constr_mat[coord[1], coord[2]] <- ifelse(direction == 'inc', 1, -1)
  }
  
  if (direction == 'inc') {
    initial_beta <- rep(1, ncol(X_big))
  } else {
    initial_beta <- rep(-1, ncol(X_big))
  }
  
  if (lik == 'linear') {
    likfun <- nlog_likelihood_linear_regression
  } else if (lik == 'logistic') {
    likfun <- nlog_likelihood_logistic_regression
  } else {
    stop('Invalid likelihood type to get_nll_big')
  }
  
  constr_out <- constrOptim(theta = initial_beta, 
                            f = likfun,
                            ui = constr_mat,
                            ci = rep(0, length(gomp_cols)), 
                            X = X_big, y = dat_sub$expr, w=w, 
                            method = "Nelder-Mead", outer.eps = 10^-5)
  
  beta_big <- constr_out$par
  beta_big <- ifelse(abs(beta_big) < 10^-3, 0, beta_big)
  nll_big <- likfun(X_big, beta_big, dat_sub$expr, w=w)
  return(nll_big)
}

get_nll_with_coef <- function(direction, lik, X_big, gomp_cols, dat_sub, w) {
  # Build constraint matrix
  constr_mat <- matrix(0, nrow = length(gomp_cols), ncol = ncol(X_big))
  coordinates <- cbind(1:length(gomp_cols), gomp_cols)
  for (i in 1:nrow(coordinates)) {
    coord <- coordinates[i,]
    constr_mat[coord[1], coord[2]] <- ifelse(direction == 'inc', 1, -1)
  }
  
  if (direction == 'inc') {
    initial_beta <- rep(1, ncol(X_big))
  } else {
    initial_beta <- rep(-1, ncol(X_big))
  }
  
  if (lik == 'linear') {
    likfun <- nlog_likelihood_linear_regression
  } else if (lik == 'logistic') {
    likfun <- nlog_likelihood_logistic_regression
  } else {
    stop('Invalid likelihood type to get_nll_big')
  }
  
  constr_out <- constrOptim(theta = initial_beta, 
                            f = likfun,
                            ui = constr_mat,
                            ci = rep(0, length(gomp_cols)), 
                            X = X_big, y = dat_sub$expr, w=w, 
                            method = "Nelder-Mead", outer.eps = 10^-5)
  
  beta_big <- constr_out$par
  beta_big <- ifelse(abs(beta_big) < 10^-3, 0, beta_big)
  nll_big <- likfun(X_big, beta_big, dat_sub$expr, w=w)
  return(list('nll' = nll_big, 'beta' = beta_big))
}

get_nll_with_coef_re <- function(formula_new, direction, lik, mm_input, gomp_cols, w=NULL) {
  if (!is.null(w)) {
    weight_sch <<- w # Apparently weights have to be a global variable
  }
  if (lik == 'linear') {
    fmod <- glFormula(formula_new, mm_input, family = gaussian, REML=FALSE)
    devfun <- do.call(mkGlmerDevfun, c(fmod, list(family=gaussian())))
  } else if (lik == 'logistic') {
    fmod <- glFormula(formula_new, mm_input, family = binomial, weights = weight_sch)
    devfun <- do.call(mkGlmerDevfun, fmod)
  } else {
    stop("Invalid lik")
  }
  opt1 <- optimizeGlmer(devfun)
  theta.lwr <- environment(devfun)$lower
  devfun <- updateGlmerDevfun(devfun, fmod$reTrms)
  rho <- environment(devfun)
  nbeta <- ncol(rho$pp$X)
  theta <- rho$pp$theta
  
  if (direction == 'inc') {
    lb <- rep(-Inf, nbeta)
    lb[gomp_cols] <- 0
    opt_big <- nloptwrap(par=c(rho$pp$theta, rep(0,nbeta)),
                         fn=devfun,
                         lower=c(theta.lwr,lb),
                         upper=c(rep(Inf,length(theta)), rep(Inf,nbeta)))
  } else {
    ub <- rep(Inf, nbeta)
    ub[gomp_cols] <- 0
    opt_big <- nloptwrap(par=c(rho$pp$theta,rep(0,nbeta)),
                         fn=devfun,
                         lower=c(theta.lwr,rep(-Inf, nbeta)),
                         upper=c(rep(Inf,length(theta)), ub))
  }
  
  if (!is.null(w)) {
    remove(weight_sch, pos = ".GlobalEnv") # Stop from being global
  }
  
  return(list('nll' = opt_big$fval, 'beta' = opt_big$par[-c(1:length(theta))]))
}

fill_left_of_ones <- function(matrix_input) {
  for (i in 1:nrow(matrix_input)) {
    # Not an augmented row
    if (max(matrix_input[i, ]) == 1) {
      one_index <- which(matrix_input[i, ] == 1)
      if (length(one_index) > 0) {
        matrix_input[i, 1:max((one_index - 1), 1)] <- 1
      }
    } else { # Augmented row
      matrix_input[i, ] = cumsum(matrix_input[i, ])
    }
  }
  
  return(matrix_input)
}

insert_zeros_at_indices <- function(values, indices) {
  result <- values
  
  for (i in seq_along(indices)) {
    index <- indices[i]
    if (index == 1) {
      result <- c(0, result)
    } else if (index >= length(result)) {
      result <- c(result, 0)
    } else {
      result <- c(result[1:(index - 1)], 0, result[index:length(result)])
    }
  }
  
  return(result)
}

correct_fval <- function(fval, compare_results, gomp_cols) {
  fval <- ifelse(abs(fval) < 0.0001, 0, fval)
  if (fval < 0) {
    if (all(compare_results[,gomp_cols] == 0)) { # If all constrained columns are 0, they're only different by convergence issues
      fval <- 0
    } else {
      fval <- NA
    }
  }
  if (all(compare_results[,gomp_cols] == 0)) {
    fval <- 0
  }
  return(fval)
}

#################
# Omp functions #
#################

omp_lm <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    omp = omps[i]
    
    # Get the full model matrix and identify the omp's columns
    X_big <- model.matrix(formula_new, dat_sub)
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(X_big)]
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      X_big[,omp_columns] <- fill_left_of_ones(X_big[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(X_big))
      
      nll_with_coef_inc <- get_nll_with_coef(direction = 'inc', 
                                             lik = 'linear', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = dat_sub, 
                                             w=rep(1, nrow(X_big)))
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef(direction = 'dec', 
                                             lik = 'linear', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = dat_sub, 
                                             w=rep(1, nrow(X_big)))
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            X_small <- X_big[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(X_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef(direction = better_direction, 
                                                     lik = 'linear', 
                                                     X_big = X_small, 
                                                     gomp_cols = omp_cols_remaining, 
                                                     dat_sub = dat_sub, 
                                                     w=rep(1, nrow(X_small)))
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(2 * (nll_with_coef_small$nll - nll_with_coef_big$nll), s=0, q=1, lower.tail = F)
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  return(results_df)
}

omp_lmer <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  # Get the full model matrix
  lmer_parsed <- lme4::lFormula(formula_new, dat_sub)
  mm_input <- data.frame(lmer_parsed$X, check.names = F)
  # Copy other columns
  mm_input[,colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm_input)]] <- 
    dat_sub[, colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm_input)]]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    mm <- mm_input
    omp = omps[i]
    
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(mm)]
    
    formula_new <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      mm[,omp_columns] <- fill_left_of_ones(mm[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(mm))
      
      nll_with_coef_inc <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'inc', 
                                                lik = 'linear', 
                                                mm = mm, 
                                                gomp_cols = omp_cols)
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'dec', 
                                                lik = 'linear', 
                                                mm = mm, 
                                                gomp_cols = omp_cols)
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        formula_less_one <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, omp_column)
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            mm_small <- mm[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(mm_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef_re(formula_new = formula_less_one,
                                                        direction = better_direction, 
                                                        lik = 'linear', 
                                                        mm = mm_small, 
                                                        gomp_cols = omp_cols_remaining)
            
            compare_results <- rbind(insert_zeros_at_indices(nll_with_coef_small$beta, omp_col),
                                     nll_with_coef_big$beta)
            
            # Test statistic
            fval <- nll_with_coef_small$nll - nll_with_coef_big$nll
            
            # Checking convergence issues v1
            fval_corrected <- FALSE
            if (fval < -0.0001 | fval > 3) {
              big_tmp <- lmer(formula_new, mm, REML = F)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(big_tmp)[omp_cols] >= 0)
              } else {
                check_condition <- all(fixef(big_tmp)[omp_cols] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(big_tmp)) < nll_with_coef_big$nll) {
                nll_with_coef_big_tmp <- list('nll' = -2 * c(logLik(big_tmp)), 'beta' = fixef(big_tmp))
              } else {
                nll_with_coef_big_tmp <- nll_with_coef_big
              }
              
              small_tmp <- lmer(formula_less_one, mm_small, REML = F)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] >= 0)
              } else {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(small_tmp)) < nll_with_coef_small$nll) {
                nll_with_coef_small_tmp <- list('nll' = -2 * c(logLik(small_tmp)), 'beta' = fixef(small_tmp))
              } else {
                nll_with_coef_small_tmp <- nll_with_coef_small
              }
              
              fval_new <- nll_with_coef_small_tmp$nll - nll_with_coef_big_tmp$nll
              
              fval_corrected <- fval_new > 0.0001 & fval_new < 3
              if (fval_corrected) {
                fval <- fval_new
              }
            }
            
            # Correct convergence issues v3
            fval <- correct_fval(fval, compare_results, omp_col)
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(fval, s=0, q=1, lower.tail = F)
            
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  return(results_df)
}

omp_glm <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    omp = omps[i]
    
    # Get the full model matrix and identify the omp's columns
    X_big <- model.matrix(formula_new, dat_sub)
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(X_big)]
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      X_big[,omp_columns] <- fill_left_of_ones(X_big[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(X_big))
      
      nll_with_coef_inc <- get_nll_with_coef(direction = 'inc', 
                                             lik = 'logistic', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = dat_sub, 
                                             w=rep(1, nrow(X_big)))
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef(direction = 'dec', 
                                             lik = 'logistic', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = dat_sub, 
                                             w=rep(1, nrow(X_big)))
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            X_small <- X_big[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(X_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef(direction = better_direction, 
                                                     lik = 'logistic', 
                                                     X_big = X_small, 
                                                     gomp_cols = omp_cols_remaining, 
                                                     dat_sub = dat_sub, 
                                                     w=rep(1, nrow(X_small)))
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(2 * (nll_with_coef_small$nll - nll_with_coef_big$nll), s=0, q=1, lower.tail = F)
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  return(results_df)
}

omp_glm_augment <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
  mm_input <- augmented_data[["mm_input"]]
  weight_scheme <- augmented_data[["weight_scheme"]]
  formula_new <- augmented_data[["new_formula"]]
  weight_scheme <- weight_scheme[rowSums(is.na(mm_input)) == 0]
  mm_input <- mm_input[rowSums(is.na(mm_input)) == 0,]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    omp = omps[i]
    
    # Get the full model matrix and identify the omp's columns
    X_big <- model.matrix(formula_new, mm_input)
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(X_big)]
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      X_big[,omp_columns] <- fill_left_of_ones(X_big[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(X_big))
      
      nll_with_coef_inc <- get_nll_with_coef(direction = 'inc', 
                                             lik = 'logistic', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = mm_input, 
                                             w=weight_scheme)
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef(direction = 'dec', 
                                             lik = 'logistic', 
                                             X_big = X_big, 
                                             gomp_cols = omp_cols, 
                                             dat_sub = mm_input, 
                                             w=weight_scheme)
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            X_small <- X_big[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(X_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef(direction = better_direction, 
                                                     lik = 'logistic', 
                                                     X_big = X_small, 
                                                     gomp_cols = omp_cols_remaining, 
                                                     dat_sub = mm_input, 
                                                     w=weight_scheme)
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(2 * (nll_with_coef_small$nll - nll_with_coef_big$nll), s=0, q=1, lower.tail = F)
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  return(results_df)
}

omp_glmer <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  dat_sub <- dat_sub[rowSums(is.na(dat_sub)) == 0,]
  # Get the full model matrix
  lmer_parsed <- lme4::lFormula(formula_new, dat_sub)
  mm_input <- data.frame(lmer_parsed$X, check.names = F)
  # Copy other columns
  mm_input[,colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm_input)]] <- 
    dat_sub[, colnames(dat_sub)[!colnames(dat_sub) %in% colnames(mm_input)]]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    mm <- mm_input
    omp = omps[i]
    
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(mm)]
    
    formula_new <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      mm[,omp_columns] <- fill_left_of_ones(mm[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(mm))
      
      nll_with_coef_inc <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'inc', 
                                                lik = 'logistic', 
                                                mm = mm, 
                                                gomp_cols = omp_cols,
                                                w = rep(1, nrow(mm)))
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'dec', 
                                                lik = 'logistic', 
                                                mm = mm, 
                                                gomp_cols = omp_cols,
                                                w = rep(1, nrow(mm)))
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        formula_less_one <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, omp_column)
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            mm_small <- mm[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(mm_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef_re(formula_new = formula_less_one,
                                                        direction = better_direction, 
                                                        lik = 'logistic', 
                                                        mm = mm_small, 
                                                        gomp_cols = omp_cols_remaining,
                                                        w = rep(1, nrow(mm)))
            
            compare_results <- rbind(insert_zeros_at_indices(nll_with_coef_small$beta, omp_col),
                                     nll_with_coef_big$beta)
            
            # Test statistic
            fval <- nll_with_coef_small$nll - nll_with_coef_big$nll
            
            # Checking convergence issues v1
            fval_corrected <- FALSE
            if (fval < -0.0001 | fval > 3) {
              big_tmp <- glmer(formula_new, mm, family = binomial)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(big_tmp)[omp_cols] >= 0)
              } else {
                check_condition <- all(fixef(big_tmp)[omp_cols] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(big_tmp)) < nll_with_coef_big$nll) {
                nll_with_coef_big_tmp <- list('nll' = -2 * c(logLik(big_tmp)), 'beta' = fixef(big_tmp))
              } else {
                nll_with_coef_big_tmp <- nll_with_coef_big
              }
              
              small_tmp <- glmer(formula_less_one, mm_small, family = binomial)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] >= 0)
              } else {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(small_tmp)) < nll_with_coef_small$nll) {
                nll_with_coef_small_tmp <- list('nll' = -2 * c(logLik(small_tmp)), 'beta' = fixef(small_tmp))
              } else {
                nll_with_coef_small_tmp <- nll_with_coef_small
              }
              
              fval_new <- nll_with_coef_small_tmp$nll - nll_with_coef_big_tmp$nll
              
              fval_corrected <- fval_new > 0.0001 & fval_new < 3
              if (fval_corrected) {
                fval <- fval_new
              }
            }
            
            # Checking convergence issues v2
            if (fval > 0.1) {
              tmp_formula_new <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps[omps != omp], "")
              
              # Model with constrained coefficients can't be worse than the model without the coefficients at all
              if (-2 * c(logLik(glmer(tmp_formula_new, mm, family = binomial))) < nll_with_coef_small$nll) {
                fval <- NA
              }
            }
            
            # Correct convergence issues v3
            fval <- correct_fval(fval, compare_results, omp_col)
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(fval, s=0, q=1, lower.tail = F)
            
            list('coef' = ifelse(fval_corrected, nll_with_coef_big_tmp$beta[omp_col], nll_with_coef_big$beta[omp_col]), 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  return(results_df)
}

omp_glmer_augment <- function(formula, random_effects_formula, dat_sub, groups, gomps, omps) {
  formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
  
  augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
  mm_input <- augmented_data[["mm_input"]]
  weight_scheme <<- augmented_data[["weight_scheme"]] # Needs to be global because of GLMER environment?!
  formula_new <- augmented_data[["new_formula"]]
  weight_scheme <<- weight_scheme[rowSums(is.na(mm_input)) == 0]
  mm_input <- mm_input[rowSums(is.na(mm_input)) == 0,]
  mm_input$`(Intercept)` <- 1
  mm_input <- mm_input[,c("(Intercept)", colnames(mm_input)[colnames(mm_input) != '(Intercept)'])]
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in seq_along(omps)) {
    mm <- mm_input
    omp = omps[i]
    
    omp_columns <- paste0(omp, levels(dat_sub[[omp]]))
    omp_columns <- omp_columns[omp_columns %in% colnames(mm)]
    
    formula_new <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, "")
    
    new_results <- tryCatch({
      # Thermometer encode the factor levels
      mm[,omp_columns] <- fill_left_of_ones(mm[,omp_columns, drop=F])
      omp_cols <- match(omp_columns, colnames(mm))
      
      nll_with_coef_inc <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'inc', 
                                                lik = 'logistic', 
                                                mm = mm, 
                                                gomp_cols = omp_cols,
                                                w = weight_scheme)
      
      # Same thing but for decreasing predictors
      nll_with_coef_dec <- get_nll_with_coef_re(formula_new = formula_new,
                                                direction = 'dec', 
                                                lik = 'logistic', 
                                                mm = mm, 
                                                gomp_cols = omp_cols,
                                                w = weight_scheme)
      
      better_direction <- ifelse(nll_with_coef_inc$nll < nll_with_coef_dec$nll, 'inc', 'dec')
      if (better_direction == 'inc') {
        nll_with_coef_big <- nll_with_coef_inc
      } else {
        nll_with_coef_big <- nll_with_coef_dec
      }
      
      new_results <- data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq_along(omp_columns)) {
        omp_column = omp_columns[j]
        omp_col = omp_cols[j]
        formula_less_one <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps, omp_column)
        
        new_result <- tryCatch({
          # If coefficient is 0, p-value must be 1
          if (nll_with_coef_big$beta[omp_col] == 0) {
            list('coef' = nll_with_coef_big$beta[omp_col], 
                 'stderr' = NA,
                 'pval' = 1,
                 'name' = omp_column)
          } else {
            mm_small <- mm[,-omp_col]
            omp_cols_remaining <- match(omp_columns, colnames(mm_small))
            omp_cols_remaining <- omp_cols_remaining[!is.na(omp_cols_remaining)]
            
            nll_with_coef_small <- get_nll_with_coef_re(formula_new = formula_less_one,
                                                        direction = better_direction, 
                                                        lik = 'logistic', 
                                                        mm = mm_small, 
                                                        gomp_cols = omp_cols_remaining,
                                                        w = weight_scheme)
            
            compare_results <- rbind(insert_zeros_at_indices(nll_with_coef_small$beta, omp_col),
                                     nll_with_coef_big$beta)
            
            # Test statistic
            fval <- nll_with_coef_small$nll - nll_with_coef_big$nll
            
            # Checking convergence issues v1
            fval_corrected <- FALSE
            if (fval < -0.0001 | fval > 3) {
              big_tmp <- glmer(formula_new, mm, family = binomial, weights = weight_scheme)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(big_tmp)[omp_cols] >= 0)
              } else {
                check_condition <- all(fixef(big_tmp)[omp_cols] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(big_tmp)) < nll_with_coef_big$nll) {
                nll_with_coef_big_tmp <- list('nll' = -2 * c(logLik(big_tmp)), 'beta' = fixef(big_tmp))
              } else {
                nll_with_coef_big_tmp <- nll_with_coef_big
              }
              
              small_tmp <- glmer(formula_less_one, mm_small, family = binomial, weights = weight_scheme)
              if (better_direction == 'inc') {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] >= 0)
              } else {
                check_condition <- all(fixef(small_tmp)[omp_cols_remaining] <= 0)
              }
              
              if (check_condition & -2 * c(logLik(small_tmp)) < nll_with_coef_small$nll) {
                nll_with_coef_small_tmp <- list('nll' = -2 * c(logLik(small_tmp)), 'beta' = fixef(small_tmp))
              } else {
                nll_with_coef_small_tmp <- nll_with_coef_small
              }
              
              fval_new <- nll_with_coef_small_tmp$nll - nll_with_coef_big_tmp$nll
              
              fval_corrected <- fval_new > 0.0001 & fval_new < 3
              if (fval_corrected) {
                fval <- fval_new
              }
            }
            
            # Checking convergence issues v2
            if (fval > 0.1) {
              tmp_formula_new <- generate_new_formula_omps(formula, random_effects_formula, dat_sub, groups, gomps, omps[omps != omp], "")
              
              # Model with constrained coefficients can't be worse than the model without the coefficients at all
              if (-2 * c(logLik(glmer(tmp_formula_new, mm, family = binomial, weights = weight_scheme))) < nll_with_coef_small$nll) {
                fval <- NA
              }
            }
            
            # Correct convergence issues v3
            fval <- correct_fval(fval, compare_results, omp_col)
            
            # Flipped order because negative log likelihoods
            pval_current_omp <- pchisqmix(fval, s=0, q=1, lower.tail = F)
            
            list('coef' = ifelse(fval_corrected, nll_with_coef_big_tmp$beta[omp_col], nll_with_coef_big$beta[omp_col]), 
                 'stderr' = NA,
                 'pval' = pval_current_omp,
                 'name' = omp_column)
          }}, error = function(err) {
            return(list('coef' = NA, 
                        'stderr' = NA,
                        'pval' = NA,
                        'name' = omp_column))
          })
        colnames(new_results) <- c('coef', 'stderr' , 'pval', 'name')
        new_results <- rbind(new_results, new_result)
      }
      new_results
    }, error = function(err) {
      n_levels <- length(omp_columns)
      new_results <- data.frame('coef' = rep(NA, n_levels),
                                'stderr' = rep(NA, n_levels),
                                'pval' = rep(NA, n_levels),
                                'name' = omp_columns)
      return(new_results)
    })
    colnames(results_df) <- c('coef', 'stderr' , 'pval', 'name')
    results_df <- rbind(results_df, new_results)
  }
  
  remove(weight_scheme, pos = ".GlobalEnv") # Stop from being global

  return(results_df)
}

# fit the data using the model selected and applying the correction
fit.model <- function(
    features,
    metadata,
    model,
    formula = NULL,
    random_effects_formula = NULL,
    correction = "BH",
    save_models = FALSE,
    augment = FALSE,
    cores = 1) {
  function_vec <- c("augment_data", "safe_deparse", "extract_special_predictor", "get_fixed_effects", "get_character_cols", "generate_new_formula_groups", 
                    "generate_new_formula_gomps", "generate_new_formula_omps", "add_joint_signif", "append_joint", "bind_and_reorder",
                    "gomp_lm", "gomp_lmer", "gomp_glm", "gomp_glm_augment", "gomp_glmer", "gomp_glmer_augment", "nlog_likelihood_linear_regression",
                    "nlog_likelihood_logistic_regression", "get_nll_big", "get_nll_with_coef", "get_nll_with_coef_re", "fill_left_of_ones",
                    "insert_zeros_at_indices", "correct_fval", "omp_lm", "omp_lmer", "omp_glm", "omp_glm_augment", "omp_glmer", "omp_glmer_augment")
      
  # Check formulas are valid
  if (is.null(random_effects_formula)) {
    if (is.null(formula)) {
      logging::logerror(
        paste("Both formula and random_effects_formula are null")
      )
      stop()
    }
  }
  
  formula <- formula(formula)
  
  # Extract group components
  extract_special_predictor_out <- extract_special_predictor(formula, 'group')
  formula <- extract_special_predictor_out[[1]]
  groups <- extract_special_predictor_out[[2]]
  
  # Extract gomp components
  extract_special_predictor_out <- extract_special_predictor(formula, 'gomp')
  formula <- extract_special_predictor_out[[1]]
  gomps <- extract_special_predictor_out[[2]]
  
  # Extract omp components
  extract_special_predictor_out <- extract_special_predictor(formula, 'omp')
  formula <- extract_special_predictor_out[[1]]
  omps <- extract_special_predictor_out[[2]]
  
  #############################################################
  # Determine the function and summary for the model selected #
  #############################################################
  
  ################
  # Linear Model #
  ################
  
  if (model == "LM") {
      if (is.null(random_effects_formula)) { # Fixed effects only
          model_function <-
              function(formula, data, na.action) {
                  return(lm(formula(formula), 
                             data = data,
                             na.action = na.action))
              }
          summary_function <- function(fit, names_to_include) {
              lm_summary <- summary(fit)$coefficients
              if (nrow(lm_summary) < length(names_to_include)) { # If deficient rank, make sure all rownames are included
                store_names <- gsub('`', '', rownames(lm_summary))
                rows_to_add = names_to_include[!(names_to_include %in% store_names)]
                lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
                rownames(lm_summary) <- c(store_names, rows_to_add)
              }
              para <- as.data.frame(lm_summary)[-1, -3]
              para$name <- rownames(lm_summary)[-1]
              return(para)
          }
          compare_function <- function(fit1, fit2) {
            return(anova(fit1, fit2, test = 'F')[2, 'Pr(>F)'])
          }
          gomp_function <- gomp_lm
          omp_function <- omp_lm
      } else { # Random effects
        ranef_function <- lme4::ranef
        model_function <-
              function(formula, data, na.action) {
                  return(lmerTest::lmer(
                      formula(formula), 
                      data = data, 
                      na.action = na.action))
              }
          summary_function <- function(fit, names_to_include) {
              lm_summary <- coef(summary(fit))
              if (nrow(lm_summary) < length(names_to_include)) { # If deficient rank, make sure all rownames are included
                store_names <- gsub('`', '', rownames(lm_summary))
                rows_to_add = names_to_include[!(names_to_include %in% store_names)]
                lm_summary <- rbind(lm_summary, matrix(rep(NaN, 5 * length(rows_to_add)), nrow=length(rows_to_add)))
                rownames(lm_summary) <- c(store_names, rows_to_add)
              }
              para <- as.data.frame(lm_summary)[-1, -c(3:4)]
              para$name <- rownames(lm_summary)[-1]
              return(para)
          }
          compare_function <- function(fit1, fit2) {
            return(anova(fit1, fit2, test = 'LRT')[2, 'Pr(>Chisq)'])
          }
          gomp_function <- gomp_lmer
          omp_function <- omp_lmer
      }
  }

  ##################
  # Logistic Model #
  ##################

  if (model == "logistic") {
    if (is.null(random_effects_formula)) { # Fixed effects only
      if (augment) {
        model_function <- function(formula, mm, weight_scheme, na.action) {
          weight_sch_current <<- weight_scheme # Needs to be global variable ?!
          
          glm_out <- glm(
            formula = formula(formula),
            family = 'binomial',
            data = mm,
            weights = weight_sch_current,
            na.action = na.action,
          )
          
          remove(weight_sch_current, pos = ".GlobalEnv") # Stop being global
          
          return(glm_out)
        }
        gomp_function <- gomp_glm_augment
        omp_function <- omp_glm_augment
      } else {
        model_function <-
          function(formula, data, na.action) {
            return(glm(
              formula(formula),
              data = data,
              family = 'binomial',
              na.action = na.action,
            ))
          }
        gomp_function <- gomp_glm
        omp_function <- omp_glm
      }
      summary_function <- function(fit, names_to_include) {
        lm_summary <- summary(fit)$coefficients
        if (nrow(lm_summary) < length(names_to_include)) { # If equal numbers of predictors and observations
          store_names <- gsub('`', '', rownames(lm_summary))
          rows_to_add = names_to_include[!(names_to_include %in% store_names)]
          lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
          rownames(lm_summary) <- c(store_names, rows_to_add)
        }
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
      compare_function <- function(fit1, fit2) {
        return(anova(fit1, fit2, test = 'Chisq')[2, 'Pr(>Chi)'])
      }
    } else { # Random effects
      ranef_function <- lme4::ranef
      if (augment) {
        model_function <-
          function(formula, mm, weight_scheme, na.action) {
            weight_sch_current <<- weight_scheme # Needs to be global variable ?!
            
            glm_out <- lme4::glmer(
              formula(formula), 
              data = mm, 
              family = 'binomial',
              na.action = na.action,
              weights = weight_sch_current,
              control = glmerControl(optimizer = 'bobyqa'))
            
            remove(weight_sch_current, pos = ".GlobalEnv") # Stop being global
            
            return(glm_out)
          }
        gomp_function <- gomp_glmer_augment
        omp_function <- omp_glmer_augment
      } else {
        model_function <-
          function(formula, data, na.action) {
            return(lme4::glmer(
              formula(formula), 
              data = data, 
              family = 'binomial',
              na.action = na.action,
              control = glmerControl(optimizer = 'bobyqa')))
          }
        gomp_function <- gomp_glmer
        omp_function <- omp_glmer
      }
      summary_function <- function(fit, names_to_include) {
        lm_summary <- coef(summary(fit))
        
        if (nrow(lm_summary) < length(names_to_include)) { # If deficient rank, make sure all rownames are included
          store_names <- gsub('`', '', rownames(lm_summary))
          rows_to_add = names_to_include[!(names_to_include %in% store_names)]
          lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
          rownames(lm_summary) <- c(store_names, rows_to_add)
        }
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
      compare_function <- function(fit1, fit2) {
        return(anova(fit1, fit2, test = 'LRT')[2, 'Pr(>Chisq)'])
      }
    }
  }
  
  #######################################
  # Init cluster for parallel computing #
  #######################################
  
  cluster <- NULL
  if (cores > 1) {
      logging::loginfo("Creating cluster of %s R processes", cores)
      cluster <- parallel::makeCluster(cores)
      parallel::clusterExport(cluster, c(ls(), function_vec),
                              envir = environment())
  }
  
  ##############################
  # Apply per-feature modeling #
  ##############################
  outputs <- pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
    # Load Required Packages
    for (lib in c(
      'dplyr',
      'pbapply',
      'lmerTest',
      'parallel',
      'lme4',
      'plyr',
      'TcGSA'
    )) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }
    
    # Extract Features One by One
    featuresVector <- features[, x]
    
    logging::loginfo(
        "Fitting model to feature number %d, %s",
        x,
        colnames(features)[x])
    
    # Make fitting matrix of features and metadata
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    
    # 0 or 1 observations
    if (length(unique(dat_sub$expr)) < 2) {
      output <- list()
      
      # List fixed effects that will be included
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, gomps, omps)
      
      # Build outputs
      output$para <- as.data.frame(matrix(NA, nrow = length(names_to_include), ncol = 3))
      output$para$name <- names_to_include
      
      output$residuals <- NA
      output$fitted <- NA
      if (!(is.null(random_effects_formula))) output$ranef <- NA
      output$fit <- NA
      
      colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
      output$para$feature <- colnames(features)[x]
      output$para$error <- ifelse(model == "logistic", "All logistic values are the same",
                                  "All LM values are the same")
      return(output)
    }
    
    # Missing first factor level
    missing_first_factor_level <- FALSE
    for (col in colnames(dat_sub)) {
      if(!is.numeric(dat_sub[,col])) {
        if (all(is.na(dat_sub$expr[dat_sub[,col] == levels(factor(dat_sub[,col]))[1]]))) {
          fixed_effects <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, gomps, omps)
          if (col %in% substr(fixed_effects, 1, nchar(col))) {
            missing_first_factor_level <- TRUE
          }
        }
      }
    }
      
    if (missing_first_factor_level) {
      output <- list()
      
      # List fixed effects that will be included
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, gomps, omps)
      
      # Build outputs
      output$para <- as.data.frame(matrix(NA, nrow = length(names_to_include), ncol = 3))
      output$para$name <- names_to_include
      
      output$residuals <- NA
      output$fitted <- NA
      if (!(is.null(random_effects_formula))) output$ranef <- NA
      output$fit <- NA
      
      colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
      output$para$feature <- colnames(features)[x]
      output$para$error <- "No data points have the baseline factor level"
      return(output)
    }

    # Augment logistic fitting
    if (augment & model == "logistic" & length(unique(featuresVector)) >= 2) {
      fit_and_message <- tryCatch({
        fit_list <- list()
        withCallingHandlers({ # Catch non-integer # successes first
          
          if (length(groups) > 0) {
            formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
            
            augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
            mm_input <- augmented_data[["mm_input"]]
            weight_scheme <- augmented_data[["weight_scheme"]]
            formula_new <- augmented_data[["new_formula"]]
            
            # Fit augmented model
            fit1 <-
              model_function(
                formula = formula_new,
                mm = mm_input, 
                weight_scheme = weight_scheme,
                na.action = na.exclude)
            
            fit_list <- list(fit1)
            
            for (group in groups) {
              # Remove terms that were added for the group
              formula_new <- generate_new_formula_groups(formula, random_effects_formula, dat_sub, groups, gomps, omps, group)

              # Fit augmented model
              fit1 <-
                model_function(
                  formula = formula_new,
                  mm = mm_input, 
                  weight_scheme = weight_scheme,
                  na.action = na.exclude)
              
              fit_list <- c(fit_list, list(fit1))
            }
            fit_list <- c(fit_list, NA)
          
            } else {
            formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
            
            # Augment data
            augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
            mm_input <- augmented_data[["mm_input"]]
            weight_scheme <- augmented_data[["weight_scheme"]]
            formula_new <- augmented_data[["new_formula"]]
            
            # Fit augmented model
            fit1 <-
              model_function(
                formula = formula_new,
                mm = mm_input, 
                weight_scheme = weight_scheme,
                na.action = na.exclude)
            
            fit_list <- list(fit1, NA)
          }
        }, warning=function(w) {
          if (w$message == "non-integer #successes in a binomial glm!") {
            # Still worked
            invokeRestart("muffleWarning")
          }
          return(fit_list)
        })
        }, warning = function(w) {
        message(paste("Feature", colnames(features)[x], ":", w))
        logging::logwarn(paste(
          "Fitting problem for feature",
          x,
          "a warning was issued"))
        
        fit_list <-
          try({
            # Data augmentation
            if (length(groups) > 0) {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
              mm_input <- augmented_data[["mm_input"]]
              weight_scheme <- augmented_data[["weight_scheme"]]
              formula_new <- augmented_data[["new_formula"]]
              
              # Fit augmented model
              fit1 <-
                model_function(
                  formula = formula_new,
                  mm = mm_input, 
                  weight_scheme = weight_scheme,
                  na.action = na.exclude)
              
              fit_list <- list(fit1)
              
              for (group in groups) {
                formula_new <- generate_new_formula_groups(formula, random_effects_formula, dat_sub, groups, gomps, omps, group)
                
                # Fit augmented model
                fit1 <-
                  model_function(
                    formula = formula_new,
                    mm = mm_input, 
                    weight_scheme = weight_scheme,
                    na.action = na.exclude)
                
                fit_list <- c(fit_list, list(fit1))
              }
              fit_list
            
              } else {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              # Augment data
              augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
              mm_input <- augmented_data[["mm_input"]]
              weight_scheme <- augmented_data[["weight_scheme"]]
              formula_new <- augmented_data[["new_formula"]]
              
              # Fit augmented model
              fit1 <-
                model_function(
                  formula = formula_new,
                  mm = mm_input, 
                  weight_scheme = weight_scheme,
                  na.action = na.exclude)
              
              list(fit1)
            }
          })
        if (inherits(fit_list, "try-error")) { # Sometimes warning process triggers error that gives unlisted result
          return(c(list(fit_list), list(w$message)))
        }
        return(c(fit_list, list(w$message)))
      
        }, error = function(err) { # Warn on augmented ataa
        fit_list <-
          try({
            # Data augmentation
            if (length(groups) > 0) {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
              mm_input <- augmented_data[["mm_input"]]
              weight_scheme <- augmented_data[["weight_scheme"]]
              formula_new <- augmented_data[["new_formula"]]
              
              # Fit augmented model
              fit1 <-
                model_function(
                  formula = formula_new,
                  mm = mm_input, 
                  weight_scheme = weight_scheme,
                  na.action = na.exclude)
              
              fit_list <- list(fit1)
              
              for (group in groups) {
                formula_new <- generate_new_formula_groups(formula, random_effects_formula, dat_sub, groups, gomps, omps, group)
                
                # Fit augmented model
                fit1 <-
                  model_function(
                    formula = formula_new,
                    mm = mm_input, 
                    weight_scheme = weight_scheme,
                    na.action = na.exclude)
                
                fit_list <- c(fit_list, list(fit1))
              }
              fit_list
            } else {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              # Augment data
              augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
              mm_input <- augmented_data[["mm_input"]]
              weight_scheme <- augmented_data[["weight_scheme"]]
              formula_new <- augmented_data[["new_formula"]]
              
              # Fit augmented model
              fit1 <-
                model_function(
                  formula = formula_new,
                  mm = mm_input, 
                  weight_scheme = weight_scheme,
                  na.action = na.exclude)
              
              list(fit1)
            }
          })
        return(c(list(fit_list), list(err$message)))
      })
    
      } else { # LM or non-augmented logistic
      fit_and_message <- tryCatch({
        if (length(groups) > 0) {
          formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
          
          fit1 <-
            model_function(
              formula_new, 
              data = dat_sub, 
              na.action = na.exclude)
          
          fit_list <- list(fit1)
          
          for (group in groups) {
            formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group], gomps, omps), collapse = " + "))

            fit1 <-
              model_function(
                formula_new, 
                data = dat_sub, 
                na.action = na.exclude)
            
            fit_list <- c(fit_list, list(fit1))
          }
          fit_list <- c(fit_list, NA)
        
          } else {
          formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
          
          fit1 <-
            model_function(
              formula_new, 
              data = dat_sub, 
              na.action = na.exclude)
          fit_list <- list(fit1, NA)
        }
        }, warning = function(w) {
        message(paste("Feature", colnames(features)[x], ":", w))
        logging::logwarn(paste(
          "Fitting problem for feature", 
          x, 
          "a warning was issued"))
        
        fit_list <-
          try({
            if (length(groups) > 0) {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))

              fit1 <-
                model_function(
                  formula_new, 
                  data = dat_sub, 
                  na.action = na.exclude)
              
              fit_list <- list(fit1)
              
              for (group in groups) {
                formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group], gomps, omps), collapse = " + "))

                fit1 <-
                  model_function(
                    formula_new, 
                    data = dat_sub, 
                    na.action = na.exclude)
                
                fit_list <- c(fit_list, list(fit1))
              }
              fit_list
            } else {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              fit1 <-
                model_function(
                  formula_new, 
                  data = dat_sub, 
                  na.action = na.exclude)
              list(fit1)
            }
          })
        if (inherits(fit_list, "try-error")) { # Sometimes warning process triggers error that gives unlisted result
          return(c(list(fit_list), list(w$message)))
        }
        return(c(fit_list, list(w$message)))
      
        }, error = function(err) {
        fit_list <-
          try({
            if (length(groups) > 0) {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))

              fit1 <-
                model_function(
                  formula_new, 
                  data = dat_sub, 
                  na.action = na.exclude)
              
              fit_list <- list(fit1)
              
              for (group in groups) {
                formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group], gomps, omps), collapse = " + "))

                fit1 <-
                  model_function(
                    formula_new, 
                    data = dat_sub, 
                    na.action = na.exclude)
                
                fit_list <- c(fit_list, list(fit1))
              }
              fit_list
            } else {
              formula_new <- formula(paste0(c(safe_deparse(formula), groups, gomps, omps), collapse = " + "))
              
              fit1 <-
                model_function(
                  formula_new, 
                  data = dat_sub, 
                  na.action = na.exclude)
              list(fit1)
            }
          })
        return(c(list(fit_list), list(err$message)))
      })
    }
    
    fit <- fit_and_message[[1]]
    
    # Gather Output
    output <- list()
    
    # Check for fitting errors and add on special predictors
    low_n_error <- FALSE
    if (all(!inherits(fit, "try-error"))) {
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, character(0), character(0), character(0))
      output$para <- summary_function(fit, c('(Intercept)', names_to_include))
      output$para <- output$para[names_to_include, , drop=F]
      
      n_uni_cols <- nrow(output$para)
      
      if (length(groups) > 0) {
        i <- 2
        for (group in groups) {
          if (all(!inherits(fit_and_message[[i]], "try-error"))) {
            output$para <- tryCatch({
              withCallingHandlers({ # Catch non-integer # successes first
                tmp_output <- rbind(output$para, list(NA, NA, compare_function(fit_and_message[[i]], fit_and_message[[1]]), group))
              }, warning=function(w) {
                if (w$message == "non-integer #successes in a binomial glm!") {
                  # Still worked
                  invokeRestart("muffleWarning")
                }
                return(tmp_output)
              })
            },
            warning = function(w) {
              rbind(output$para, list(NA, NA, NA, group))
            },
            error = function(err) {
              rbind(output$para, list(NA, NA, NA, group))
            })
            rownames(output$para) <- c(rownames(output$para)[-nrow(output$para)], group)
            
          } else {
            output$para <- rbind(output$para, list(NA, NA, NA, NA))
          }
          
          i <- i + 1
        }
      }
      
      if (length(gomps) > 0) {
        output$para <- tryCatch({
          withCallingHandlers({ # Catch non-integer # successes first
            gomp_output <- gomp_function(formula, random_effects_formula, dat_sub, groups, gomps, omps)
            rownames(gomp_output) <- gomp_output[,4]
            colnames(gomp_output) <- colnames(output$para)
            tmp_output <- rbind(output$para, gomp_output) 
            tmp_output
          }, warning=function(w) {
            if (w$message == "non-integer #successes in a binomial glm!") {
              # Still worked
              invokeRestart("muffleWarning")
            }
            return(tmp_output)
          })
        },
        warning = function(w) {
          org_row_num <- nrow(output$para)
          output$para[(org_row_num + 1) : (org_row_num + length(gomps)), 1:3] <- NA
          output$para[(org_row_num + 1) : (org_row_num + length(gomps)), 4] <- gomps
          output$para
        },
        error = function(err) {
          org_row_num <- nrow(output$para)
          output$para[(org_row_num + 1) : (org_row_num + length(gomps)), 1:3] <- NA
          output$para[(org_row_num + 1) : (org_row_num + length(gomps)), 4] <- gomps
          output$para
        })
        rownames(output$para) <- c(rownames(output$para)[-((nrow(output$para) - length(gomps) + 1) : 
                                                             nrow(output$para))], gomps)
      }
      
      if (length(omps) > 0) {
        for (omp in omps) {
          omp_levels <- paste0(omp, levels(dat_sub[[omp]]))
          output$para <- tryCatch({
            withCallingHandlers({ # Catch non-integer # successes first
              omp_output <- omp_function(formula, random_effects_formula, dat_sub, groups, gomps, omps)
              rownames(omp_output) <- omp_output[,4]
              colnames(omp_output) <- colnames(output$para)
              tmp_output <- rbind(output$para, omp_output) 
            }, warning=function(w) {
              if (w$message == "non-integer #successes in a binomial glm!") {
                # Still worked
                invokeRestart("muffleWarning")
              }
              return(tmp_output)
            })
          },
          warning = function(w) {
            org_row_num <- nrow(output$para)
            output$para[(org_row_num + 1) : (org_row_num + length(omp_levels)), 1:3] <- NA
            output$para[(org_row_num + 1) : (org_row_num + length(omp_levels)), 4] <- omp_levels
            output$para
          },
          error = function(err) {
            org_row_num <- nrow(output$para)
            output$para[(org_row_num + 1) : (org_row_num + length(omp_levels)), 1:3] <- NA
            output$para[(org_row_num + 1) : (org_row_num + length(omp_levels)), 4] <- omp_levels
            output$para
          })
        }
      }
      
      if (sum(!is.na(dat_sub$expr)) < 20 & nrow(output$para) > n_uni_cols) {
        output$para[(n_uni_cols + 1):nrow(output$para),1:3] <- NA
        low_n_error <- TRUE
      }

      # Check whether summaries are correct
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, gomps, omps)
      if (any(!(names_to_include %in% rownames(output$para)))) {
        
        # Don't worry about dropped factor levels
        missing_names <- names_to_include[!(names_to_include %in% rownames(output$para))]
        character_cols <- get_character_cols(dat_sub)
        if (!all(missing_names %in% character_cols)) {
          fit_properly <- FALSE
          fit_and_message[[length(fit_and_message)]] <- "Metadata dropped during fitting (rank deficient)"
        } else {
          fit_properly <- TRUE
        }
      } else { # No errors, summaries are correct
        fit_properly <- TRUE
      }
      
      # Prevents individual group levels from ending up in results (I think?)
      output$para <- output$para[rownames(output$para) %in% names_to_include,]
    } else { # Fit issue occurred
        fit_properly <- FALSE
    }
    
    if (fit_properly) {
      output$residuals <- residuals(fit)
      output$fitted <- fitted(fit)
      if (!(is.null(random_effects_formula))) {
        # Returns a list with a table for each random effect
        l <- ranef_function(fit)
        
        # Rename rows as random effect labels if only random intercepts
        if (length(l) == 1 & ncol(l[[1]]) == 1 & colnames(l[[1]])[1] == "(Intercept)") {
          d<-as.vector(unlist(l))
          names(d)<-unlist(lapply(l, row.names))
          d[setdiff(unique(metadata[,names(l)]), names(d))] <- NA
          d <- d[order(unique(metadata[,names(l)]))]
          output$ranef<-d
        } else {
          # Otherwise return the random effects list
          output$ranef <- l
        }
      }
      if (save_models) {
        output$fit <- fit
      } else {
        output$fit <- NA
      }
    
      } else { # Fitting issue
      logging::logwarn(paste(
        "Fitting problem for feature", 
        x, 
        "returning NA"))
      
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, gomps, omps)
      
      # Store NA values for missing outputs
      output$para <- as.data.frame(matrix(NA, nrow = length(names_to_include), ncol = 3))
      output$para$name <- names_to_include
      
      output$residuals <- NA
      output$fitted <- NA
      if (!(is.null(random_effects_formula))) output$ranef <- NA
      output$fit <- NA
    }

    colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
    output$para$feature <- colnames(features)[x]
    output$para$error <- fit_and_message[[length(fit_and_message)]]
    output$para$error[!is.na(output$para$error) & is.na(output$para$pval)] <- 'Fitting error'
    if (low_n_error) {
      output$para[(n_uni_cols + 1):nrow(output$para),"error"] <- "Too few observations for asymptotic test"
    }
    
    return(output)
    })
  
  # stop the cluster
  if (!is.null(cluster))
      parallel::stopCluster(cluster)
  
  # bind the results for each feature
  paras <-
      do.call(rbind, lapply(outputs, function(x) {
          return(x$para)
      }))
  residuals <-
      do.call(rbind, lapply(outputs, function(x) {
          return(x$residuals)
      }))
  
  row.names(residuals) <- colnames(features)

  fitted <-
    do.call(rbind, lapply(outputs, function(x) {
      return(x$fitted)
    }))
  row.names(fitted) <- colnames(features)   
  
  fits <-
    lapply(outputs, function(x) {
      return(x$fit)
    })
  names(fits) <- colnames(features)  
  
  # Return NULL rather than empty object if fits aren't saved
  if (all(is.na(fits))) {
    fits <- NULL
  }
  
  if (!(is.null(random_effects_formula))) {
    ranef <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$ranef)
      }))
    row.names(ranef) <- colnames(features) 
  }
    
  ################################
  # Apply correction to p-values #
  ################################
  
  paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
  
  #####################################################
  # Determine the metadata names from the model names #
  #####################################################
  
  metadata_names <- colnames(metadata)
  # order the metadata names by decreasing length
  metadata_names_ordered <-
      metadata_names[order(
          nchar(metadata_names), decreasing = TRUE)]
  # find the metadata name based on the match 
  # to the beginning of the string
  extract_metadata_name <- function(name) {
    tmp_val <- metadata_names_ordered[mapply(
      startsWith, 
      name, 
      metadata_names_ordered)][1]
    if (is.na(tmp_val)) {
      metadata_names_ordered[mapply(
        grepl, 
        metadata_names_ordered,
        name)][1]
    } else {
      return(tmp_val)
    }
  }
  paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
  # compute the value as the model contrast minus metadata
  paras$value <-
      mapply(function(x, y) {
          if (x == y)
              x
          else
              gsub("^\\:", "", sub(x, "", y))
      }, paras$metadata, paras$name)
  
  ##############################
  # Sort by decreasing q-value #
  ##############################
  
  paras <- paras[order(paras$qval, decreasing = FALSE), ]
  paras <-
      dplyr::select(
          paras,
          c('feature', 'metadata', 'value'),
          dplyr::everything())
  paras$model <- model
  rownames(paras)<-NULL
  
  if (!(is.null(random_effects_formula))) {
    return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = ranef, "fits" = fits))
  } else {
    return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = NULL, "fits" = fits))
  }
}





