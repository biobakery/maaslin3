#!/usr/bin/env Rscript
# Load Required Packages
for (lib in c(
  'dplyr',
  'pbapply',
  'lmerTest',
  'parallel',
  'lme4',
  'multcomp'
)) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# Function to augment data for logistic fitting
augment_data <- function(formula, random_effects_formula, dat_sub){
  dat_sub_new <- rbind(dat_sub, dat_sub, dat_sub)
  dat_sub_new$expr[(nrow(dat_sub) + 1) : (2 * nrow(dat_sub))] <- 1
  dat_sub_new$expr[(nrow(dat_sub) * 2 + 1) : (3 * nrow(dat_sub))] <- 0
  
  formula <- formula(formula)
  
  if (is.null(random_effects_formula)) { # No random effects
    p <- ncol(model.matrix(formula, dat_sub)) - 1
  } else { # With random effects
    p <- ncol(lme4::lFormula(formula, dat_sub)$X) - 1
  }
  
  # Calculate the weights
  weight_scheme <- c(rep(1, nrow(dat_sub)), rep(p / (2 * nrow(dat_sub)), 2 * nrow(dat_sub)))
  
  return(list(mm_input = dat_sub_new, weight_scheme = weight_scheme, new_formula = formula))
}

safe_deparse <- function(formula) {
  paste0(trimws(deparse(formula)), collapse = " ")
}

# Extract a predictor of the form: predictor_type(variable)
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

# Get all fixed effects, not assuming a model has or can be fit
get_fixed_effects <- function(formula, random_effects_formula, dat_sub, groups, ordereds) {
  names_to_include <- c()
  if (is.null(random_effects_formula)) { # Fixed and group effects only
    names_to_include <- colnames(model.matrix(formula(gsub("^expr ", "", safe_deparse(formula))), dat_sub))
    names_to_include <- names_to_include[names_to_include != "(Intercept)"]
  } else { # Random effects
    patterns <- paste0("(", unlist(lme4::findbars(formula(gsub("^expr ", "", safe_deparse(formula))))), ")")
    
    for (pattern in patterns) {
      fixed_effects_only <- gsub(pattern, "", 
                       paste0(trimws(safe_deparse(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = " "), 
                       fixed = TRUE)
      fixed_effects_only <- gsub("[+ ]+$", "", fixed_effects_only)
      fixed_effects_only <- gsub("\\+\\s*\\++", "+", fixed_effects_only)
      formula <- formula(fixed_effects_only)
    }
    
    names_to_include <- colnames(model.matrix(formula, dat_sub))
    names_to_include <- names_to_include[names_to_include != "(Intercept)"]
  }
  
  ordered_levels <- c()
  for (ordered in ordereds) {
    ordered_levels <- c(ordered_levels, paste0(ordered, levels(dat_sub[[ordered]])[-1]))
  }
  
  names_to_include <- c(names_to_include, groups, ordered_levels)
  return(names_to_include)
}

# Get non-baseline names for all non-numeric columns
get_character_cols <- function(dat_sub) {
  all_factors <- c()
  for (col in colnames(dat_sub)) {
    if (!is.numeric(dat_sub[,col])) {
      if (is.factor(dat_sub[,col])) {
        factor_levels <- levels(dat_sub[,col])
      } else {
        factor_levels <- levels(factor(dat_sub[,col]))
      }
      # All factor levels except baseline
      all_factors <- c(all_factors, paste0(col, unique(factor_levels[-1])))
    }
  }
  return(all_factors)
}

# Get joint significance for zeros and non-zeros
add_joint_signif <- function(fit_data_abundance, fit_data_prevalence, analysis_method, correction) {
  # Subset to shared columns
  fit_data_prevalence_signif <- fit_data_prevalence$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_prevalence_signif) <- c("feature", "metadata", "value", "name", "logistic", "logistic_error")
  fit_data_abundance_signif <- fit_data_abundance$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_abundance_signif) <- c("feature", "metadata", "value", "name", analysis_method, "LM_error")
  
  # Join and check linear and logistic pieces
  merged_signif <- dplyr::full_join(unique(fit_data_prevalence_signif), unique(fit_data_abundance_signif), 
                             by=c("feature", "metadata", "value", "name"))
  
  # Stop everything and show the difference between the logistic and linear models fit
  if (nrow(merged_signif) != nrow(unique(fit_data_prevalence_signif)) | nrow(merged_signif) != nrow(unique(fit_data_abundance_signif))) {
    print(nrow(unique(merged_signif)))
    print(nrow(unique(fit_data_prevalence_signif)))
    print(nrow(unique(fit_data_abundance_signif)))
    print(dplyr::anti_join(unique(fit_data_prevalence_signif), unique(fit_data_abundance_signif), by=c("feature", "metadata", "value", "name")))
    print(dplyr::anti_join(unique(fit_data_abundance_signif), unique(fit_data_prevalence_signif), by=c("feature", "metadata", "value", "name")))
    stop("Merged significance tables have different associations. This is likely a package error due to unexpected data or models.")
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
  
  return(list(append_joint(fit_data_abundance, merged_signif), append_joint(fit_data_prevalence, merged_signif)))
}

# Take logistic or LM component and add on the merged significance pieces
append_joint <- function(outputs, merged_signif) {
  merged_signif <- merged_signif[,c("feature", "metadata", "value", "name", "pval_joint", "qval_joint")]
  tmp_colnames <- colnames(outputs$results)
  tmp_colnames <- dplyr::case_when(tmp_colnames == "pval" ~ "pval_individual",
                            tmp_colnames == "qval" ~ "qval_individual",
                            TRUE ~ tmp_colnames)
  colnames(outputs$results) <- tmp_colnames
  
  merged_signif <- merge(outputs$results, merged_signif, 
                         by=c("feature", "metadata", "value", "name"))
  
  merged_signif <- merged_signif[order(merged_signif$qval_joint),]
  
  return(merged_signif)
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
    cores = 1,
    median_comparison = FALSE,
    median_comparison_threshold = 0.1,
    feature_specific_covariate = NULL,
    feature_specific_covariate_name = NULL,
    feature_specific_covariate_record = NULL) {
  function_vec <- c("augment_data", "safe_deparse", "extract_special_predictor", 
                    "get_fixed_effects", "get_character_cols", 
                    "add_joint_signif", "append_joint")
      
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
  
  # Extract ordered components
  extract_special_predictor_out <- extract_special_predictor(formula, 'ordered')
  formula <- extract_special_predictor_out[[1]]
  ordereds <- extract_special_predictor_out[[2]]
  
  #############################################################
  # Determine the function and summary for the model selected #
  #############################################################
  
  optimizers <- c('nloptwrap', 'nlminbwrap', 'bobyqa', 'Nelder_Mead')
  optCtrlList <- list(list(maxeval = 100000), list(maxit = 1500), 
                      list(maxfun = 100000), list(maxfun = 100000))
  
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
      } else { # Random effects
        ranef_function <- lme4::ranef
        model_function <-
              function(formula, data, na.action) {
                index <- 1
                
                while (index < length(optimizers)) {
                  tryCatch({
                    return(lmerTest::lmer(
                      formula(formula), 
                      data = data, 
                      na.action = na.action, 
                      control = lme4::lmerControl(optimizer = optimizers[index],
                                            optCtrl = optCtrlList[[index]])))
                  }, warning = function(w) {
                    'warning'
                  }, error = function(e) {
                    'error'
                  })
                  
                  # Something warned or errored if here
                  index <- index + 1
                }
                
                return(lmerTest::lmer(
                  formula(formula), 
                  data = data, 
                  na.action = na.action, 
                  control = lme4::lmerControl(optimizer = optimizers[index],
                                        optCtrl = optCtrlList[[index]])))
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
      }
  }

  ##################
  # Logistic Model #
  ##################

  if (model == "logistic") {
    if (is.null(random_effects_formula)) { # Fixed effects only
      if (augment) {
        model_function <- function(formula, mm, weight_scheme, na.action) {
          assign("weight_sch_current", weight_scheme, envir = environment(formula))
          
          glm_out <- glm(
            formula = formula(formula),
            family = 'binomial',
            data = mm,
            weights = weight_sch_current,
            na.action = na.action,
          )
          
          return(glm_out)
        }
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
    } else { # Random effects
      ranef_function <- lme4::ranef
      if (augment) {
        model_function <-
          function(formula, mm, weight_scheme, na.action) {
            assign("weight_sch_current", weight_scheme, envir = environment(formula))

            index <- 1
            
            while (index < length(optimizers)) {
              glm_out <- tryCatch({
                withCallingHandlers({ # Catch non-integer # successes first
                  fit1 <- lme4::glmer(
                    formula(formula), 
                    data = mm, 
                    family = 'binomial',
                    na.action = na.action,
                    weights = weight_sch_current,
                    control = lme4::glmerControl(optimizer = optimizers[index],
                                           optCtrl = optCtrlList[[index]]))
                }, warning=function(w) {
                  if (w$message == "non-integer #successes in a binomial glm!") {
                    # Still worked
                    invokeRestart("muffleWarning")
                  }
                })
              }, warning = function(w) {
                'warning'
              }, error = function(e) {
                'error'
              })
              
              # Something warned or errored if here
              if (is.character(glm_out)) {
                index <- index + 1
              } else  {
                break
              }
            }
            
            if (is.character(glm_out)) {
              withCallingHandlers({ # Catch non-integer # successes first
                  fit1 <- lme4::glmer(
                    formula(formula), 
                    data = mm, 
                    family = 'binomial',
                    na.action = na.action,
                    weights = weight_sch_current,
                    control = lme4::glmerControl(optimizer = optimizers[index],
                                           optCtrl = optCtrlList[[index]]))
              }, warning=function(w) {
                if (w$message == "non-integer #successes in a binomial glm!") {
                  # Still worked
                  invokeRestart("muffleWarning")
                }
              })
              return(fit1)
            } else  {
              return(glm_out)
            }
          }
      } else {
        model_function <-
          function(formula, data, na.action) {
            index <- 1
            
            while (index < length(optimizers)) {
              glm_out <- tryCatch({
                lme4::glmer(
                  formula(formula), 
                  data = data, 
                  family = 'binomial',
                  na.action = na.action,
                  control = lme4::glmerControl(optimizer = optimizers[index],
                                         optCtrl = optCtrlList[[index]]))
              }, warning = function(w) {
                'warning'
              }, error = function(e) {
                'error'
              })
              
              # Something warned or errored if here
              if (is.character(glm_out)) {
                index <- index + 1
              } else  {
                break
              }
            }
            
            if (is.character(glm_out)) {
              return(lme4::glmer(
                formula(formula), 
                data = data, 
                family = 'binomial',
                na.action = na.action,
                control = lme4::glmerControl(optimizer = optimizers[index],
                                       optCtrl = optCtrlList[[index]])))
            } else  {
              return(glm_out)
            }
          }
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
      'multcomp'
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
    if (!is.null(feature_specific_covariate)) {
      covariateVector <- feature_specific_covariate[, x]
      
      dat_sub <- data.frame(expr = as.numeric(featuresVector), 
                            feature_specific_covariate = covariateVector,
                            metadata)
      
      colnames(dat_sub)[colnames(dat_sub) == 'feature_specific_covariate'] <- 
        feature_specific_covariate_name
    } else {
      dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    }

    # 0 or 1 observations
    if (length(unique(dat_sub$expr)) < 2) {
      output <- list()
      
      # List fixed effects that will be included
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, ordereds)
      
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
      if(is.factor(dat_sub[,col])) {
        if (all(is.na(dat_sub$expr[dat_sub[,col] == levels(dat_sub[,col])[1]]))) {
          fixed_effects <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, ordereds)
          if (col %in% substr(fixed_effects, 1, nchar(col))) {
            missing_first_factor_level <- TRUE
          }
        }
      }
    }
    
    if (missing_first_factor_level) {
      output <- list()
      
      # List fixed effects that will be included
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, ordereds)
      
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
      warning_message <- NA
      error_message <- NA
      fit1 <- tryCatch({
        withCallingHandlers({
          withCallingHandlers({ # Catch non-integer # successes first
            formula_new <- formula(paste0(c(safe_deparse(formula), groups, ordereds), collapse = " + "))
            
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
            fit1
          }, warning=function(w) {
            if (w$message == "non-integer #successes in a binomial glm!") {
              # Still worked
              invokeRestart("muffleWarning")
            }
          })
        }, warning = function(w) {
          message(paste("Feature", colnames(features)[x], ":", w))
          logging::logwarn(paste(
            "Fitting problem for feature", 
            x, 
            "a warning was issued"))
          
          warning_message <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })}, 
        error = function(err) {
          error_message <<- err$message
          error_obj <- structure(list(message = conditionMessage(err)), class = "try-error")
          return(error_obj)
        })
        if (!is.na(error_message)) {
          fit_and_message <- c(list(fit1), list(error_message))
          suppressWarnings(remove(error_message, pos = ".GlobalEnv"))
        } else if (!is.na(warning_message)) {
          fit_and_message <- c(list(fit1), list(warning_message))
          suppressWarnings(remove(warning_message, pos = ".GlobalEnv"))
        } else {
          fit_and_message <- c(list(fit1), NA)
        }
      } else { # LM or non-augmented logistic
        warning_message <- NA
        error_message <- NA
        fit1 <- tryCatch({
          withCallingHandlers({
            formula_new <- formula(paste0(c(safe_deparse(formula), groups, ordereds), collapse = " + "))
            
            fit1 <-
              model_function(
                formula_new, 
                data = dat_sub, 
                na.action = na.exclude)
            fit1
          }, warning = function(w) {
            message(paste("Feature", colnames(features)[x], ":", w))
            logging::logwarn(paste(
              "Fitting problem for feature", 
              x, 
              "a warning was issued"))
            
            warning_message <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          })}, 
          error = function(err) {
            error_message <<- err$message
            error_obj <- structure(list(message = conditionMessage(err)), class = "try-error")
            return(error_obj)
          })
        if (!is.na(error_message)) {
          fit_and_message <- c(list(fit1), list(error_message))
          suppressWarnings(remove(error_message, pos = ".GlobalEnv"))
        } else if (!is.na(warning_message)) {
          fit_and_message <- c(list(fit1), list(warning_message))
          suppressWarnings(remove(warning_message, pos = ".GlobalEnv"))
        } else {
          fit_and_message <- c(list(fit1), NA)
        }
        
        
      }
    
    fit <- fit_and_message[[1]]
    
    # Gather Output
    output <- list()
    
    # Check for fitting errors and add on special predictors
    low_n_error <- FALSE
    if (all(!inherits(fit, "try-error"))) {
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, character(0), character(0))
      output$para <- suppressWarnings(summary_function(fit, c('(Intercept)', names_to_include))) # Suppress warnings about variance-covariance matrix calculation
      output$para <- output$para[names_to_include, , drop=F]
      
      n_uni_cols <- nrow(output$para)
      
      if (length(groups) > 0) {
        for (group in groups) {
          output$para <- tryCatch({
            withCallingHandlers({ # Catch non-integer # successes first
              if(is.null(random_effects_formula)) { # Fixed effects
                if (model == "logistic") {
                  pval_new <- tryCatch({anova(fit, test = 'LRT')[group,'Pr(>Chi)']},
                                       error = function(err) { NA })
                } else {
                  pval_new <- tryCatch({anova(fit)[group,'Pr(>F)']},
                                       error = function(err) { NA })
                }
              } else { # Random effects
                if (model == "logistic") {
                  if (augment) {
                    assign("weight_sch_current", weight_scheme, envir = environment(formula))

                    fit_new <-
                      model_function(
                        formula = update.formula(formula(fit), formula(paste0('~.-', group))),
                        mm = mm_input, 
                        weight_scheme = weight_sch_current,
                        na.action = na.exclude)
                  } else {
                    fit_new <-
                      model_function(
                        update.formula(formula(fit), formula(paste0('~.-', group))), 
                        data = dat_sub, 
                        na.action = na.exclude)
                  }
                  
                  pval_new <- tryCatch({anova(fit_new, fit)[2, 'Pr(>Chisq)']},
                                       error = function(err) { NA })
                } else {
                  contrast_mat <- matrix(0, ncol = length(lme4::fixef(fit)), nrow = length(levels(dat_sub[[group]])[-1]))
                  contrast_mat[1:length(levels(dat_sub[[group]])[-1]),
                               which(names(lme4::fixef(fit)) %in% paste0(group, levels(dat_sub[[group]])[-1]))] <- 
                    diag(1, nrow = length(levels(dat_sub[[group]])[-1]))
                  
                  pval_new <- tryCatch({lmerTest::contest(fit, contrast_mat, rhs = rep(0, length(levels(dat_sub[[group]])[-1])))[['Pr(>F)']]},
                                          error = function(err) { NA })
                }
              }
              
              tmp_output <- rbind(output$para, list(NA, NA, pval_new, group))
            }, warning=function(w) {
              if (w$message == "non-integer #successes in a binomial glm!") {
                # Still worked
                invokeRestart("muffleWarning")
              }
            })
          },
          warning = function(w) {
            rbind(output$para, list(NA, NA, NA, group))
          },
          error = function(err) {
            rbind(output$para, list(NA, NA, NA, group))
          })
          rownames(output$para) <- c(rownames(output$para)[-nrow(output$para)], group)
          
        }
      }
      
      if (length(ordereds) > 0) {
        for (ordered in ordereds) {
          ordered_levels <- paste0(ordered, levels(dat_sub[[ordered]])[-1])
          output$para <- tryCatch({
            withCallingHandlers({ # Catch non-integer # successes first
              if(is.null(random_effects_formula)) { # Fixed effects
                if (any(!ordered_levels %in% names(coef(fit)))) {
                  fit_and_message[[length(fit_and_message)]] <- "Error: Some ordered levels are missing"
                  stop("Some ordered levels are missing")
                }
                
                contrast_mat <- matrix(0, ncol = length(coef(fit)), nrow = length(levels(dat_sub[[ordered]])[-1]))
                
                cols_to_add_1s <- which(names(coef(fit)) %in% ordered_levels)
                contrast_mat[1, cols_to_add_1s[1]] <- 1
                for (i in seq_along(cols_to_add_1s[-1])) {
                  contrast_mat[i + 1, cols_to_add_1s[-1][i]] <- 1
                  contrast_mat[i + 1, cols_to_add_1s[i]] <- -1
                }
                
                pvals_new <- c()
                coefs_new <- c()
                sigmas_new <- c()
                for (row_num in 1:nrow(contrast_mat)) {
                  contrast_vec <- t(matrix(contrast_mat[row_num,]))
                  pvals_new <- c(pvals_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                      rhs = 0))$test$pvalues},
                                        error = function(err) { NA }))
                  coefs_new <- c(coefs_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                                   rhs = 0))$test$coefficients},
                                                     error = function(err) { NA }))
                  sigmas_new <- c(sigmas_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                                   rhs = 0))$test$sigma},
                                                     error = function(err) { NA }))
                }
              } else { # Random effects linear
                if (any(!ordered_levels %in% names(lme4::fixef(fit)))) {
                  fit_and_message[[length(fit_and_message)]] <- "Error: Some ordered levels are missing"
                  stop("Some ordered levels are missing")
                }
                
                contrast_mat <- matrix(0, ncol = length(lme4::fixef(fit)), nrow = length(levels(dat_sub[[ordered]])[-1]))
                
                cols_to_add_1s <- which(names(lme4::fixef(fit)) %in% ordered_levels)
                contrast_mat[1, cols_to_add_1s[1]] <- 1
                for (i in seq_along(cols_to_add_1s[-1])) {
                  contrast_mat[i + 1, cols_to_add_1s[-1][i]] <- 1
                  contrast_mat[i + 1, cols_to_add_1s[i]] <- -1
                }
                
                if (model == "logistic") {
                  pvals_new <- c()
                  coefs_new <- c()
                  sigmas_new <- c()
                  for (row_num in 1:nrow(contrast_mat)) {
                    contrast_vec <- t(matrix(contrast_mat[row_num,]))
                    pvals_new <- c(pvals_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                                     rhs = 0))$test$pvalues},
                                                       error = function(err) { NA }))
                    coefs_new <- c(coefs_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                                     rhs = 0))$test$coefficients},
                                                       error = function(err) { NA }))
                    sigmas_new <- c(sigmas_new, tryCatch({summary(multcomp::glht(fit, linfct = contrast_vec, 
                                                                       rhs = 0))$test$sigma},
                                                         error = function(err) { NA }))
                  }
                } else {
                  pvals_new <- c()
                  coefs_new <- c()
                  sigmas_new <- c()
                  for (row_num in 1:nrow(contrast_mat)) {
                    contrast_vec <- t(matrix(contrast_mat[row_num,]))
                    pvals_new <- c(pvals_new, tryCatch({lmerTest::contest(fit, matrix(contrast_vec, T), rhs = 0)[['Pr(>F)']]},
                                                       error = function(err) { NA }))
                    coefs_new <- c(coefs_new, tryCatch({contrast_vec %*% lme4::fixef(fit)},
                                                       error = function(err) { NA }))
                    sigmas_new <- c(sigmas_new, tryCatch({sqrt((contrast_vec %*% vcov(fit) %*% t(contrast_vec))[1,1])},
                                                         error = function(err) { NA }))
                  }
                }
              }
              
              tmp_output <- rbind(output$para, list(coefs_new, 
                                                    sigmas_new,
                                                    pvals_new, ordered_levels))
              rownames(tmp_output)[(nrow(tmp_output) - length(ordered_levels) + 1) : nrow(tmp_output)] <- ordered_levels
              tmp_output
            }, warning=function(w) {
              if (w$message == "non-integer #successes in a binomial glm!") {
                # Still worked
                invokeRestart("muffleWarning")
              }
            })
          },
          warning = function(w) {
            org_row_num <- nrow(output$para)
            output$para[(org_row_num + 1) : (org_row_num + length(ordered_levels)), 1:3] <- NA
            output$para[(org_row_num + 1) : (org_row_num + length(ordered_levels)), 4] <- ordered_levels
            rownames(output$para)[(org_row_num + 1) : (org_row_num + length(ordered_levels))] <- ordered_levels
            output$para
          },
          error = function(err) {
            org_row_num <- nrow(output$para)
            output$para[(org_row_num + 1) : (org_row_num + length(ordered_levels)), 1:3] <- NA
            output$para[(org_row_num + 1) : (org_row_num + length(ordered_levels)), 4] <- ordered_levels
            rownames(output$para)[(org_row_num + 1) : (org_row_num + length(ordered_levels))] <- ordered_levels
            output$para
          })
        }
      }
      
      # If any groups, gomps, or omps are added, make sure the sample size is big enough for an asymptotic test
      # if (sum(!is.na(dat_sub$expr)) < 20 & nrow(output$para) > n_uni_cols) {
      #   output$para[(n_uni_cols + 1):nrow(output$para),1:3] <- NA
      #   low_n_error <- TRUE
      # }

      # Check whether summaries are correct
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, ordereds)
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
      if (median_comparison) {
        output$fit <- fit
      } else {
        if (save_models) {
          output$fit <- fit
        } else {
          output$fit <- NA
        }
      }
    } else { # Fitting issue
      logging::logwarn(paste(
        "Fitting problem for feature", 
        x, 
        "returning NA"))
      
      names_to_include <- get_fixed_effects(formula, random_effects_formula, dat_sub, groups, ordereds)
      
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
    output$para$error[is.na(output$para$error) & is.na(output$para$pval)] <- 'Fitting error (NA p-value returned from fitting procedure)'
    
    # Checks for poor model fits or possible linear separability
    # if (fit_properly & augment & model == 'logistic') {
    #   cooks_distances <- suppressWarnings(cooks.distance(fit))
    #   weights_large_cooks_distances <- weight_scheme[!is.na(cooks_distances) & cooks_distances > 1]
    #   if (mean(weights_large_cooks_distances != 1) > 0.95 & 
    #       length(weights_large_cooks_distances) > 0.01 * nrow(dat_sub)) {
    #     output$para$error <- ifelse(!is.na(output$para$error), output$para$error,
    #                                 paste0(100 * round(mean(weights_large_cooks_distances != 1), 3), 
    #                                        "% of Cooks distances >1 came from the augmented values, suggesting the model fit is heavily influenced by the data augmentation. Check the model before using resulting small p-values."))
    #   }
    # 
    #   kept_stderr <- output$para$name[!is.na(output$para$stderr)]
    #   if (any(output$para$stderr[!is.na(output$para$stderr)] < 
    #           0.01 / apply(model.matrix(fit)[, -1, drop = F], 2, sd)[kept_stderr])) {
    #     output$para$error <- ifelse(!is.na(output$para$error), output$para$error,
    #                                 'StdErr is below 0.01/SD(X) for at least one coefficient, suggesting a silent poor model fit. Check the model before using resulting small p-values.')
    #   }
    # }
      
    # if (low_n_error) {
    #   output$para[(n_uni_cols + 1):nrow(output$para),"error"] <- "Too few observations for asymptotic test"
    # }
    
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
  
  # Adjust the significance levels based on tests against the median
  if (median_comparison) {
    logging::loginfo(
      "Performing tests against medians")
    
    if (length(ordereds) > 0) {
      ordered_levels <- c()
      for (ordered in ordereds) {
        ordered_levels <- c(ordered_levels, paste0(ordered, levels(metadata[[ordered]])[-1]))
      }
    } else {
      ordered_levels <- c()
    }

    final_paras <- paras[!is.na(paras$error) | paras$name %in% groups,]
    paras <- paras[is.na(paras$error) & !paras$name %in% groups,]
    for (metadata_variable in unique(paras$name)) {
      paras_sub <- paras[paras$name == metadata_variable,]
      
      # Get the current median, excluding coefficients with NA p-values or 
      # p-values close to 1 since these could be model misfits
      cur_median <- median(paras_sub$coef[!is.na(paras_sub$pval) & paras_sub$pval < 0.95], na.rm=T)
      
      pvals_new <- vector()
      if (metadata_variable %in% ordered_levels) {
        ordered <- ordereds[which(startsWith(metadata_variable, ordereds))]
        
        for (feature in paras_sub$feature) {
          if(is.null(random_effects_formula)) { # Fixed effects
            cur_fit <- fits[[feature]]
            
            if (!metadata_variable %in% names(coef(cur_fit))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            if (any(!paste0(ordered, levels(metadata[[ordered]])[-1]) %in% names(coef(cur_fit)))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            mm_variable <- model.matrix(cur_fit)[, metadata_variable]
            if (any(!unique(mm_variable[!is.na(mm_variable)]) %in% c(0,1))) {
              median_comparison_threshold_updated <- median_comparison_threshold / sd(mm_variable)
            } else {
              median_comparison_threshold_updated <- median_comparison_threshold
            }
            
            if (is.na(coef(cur_fit)[which(names(coef(cur_fit)) == metadata_variable)])) {
              pval_new_current <- NA
            } else if (abs(coef(cur_fit)[which(names(coef(cur_fit)) == metadata_variable)] - 
                    cur_median) < median_comparison_threshold_updated) {
              pval_new_current <- 1
            } else {
              contrast_mat <- matrix(0, ncol = length(coef(cur_fit)), nrow = length(levels(metadata[[ordered]])[-1]))
              
              cols_to_add_1s <- which(names(coef(cur_fit)) %in% paste0(ordered, levels(metadata[[ordered]])[-1]))
              contrast_mat[1, cols_to_add_1s[1]] <- 1
              for (i in seq_along(cols_to_add_1s[-1])) {
                contrast_mat[i + 1, cols_to_add_1s[-1][i]] <- 1
                contrast_mat[i + 1, cols_to_add_1s[i]] <- -1
              }
  
              contrast_vec <- t(matrix(contrast_mat[which(paste0(ordered, levels(metadata[[ordered]])[-1]) == metadata_variable),]))
              pval_new_current <- tryCatch({summary(multcomp::glht(cur_fit, linfct = contrast_vec, 
                                                               rhs = cur_median))$test$pvalues},
                                                 error = function(err) { NA })
            }

            pvals_new <- c(pvals_new, pval_new_current)
          } else { # Random effects linear
            cur_fit <- fits[[feature]]
            
            if (!metadata_variable %in% names(lme4::fixef(cur_fit))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            if (any(!paste0(ordered, levels(metadata[[ordered]])[-1]) %in% names(lme4::fixef(cur_fit)))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            mm_variable <- model.matrix(cur_fit)[, metadata_variable]
            if (any(!unique(mm_variable[!is.na(mm_variable)]) %in% c(0,1))) {
              median_comparison_threshold_updated <- median_comparison_threshold / sd(mm_variable)
            } else {
              median_comparison_threshold_updated <- median_comparison_threshold
            }
            
            contrast_mat <- matrix(0, ncol = length(lme4::fixef(cur_fit)), nrow = length(levels(metadata[[ordered]])[-1]))
            cols_to_add_1s <- which(names(lme4::fixef(cur_fit)) %in% paste0(ordered, levels(metadata[[ordered]])[-1]))
            contrast_mat[1, cols_to_add_1s[1]] <- 1
            for (i in seq_along(cols_to_add_1s[-1])) {
              contrast_mat[i + 1, cols_to_add_1s[-1][i]] <- 1
              contrast_mat[i + 1, cols_to_add_1s[i]] <- -1
            }
            contrast_vec <- t(matrix(contrast_mat[which(paste0(ordered, levels(metadata[[ordered]])[-1]) == metadata_variable),]))
            
            if (model == "logistic") {
              if (is.na(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)])) {
                pval_new_current <- NA
              } else if (abs(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)] - 
                             cur_median) < median_comparison_threshold_updated) {
                pval_new_current <- 1
              } else {
                pval_new_current <- tryCatch({summary(multcomp::glht(cur_fit, linfct = contrast_vec, 
                                                           rhs = cur_median))$test$pvalues},
                                             error = function(err) { NA })
              }
              pvals_new <- c(pvals_new, pval_new_current)
            } else {
              if (is.na(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)])) {
                pval_new_current <- NA
              } else if (abs(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)] - 
                             cur_median) < median_comparison_threshold_updated) {
                pval_new_current <- 1
              } else {
                pval_new_current <- tryCatch({lmerTest::contest(cur_fit, matrix(contrast_vec, T), rhs = cur_median)[['Pr(>F)']]},
                                                   error = function(err) { NA })
              }
              
              pvals_new <- c(pvals_new, pval_new_current)
            }
          }
        }
      } else {
        for (feature in paras_sub$feature) {
          if(is.null(random_effects_formula)) { # Fixed effects
            cur_fit <- fits[[feature]]
            
            if (!metadata_variable %in% names(coef(cur_fit))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            mm_variable <- model.matrix(cur_fit)[, metadata_variable]
            if (any(!unique(mm_variable[!is.na(mm_variable)]) %in% c(0,1))) {
              median_comparison_threshold_updated <- median_comparison_threshold / sd(mm_variable)
            } else {
              median_comparison_threshold_updated <- median_comparison_threshold
            }
            
            contrast_vec <- rep(0, length(coef(cur_fit)))
            contrast_vec[which(names(coef(cur_fit)) == metadata_variable)] <- 1
            
            if (is.na(coef(cur_fit)[which(names(coef(cur_fit)) == metadata_variable)])) {
              pval_new_current <- NA
            } else if (abs(coef(cur_fit)[which(names(coef(cur_fit)) == metadata_variable)] - 
                           cur_median) < median_comparison_threshold_updated) {
              pval_new_current <- 1
            } else {
              pval_new_current <- tryCatch({summary(multcomp::glht(fits[[feature]], linfct = matrix(contrast_vec, T), 
                                                         rhs = cur_median))$test$pvalues[1]},
                                           error = function(err) { NA })
            }
            
            pvals_new <- c(pvals_new, pval_new_current)
          } else { # Random effects linear
            cur_fit <- fits[[feature]]
            
            if (!metadata_variable %in% names(lme4::fixef(cur_fit))) {
              pvals_new <- c(pvals_new, NA)
              next
            }
            
            mm_variable <- model.matrix(cur_fit)[, metadata_variable]
            if (any(!unique(mm_variable[!is.na(mm_variable)]) %in% c(0,1))) {
              median_comparison_threshold_updated <- median_comparison_threshold / sd(mm_variable)
            } else {
              median_comparison_threshold_updated <- median_comparison_threshold
            }
            
            contrast_vec <- rep(0, length(lme4::fixef(cur_fit)))
            contrast_vec[which(names(lme4::fixef(cur_fit)) == metadata_variable)] <- 1
            
            if (model == "logistic") {
              if (is.na(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)])) {
                pval_new_current <- NA
              } else if (abs(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)] - 
                             cur_median) < median_comparison_threshold_updated) {
                pval_new_current <- 1
              } else {
                pval_new_current <- tryCatch({summary(multcomp::glht(fits[[feature]], linfct = matrix(contrast_vec, T), 
                                                           rhs = cur_median))$test$pvalues[1]},
                                             error = function(err) { NA })
              }
              
              pvals_new <- c(pvals_new, pval_new_current)
            } else {
              if (is.na(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)])) {
                pval_new_current <- NA
              } else if (abs(lme4::fixef(cur_fit)[which(names(lme4::fixef(cur_fit)) == metadata_variable)] - 
                             cur_median) < median_comparison_threshold_updated) {
                pval_new_current <- 1
              } else {
                pval_new_current <- tryCatch({lmerTest::contest(cur_fit, matrix(contrast_vec, T), rhs = cur_median)[['Pr(>F)']]},
                                             error = function(err) { NA })
              }
              
              pvals_new <- c(pvals_new, pval_new_current)
            }
          }
        }
      }
      
      paras_sub$pval <- pvals_new
      final_paras <- rbind(final_paras, paras_sub)
    }
    
    paras <- final_paras
  }
  
  # Return NULL rather than empty object if fits aren't saved
  if (all(is.na(fits)) | !save_models) {
    fits <- NULL
  }
  
  if (!(is.null(random_effects_formula))) {
    ranef <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$ranef)
      }))
    row.names(ranef) <- colnames(features) 
  }
    
  #####################################################
  # Determine the metadata names from the model names #
  #####################################################
  
  metadata_names <- colnames(metadata)
  if (!is.null(feature_specific_covariate)) {
    metadata_names <- union(metadata_names, feature_specific_covariate_name)
  }
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
  
  if (!is.null(feature_specific_covariate_record)) {
    if (!feature_specific_covariate_record) {
      paras <- paras[!grepl(feature_specific_covariate_name, paras$name),]
    }
  }
  
  ################################
  # Apply correction to p-values #
  ################################
  
  paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
  
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





