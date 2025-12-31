n = 30

dat_sub = data.frame(expr = rnorm(n),
               g = sample(letters[1:3], size = n, replace = TRUE),
               y = factor(sample(LETTERS[1:3], size = n, replace = TRUE),
                          levels = LETTERS[1:3]),
               x = rnorm(n)
)

dat_sub$expr[dat_sub$y == "A"] = NA

formula = ~ x + y + (1|g)
random_effects_formula = expr ~ (1|g)
groups = NULL
ordereds = NULL

features = matrix(rnorm(n*3), nrow = n)
x = 1
feature_specific_covariate_name = NULL

# copied from the original function
expect_true(any(c(vapply(colnames(dat_sub), 
                         function(col) {
                             if (is.factor(dat_sub[, col])) {
                                 if (all(is.na(dat_sub$expr[dat_sub[, col] == 
                                                            levels(dat_sub[, col])[1]]))) {
                                     fixed_effects <-
                                         get_fixed_effects(formula,
                                                           random_effects_formula,
                                                           dat_sub,
                                                           groups,
                                                           ordereds,
                                                           feature_specific_covariate_name)
                                     if (col %in% substr(fixed_effects, 1, nchar(col))) {
                                         return(TRUE)
                                     }
                                 }
                             }
                             return(FALSE)
                         }, logical(1)))))
fxf <-
        get_fixed_effects(formula,
                          random_effects_formula,
                          dat_sub,
                          groups,
                          ordereds,
                          feature_specific_covariate_name)

# Check my simplified function:
expect_equal(any(c(vapply(colnames(dat_sub), 
                          function(col) {
                              if (is.factor(dat_sub[, col])) {
                                  if (all(is.na(dat_sub$expr[dat_sub[, col] == 
                                                             levels(dat_sub[, col])[1]]))) {
                                      fixed_effects <-
                                          get_fixed_effects(formula,
                                                            random_effects_formula,
                                                            dat_sub,
                                                            groups,
                                                            ordereds,
                                                            feature_specific_covariate_name)
                                      if (col %in% substr(fixed_effects, 1, nchar(col))) {
                                          return(TRUE)
                                      }
                                  }
                              }
                              return(FALSE)
                          }, logical(1)))),
             any(mapply(FUN = check_mffl_one, 
                        as.list(dat_sub),
                        colnames(dat_sub),
                        MoreArgs = list(ex = dat_sub$expr,
                                        fxf = fxf))))
             
goal = list(para = structure(list(coef = c(NA, NA, NA), 
                                  stderr = c(NA, NA, NA),
                                  pval = c(NA, NA, NA),
                                  name = c("x", "yB", "yC"), 
                                  error = c("No data points have the baseline factor level", 
                                            "No data points have the baseline factor level",
                                            "No data points have the baseline factor level")), 
                             row.names = c(NA, -3L), 
                             class = "data.frame"), 
            residuals = NA, 
            fitted = NA, 
            ranef = NA, 
            fit = NA)

expect_equal(goal,
             check_missing_first_factor_level(formula, random_effects_formula,
                                              dat_sub, groups, ordereds, features, x, 
                                              feature_specific_covariate_name))
